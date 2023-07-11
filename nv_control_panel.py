# -*- coding: utf-8 -*-
"""
This file contains functions to control the CFM. Just change the function call
in the main section at the bottom of this file and run the file. Shared or
frequently changed parameters are in the __main__ body and relatively static
parameters are in the function definitions.

Created on Sun Nov 25 14:00:28 2018

@author: mccambria

use file 5/5/2022
"""


# %% Imports

import labrad
import utils.tool_belt as tool_belt
import utils.positioning as positioning
import utils.common as common
import majorroutines.image_sample as image_sample
import majorroutines.optimize as optimize
import majorroutines.stationary_count as stationary_count
import majorroutines.resonance as resonance
import majorroutines.pulsed_resonance as pulsed_resonance
import majorroutines.optimize_magnet_angle as optimize_magnet_angle
import majorroutines.rabi as rabi
import majorroutines.ramsey as ramsey
import majorroutines.spin_echo as spin_echo
from utils.tool_belt import States, NormStyle 
import utils.auto_tracker as auto_tracker
import time
import numpy as np
import copy
import csv

# %% Major Routines

def do_auto_check_location(nv_sig=None,close_plot=False, haystack_fname = None):
    tool_belt.check_exp_lock()
    tool_belt.set_exp_lock()
    
    # set up a temperorary nv_sig
    green_power =10
    sample_name = "WiDQ"
    green_laser = "cobolt_515"
    nv_sig_basic =  {
        "coords":[5, 5, 5],  
        "name": "{}-autotracking".format(sample_name,),
        "expected_count_rate":None,        
        "imaging_laser":green_laser,
        "imaging_laser_power": green_power,
        "imaging_readout_dur": 1e7,
        "collection_filter": "630_lp",
        "magnet_angle": 0, 
        'norm_style':NormStyle.SINGLE_VALUED}  # 14.5 max. units is dBm

    # get the info of the haystack file and nv coord that we use for autotracking
    with labrad.connect() as cxn:
        haystack_fname = common.get_registry_entry(cxn,"haystack_fname" , ["", "Config", "AutoTracking"])
        haystack_coord_x_um = common.get_registry_entry(cxn,"haystack_nv_coord_x" , ["", "Config", "AutoTracking"])
        haystack_coord_y_um = common.get_registry_entry(cxn,"haystack_nv_coord_y" , ["", "Config", "AutoTracking"])
        haystack_coord_z_um = common.get_registry_entry(cxn,"haystack_nv_coord_z" , ["", "Config", "AutoTracking"])
    
    # The NV coords saved in the registry are in um, so convert back to V
    haystack_coord_x = haystack_coord_x_um/20
    haystack_coord_y = haystack_coord_y_um/20
    haystack_coord_z = haystack_coord_z_um/20
    nv_sig_basic['coords'] = [haystack_coord_x, haystack_coord_y, haystack_coord_z]
    
    # collect an image that is smaller than haystack image and has same resolution
    needle_fname = do_image_sample(nv_sig_basic,scan_size='needle',close_plot=close_plot, standalone_exp = False)
    
    #run the auto tracker image processing to locate needle image in haystack image
    x_shift, y_shift = auto_tracker.get_shift(nv_sig_basic, haystack_fname, needle_fname,close_plot=close_plot)
    
    # add the shift of the two images to the current drift
    with labrad.connect() as cxn:
        drift = positioning.get_drift(cxn)
        new_drift = [drift[0] - x_shift, drift[1] - y_shift, drift[2]]
        positioning.set_drift(cxn, new_drift)
        
    # With updated drift, optimize on NV to accurately find drift in all three dimensions
    nv_sig_copy = copy.deepcopy(nv_sig_basic)
    nv_sig_copy['expected_count_rate'] = None
    opti_coords, opti_count_rate = do_optimize(nv_sig_copy,plot_data=False,save_data = False, close_plot=close_plot,
                                               standalone_exp = False)      
        
    if opti_count_rate > 8:
        return
    else:
        raise RuntimeError('counts too low at opti coords')
     
    nv_sig['expected_count_rate'] = opti_count_rate
    
    tool_belt.set_exp_unlock()
    
def do_update_haystack_file(nv_sig):
    ''' This function takes an image wit hthe "haystack" setting, and will
    update the filename for the current haystack file that is used for autotracking.
    
    It is suggested that you perform this centered on an NV, right after optimizing on it.
    '''
    
    haystack_fname = do_image_sample(nv_sig, 'haystack')
    coords = nv_sig['coords']
    coords_um = [coords[0]*20, coords[1]*20, coords[2]*20]
    
    with labrad.connect() as cxn:
        p = cxn.registry()
        p.cd("", "Config", "AutoTracking")
        p.set("haystack_fname", haystack_fname)
        p.set("haystack_nv_coord_x", coords_um[0])
        p.set("haystack_nv_coord_y", coords_um[1])
        p.set("haystack_nv_coord_z", coords_um[2])
        
    admin_webpage = 'https://pub.physics.wisc.edu/qubit/admin/?s=edit_station&name=NV-Center+Station+1'
    print('\n***HAYSTACK FILE UPDATED***\nCenter coords of haystack image (in um): [{:.2f}, {:.2f}, {:.2f}]\nUpdate um coords at {}'.format(coords_um[0],
                                                                                         coords_um[1], coords_um[2], admin_webpage))
    
    

def do_image_sample(nv_sig, scan_size='medium', um_plot = False, close_plot=False, 
                    widqol = False, standalone_exp = True):
    scan_options=['huge','medium','big-ish','small','small-ish','needle','haystack','big','test','bigger-highres']
    if scan_size not in scan_options:
    #     raise Exception():
        print('scan_size must be in: ', scan_options)
        return 
    if scan_size == 'huge':
        scan_range = 1.5 # large scan
        num_steps = 151
    elif scan_size == 'big':
        scan_range = .75 # large scan
        num_steps = 76
    elif scan_size == 'bigger-highres':
        scan_range = 1.# large scan
        num_steps = 100
    elif scan_size == 'medium':
        scan_range = 0.4 # large scan
        num_steps = 81
    elif scan_size == 'small-ish':
        scan_range = 0.3 # large scan
        num_steps = 30
    elif scan_size == 'small':
        scan_range = 0.1 # large scan
        num_steps = 21
    elif scan_size == 'big-ish':
        scan_range = 0.8
        num_steps = 60
    elif scan_size == 'test':
        scan_range = .3
        num_steps = 10
    elif scan_size == 'needle':
        scan_range = 0.4 
        num_steps = 40
    elif scan_size == 'haystack':
        scan_range = 1.2 # large scan
        num_steps = 120
        
    # For now we only support square scans so pass scan_range twice
    fname = image_sample.main(nv_sig, scan_range, scan_range, num_steps,um_plot, 
                              close_plot=close_plot,
                              widqol = widqol,
                              standalone_exp = standalone_exp)
    return fname



def do_optimize(nv_sig,set_to_opti_coords=False,save_data=True,plot_data=True,close_plot=False, 
                standalone_exp = True):

    opti_coords, opti_count_rate = optimize.main(
        nv_sig,
        set_to_opti_coords,
        save_data,
        plot_data, 
        close_plot=close_plot,
        standalone_exp = standalone_exp
    )
    
    return opti_coords, opti_count_rate



def do_stationary_count(nv_sig):

    run_time = 3 * 60 * 10 ** 9  # ns

    stationary_count.main(nv_sig, run_time)



def do_resonance(nv_sig,  freq_center=2.87, freq_range=0.2, uwave_power=-5.0, 
                 num_steps = 101, num_runs = 40,close_plot=False, widqol = False):
    
    resonance.main(
        nv_sig,
        freq_center,
        freq_range,
        num_steps,
        num_runs,
        uwave_power,
        state=States.HIGH,
        opti_nv_sig = nv_sig,
        close_plot=close_plot,
        widqol = widqol
    )


def do_pulsed_resonance(nv_sig, freq_center=2.87, freq_range=0.2,num_runs=30,close_plot=False):

    num_steps = 51
    num_reps = 2e4
    # num_runs = runs
    uwave_power = 14
    uwave_pulse_dur = int(nv_sig["rabi_LOW"]/2)

    pulsed_resonance.main(
        nv_sig,
        freq_center,
        freq_range,
        num_steps,
        num_reps,
        num_runs,
        uwave_power,
        uwave_pulse_dur,
        opti_nv_sig = nv_sig,
        close_plot=close_plot
    )
    
def do_optimize_magnet_angle(nv_sig):
    num_angle_steps = 6
    angle_range = [0, 150]
    freq_center = 2.87
    freq_range= 0.2
    num_freq_steps = 75
    num_freq_runs = 5
    
    # Pulsed
    # uwave_power = 14
    # uwave_pulse_dur = 64/2
    # num_freq_reps = int(1e4)
    
    # CW
    uwave_power = -5
    uwave_pulse_dur = None
    num_freq_reps = None
    
    
    angle = optimize_magnet_angle.main(
        nv_sig,
        angle_range, 
        num_angle_steps,
        freq_center, 
        freq_range,
        num_freq_steps,
        num_freq_reps, 
        num_freq_runs,
        uwave_power, 
        uwave_pulse_dur)

    return angle
def do_rabi(nv_sig,  state, uwave_time_range=[0, 200], num_steps = 51, num_reps = 2e4, num_runs=20,close_plot=False, widqol = False):

    num_reps = int(num_reps)

    period = rabi.main(
        nv_sig,
        uwave_time_range,
        state,
        num_steps,
        num_reps,
        num_runs,
        opti_nv_sig = nv_sig,
        close_plot=close_plot,
        widqol = widqol
    )
    nv_sig["rabi_{}".format(state.name)] = period


def do_ramsey(nv_sig,  precession_time_range = [0, 0.2 * 10 ** 4], num_steps = 101,
                      set_detuning=0,num_reps = 2e4, num_runs=10, state=States.LOW,close_plot=False, do_fft  =True, widqol = False):

    detuning = set_detuning  # MHz
    # precession_time_range = [0, 1 * 10 ** 4]
    # precession_time_range = [0, .6 * 10 ** 3]
    # num_steps = 101
    num_reps = int(num_reps)
    # num_runs = runs

    ramsey.main(
        nv_sig,
        detuning,
        precession_time_range,
        num_steps,
        num_reps,
        num_runs,
        state,
        opti_nv_sig = nv_sig,
        close_plot=close_plot,
        do_fft = do_fft,
        widqol = widqol
    )


def do_spin_echo(nv_sig, echo_time_range = [0, 80 * 10 ** 3],num_steps = 81,
                 num_reps = 1e4, num_runs=40, state = States.LOW,close_plot=False, 
                 calc_theta_B = False, widqol = False):

    # T2* in nanodiamond NVs is just a couple us at 300 K
    # In bulk it's more like 100 us at 300 K

    num_reps = int(num_reps)


    # state = States.LOW

    angle = spin_echo.main(
        nv_sig,
        echo_time_range,
        num_steps,
        num_reps,
        num_runs,
        state,
        close_plot=close_plot,
        calc_theta_B = calc_theta_B,
        widqol = widqol
    )
    return angle

def get_drift():
    with labrad.connect() as cxn:
        drift = positioning.get_drift(cxn)
    return drift


def set_drift(drift_to_set):
    with labrad.connect() as cxn:
        positioning.set_drift(cxn,np.asarray(drift_to_set))

def reset_xy_drift():
    cur_drift = get_drift()
    
    with labrad.connect() as cxn:
        positioning.set_drift(cxn,np.array([0,0,cur_drift[2]]))

def reset_xyz_drift():
    
    with labrad.connect() as cxn:
        positioning.set_drift(cxn,np.array([0,0,0]))



# %% Run the file


if __name__ == '__main__':


    # %% Shared parameters
    

    green_power =10
    sample_name = "WiDQ"
    green_laser = "cobolt_515"
    
        
    nv_sig = {
        "coords":[4.832, 4.807, 5.09],  
        "name": "{}-nv1".format(sample_name,),
        "expected_count_rate":None,
        "disable_opt":False,
        "ramp_voltages": False,
        
        "spin_laser": green_laser,
        "spin_laser_power": green_power,
        "spin_pol_dur": 1e4,
        "spin_readout_laser_power": green_power,
        "spin_readout_dur": 350,
        
        "imaging_laser":green_laser,
        "imaging_laser_power": green_power,
        "imaging_readout_dur": 1e7,
        "collection_filter": "630_lp",
        
        "magnet_angle": 0, 
        "resonance_LOW":2.833 ,"rabi_LOW": 87.9, "uwave_power_LOW": 14,  # 15.5 max. units is dBm
        "resonance_HIGH": 2.904, "rabi_HIGH": 60, "uwave_power_HIGH": 14, 
        'norm_style':NormStyle.SINGLE_VALUED}  # 14.5 max. units is dBm
    
    nv_sig = nv_sig

    # %% Functions to run
    tool_belt.check_exp_lock()
    
    try:

        # reset_xy_drift()
        # reset_xyz_drift()
        # positioning.set_xyz (labrad.connect(), [5,5,5])
        
        # with labrad.connect() as cxn:
        #     tool_belt.laser_on(cxn, 'cobolt_515') # turn the laser on
            # tool_belt.laser_off(cxn, 'cobolt_515') # turn the laser on
        
        do_auto_check_location(nv_sig,close_plot=False)
        do_update_haystack_file(nv_sig)

        
        # do_optimize(nv_sig)
        # do_image_sample(nv_sig, scan_size='small')
        # do_image_sample(nv_sig,  scan_size='needle')
        # do_image_sample(nv_sig,  scan_size='medium', um_plot = False)
        # do_image_sample(nv_sig,  scan_size='haystack')
        # do_image_sample(nv_sig,  scan_size='big')
        # do_image_sample(nv_sig,  scan_size='small-ish')
        # do_image_sample(nv_sig,  scan_size='bigger-highres')
        
        # do_image_sample(nv_sig,  scan_size='test')
        
        # do_image_sample(nv_sig,  scan_size='auto-tracker')
        # do_image_sample(nv_sig, scan_size='huge')
        # z_list = np.arange(3.47,4.87,0.15)
        # for z in z_list:
        #     nv_sig['coords'][2]=z
        #     do_image_sample(nv_sig, scan_size='big-ish')
        # do_optimize(nv_sig)
        # nv_sig['disable_opt']=True
        # do_stationary_count(nv_sig, )
        
        # do_optimize_magnet_angle(nv_sig)
        # do_pulsed_resonance(nv_sig, freq_center=2.87, freq_range=0.2,num_runs=5)
        # mangles = [0,30,60,90,120,150]
        # for m in mangles:
        #     nv_sig['magnet_angle'] = m
        #     do_resonance(nv_sig, 2.87, 0.25, num_runs = 15)
        # nv_sig['disable_opt']=True
        # do_resonance(nv_sig, 2.87, 0.2,num_steps=51,num_runs=5)
        # do_resonance_state(nv_sig , States.LOW)
                
        do_rabi(nv_sig,  States.LOW, uwave_time_range=[0, 200],num_steps = 51, num_runs=5)
        # do_rabi(nv_sig,  States.HIGH, uwave_time_range=[0, 250],num_runs=30)
        
        # do_ramsey(nv_sig, set_detuning=0,num_runs=25, precession_time_range = [0, 1.75 * 10 ** 3],num_steps = 75)  
       
        # do_spin_echo(nv_sig,echo_time_range = [0, 110 * 10 ** 3], num_steps=71, num_runs=50, calc_theta_B = True) 
        pass
    finally:

        # Make sure everything is reset
        tool_belt.reset_cfm()
        tool_belt.set_exp_unlock()
        tool_belt.reset_safe_stop()
