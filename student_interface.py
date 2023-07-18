# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 15:18:03 2022

@author: Carter Fox

This will be the interface where students can run all the microscope commands/experiments. 
It will run nv_control_panel.py with the inputted parameters

"""
import utils.positioning as positioning
import utils.tool_belt as tool_belt
import utils.common as common
from utils.tool_belt import States, NormStyle 
import time
import numpy as np
import nv_control_panel as nv

# %%

if __name__ == "__main__":
    tool_belt.check_exp_lock()
    
    # %%%%%%%%%%%%%%% NV Parameters %%%%%%%%%%%%%%%
    
    nv_coords = [5.042, 5.042, 4.00]# V  #
    # expected_count_rate = None
    expected_count_rate = 18
    # kps
    # magnet_angle = 85 # deg
    magnet_angle = 153+(1/3) # deg
    
    resonance_LOW =  2.841     # GHz
    rabi_LOW = 105.7             # ns   
    uwave_power_LOW = 2    # dBm  15.5 max
    
    resonance_HIGH = 2.900     # GHz
    rabi_HIGH = 94.1            # ns 
    uwave_power_HIGH = 14.5     # dBm  14.5 max 
    
    #%%  Prepare nv_sig with nv parameters  (do not alter nv_sig)
    
    green_power = 10
    sample_name = "E6"
    green_laser = "cobolt_515"
    
        
    nv_sig = {
        "coords": nv_coords,
        
        "name": "{}-nv1".format(sample_name,),"disable_opt":False,"ramp_voltages": False,
        "spin_laser": green_laser, 
        "spin_laser_power": green_power,
        "spin_pol_dur": 1e4, 
        "spin_readout_laser_power": green_power, 
        "spin_readout_dur": 350,
        'norm_style':NormStyle.SINGLE_VALUED,
        
        "imaging_laser":green_laser, 
        "imaging_readout_dur": 1e7, "collection_filter": "630_lp",
        
        "expected_count_rate":expected_count_rate,
        "magnet_angle": magnet_angle, 
        
        "resonance_LOW":resonance_LOW ,"rabi_LOW": rabi_LOW, "uwave_power_LOW": uwave_power_LOW,  
        "resonance_HIGH": resonance_HIGH , "rabi_HIGH": rabi_HIGH, "uwave_power_HIGH": uwave_power_HIGH,
        }   
    
    nv_sig = nv_sig

    
    # %% %%%%%%%%%%%%%%% Experimental section %%%%%%%%%%%%%%%
    
    try:

        ####### Useful global functions #######
        ### Get/Set drift
        # nv.set_drift([0,0,0])
        # nv.reset_xy_drift() #Check that this is noted in lab manual
        # nv.reset_xyz_drift()
        # print(nv.get_drift())
        # nv_sig['disable_opt']=True
        # nv.do_stationary_count(nv_sig)
        
        ### Autotracking functions
        # nv.do_auto_check_location(nv_sig)
        # nv.do_update_haystack_file(nv_sig)
        
        ### Turn laser on 
        # tool_belt.laser_on('cobolt_515') # turn the laser
    
        ####### EXPERIMENT 0: Finding an nv #######
        ### Take confocal image
        ### xy scans can be ['small', 'medium', 'big-ish', 'big', 'huge']
        # nv.do_image_sample(nv_sig,  scan_size='small')  
        # nv.do_image_sample(nv_sig, scan_size='medium') 
        # nv.do_image_sample(nv_sig, scan_size='big')
        # nv.do_image_sample(nv_sig,  scan_size='big-ish')
        # nv.do_image_sample(nv_sig, scan_size='huge')
        # nv.do_image_sample(nv_sig, scan_size='needle')
                
        
        # Optimize on NV
        # nv.do_optimize(nv_sig)
            
        
        ####### EXPERIMENT 1: CW electron spin resonance #######
        ### Measure CW resonance
        # mangles = [0,30,60,90,120,150]
        # nv.do_resonance(nv_sig, freq_center=2.87, freq_range=0.2, uwave_power=-5.0, num_runs=15, num_steps=61)
    
        ####### EXPERIMENT 2: Rabi oscillations #######
        # mpowers = [-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,15]
        # for i in mpowers:
        #     nv_sig["uwave_power_LOW"]=i
        # nv.do_rabi(nv_sig,  States.LOW , uwave_time_range=[0, 200], num_runs=15, num_steps=51, num_reps=1e4)
         # nv.do_rabi(nv_sig,  States.HIGH, uwave_time_range=[0, 500], num_runs=20, num_steps=51, num_reps=2e4)
        
        
        ####### EXPERIMENT 3: Ramsey experiment #######
        # nv.do_ramsey(nv_sig, state=States.LOW, precession_time_range = [0, 2000], set_detuning=4, num_runs=2, num_steps = 101, num_reps=2e4)  

        # ####### EXPERIMENT 4: Spim echo #######
        # nv.do_spin_echo(nv_sig, state=States.LOW, echo_time_range = [0, 50000], 
        #                 num_runs=2, num_steps=41, num_reps=2e4) 
        
    finally:

        # Make sure everything is reset
        tool_belt.set_exp_unlock()
        tool_belt.reset_cfm()
        tool_belt.reset_safe_stop()

