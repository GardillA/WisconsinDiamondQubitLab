# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 11:16:50 2024

@author: student
"""
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
import labrad

# %%

if __name__ == "__main__":
    tool_belt.check_exp_lock()

    # %%%%%%%%%%%%%%% NV Parameters %%%%%%%%%%%%%%%

    nv_coords = [5, 5, 5]  # V  #
    # nv_coords = [4.92, 6.6, 1.6]  # V  4, 4.4, 53.977, 4.199, 3.915, 4.393, 1.60   
    # nv_coords = [5.103, 6.665, 1.60]
    # nv_coords = [4.841, 6.728, 1.60]
    # nv_coords = [4.586, 6.942, 1.60]
    # nv_coords = [4.651, 7.123, 1.60]
    # nv_coords = [4.649, 7.309, 1.60]    
    # nv_coords = [4.542, 7.457, 1.60]   
    # nv_coords = [2.445, 9.01, 1.60]
    # nv_coords = [3.57, 8.865, 1.60]      
    # nv_coords = [2.256, 8.314, 1.60]
    
    # nv_coords = [4.938, 6.566, 1.60]
    # nv_coords = [4.434, 7.501, 1.60]     
    # 
    # nv_coords = [3.705, 4.164, 1.12]
    # nv_coords = [3.824, 3.859, 2.00]
    
    expected_count_rate = None
    # kps
    magnet_angle = 85 # deg
    magnet_angle = 60  # deg   

    resonance_LOW = 2.831   # GHz
    rabi_LOW = 78.8            # ns
    uwave_power_LOW = 14    # dBm  15.5 max

    resonance_HIGH = 2.909     # GHz
    rabi_HIGH = 101.4            # ns
    uwave_power_HIGH = 14     # dBm  14.5 max

    # %%  Prepare nv_sig with nv parameters  (do not alter nv_sig)

    green_power = 10
    sample_name = "E6"
    green_laser = "cobolt_515"

    nv_sig = {
        "coords": nv_coords,

        "name": "{}-nv1".format(sample_name,), "disable_opt": False, "ramp_voltages": False,
        "spin_laser": green_laser,
        "spin_laser_power": green_power,
        "spin_pol_dur": 1e4,
        "spin_readout_laser_power": green_power,
        "spin_readout_dur": 350,
        'norm_style': NormStyle.SINGLE_VALUED,

        "imaging_laser": green_laser,
        "imaging_readout_dur": 1e7, "collection_filter": "630_lp",

        "expected_count_rate": expected_count_rate,
        "magnet_angle": magnet_angle,

        "resonance_LOW": resonance_LOW, "rabi_LOW": rabi_LOW, "uwave_power_LOW": uwave_power_LOW,
        "resonance_HIGH": resonance_HIGH, "rabi_HIGH": rabi_HIGH, "uwave_power_HIGH": uwave_power_HIGH,
    }

    # %% %%%%%%%%%%%%%%% Experimental section %%%%%%%%%%%%%%%

    try:
        
        ## MCC testing
        
        # with labrad.connect() as cxn:
        #     positioning.set_xyz_on_nv(cxn, nv_sig,drift_adjust=False)
        
        # tool_belt.init_safe_stop()
        # # for x in np.linspace(3, 7, 4):
        # #     for y in np.linspace(3, 7, 4):
        # for z in np.linspace(1.11, 1.15, 5):
        #     print(z)
        #     if tool_belt.safe_stop():
        #         break
        #     nv_sig["coords"][2] = z
        #     # nv_sig["coords"] = [x, y, z]
        #     nv.do_image_sample(nv_sig, scan_size='small')

        ####### Useful global functions #######
        # Get/Set drift
        # nv.set_drift([0, 0, 0])
        # nv.reset_xy_drift() #Check that this is noted in lab manual
        # nv.reset_xyz_drift()
        # print(nv.get_drift())
        # nv_sig['disable_opt']=True

        # nv.do_stationary_count(nv_sig)
  
        # Autotracking functions
        # nv.do_auto_check_location(nv_sig)
        # nv.do_update_haystack_file(nv_sig)

        # Turn laser on
        # tool_belt.laser_on('cobolt_515') # turn the laser
        

        ####### EXPERIMENT 0: Finding an nv #######
        # Take confocal image
        # xy scans can be ['small', 'medium', 'big-ish', 'big', 'huge']
        nv.do_image_sample(nv_sig,  scan_size='test')
        # nv.do_image_sample(nv_sig,  scan_size='small-ish')
        # nv.do_image_sample(nv_sig, scan_size='medium')
        # nv.do_image_sample(nv_sig, scan_size='big')
        # nv.do_image_sample(nv_sig,  scan_size='big-ish')
        # nv.do_image_sample(nv_sig, scan_size='huge')
        # nv.do_image_sample(nv_sig, scan_size='needle')
        # nv.do_image_sample(nv_sig,  scan_size='small')
        
        ###Ishita###
        ###Z scans###
        startz = 1.54
        stopz = 1.64
        stepz = .01
        nofzscans = (stopz - startz)/stepz
        # print('number of z scans=',nofzscans )
        def scanzs(startz, stopz, stepz, nv_coords, nv_sig, scan_size):
            zzz = np.arange(startz, stopz, stepz)
            nofzscans = (stopz - startz)/stepz
            print('number of z scans=',nofzscans )
            # zzz = np.arange(1, 2, 0.3)
            count = 1
            for i in zzz:
                print(count, ', z= ', i)                
                nv_coords[2] = i
                nv.do_image_sample(nv_sig,  scan_size)
                count +=1
        scan_size = 'test'
        # scanzs(startz, stopz, stepz, nv_coords, nv_sig, scan_size)
        
        ###XY scans sets###
        coordset = [[4.973, 5.488],[5.401, 4.441],[4.544, 4.584], [1.76, 5.607],[4.942, 5.103],[4.925, 5.298],[7.09, 2.657],[2.545, 8.177]]
        arrcoord = np.array(coordset)
        # print(len(arrcoord))
        def scanxys(arrcoord, nv_coords, nv_sig, scan_size):
            nofscans = len(arrcoord)
            print('number of scans=',nofscans )
            # zzz = np.arange(1, 2, 0.3)
            count = 1
            for i in range(0, nofscans):
                print('scan number = ', count, ', x y= ', arrcoord[i,:] )     
                # nv_coords = arrcoord[i,:]
                nv_coords[0] = arrcoord[i,0]             
                nv_coords[1] = arrcoord[i,1]                
                nv.do_image_sample(nv_sig,  scan_size)
                count +=1
        scan_size = 'test'
        # scanxys(arrcoord, nv_coords, nv_sig, scan_size)
# 
        # Optimize on NV
        # nv.do_optimize(nv_sig)

        ####### EXPERIMENT 1: CW electron spin resonance #######
        # Measure CW resonance
        # mangles = [0,30,60,90,120,150]
        # nv.do_resonance(nv_sig, freq_center=2.87, freq_range=0.2, uwave_power=-15.0, num_runs=25, num_steps=101)

        ####### EXPERIMENT 2: Rabi oscillations #######
        # mpowers = [-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,15]
        # for i in mpowers:
        #     nv_sig["uwave_power_LOW"]=i
        # nv.do_rabi(nv_sig,  States.LOW , uwave_time_range=[0, 200], num_runs=30, num_steps=101, num_reps=1e4)
        # nv.do_rabi(nv_sig,  States.HIGH, uwave_time_range=[0, 200], num_runs=15, num_steps=51, num_reps=1e4)

        ####### EXPERIMENT 3: Ramsey experiment #######
        # nv.do_ramsey(nv_sig, state=States.LOW, precession_time_range = [0, 6000], set_detuning=4, num_runs=80, num_steps = 76, num_reps=1e4)

        # ####### EXPERIMENT 4: Spim echo #######

        # nv.do_spin_echo(nv_sig, state=States.LOW, echo_time_range = [0, 40000],num_runs=50, num_steps=51, num_reps=1e4)

    finally:

        # Make sure everything is reset
        tool_belt.set_exp_unlock()
        tool_belt.reset_cfm()
        tool_belt.reset_safe_stop()
