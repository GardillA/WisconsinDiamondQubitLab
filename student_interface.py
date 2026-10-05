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
    #4.251, 4.177
    nv_coords = [3.19, 2.154, 5.51]
  
    
    expected_count_rate = None
 
    magnet_angle = 90 #30 #deg

    resonance_LOW = 2.8108#2.8179   # GHz
    rabi_LOW = 152.9       # 120.30  # ns
    uwave_power_LOW = -5      # dBm  15.5 max

    resonance_HIGH = 2.9317 #2.9268  # GHz
    rabi_HIGH = 211.6      # 126.11  # ns
    uwave_power_HIGH = -4     # dBm  14.5 max
    # %%  Prepare nv_sig with nv parameters  (do not alter nv_sig)
    green_power = 15 # mW
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
        # # for x in np.linspace(2, 8, 3):
        # #     for y in np.linspace(2, 8, 3):
        # for z in np.linspace(2, 8, 5):
        #     print(z)
        #     if tool_belt.safe_stop():
        #         break
        #     # nv_sig["coords"][0] = x
        #     # nv_sig["coords"][1] = y
        #     nv_sig["coords"][2] = z
        #     # nv_sig["coords"] = [x, y, z]
        #     nv.do_image_sample(nv_sig, scan_size='big')

        ####### Useful global functions #######
        # Get/Set drift
        nv.set_drift([0, 0, 0])
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
        
        #### Imaging for Alignment #######
        # zsteps = np.flip(np.linspace(4.0,5.0,5))
        # for n in zsteps:
        #     nv_coords = [4.25, 4.25, n]
        #     nv_sig = {
        #             "coords": nv_coords,
    
        #             "name": "{}-nv1".format(sample_name,), "disable_opt": False, "ramp_voltages": False,
        #             "spin_laser": green_laser,
        #             "spin_laser_power": green_power,
        #             "spin_pol_dur": 1e4,
        #             "spin_readout_laser_power": green_power,
        #             "spin_readout_dur": 350,
        #             'norm_style': NormStyle.SINGLE_VALUED,
    
        #             "imaging_laser": green_laser,
        #             "imaging_readout_dur": 1e7, "collection_filter": "630_lp",
    
        #             "expected_count_rate": expected_count_rate,
        #             "magnet_angle": magnet_angle,
    
        #             "resonance_LOW": resonance_LOW, "rabi_LOW": rabi_LOW, "uwave_power_LOW": uwave_power_LOW,
        #             "resonance_HIGH": resonance_HIGH, "rabi_HIGH": rabi_HIGH, "uwave_power_HIGH": uwave_power_HIGH,
        #         }
        # nv.set_drift([0.003, -0.004, -0.02])
        # nv.do_image_sample(nv_sig, scan_size='big')
        
            
        # ####### EXPERIMENT 0: Finding an nv #######
        # # Take confocal image
        # # xy scans can be ['small', 'medium', 'big-ish', 'big', 'huge']
        # nv.do_image_sample(nv_sig,  scan_size='test')
        # nv.do_image_sample(nv_sig,  scan_size='small')
        # 
        # nv.do_image_sample(nv_sig, scan_size='medium')
        #nv.do_image_sample(nv_sig, scan_size='big')
        # nv.do_image_sample(nv_sig,  scan_size='big-ish')
        # nv.do_image_sample(nv_sig, scan_size='huge')
        # nv.do_image_sample(nv_sig, scan_size='needle')
        # nv.do_image_sample(nv_sig, scan_size='haystack')
        
        for z in [2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8, 5.2]:#[2, 2.5, 3, 3.5, 4, 4.5, 5]::
            nv_coords[2] = z
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
            nv.do_image_sample(nv_sig,  scan_size='big-ish')



        # Optimize on NV
        # nv.do_optimize(nv_sig)

        ####### EXPERIMENT 1: CW electron spin resonance #######
        # Measure CW resonance

        # nv.do_resonance(nv_sig, freq_center=2.87, freq_range=0.2, uwave_power=-15, num_runs=25, num_steps=101)
        
        # A loop to check multiple magnet angles
        # mag_list = [0, 30, 60, 90, 120, 150, 180]
        # for n in mag_list:

        #     mag_ang = n
            
            
        #     nv_sig_temp = {
        #         "coords": nv_coords,
    
        #         "name": "{}_{}-nv1".format(sample_name,magnet_angle), "disable_opt": False, "ramp_voltages": False,
        #         "spin_laser": green_laser,
        #         "spin_laser_power": green_power,
        #         "spin_pol_dur": 1e4,
        #         "spin_readout_laser_power": green_power,
        #         "spin_readout_dur": 350,
        #         'norm_style': NormStyle.SINGLE_VALUED,
    
        #         "imaging_laser": green_laser,
        #         "imaging_readout_dur": 1e7, "collection_filter": "630_lp",
    
        #         "expected_count_rate": expected_count_rate,
        #         "magnet_angle": mag_ang,
    
        #         "resonance_LOW": resonance_LOW, "rabi_LOW": rabi_LOW, "uwave_power_LOW": uwave_power_LOW,
        #         "resonance_HIGH": resonance_HIGH, "rabi_HIGH": rabi_HIGH, "uwave_power_HIGH": uwave_power_HIGH}
        #     print("Angle: ", n)    
        #     nv.do_resonance(nv_sig_temp, freq_center=2.87, freq_range=0.2, uwave_power=-15.0, num_runs=25, num_steps=101)

        ####### EXPERIMENT 2: Rabi oscillations #######
        
        ##Note: make sure angle is set to corresponding LOW and HIGH microwave frequency pair
        
        # mpowers = [15, 4,-4,-8,-12]
        # for i in mpowers:
        #     nv_sig["uwave_power_LOW"]=i
        #     nv.do_rabi(nv_sig,  States.LOW , uwave_time_range=[0, 200], num_runs=15, num_steps=51, num_reps=1e4)
        
        # nv.do_rabi(nv_sig,  States.LOW, uwave_time_range=[0, 200], num_runs=15, num_steps=51, num_reps=1e4)
        # nv.do_rabi(nv_sig,  States.HIGH, uwave_time_range=[0, 200], num_runs=5, num_steps=51, num_reps=1e4)
        # print("mpowers done")
        # powers = [11, 0, -11]
        # for p in powers:
        #       print("Power:", p)
        #       nv_sig["uwave_power_LOW"] = p
        #       nv_sig["uwave_power_HIGH"] = p
        #       nv.do_rabi(nv_sig,  States.LOW, uwave_time_range=[0, 200], num_runs=15, num_steps=51, num_reps=1e4)

        # for p in powers:
        #         print("Power:", p)
        #         nv_sig["uwave_power_LOW"] = p
        #         nv_sig["uwave_power_HIGH"] = p
        #         nv.do_rabi(nv_sig,  States.HIGH, uwave_time_range=[0, 200], num_runs=15, num_steps=51, num_reps=1e4)


        ####### EXPERIMENT 3: Ramsey experiment #######
        # Run whichever has lower period
        # nv.do_ramsey(nv_sig, state=States.LOW, precession_time_range = [0, 8000], set_detuning=4, num_runs=50, num_steps = 101, num_reps=1e4)
        # nv.do_ramsey(nv_sig, state=States.HIGH, precession_time_range = [0, 8000], set_detuning=4, num_runs=50, num_steps = 101, num_reps=1e4)

        # ####### EXPERIMENT 4: Spin echo #######

        # nv.do_spin_echo(nv_sig, state=States.LOW, echo_time_range = [0, 150000],num_runs=50, num_steps=151, num_reps=1e4)
        # nv.do_spin_echo(nv_sig, state=States.HIGH, echo_time_range = [0, 150000],num_runs=50, num_steps=151, num_reps=1e4)

        # nv.do_ramsey(nv_sig, state=States.LOW, precession_time_range = [0, 8000], set_detuning=4, num_runs=50, num_steps = 101, num_reps=1e4)
        #Scott's Benchmark tests for daily drift
        # nv.do_image_sample(nv_sig,  scan_size='small')
        #nv.do_image_sample(nv_sig, scan_size='medium')
        # nv.do_image_sample(nv_sig, scan_size='huge')
        
        # mag_list = [30,60,90,120,150,180]
        # for n in mag_list:
        #     magnet_angle = n
        #     nv_sig = {
        #         "coords": nv_coords,
    
        #         "name": "{}-nv1".format(sample_name,), "disable_opt": False, "ramp_voltages": False,
        #         "spin_laser": green_laser,
        #         "spin_laser_power": green_power,
        #         "spin_pol_dur": 1e4,
        #         "spin_readout_laser_power": green_power,
        #         "spin_readout_dur": 350,
        #         'norm_style': NormStyle.SINGLE_VALUED,
    
        #         "imaging_laser": green_laser,
        #         "imaging_readout_dur": 1e7, "collection_filter": "630_lp",
    
        #         "expected_count_rate": expected_count_rate,
        #         "magnet_angle": magnet_angle,
    
        #         "resonance_LOW": resonance_LOW, "rabi_LOW": rabi_LOW, "uwave_power_LOW": uwave_power_LOW,
        #         "resonance_HIGH": resonance_HIGH, "rabi_HIGH": rabi_HIGH, "uwave_power_HIGH": uwave_power_HIGH}
        
        #nv.do_resonance(nv_sig, freq_center=2.87, freq_range=0.2, uwave_power=-15.0, num_runs=25, num_steps=101)

        
        # nv.do_resonance(nv_sig, freq_center=2.87, freq_range=0.2, uwave_power=-15.0, num_runs=15, num_steps=51)
        
        # nv.do_rabi(nv_sig,  States.HIGH , uwave_time_range=[0, 200], num_runs=15, num_steps=51, num_reps=1e4)
        # nv.do_rabi(nv_sig,  States.LOW, uwave_time_range=[0, 200], num_runs=15, num_steps=51, num_reps=1e4)

        #nv.do_spin_echo(nv_sig, state=States.LOW, echo_time_range = [0, 150000],num_runs=50, num_steps=151, num_reps=1e4)
        
    finally:
        # Make sure everything is reset
        #tool_belt.set_exp_unlock()
        # tool_belt.reset_cfm()
        #tool_belt.reset_safe_stop()
        # tool_belt.reset_xy_drift()
        pass
    