# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 15:18:03 2022

@author: Carter Fox and Dan Bradley

This will be the interface where the web user interface can run all the microscope commands/experiments. 
It will run nv_control_panel with the inputted parameters

"""

# avoid login prompts
import os
os.environ['LABRADUSER'] = ""
os.environ['LABRADPASSWORD'] = ""

import nv_control_panel as nv
import utils.tool_belt as tool_belt
from utils.tool_belt import States, NormStyle
import labrad
import time
import sys
import argparse
from pathlib import Path
import copy
# from utils.tool_belt import reset_xy_drift as reset_xy_drift
# import utils.tool_belt.reset_drift as reset_xyz_drift
# import utils.tool_belt.laser_on_no_cxn as laser_on
# import utils.tool_belt.laser_off_no_cxn as laser_off

# %%
# Routines that require mw args
routines_mw = ["ESR", "rabi", "ramsey", "spin-echo"]

#%%

def cwd_get_file_path(file,timestamp,fname,subfolder=None):
    """write output files into the current working directory
       rather than the default location
    """
    return Path(os.getcwd()) / fname

def cwd_get_raw_data_path(
    file_name,
    path_from_nvdata=None,
    nvdata_dir=None,
):
    """retrieve files from the current working directory if they exist there"""
    file_name_ext = "{}.txt".format(file_name)
    if os.path.exists(file_name_ext):
        return Path(file_name_ext)
    return tool_belt.orig_get_raw_data_path(file_name,path_from_nvdata,nvdata_dir)

def safeFileName(fname):
    fname = fname.replace("/","_")
    fname = fname.replace("\\","_")
    fname = fname.replace(" ","_")
    return fname

if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog = sys.argv[0],description = 'Run nvcenter experiments executed by the web worker process.')
    parser.add_argument('--experiment-type',action='store',required=True)
    parser.add_argument('--name',action='store')
    parser.add_argument('--test-frequency-min',action='store',type=float)
    parser.add_argument('--test-frequency-max',action='store',type=float)
    parser.add_argument('--num-test-points',action='store',type=int) 
    parser.add_argument('--num-test-averages',action='store',type=int) 
    parser.add_argument('--image-size',action='store',type=str)
    parser.add_argument('--set-frequency',action='store',type=float)
    parser.add_argument('--test-uwave-time-max',action='store',type=float) 
    parser.add_argument('--test-precession-time-max',action='store',type=float) 
    parser.add_argument('--test-echo-time-max',action='store',type=float)
    parser.add_argument('--x',action='store',type=float,required=True)
    parser.add_argument('--y',action='store',type=float,required=True)
    parser.add_argument('--z',action='store',type=float,required=True)
    parser.add_argument('--set-magnet_angle', action='store',type=float)  
    parser.add_argument('--set-uwave-power', action='store',type=float) 
    parser.add_argument('--set-uwave-frequency', action='store',type=float)  
    parser.add_argument('--set-pi-pulse-period', action='store',type=float)  
    args = parser.parse_args()

    tool_belt.get_file_path = cwd_get_file_path
    tool_belt.orig_get_raw_data_path = tool_belt.get_raw_data_path
    tool_belt.get_raw_data_path = cwd_get_raw_data_path

    # %%%%%%%%%%%%%%% NV Parameters %%%%%%%%%%%%%%%
    
    sample_name = args.name
    if not sample_name:
        sample_name = args.experiment_type
    sample_name = safeFileName(sample_name)

    nv_coords = [args.x,args.y,args.z] # V
    expected_count_rate = 16   # kcps
    
    
    #%%  Prepare nv_sig with nv parameters  (do not alter nv_sig)

    green_power = 10
    green_laser = "cobolt_515"

    nv_sig = { 
        "coords":nv_coords,
        "name": sample_name,
        "disable_opt":False, "ramp_voltages": False,
        "expected_count_rate":expected_count_rate,
        
        "spin_laser": green_laser, "spin_laser_power": green_power, 
        "spin_readout_laser_power": green_power,
        "spin_pol_dur": 1e4, "spin_readout_dur": 350,
        "imaging_laser":green_laser, "imaging_laser_power": green_power,
        "imaging_readout_dur": 1e7,
        'norm_style':NormStyle.SINGLE_VALUED,

        "collection_filter": "630_lp",
        "magnet_angle": 0,
        "resonance_LOW":2.87,"rabi_LOW": 60,
        "uwave_power_LOW": 14, 
        "resonance_HIGH": 2.87,
        "rabi_HIGH": 60, 
        "uwave_power_HIGH": 14 }  
    
    
    # %% %%%%%%%%%%%%%%% Experimental section %%%%%%%%%%%%%%%
    
    try:

        ####### Useful global functions #######
        ### Reset drift
        # tool_belt.reset_xy_drift()
        # print(tool_belt.get_drift())
        # tool_belt.set_drift([0.0, 0.0, 0.0]) 
        
        # Perform optimize, if counts are below a certain value, run autotracker to try to find the NV.
        _, opti_count_rate = nv.do_optimize(nv_sig, close_plot=True)
        if opti_count_rate < 8:
            nv.do_auto_check_location(nv_sig,close_plot=True)
        
        
        if args.experiment_type == 'auto-tracker':
            #haystack_fname='2023_05_31-11_21_50-WiQD-nv1_XY' add the file into this file
            nv.do_auto_check_location(nv_sig)
        
        elif args.experiment_type == "image":
            fname = nv.do_image_sample(nv_sig, scan_size=args.image_size, um_plot = True, close_plot=True)
            # print("Image fname: ",fname)
        elif args.experiment_type == "optimize":
            nv.do_optimize(nv_sig, close_plot=True)
 
        elif args.experiment_type in routines_mw:
            
            # for the experiments that use MWs, make a copy of the nv_sig and put in the arguments
            nv_sig_run = copy.deepcopy(nv_sig)
            
            # for all these routines, the magnet angle will need to be set. 
            #And we will set the state to LOW automatically
            nv_sig_run['magnet_angle'] = args.set_magnet_angle
            state_input = States.LOW
            
            #For MW experiments other than ESR, they will need the MW power 
            #and the frequency
            if args.experiment_type != "ESR":
                set_uwave_power_dbm = tool_belt.mW_to_dBm(args.set_uwave_power)
                # nv_sig_run = copy.deepcopy(nv_sig)
                nv_sig_run['resonance_LOW'] = args.set_uwave_frequency
                nv_sig_run['uwave_power_LOW'] = set_uwave_power_dbm
                
                #If the MW experiment is not ESR or Rabi, they will need the pi pulse time
                if args.experiment_type != "rabi":
                    set_pi_pulse = args.set_pi_pulse_period
                    #our code expects the rabi period (or twice the pi pulse duration), so we will provide that
                    nv_sig_run['rabi_LOW'] = 2*set_pi_pulse
            
            
            if args.experiment_type == "ESR":
                # take the min and max frequency inputs, and convert to the 
                # center freq and the range, which is what our code uses
                freq_range_ = abs(args.test_frequency_max - args.test_frequency_min)
                # make sure we add the half-range to the lower frequency input
                lower_freq_min = min(args.test_frequency_max, args.test_frequency_min)
                freq_center_ = freq_range_/2 + lower_freq_min
                
                nv.do_resonance(nv_sig_run, freq_center=freq_center_, freq_range=freq_range_, 
                                num_steps=args.num_test_points, num_runs=args.num_test_averages, close_plot=True)
            
            elif args.experiment_type == "rabi":
                nv.do_rabi(nv_sig_run,state=state_input,uwave_time_range=[0,args.test_uwave_time_max], 
                           num_steps=args.num_test_points, num_runs=args.num_test_averages, close_plot=True)
            
            elif args.experiment_type == "ramsey":
                nv.do_ramsey(nv_sig_run,state=state_input,set_detuning=0,precession_time_range=[0,1000*args.test_precession_time_max], 
                             num_steps=args.num_test_points, num_runs=args.num_test_averages, close_plot=True)
            
            elif args.experiment_type == "spin-echo":
                # users will input the total precession time, and our code expects
                # the time to be half the total precession time, T = 2*tau.
                
                # So divide the user's input by two to use as the input to our function
                
                spin_echo_tau_max = args.test_echo_time_max / 2
                nv.do_spin_echo(nv_sig_run,state=state_input,echo_time_range=[0,spin_echo_tau_max], 
                                num_steps=args.num_test_points, num_runs=args.num_test_averages, close_plot=True)
        
        else:
            raise Exception("Unsupported experiment type: " + repr(args.experiment_type))
    
    except Exception as exc:
        print("Code crashed. Press enter to see error")
        raise exc
    
    finally:
        # Reset our hardware - this should be done in each routine, but
        # let's double check here
        tool_belt.reset_cfm()
        tool_belt.reset_safe_stop()
