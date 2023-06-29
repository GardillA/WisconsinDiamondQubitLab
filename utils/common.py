# -*- coding: utf-8 -*-
"""
Functions, etc to be referenced only by other utils. If you're running into
a circular reference in utils, put the function or whatever here. 

Created September 10th, 2021

@author: mccambria
"""

import platform
from pathlib import Path
from utils.tool_belt import ExpLock 
### Lab-specific stuff here

shared_email = "kolkowitznvlab@gmail.com"
windows_nvdata_dir = Path("E:/Shared drives/Kolkowitz Lab Group/nvdata")
linux_nvdata_dir = Path.home() / "E/nvdata"

###

def check_exp_lock(cxn, timestamp=None, runtime=None): #potentially pst the runtime and the timestamp?
    '''
    When an experiment runs, it will update a variable in the labrad registry.
    
    This will check what that entry is, and if it is locked, it will throw an error
    '''
    locked_state = ExpLock.LOCK
    lock_val = get_registry_entry(cxn, 'ExperimentalLock', ("Config", "ExperimentalLock"))
    
    if lock_val == locked_state.value:
        raise Exception('Experiment is currently being used, please try to run again later')
        
    return lock_val

def set_exp_lock(cxn):
    '''
    When an experiment runs, it will update a variable in the labrad registry.
    
    This will set the variable to the "lock" state an not allow other users to run experiments on the same set up
    '''
    
    locked_state = ExpLock.LOCK
    p = cxn.registry()
    p.cd("", "Config", "ExperimentalLock")
    p.set("ExperimentalLock", locked_state.value)
    # return value
    
def set_exp_unlock(cxn):
    '''
    When an experiment runs, it will update a variable in the labrad registry.
    
    This will set the variable to an unlock state.
    '''
    
    unlocked_state = ExpLock.UNLOCK
    p = cxn.registry()
    p.cd("", "Config", "ExperimentalLock")
    p.set("ExperimentalLock", unlocked_state.value)
    # return value
    
def get_nvdata_dir():
    """Returns an OS-dependent Path to the nvdata directory (configured above)"""
    
    check_if_instructional_lab_pc_name = True
    
    if check_if_instructional_lab_pc_name != None:
        nvdata_dir = Path("C:/Users/student/Documents/LAB_DATA")
    else:
        os_name = platform.system()
        if os_name == "Windows":
            nvdata_dir = windows_nvdata_dir
        elif os_name == "Linux":
            nvdata_dir = linux_nvdata_dir

    return nvdata_dir


def get_registry_entry(cxn, key, directory):
    """Return the value for the specified key. Directory as a list,
    where an empty string indicates the top of the registry
    """

    p = cxn.registry.packet()
    p.cd("", *directory)
    p.get(key)
    return p.send()["get"]


def get_server(cxn, server_type):
    """Helper function for server getters in tool_belt. Return None if we can't
    make the connection for whatever reason (e.g. the key does not exist in 
    the registry
    """
    try:
        server_name = get_registry_entry(cxn, server_type, ["", "Config", "Servers"])
        server = getattr(cxn, server_name)
    except Exception as exc:
        # print(f"Could not get server type {server_type}")
        server = None
    return server
