# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:30:52 2024

@author: jmbaer

Testing I/O with Swabian 8/2 Pulse Streamer on unused channels


"""

from pulsestreamer import Sequence, PulseStreamer, findPulseStreamers, TriggerStart, TriggerRearm

'''
load the pulse pattern as a list of tuples

params:
    [(time in ns, high=1 or low=0),...(),...()]

'''
pattern = [(50000, 1), (10000, 0)]

'''
static IP from labrad config
'''
ip = '10.130.33.247'
ps = PulseStreamer(ip)


'''
create an instance of a sequence, load the pattern to channel 3
and create a stream for 1000 repititions of the pattern
'''
sequence = ps.createSequence()
sequence.setDigital(4,pattern)
ps.stream(sequence,1000)

# %%
'''
triggering with software

if the start is SOFTWARE then you need to initiate the pulse with startNow()
'''
# ps.setTrigger(start = TriggerStart.SOFTWARE)
# ps.startNow()

# %%
''' 
if the start is HARDWARE then the swabian is auto set to start the pulse train
when it gets the trigger input. the swabian that we have is hardware 2.3 so
the trigger has to be above 2.0 volts. in the documentation make sure that
you are looking at the correct hardware version for restraints
'''
ps.setTrigger(start = TriggerStart.HARDWARE_RISING_AND_FALLING,rearm=TriggerRearm.AUTO)


# %% troubleshooting
print(f'hasSequence = {ps.hasSequence()}\n')
print(f'isStreaming = {ps.isStreaming()}\n')
print(f'hasFinished = {ps.hasFinished()}\n')
# print(f'triggerStart = {ps.getTriggerStart()}\n')
