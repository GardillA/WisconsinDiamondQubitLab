% This script demonstrates changes needed for making code written for old 
% API v0.9 work with new API v1.1.
clear all

% All PulseStreamer API functions and classes are placed in the
% "PulseStreamer" namespace. In order to shorter the names we import the
% namespace contents to the global namespace using "import" command. Then
% instead of "PulseStreamer.Pulsestreamer('ip')" we can write shorter
% version "PulseStreamer('ip')" 
import PulseStreamer.*

% Compatibility functions and classes can be found in the
% "PulseStreamer.compat" namespace so we import them too.
import PulseStreamer.compat.*

ip = 'pulsestreamer';

% connect to the pulse streamer
ps = PulseStreamer(ip);

%% basic settings
% MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW %
% OutputState class now accepts list of channels which are to be set high
% instead of bitmask.
%OLD:
% Set all channels to 0
% outputZero = OutputState(0,0,0); %0 or -1 means that the sequence is repeated until the power is turned off
%NEW: 
% Set all channels to 0
outputZero = OutputState([],0,0); % or OutputState.ZERO()

% initialOutputState = outputZero;   %OLD
finalOutputState = outputZero;
% underflowOutputState = outputZero; %OLD

% MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW %
%OLD: start = PSStart.Immediate;
start = TriggerStart.IMMEDIATE;
ps.setTrigger(start, TriggerRearm.MANUAL);

% settings for sequence generation
numberOfSequences = 10;
pulsesPerSequence = 100;
n_runs = 100;

disp(['Test performance for ' num2str(numberOfSequences * pulsesPerSequence * n_runs) ' of pulses in total.']);

fprintf('\n')

%% Generate sequences
%
% MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW %
%
% The way how sequences are generated has changed in the v1.0. P and PH
% classes are deprecated. For easier migration we provide modified 
% versions of P and PH classes and a function that converts P and PH
% objects to a Sequence object. Further use of P and PH classes is 
% discouraged as they will be removed in the future. 

disp(['Generating ' num2str(numberOfSequences) ' sequences with ' num2str(pulsesPerSequence) ' PHs each, that means in total ' num2str(numberOfSequences * pulsesPerSequence) ' PHs.']);
tic 
% we first create "numberOfSequences" different sequence groups (S1, S2, ...)"
sequences = cell(1,numberOfSequences);        
for iSeq=1:numberOfSequences
    sequences{iSeq} = [];        
    for iPulse = 1:pulsesPerSequence
        % the content of the PHs is more or less arbitary
        sequences{iSeq} = sequences{iSeq} + PH(1000, mod(iPulse * iSeq, 256), 0, 0);
    end
end
toc

fprintf('\n')
% case one - output one sequences after another and loop this "n_runs" times
% (S1, S2, ...) * n_runs
disp(['a) Output one sequences after another and loop this ' num2str(n_runs) ' times.']);
disp(['Total number of pulses: ' num2str(numberOfSequences * pulsesPerSequence * n_runs)]);
disp(['The number of times all sequences are repeated are passed to the "stream" method, which is the most efficient way for repeating the whole sequence.']);
tic
runs = n_runs;
pSeq = [];
for iSeq=1:numberOfSequences
    pSeq = pSeq + sequences{iSeq};
end

% MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW %
% OLD: 
%   ps.stream(pSeq, runs, initialOutputState, finalOutputState, underflowOutputState, start);
% NEW:
% method "stream" has new signature: stream(sequence, n_runs, finalState)
% where "sequence" MUST BE an object of class "Sequence".
%
% You can convert an array of P/PH objects created above into 
% Sequence object using the following compatibility function.

pSeqNew = convert_PPH_to_PSSequence(pSeq); % This is REQUIRED! % MIGRATION: CONVERSION FUNCTION

ps.stream(pSeqNew, runs, finalOutputState); %NEW
toc

%%
% case two - output each sequence "n_runs" times. Then continue with the
% next sequence for "n_runs" times and so on...
% (S1*n_runs, S2*n_runs, ...)
fprintf('\n')
disp(['b) Output each sequence ' num2str(n_runs) ' times and continue with the next sequence the same way.']);
disp(['Total number of pulses: ' num2str(numberOfSequences * pulsesPerSequence * n_runs)]);
disp(['The repetition is build up with "*" which is slower compared to the method of a).'])
disp(['The advantage is that not only the whole sequence in total can be repeated, but also sub-sequences as shown here.']);
tic
pSeq = [];
for iSeq=1:numberOfSequences
    pSeq = pSeq + sequences{iSeq} * n_runs;
end
%
runs = 1;

% MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW % MIGRATION: NEW %
% OLD: 
%   ps.stream(pSeq, runs, initialOutputState, finalOutputState, underflowOutputState, start);
% NEW:
% method "stream" has new signature: stream(sequence, n_runs, finalState)
% where "sequence" MUST BE an object of class "Sequence".
%
% You can convert an array of P/PH objects created above into 
% Sequence object using the following compatibility function.

pSeqNew = convert_PPH_to_PSSequence(pSeq); % This is REQUIRED! % MIGRATION: CONVERSION FUNCTION

ps.stream(pSeqNew, runs, finalOutputState); %NEW
toc

