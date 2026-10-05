%% Example 4: Creating complex sequences
% This example provides step by step guidance on how very complex sequences
% can be created using Pulse Streamer API classes.
% We start from creating signal patterns as simple elements of the desired 
% sequence. Then, we assign these patterns to channels and build sub-sequences. 
% Finally, we use these elements to compose one complex sequence.

clear all;

% All PulseStreamer API functions and classes are placed in the
% "PulseStreamer" namespace. In order to shorter the names we import the
% namespace contents to the global namespace using "import" command. Then
% instead of "PulseStreamer.Pulsestreamer('ip')" we can write shorter
% version "PulseStreamer('ip')" 
import PulseStreamer.*

% DHCP is activated in factory settings
% Hostname or IP address of the pulse streamer
ip = 'pulsestreamer';

% Connect to the Pulse Streamer
ps = PulseStreamer(ip);

% Make sure that Pulse Streamer is in default state
ps.reset();

%% Create necessary signal patterns

% A pulse of 100us period (50us high, 50us low)
t_p1 = 100000; % ns
pattern1 = {fix(t_p1/2), true; fix(t_p1/2), false};

% Here we produce pulses with 1us period such that the total duration 
% is the same as high level of the pattern1
t_p2 = 1000; % ns, period
nPulses2 = ceil(fix(t_p1/2)/t_p2); % number of pulses that fit into high level of pattern1
pattern2 = cell(nPulses2, 2);
for ii = 1:nPulses2
    state = mod(ii,2); % alternating high/low starting from high
    pattern2(ii,:) = {fix(t_p2), state}; 
end

% Pulse of high level for 8us pulse followed by low level for 2us  
pattern3 = {8000,true; 2000, false};

% One period of a pattern that we want to be repeated over duration of the
% final sequence. 
pattern_clk = {25000,true; 25000, false};

%% Sub-sequence 1. Start marker. Pattern3
% Short pulse on channel 0 followed by long low level while 
% all other channels kept low.

% Initialize sequence
subseq1 = PulseStreamer.Sequence();

% Assign pattern to channel
subseq1.setDigital(0, pattern3);

%% Sub-sequence 2. Slow pulses on channel 1 only

% Discard existing sequence builder and create new
subseq2 = PulseStreamer.Sequence();

% Assign patterns to channels
subseq2.setDigital(1, pattern1);

%% Sub-sequence 3. Slow and fast pulses on channels 1 and 2
% Pattern1 and Pattern2

% Discard existing sequence builder and create new
subseq3 = PulseStreamer.Sequence();

% Assign patterns to channels
subseq3.setDigital(1, pattern1);
subseq3.setDigital(2, pattern2);

%% Create final sequence as a combination of the created sub-sequences

% At the beginning, we want to have to repeat twice a short pulse as trigger 
% so, we use subsequence1 as a starting point
sequence = subseq1;

% Now we add sub-sequence 2 repeated 2 times
sequence = [sequence, subseq2.repeat(2)];

% Here we add a single trigger pulse (sub-sequence 1), then 
% sub-sequence 3 repeated 2 times and followed by sub-sequence 2
sequence = [sequence, subseq1, subseq3.repeat(2), subseq2];

% Here we add more combinations
sequence = [sequence, subseq1, subseq3.repeat(3), subseq2];

%% Create clock pattern on channel 3 for the total sequence duration
% Calculate how many periods fint into sequence duration
seq_duration = sequence.getDuration();
pattern_clk_dur = sum(cell2mat(pattern_clk(:,1)));
n_periods = floor(seq_duration / pattern_clk_dur);

% Repeat pattern_clk to fill whole sequence and set it to channel 3.
% By doing this, we overwrite whatever was mapped to channel 3 with the new
% pattern.
sequence.setDigital(3, repmat(pattern_clk, n_periods, 1));

%% Send the sequence to the Pulse Streamer
finalState = PulseStreamer.OutputState([],0,0); % After sequence has finished this values will be set.
n_runs = 1; % Number of times to replay the sequence that was uploaded to the Pulse Streamer
ps.stream(sequence, n_runs, finalState)


%% Visualize the sequence 
% The Sequence object has a "plot()" method that plots the sequence data in
% a MATLAB figure. This method does not consider repetitions requested  
% from the hardware by the use of "nRuns" parameter.  
plot(sequence); % or sequence.plot()

%% How to replay this sequence without uploading it again
% As you may have noticed, it took a while to upload the sequence 
% we have created. This sequence stays in the memory of the
% Pulse Streamer until another sequence is sent, constant output is set or
% the Pulse Streamer is reset.

% First, we check that the Pulse Streamer still has the sequence
% this will print some status information to MATLAB's console
boolToYesNo = containers.Map({false, true},{'No', 'Yes'});
fprintf('hasSequence:\t %s\n', boolToYesNo(ps.hasSequence));
fprintf('isStreaming:\t %s\n', boolToYesNo(ps.isStreaming));
fprintf('hasFinished:\t %s\n', boolToYesNo(ps.hasFinished));

% You can tell the Pulse Streamer to replay this sequence again.
% Note that nRuns specified in ps.stream() function is still valid
ps.startNow();

% This will work as long as the trigger rearm is set to AUTO as
ps.setTrigger(TriggerStart.IMMEDIATE, TriggerRearm.AUTO); % This is the default setting

% In order to make sure that the sequence is streamed only once, either 
% with software or hardware trigger, we have to set trigger rearm to
% "MANUAL".
ps.setTrigger(TriggerStart.IMMEDIATE, TriggerRearm.MANUAL);

% From now on, the sequence will be streamed only once, immediately after 
% it is uploaded. However, it is possible to rearm the trigger without 
% uploading the sequence again. For that we use "rearm()" method.
if ps.rearm()
    disp('Rearmed.');
else
    disp('Not rearmed. Because not yet triggered or trigger mode is NORMAL');
end


