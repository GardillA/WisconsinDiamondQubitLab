%% Example 3: This script shows various aspects of the API for the Pulse Streamer 
clear all

% All PulseStreamer API functions and classes are placed in the
% "PulseStreamer" namespace. In order to shorter the names we import the
% namespace contents to the global namespace using "import" command. Then
% instead of "PulseStreamer.Pulsestreamer('ip')" we can write shorter
% version "PulseStreamer('ip')" 
import PulseStreamer.*

% DHCP is activated in factory settings
% IP address of the pulse streamer (default hostname is PulseStreamer)
ip = 'pulsestreamer';

% Connect to the pulse streamer
ps = PulseStreamer(ip);

% Make sure that Pulse Streamer is in the default state
ps.reset();

%% Choose the active example here
exampleNo = 1;

% 1: all zero
% 2: channels [0,2,4,6] max
% 3: jumping between zero and max
% 4: callback function after sequence streaming has finished
% 5: sequence with alternating pulses (random length)
% 6: sequence duration of two seconds with callback when finished

switch exampleNo
    case 1
        % Set all outputs to zero
        
        % First, create the OutputState object which stores the output
        % values.
        % channels list to be set HIGH (first parameter) and the
        % analog channels second and third parameter
        % OutputState(channels, analog0, analog1)
        outputState = OutputState([],0,0);
        
        % The following method applies the state to the outputs
        % immediately. Use this method for ocassional state changes only.
        ps.constant(outputState);
        
    case 2
        % Set outputs to max value
        % digital channels
        channels = [0,2,4,6]; % Channels to be set HIGH (+3V)
        
        % Analog channel values shall be specified in Volts.
        % Analog channels have 16 bit resolution
        analog0 = 1;  % +1V
        analog1 = -1; % -1V
        
        outputMax = OutputState(channels, analog0, analog1);
        ps.constant(outputMax);
        
    case 3
        % This case shows how to stream a sequence of values 
        % (switch between the two values from case 1 and 2)
        % In the previous 2 cases you have learned how to set a constant
        % value at the Pulse Streamer output. When you want to output
        % arbitrary sequence of values with hardware controlled timing, 
        % the "stream()" method shall be used instead of "constant()".
       
        %duration of the pulses
        duration = 1000; %ns
        
        pattern = {duration, true; duration, false};
        
        sequence = ps.createSequence();
        
        % Set the pattern to all digital outputs
        for channel = 0:7
            sequence.setDigital(channel, pattern);
        end
        
        
        % nRuns parameter specifies the number of times the sequence 
		% should be repeated by hardware.
        % When n_runs<0 the sequence will be repeated indefinitely.
        n_runs = PulseStreamer.REPEAT_INFINITELY;
        
        % We also have to define the output state that shall be set
        % after the sequence streaming has completed.
        final_state = OutputState([],0,0);
        
        % Finally, the way to start the sequence must be given. Here we
        % start the sequence as soon as it sent to the Pulse Streamer
        
        % Sent the data to the Pulse Streamer.
        ps.stream(sequence, n_runs, final_state);
        
    case 4
        % Use a callback method as a way to find out whether the
        % Pulse Streamer output has finished
        
        % We take the sequence from case 3 but with a duration of 1e9 ns for
        % each pulse so that the duration of the sequence is 2s.
        % In total, we run it 3 times that means 3x2 seconds.
        
        n_runs = 3;
        final_state = OutputState([],0,0);
        
        duration = 1e9; %ns
        
        pattern = {duration,true; duration, false};
        
        sequence = ps.createSequence();
        % Set the pattern to first 4 digital outputs
        sequence.setDigital(0:3, pattern);
        
        % Here we add the callback function and start streaming the sequence.
        ps.setCallbackFinished(@Example3_VariousCallbackMethod);
        ps.stream(sequence, n_runs, final_state);
        
        % The callback function will be called after 6s and will print
        % a message to console on completion
        
    case 5
        %
        %   Generate a sequence of alternating high low pulses with random pulse lengths on the digital
        %   channels 1-7 and the two analog channels.
        %
        %   Digital channel 0 is used to produce trigger signal.
        %
        %   The generated sequence runs in an infinite loop.
        n_pulses = 1000;
        minPulseLength = 100;
        maxPulseLength = 1000;
        n_runs = PulseStreamer.REPEAT_INFINITELY;
        
        times = randi([minPulseLength, maxPulseLength], 1, n_pulses); % random pulse durations
        
        sequence = ps.createSequence();
        
        % Trigger signal high during random sequence 
        % and low for 10000 ns after it.
        sequence.setDigital(0, {sum(times), true; 10000, false});
        
        for channel = 1:7
            % Digital patterns
            levels = randi([0,1],1,n_pulses);
            pattern = num2cell([times(:), levels(:)]);
            sequence.setDigital(channel, pattern);
        end
        
        for channel = 0:1
           % Analog patterns
           levels = rand(1,n_pulses)*2-1;
           pattern = num2cell([times(:), levels(:)]);
           sequence.setAnalog(channel, pattern);
        end
        
        final_state = OutputState([],0,0);
        ps.stream(sequence, n_runs, final_state);
        
    case 6
        %   
        %   Generate a sequence of 100 alternating high low pulses (low: 10us,
        %   high: 10us) which are repeated 1000 times (total duration: 2s).
        %   A callback function is registered so that at the end the 
        %   Example2_VariousCallbackMethod
        %   is executed showing 
        %   "hasFinishedCallback - Pulse Streamer finished." in the console
        n_pulses = 100;
        n_runs = 1000;
        pattern = cell(n_pulses, 2);
        for i = 1:n_pulses
            state = mod(i,2);
            pattern(i,:) = {10000, state};
        end
        sequence = ps.createSequence();
        sequence.setDigital(0, pattern);
        
        final_state = OutputState([],0,0);
        ps.setCallbackFinished(@Example3_VariousCallbackMethod)
        ps.stream(sequence, n_runs, final_state);
end
