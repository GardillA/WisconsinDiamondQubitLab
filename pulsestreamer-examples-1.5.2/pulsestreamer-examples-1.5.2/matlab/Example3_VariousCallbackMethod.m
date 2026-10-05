%% Callback function for example 3.

function Example3_VariousCallbackMethod(pstreamer)
    % this is the test callback function for Example3_Various.m - case 4
    
    disp('hasFinishedCallback - Pulse Streamer finished.');
    
    % With "pstreamer" we have access the PulseStreamer object that called
    % this callback. Here we simply print sequence streaming status.
    YesNo = containers.Map({false, true}, {'No', 'Yes'});
    fprintf('hasSequence:\t %s\n', YesNo(pstreamer.hasSequence()));
    fprintf('isStreaming:\t %s\n', YesNo(pstreamer.isStreaming()));
    fprintf('hasFinished:\t %s\n', YesNo(pstreamer.hasFinished()));
end
