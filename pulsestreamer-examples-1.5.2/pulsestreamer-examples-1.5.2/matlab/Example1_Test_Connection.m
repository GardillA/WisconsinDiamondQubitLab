%% Example 1: Test connection and get hardware info.
% This script connects to PulseStreamer and prints: 
%  - MAC address of the ethernet interface 
%  - Serial number of the Pulse Streamer
%  - Firmware version
%  - Version of the API class
%

clear all;

% All PulseStreamer API functions and classes are placed in the
% "PulseStreamer" namespace. In order to shorter the names we import the
% namespace contents to the global namespace using "import" command. Then
% instead of "PulseStreamer.Pulsestreamer('ip')" we can write shorter
% version "PulseStreamer('ip')" 
import PulseStreamer.*


devices = findPulseStreamers();

disp(devices);

% DHCP is activated in factory settings
if numel(devices) > 0
    % Connect to the first discovered Pulse Streamer
    ip = devices(1).ip;
else
    % if discovery failed try to connect by the default hostname
    % IP address of the pulse streamer (default hostname is 'pulsestreamer')
    ip = 'pulsestreamer';
end

% connect to the pulse streamer
ps = PulseStreamer(ip);

% Print summary to MATLAB console
fprintf('*** Pulse Streamer ***\n');
fprintf('Serial:         %s\n', ps.getSerial());
fprintf('FPGA ID:        %s\n', ps.getFPGAID());
fprintf('Firmware ver:   %s\n', ps.getFirmwareVersion());
fprintf('Client version: %s\n', PulseStreamer.client_version());

