function [signal_noise] = NoiseAddition(signal,fs,Nbit)

% INPUTS:
% - signal : vector of input signal without noise
%
% OUTPUTS:
% - noiseSignal : vector of ouput signal with AWGN

%% Parameters

EbNo = 50; %Energy of one by over the PSD of the noise ratio (in dB)

%% Calculs

signal_energy = (trapz(abs(signal).^2))*(1/fs);     %energy of the signal
bit_energy = (signal_energy / Nbit)/2;   %divided by 2 bcs it's a baseband signal
No = bit_energy / EbNo;     % noise PSD
noise_power = 2*No*fs;
noise = sqrt(noise_power/2)*(randn(1,length(signal))+1i*randn(1,length(signal)));

signal_noise = signal + noise;

end

