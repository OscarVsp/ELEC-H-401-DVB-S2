function [signal_noise] = NoiseAddition(signal,EbNo,fs,Nbit)

% INPUTS:
% - signal : vector of input signal without noise
%
% OUTPUTS:
% - noiseSignal : vector of ouput signal with AWGN
signal_energy = (trapz(abs(signal).^2))*(1/fs);     %energy of the signal
eb = (signal_energy) /(2* Nbit);   %divided by 2 bcs it's a baseband signal
No = eb /(10.^(EbNo/10));     % noise PSD
noise_power = 2*No*fs; %psd 2 bp
noise = sqrt(noise_power/2)*(randn(1,length(signal))+1i*randn(1,length(signal)));

signal_noise = signal + noise;

end

