function [corrupted_signal] = Add_CFO_and_phase_shift(r,CFO,fsamp)
%This function adds CFO and phase shift errors. It has to be implemented
%just before the matched filter in the communication channel.
%   r is the signal affected by the errors
%   CFO is the carrier frequency offset in Hz of the carrier frequency
%   fc is the carrier frequency
%   phase_shift is the phase shift between 0 and 2*pi

%phase_shift = rand(1)*2*pi;
phase_shift = 0;
N = length(r); n = 0:N-1;
T = 1/fsamp;
corrupted_signal = r.*exp(1j*(2*pi*CFO*T.*n +phase_shift));

end

