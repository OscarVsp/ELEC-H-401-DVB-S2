function [corrupted_signal,phase_shift,time_shift] = Add_sync_errors(r,CFO,fsamp,M,Nbps)
%This function adds CFO, time and phase shift errors. It has to be implemented
%just before the matched filter in the communication channel.

%   r is the signal affected by the errors
%   CFO is the carrier frequency offset in Hz of the carrier frequency
%   fc is the carrier frequency
%   time_shift is the time shift of the signal which will be added with
%   circshift function
%   phase_shift is the phase shift between 0 and 2*pi

eps = (rand(1)-0.5); % values of eps in [-0.5,0.5] time shift is eps*T (symbol duration)
time_shift= round(eps*(M-1)); % We can have a time shift up to 1 symbol in indices
%time_shift=-4;
phase_shift = rand(1)*2*pi;
%phase_shift = 0;
N = length(r); n = 0:N-1;
T = 1/fsamp;
corrupted_signal = r;
%corrupted_signal = r.*exp(1j*(2*pi*CFO*T.*n +phase_shift)); %multiplication by the complex exponential (CFO and phase shift)
corrupted_signal = circshift(corrupted_signal,time_shift); %Time shift represented by a circular shift
end

