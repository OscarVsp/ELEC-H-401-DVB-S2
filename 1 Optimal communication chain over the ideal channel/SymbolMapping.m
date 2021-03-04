clear all; clc; close all;


%% Parameters

Nbit = 2000;             %Nombre of bits
Nbps = 4;               %Nombre of bits per symbol
M = 4;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
fs = 10*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut)



%% Bit Generator

bit_tx = randi(2,1,Nbit)-1; 


%% Mapping

if (Nbps > 1)
    symb_tx = mapping(bit_tx', Nbps, 'qam')';
else
    symb_tx = mapping(bit_tx', Nbps, 'pam')';
end

%% Upsampling

upsampled_symb_tx = UpSampling(symb_tx,Nbit,Nbps,M);


%% Transmitter Filter

signal_tx = HalfrootNyquistFilter(upsampled_symb_tx,fs);


%% Transmission Channel

%signal_rx = NoiseAddition(signal_tx,fs,Nbit);      %With Noise
signal_rx = signal_tx;                              %Without Nosie

%fig_signal_tx = figure('Name','signal_tx','NumberTitle','off');plot(signal_rx,'b.');grid on;hold on;plot(signal_tx,'ro');



%% Receiver Filter

upsampled_symb_rx = HalfrootNyquistFilter(signal_rx,fs);

%% Downsampling

symb_rx = DownSampling(upsampled_symb_rx,Nbit,Nbps,M);


%% Demapping
 
if (Nbps > 1)
    bit_rx = demapping(symb_rx', Nbps, 'qam');
else
    bit_rx = demapping(real(symb_rx)', Nbps, 'pam');
end


%% Bits check

check_bits = norm(bit_tx - bit_rx)