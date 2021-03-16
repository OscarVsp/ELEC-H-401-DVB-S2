clear all; clc; close all;


%% Parameters

Nbit = 600;
Nbps = 4;               %Nombre of bits per symbol (1 = BPSK, 2 = 4QAM, 4 = 16QAM, 6 = 64QAM)
M = 4;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
T_symb = 1/(2*f_cut);
fs = M*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
EbNo = 20; %Energy of one by over the PSD of the noise ratio (in dB)
EbNo_array = [32 24 16 12 8 4 2];
%% Bit Generator

bit_tx = randi(2,1,Nbit)-1; 

%% Mapping
%Maps the bits into desired symbols in the complex plane 
if (Nbps > 1)
    symb_tx = mapping(bit_tx', Nbps, 'qam')';
else
    symb_tx = mapping(bit_tx', Nbps, 'pam')';
end

figure(1);
title('symbol with noise')
plot(symb_tx,'o'); hold on;
legend_array = zeros(1,length(EbNo_array)+1);
legend_array(1) = "tx";
%% Upsampling

upsampled_symb_tx = UpSampling(symb_tx,Nbit,Nbps,M);


%% Transmitter Filter

signal_tx = upsampled_symb_tx;


%% Transmission Channel

signal_rx = zeros(length(EbNo_array),length(signal_tx));
for i = length(EbNo_array)
    signal_rx(i) = NoiseAddition(signal_tx,EbNo_array(i),fs,Nbit);
    legend_array(i+1) = "EbNo = "+int2str(EbNo(i))+" ,";
end


%% Receiver Filter

upsampled_symb_rx = zeros(length(EbNo_array),length(upsampled_signal_tx));
for i = length(EbNoArray)
    upsampled_symb_rx(i) = signal_rx(i);
end



%% Downsampling

symb_rx = zeros(length(EbNo_array),length(symb_tx));
for i = length(EbNoArray)
    symb_rx(i) = DownSampling(upsampled_symb_rx(i),Nbit,Nbps,M);
    plot(symb_rx(i),'.'); hold on;
end



%% Demapping

 
bit_rx = zeros(length(EbNo_array),length(bit_tx));
ErrorRatio = zeros(1,length(EbNo_array));
for i = length(EbNoArray)
    if (Nbps > 1)
        bit_rx(i) = demapping(symb_rx(i)', Nbps, 'qam')';
    else
        bit_rx(i) = demapping(real(symb_rx(i))', Nbps, 'pam')';
    end
    legend_array(i+1) = legend_array(i+1) + "BER = "+int2str(ErrorCalculator(bit_rx(i),bit_tx));
end
legend(legend_array);