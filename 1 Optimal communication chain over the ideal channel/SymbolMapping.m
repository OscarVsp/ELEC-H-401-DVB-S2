clear all; clc; close all;


%% Parameters

% image_tx = 'cp.png'; %Only 3x8bits images
% [bit_tx,nbl,nbc,nbd] = ImageToBits(image_tx);

Nbit = 1000;
Nbps = 4;               %Nombre of bits per symbol
M = 10;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
fs = M*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
EbNo =10; %Energy of one by over the PSD of the noise ratio (in dB)

%% Bit Generator

bit_tx = randi(2,1,Nbit)-1; 
Nbit = length(bit_tx);

%% Mapping

if (Nbps > 1)
    symb_tx = mapping(bit_tx', Nbps, 'qam')';
else
    symb_tx = mapping(bit_tx', Nbps, 'pam')';
end

%figure(12); grid on; title("Mapped signal in absolute value"); stem(abs(symb_tx));
%% Upsampling

upsampled_symb_tx = UpSampling(symb_tx,Nbit,Nbps,M);


%% Transmitter Filter

filter = HalfrootNyquistFilter(fs); %So we compute de filter only 1 time as this is the same
N_filter = length(filter);

signal_tx = conv(upsampled_symb_tx,filter);%Convolution of the signal with the filter



%figure(7); stem(abs(signal_tx)); title("transmited signal")
%signal_tx = upsampled_symb_tx;


%% Transmission Channel

%signal_rx = NoiseAddition(signal_tx,EbNo,fs,Nbit);      %With Noise
signal_rx = signal_tx;                              %Without Nosie

%fig_signal_tx = figure('Name','signal_tx','NumberTitle','off');plot(signal_rx,'b.');grid on;hold on;plot(signal_tx,'ro');



%% Receiver Filter
matched_filter = flip(filter); %Get the time reversal of the filter g(-t) 
upsampled_symb_rx = conv(signal_rx,matched_filter); %Matched filter convolved with the signal


%figure(8);stem(abs(upsampled_symb_rx)); title("upsampled received signal"); grid on;
upsampled_symb_rx = upsampled_symb_rx( (N_filter-1):(length(upsampled_symb_rx)-N_filter-1) ); %Removing unecessary parts due to convolution (conv length = N + M -1 and we need to stay at M) --> to be sure to start at t=0
%figure(9);stem(abs(upsampled_symb_rx)); title("Received signal truncated to start at first tap"); grid on;
%upsampled_symb_rx = signal_rx;

%% Downsampling

symb_rx = DownSampling(upsampled_symb_rx,Nbit,Nbps,M);

%figure(10); grid on; stem(abs(symb_rx)); title("Received downsampled signals")
%% Demapping
 
if (Nbps > 1)
    bit_rx = demapping(symb_rx', Nbps, 'qam')';
else
    bit_rx = demapping(real(symb_rx)', Nbps, 'pam')';
end


%% Bits check

ErrorRatio = ErrorCalculator(bit_rx,bit_tx)

%image_rx = BitsToImage(bit_rx,nbl,nbc,nbd);
%figure('name',"BER = "+num2str(ErrorRatio,4));image(image_rx);