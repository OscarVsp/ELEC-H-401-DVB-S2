clear all; clc; close all;


%% Parameters

% image_tx = 'cp.png'; %Only 3x8bits images
% [bit_tx,nbl,nbc,nbd] = ImageToBits(image_tx);

Nbit = 20;
Nbps = 4;               %Nombre of bits per symbol (1 = BPSK, 2 = 4QAM, 4 = 16QAM, 6 = 64QAM)
M = 16;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
T_symb = 1/(2*f_cut);
fs = M*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
EbNo = 15; %Energy of one by over the PSD of the noise ratio (in dB)

%% Bit Generator

bit_tx = randi(2,1,Nbit)-1; 
Nbit = length(bit_tx);

%% Mapping
%Maps the bits into desired symbols in the complex plane 
if (Nbps > 1)
    symb_tx = mapping(bit_tx', Nbps, 'qam')';
else
    symb_tx = mapping(bit_tx', Nbps, 'pam')';
end

%% Upsampling


upsampled_symb_tx = UpSampling(symb_tx,Nbit,Nbps,M);

%% Transmitter Filter

filter = HalfrootNyquistFilter(fs,T_symb); %So we compute de filter only 1 time as this is the same

signal_tx = conv(upsampled_symb_tx,filter);%Convolution of the signal with the filter
signal_tx_2 = upsampled_symb_tx;





%% Transmission Channel

signal_rx = NoiseAddition(signal_tx,EbNo,fs,Nbit);      %With Noise
signal_rx_2 = NoiseAddition(signal_tx_2,EbNo,fs,Nbit);                          %Without Nosie





%% Receiver Filter
matched_filter = flip(filter); %Get the time reversal of the filter g(-t) 
upsampled_symb_rx = conv(signal_rx,matched_filter,'valid'); %Matched filter convolved with the signal; 'valid' to remove padded zeroes added by the convolution
upsampled_symb_rx_2 = signal_rx_2; %Matched filter convolved with the signal; 'valid' to remove padded zeroes added by the convolution


%% Downsampling

symb_rx = DownSampling(upsampled_symb_rx,Nbit,Nbps,M);
symb_rx_2 = DownSampling(upsampled_symb_rx_2,Nbit,Nbps,M);

%figure(10); grid on; stem(abs(symb_rx)); title("Received downsampled signals")
%% Demapping
 
if (Nbps > 1)
    bit_rx = demapping(symb_rx', Nbps, 'qam')';
else
    bit_rx = demapping(real(symb_rx)', Nbps, 'pam')';
end

if (Nbps > 1)
    bit_rx_2 = demapping(symb_rx_2', Nbps, 'qam')';
else
    bit_rx_2 = demapping(real(symb_rx_2)', Nbps, 'pam')';
end


%% Bits check

%ErrorRatio = ErrorCalculator(bit_rx,bit_tx)

figure(1); stem(real(bit_tx)); hold on; stem(real(bit_rx)); hold on; stem(real(bit_rx_2)); hold on;legend('tx','rx filter','rx no filter'); title("Bits");
figure(2); stem(real(symb_tx)); hold on; stem(real(symb_rx)); hold on; stem(real(symb_rx_2)); hold on; legend('tx','rx filter','rx no filter');title("Symbol");
figure(3); stem(real(upsampled_symb_tx)); hold on; stem(real(upsampled_symb_rx));hold on; stem(real(upsampled_symb_rx_2));hold on;legend('tx','rx filter','rx no filter'); title("Symbol upsampled");
figure(4); stem(real(signal_tx)); hold on; stem(real(signal_rx));hold on; legend('signal filter','signal filter with noise');title("Signal tx after fitler");
figure(5); stem(real(signal_tx_2)); hold on; stem(real(signal_rx_2));hold on; legend('signal no filter','signal no filter with noise');title("Signal rx after fitler");
%image_rx = BitsToImage(bit_rx,nbl,nbc,nbd);
%figure('name',"BER = "+num2str(ErrorRatio,4));image(image_rx);