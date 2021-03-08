clear all; clc; close all;
%% Parameters

Ns = 500; %Number of sample to send for each test

Nbit = 48000;
Nbps = 4;               %Nombre of bits per symbol
M = 4;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
fs = 10*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut)
EbNo = 1; %Energy of one by over the PSD of the noise ratio (in dB)

test_value_array = (0:0.1:30); %the values to test for the parameters
BER_array = zeros(1,length(test_value_array));

for i = 1:length(test_value_array)
    compteur = i*100/length(test_value_array) %Just a percent compteur for us
    EbNo = test_value_array(i); %HERE REMPLACE THE PARAMETER
    BER_samples = zeros(1,Ns);
    for n = 1:Ns
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

        filter = HalfrootNyquistFilter(fs); %So we compute de filter only 1 time as this is the same

        %signal_tx = conv(upsampled_symb_tx,filter); %doesn't work for now, bcs
        %filter is not a vector, it's a 51x51 complexe matrix

        signal_tx = upsampled_symb_tx;


        %% Transmission Channel

        signal_rx = NoiseAddition(signal_tx,EbNo,fs,Nbit);      %With Noise
        %signal_rx = signal_tx;                              %Without Nosie

        %fig_signal_tx = figure('Name','signal_tx','NumberTitle','off');plot(signal_rx,'b.');grid on;hold on;plot(signal_tx,'ro');



        %% Receiver Filter

        %upsampled_symb_rx = conv(signal_rx,filter); %same as transmitter

        upsampled_symb_rx = signal_rx;

        %% Downsampling

        symb_rx = DownSampling(upsampled_symb_rx,Nbit,Nbps,M);


        %% Demapping

        if (Nbps > 1)
            bit_rx = demapping(symb_rx', Nbps, 'qam')';
        else
            bit_rx = demapping(real(symb_rx)', Nbps, 'pam')';
        end


        %% Bits check

        BER_samples(n) = ErrorCalculator(bit_rx,bit_tx);
    end
    BER_array(i) = mean(BER_samples);
end


fig1 = figure(1);
semilogy(test_value_array,BER_array,'b');
hold on;
limite = ones(1,length(test_value_array))./Nbit;
semilogy(test_value_array,limite,'r');
hold on;
title('BER fonction of Eb/No.');
xlabel('Eb/No (dB)');
ylabel('BER (dB)');
