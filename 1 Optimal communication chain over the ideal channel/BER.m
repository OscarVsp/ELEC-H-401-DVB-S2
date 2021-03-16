clear all; clc; close all;
%% Parameters

N_average = 200; %Number of sample to send for each test

Nbit = 12000;
Nbps = 4;               %Nombre of bits per symbol
M = 4;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
T_symb = 1/(2*f_cut);
fs = M*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
EbNo = 10; %Energy of one by over the PSD of the noise ratio (in dB)

first_value_array = [4 8 16 32 48 64 96 128]; %M
second_value_array = [1 2 4 6]; %NBPS
third_value_array = (1:0.5:30); %EbNo

for k = 1:length(first_value_array)
    
    M = first_value_array(k);
    figure(k);
    BER_array = zeros(4,length(third_value_array));
    for j = 1:length(second_value_array)

        Nbps = second_value_array(j); %Nbps

        for i = 1:length(third_value_array)

            compteur = ((i-1) + (j-1)*length(third_value_array) + (k-1)*(length(third_value_array)*length(second_value_array)))*(100)/(length(first_value_array)*length(second_value_array)*length(third_value_array)) %Percent compteur for us
            EbNo = third_value_array(i); %EbNo
            BER_samples = zeros(1,N_average);
            for n = 1:N_average
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


                %% Transmission Channel

                signal_rx = NoiseAddition(signal_tx,EbNo,fs,Nbit);      %With Noise


                %% Receiver Filter

                matched_filter = flip(filter); %Get the time reversal of the filter g(-t) 
                upsampled_symb_rx = conv(signal_rx,matched_filter,'valid'); %Matched filter convolved with the signal; 'valid' to remove padded zeroes added by the convolution


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
            BER_array(j,i) = mean(BER_samples); 
        end
        semilogy(third_value_array,BER_array(j,:));
        hold on;
    end

    title("BER fonction of Eb/No M =" + int2str(M)+", Nbit ="+int2str(Nbit)+", N average ="+int2str(N_average));
    xlabel('Eb/No (dB)');
    ylabel('BER');
    %legend('M = 2','M = 4','M = 8','M = 16','M = 32','M = 48','M = 64');
    legend('BPSK','QPSK','16QAM','64QAM');

end
