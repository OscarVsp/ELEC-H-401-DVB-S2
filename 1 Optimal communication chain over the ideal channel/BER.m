clear all; clc; close all;
%% Parameters

N_average = 10; %Number of sample to send for each test

Nbit = 12000;
f_cut = 1e6; %Hz Cutoff frequency
T_symb = 1/(2*f_cut);
M = 16;
N_taps = 201; %number of taps of the filter

M_array = [1 2 3 4 5 6 8 12 16]; %M
N_taps_array = [3 9 15 29 51 101 201]; %N_taps
Nbps_array = [6]; %NBPS
EbNo_array = (1:1:30); %EbNo

for k = 1:length(Nbps_array)
    
    Nbps = Nbps_array(k);
    figure(k);
    BER_array = zeros(4,length(EbNo_array));
    for j = 1:length(N_taps_array)

        N_taps = N_taps_array(j);

        for i = 1:length(EbNo_array)

            compteur = ((i-1) + (j-1)*length(EbNo_array) + (k-1)*(length(EbNo_array)*length(N_taps_array)))*(100)/(length(N_taps_array)*length(Nbps_array)*length(EbNo_array)) %Percent compteur for us
            EbNo = EbNo_array(i); %EbNo
            fs = M*f_cut;
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

                filter = HalfrootNyquistFilter(fs,T_symb,N_taps); %So we compute de filter only 1 time as this is the same
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
        semilogy(EbNo_array,BER_array(j,:));grid on;
        hold on;
    end

    title("BER fonction of Eb/No M =" + int2str(M)+", Nbps =" + int2str(Nbps)+", Nbit ="+int2str(Nbit)+", N average ="+int2str(N_average));
    %title("BER fonction of Eb/No M =" + int2str(M)+", Nbit ="+int2str(Nbit)+", N average ="+int2str(N_average));
    xlabel('E_b/N_o (dB)');
    ylabel('BER');
    %legend('M = 1','M = 2','M = 3','M = 4','M = 5','M = 6','M = 8','M = 12','M = 16');
    %legend('BPSK','QPSK','16QAM','64QAM');
    legend('N = 3','N = 9','N = 15','N = 29','N = 51','N = 101','N = 201');
end