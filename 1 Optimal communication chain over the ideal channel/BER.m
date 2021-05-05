clear all; clc; close all;
%% Parameters

N_average = 50; %Number of sample to send for each test
Nbit = 40000;
Nbps = 4;    %Nombre of bits per symbol (1 = BPSK, 2 = 4QAM, 4 = 16QAM, 6 = 64QAM)
f_cut = 1e6; %Hz Cutoff frequency
fsymb = 2*f_cut;
T_symb = 1/fsymb;
fsamp = 16*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
M = 32; %Upsampling factor (au moins plus grand que 4)
EbNo = 8; %Energy of one by over the PSD of the noise ratio (in dB)
N_taps = 101; %number of taps of the filter
beta = 0.3; %Makes the window smoother as beta increases // roll-off factor given in the specifications

%% BER parameters arrays

M_array = [4 8 12 16 20 24 28 32 40 48]; %M
N_taps_array = [11 31 61 101 151 201]; %M
Nbps_array = [1 2 4 6]; %NBPS
EbNo_array = (1:1:25); %EbNo

%% Loops

for k = 1:length(Nbps_array)
    
    Nbps = Nbps_array(k);
    figure(k);
    BER_array = zeros(4,length(EbNo_array));
    BER_array_simple = zeros(4,length(EbNo_array));
        
    for j = 1:length(N_taps_array)

        N_taps = N_taps_array(j);6

        for i = 1:length(EbNo_array)

            compteur = ((i-1) + (j-1)*length(EbNo_array) + (k-1)*(length(EbNo_array)*length(N_taps_array)))*(100)/(length(Nbps_array)*length(N_taps_array)*length(EbNo_array)) %Percent compteur for us
            EbNo = EbNo_array(i); %EbNo
            BER_samples = zeros(1,N_average);
            BER_samples_simple = zeros(1,N_average);
            for n = 1:N_average
                %% Bit Generator
 
                bits = randi(2,1,Nbit)-1;
                check_mult = mod(Nbit,Nbps);
                if check_mult == 0
                    bits_tx = bits;
                else
                    bits_tx = [bits  zeros(1,Nbps - check_mult)];
                end
                Nbit_tx = length(bits_tx);
                

                %% Mapping
                %Maps the bits into desired symbols in the complex plane 
                if (Nbps > 1)
                    symb_tx = mapping(bits_tx', Nbps, 'qam')';
                else
                    symb_tx = mapping(bits_tx', Nbps, 'pam')';
                end

                %% Upsampling

                upsampled_symb_tx = UpSampling(symb_tx,Nbit_tx,Nbps,M);


                %% Transmitter Filter

                filter = HalfrootNyquistFilter(fsamp,T_symb,N_taps); %So we compute de filter only 1 time as this is the same
                signal_tx = conv(upsampled_symb_tx,filter);%Convolution of the signal with the filter

                %% Transmission Channel

                signal_rx = NoiseAddition(signal_tx,EbNo,fsamp,Nbit_tx);      %With Noise


                %% Receiver Filter

                matched_filter = flip(filter); %Get the time reversal of the filter g(-t) 
                upsampled_symb_rx = conv(signal_rx,matched_filter,'valid'); %Matched filter convolved with the signal; 'valid' to remove padded zeroes added by the convolution


                %% Downsampling

                symb_rx = DownSampling(upsampled_symb_rx,Nbit_tx,Nbps,M);
        
                %% Demapping

                if (Nbps > 1)
                    bits_rx = demapping(symb_rx', Nbps, 'qam')';
                else
                    bits_rx = demapping(real(symb_rx)', Nbps, 'pam')';
                end
                bits_down_scaled = bits_rx(1:Nbit);
                %% Bits check

                BER_samples(n) = ErrorCalculator(bits_rx,bits_tx);
            end
            BER_array(j,i) = mean(BER_samples);
        end
        semilogy(EbNo_array,BER_array(j,:),'-*');grid on;
        hold on;
    end

%     title("BER fonction of E_b/N_o M =" + int2str(M)+", Nbit ="+int2str(Nbit)+", N average ="+int2str(N_average));
%     title("BER fonction of E_b/N_o Nbps =" + int2str(Nbps)+", Nbit ="+int2str(Nbit)+", N average ="+int2str(N_average));
    title("BER fonction of E_b/N_o M =" + int2str(M)+", Nbps =" + int2str(Nbps)+", Nbit ="+int2str(Nbit)+", N average ="+int2str(N_average));
    xlabel('E_b/N_o (dB)');
    ylabel('BER');
%     legend('BPSK','QPSK','16-QAM','64-QAM');
%     legend('M = 4','M = 8','M = 12','M = 16','M = 20','M = 24','M = 28','M = 32','M = 40','M = 48');
    legend('#taps = 11','#taps = 21','#taps = 61','#taps = 101','#taps = 151','#taps = 201');
end