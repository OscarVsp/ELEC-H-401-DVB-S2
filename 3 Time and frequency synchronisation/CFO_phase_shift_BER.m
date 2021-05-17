clc;clear;close all;
addpath('../1 Optimal communication chain over the ideal channel');

%% Parameters

N_packet = 96*1;
N_bit_per_pack = 128;
CodeRate = 1/2;
N_bits=N_bit_per_pack*N_packet;
Nbps = 4;    %Nombre of bits per symbol (1 = BPSK, 2 = 4QAM, 4 = 16QAM, 6 = 64QAM)
f_cut = 1e6; %Hz Cutoff frequency
fsymb = 2*f_cut; 
T_symb = 1/fsymb;
fsamp = 16*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
M = 24; %Upsampling factor (link to fsamp/fsymb)
EbNo = 8; %Energy of one by over the PSD of the noise ratio (in dB)
N_taps = 101; %number of taps of the filter
beta = 0.3; %Makes the window smoother as beta increases // roll-off factor given in the specifications
fs_passband= 2e9; %Hz
CFO = 2; % CFO in ppm (2 to see something)

%% BER parameters
Average = 10;
EbNoArray = -4:2:20; %Energy of one by over the PSD of the noise ratio (in dB)


for m=[1]
    
    Nbps = m;

    BER_NO_CFO = zeros(1,length(EbNoArray));
    BER_CFO_1 = zeros(1,length(EbNoArray));
    BER_CFO_2 = zeros(1,length(EbNoArray));
    BER_CFO_3 = zeros(1,length(EbNoArray));
    BER_CFO_4 = zeros(1,length(EbNoArray));
    
    for k=1:length(EbNoArray)
        k
        EbNo=EbNoArray(k);

        BER_NO_CFO_temp = zeros(1,length(Average));
        BER_CFO_1_temp = zeros(1,length(EbNoArray));
        BER_CFO_2_temp = zeros(1,length(EbNoArray));
        BER_CFO_3_temp = zeros(1,length(EbNoArray));
        BER_CFO_4_temp = zeros(1,length(EbNoArray));

        for n=1:Average

            %% Bit Generator

            bits = randi(2,1,N_bits)-1;
            check_mult = mod(N_bits,Nbps);
            if check_mult == 0
                bits_tx = bits;
            else
                bits_tx = [bits  zeros(1,Nbps - check_mult)];
            end
            Nbit_tx = length(bits_tx);

            %% Mapping

            if (Nbps > 1)
                symb_tx = mapping(bits_tx', Nbps, 'qam')';
            else
                symb_tx = mapping(bits_tx', Nbps, 'pam')';
            end

            %% Upsampling

            upsampled_symb_tx = UpSampling(symb_tx,Nbit_tx,Nbps,M);
            

            %% Transmitter Filter

            filter = HalfrootNyquistFilter(fsamp,T_symb,N_taps); 
            signal_tx = conv(upsampled_symb_tx,filter);   

            %% Transmission Channel   
            signal_rx = NoiseAddition(signal_tx,EbNo,fsamp,Nbit_tx); 
            
            %% CFO          
            CFO_1 = fs_passband*0*1e-6; %CFO in Hz
            CFO_2 = fs_passband*0.5*1e-6; %CFO in Hz
            CFO_3 = fs_passband*1*1e-6; %CFO in Hz
            CFO_4 = fs_passband*12*1e-6; %CFO in Hz
            [signal_rx_CFO_1, phase_shift] = Add_CFO_and_phase_shift(signal_rx,CFO_1,fsamp);
            [signal_rx_CFO_2, phase_shift] = Add_CFO_and_phase_shift(signal_rx,CFO_2,fsamp);
            [signal_rx_CFO_3, phase_shift] = Add_CFO_and_phase_shift(signal_rx,CFO_3,fsamp);
            [signal_rx_CFO_4, phase_shift] = Add_CFO_and_phase_shift(signal_rx,CFO_4,fsamp);

            %% Receiver Filter
            matched_filter = flip(filter); 
            upsampled_symb_rx = conv(signal_rx,matched_filter,'valid'); 
            upsampled_symb_rx_CFO_1 = conv(signal_rx_CFO_1,matched_filter,'valid');
            upsampled_symb_rx_CFO_2 = conv(signal_rx_CFO_2,matched_filter,'valid'); 
            upsampled_symb_rx_CFO_3 = conv(signal_rx_CFO_3,matched_filter,'valid'); 
            upsampled_symb_rx_CFO_4 = conv(signal_rx_CFO_4,matched_filter,'valid'); 


            %% Downsampling
            symb_rx = DownSampling(upsampled_symb_rx,Nbit_tx,Nbps,M);
            
            sampling_time_shift = 0;
            symb_rx_CFO_1 = DownSampling(upsampled_symb_rx_CFO_1,Nbit_tx,Nbps,M,sampling_time_shift);
            symb_rx_CFO_2 = DownSampling(upsampled_symb_rx_CFO_2,Nbit_tx,Nbps,M,sampling_time_shift);
            symb_rx_CFO_3 = DownSampling(upsampled_symb_rx_CFO_3,Nbit_tx,Nbps,M,sampling_time_shift);
            symb_rx_CFO_4 = DownSampling(upsampled_symb_rx_CFO_4,Nbit_tx,Nbps,M,sampling_time_shift);
            

            %% Demapping
            %Only for the uncoded and hard decoder symb

            if (Nbps > 1)
                bit_rx = demapping(symb_rx', Nbps, 'qam')';
                bit_rx_CFO_1 = demapping(symb_rx_CFO_1', Nbps, 'qam')';
                bit_rx_CFO_2 = demapping(symb_rx_CFO_2', Nbps, 'qam')';
                bit_rx_CFO_3 = demapping(symb_rx_CFO_3', Nbps, 'qam')';
                bit_rx_CFO_4 = demapping(symb_rx_CFO_4', Nbps, 'qam')';                
            else
                bit_rx = demapping(real(symb_rx)', Nbps, 'pam')';
                bit_rx_CFO_1 = demapping(real(symb_rx_CFO_1)', Nbps, 'pam')';
                bit_rx_CFO_2 = demapping(real(symb_rx_CFO_2)', Nbps, 'pam')';
                bit_rx_CFO_3 = demapping(real(symb_rx_CFO_3)', Nbps, 'pam')';
                bit_rx_CFO_4 = demapping(real(symb_rx_CFO_4)', Nbps, 'pam')';
            end
            
            
            bits_downscale = bit_rx(1:N_bits);
            bits_downscale_CFO_1 = bit_rx_CFO_1(1:N_bits);
            bits_downscale_CFO_2 = bit_rx_CFO_2(1:N_bits);
            bits_downscale_CFO_3 = bit_rx_CFO_3(1:N_bits);
            bits_downscale_CFO_4 = bit_rx_CFO_4(1:N_bits);
            

            BER_NO_CFO_temp(1,n) = ErrorCalculator(bits_downscale,bits);
            BER_CFO_1_temp(1,n) = ErrorCalculator(bits_downscale_CFO_1,bits);
            BER_CFO_2_temp(1,n) = ErrorCalculator(bits_downscale_CFO_2,bits);
            BER_CFO_3_temp(1,n) = ErrorCalculator(bits_downscale_CFO_3,bits);
            BER_CFO_4_temp(1,n) = ErrorCalculator(bits_downscale_CFO_4,bits);

        end

        BER_NO_CFO(1,k) = mean(BER_NO_CFO_temp);
        BER_CFO_1(1,k) = mean(BER_CFO_1_temp);
        BER_CFO_2(1,k) = mean(BER_CFO_2_temp);
        BER_CFO_3(1,k) = mean(BER_CFO_3_temp);
        BER_CFO_4(1,k) = mean(BER_CFO_4_temp);
        
    end

    figure(m);
    semilogy(EbNoArray,BER_NO_CFO,'-*');hold on;
    semilogy(EbNoArray,BER_CFO_1,'-*');hold on;
    semilogy(EbNoArray,BER_CFO_2,'-*');hold on;
    semilogy(EbNoArray,BER_CFO_3,'-*');hold on;
    semilogy(EbNoArray,BER_CFO_4,'-*');hold on;
    grid on; title("BER Nbps = " + int2str(m));
    legend('CFO = 0','CFO = 0.1','CFO = 0.5','CFO = 1','CFO = 2');
    xlabel("E_b/N_0 [dB]");
    ylabel("BER");

end

