clc;clear;close all;
addpath('../1 Optimal communication chain over the ideal channel');

%% Parameters

N_packet = 96*4;
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

%% BER parameters
Average = 50;
EbNoArray = -4:1:20; %Energy of one by over the PSD of the noise ratio (in dB)
max_iter_array = [1 2 4 8];

H0 = makeLdpc(N_bit_per_pack, N_bit_per_pack/CodeRate,0,1,3); % Create initial parity check matrix

for m=[1 2 4 6]
    
    Nbps = m;

    BER_uncoded = zeros(1,length(EbNoArray));
    BER_HardDecoded_1 = zeros(1,length(EbNoArray));
    BER_HardDecoded_2 = zeros(1,length(EbNoArray));
    BER_HardDecoded_4 = zeros(1,length(EbNoArray));
    BER_HardDecoded_8 = zeros(1,length(EbNoArray));
    
    for k=1:length(EbNoArray)
        k
        EbNo=EbNoArray(k);

        BER_uncoded_temp = zeros(1,length(Average));
        BER_HardDecoded_1_temp = zeros(1,length(EbNoArray));
        BER_HardDecoded_2_temp = zeros(1,length(EbNoArray));
        BER_HardDecoded_4_temp = zeros(1,length(EbNoArray));
        BER_HardDecoded_8_temp = zeros(1,length(EbNoArray));

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

            %% Encoder

            blocks=(reshape(bits_tx,N_bit_per_pack,N_packet))';
            [checkbits, H] = makeParityChk(blocks', H0, 0);
            checkbits = (checkbits)';
            codedbits = (horzcat(checkbits,blocks));
            codedbits_tx=(reshape(codedbits',[],1))';

            %% Mapping

            if (Nbps > 1)
                symb_coded_tx = mapping(codedbits_tx', Nbps, 'qam')';
                symb_uncoded_tx = mapping(bits_tx',Nbps,'qam')';
            else
                symb_coded_tx = mapping(codedbits_tx', Nbps, 'pam')';
                symb_uncoded_tx = mapping(bits_tx',Nbps,'pam')';
            end

            %% Upsampling

            upsampled_symb_coded_tx = UpSampling(symb_coded_tx,Nbit_tx*2,Nbps,M);
            upsampled_symb_uncoded_tx = UpSampling(symb_uncoded_tx,Nbit_tx,Nbps,M);

            %% Transmitter Filter

            filter = HalfrootNyquistFilter(fsamp,T_symb,N_taps); 
            signal_coded_tx = conv(upsampled_symb_coded_tx,filter);
            signal_uncoded_tx = conv(upsampled_symb_uncoded_tx,filter);   

            %% Transmission Channel
            [signal_coded_rx,No] = NoiseAddition(signal_coded_tx,EbNo,fsamp,Nbit_tx/CodeRate);    
            signal_uncoded_rx = NoiseAddition(signal_uncoded_tx,EbNo,fsamp,Nbit_tx); 

            %% Receiver Filter
            matched_filter = flip(filter); 
            upsampled_symb_coded_rx = conv(signal_coded_rx,matched_filter,'valid'); 
            upsampled_symb_uncoded_rx = conv(signal_uncoded_rx,matched_filter,'valid'); 

            %% Downsampling

            symb_coded_rx = DownSampling(upsampled_symb_coded_rx,Nbit_tx*2,Nbps,M);
            symb_uncoded_rx = DownSampling(upsampled_symb_uncoded_rx,Nbit_tx,Nbps,M);

            %% Demapping
            %Only for the uncoded and hard decoder symb

            if (Nbps > 1)
                bit_uncoded_rx = demapping(symb_uncoded_rx', Nbps, 'qam')';
                bit_coded_rx_hard = demapping(symb_coded_rx', Nbps, 'qam')';
            else
                bit_uncoded_rx = demapping(real(symb_uncoded_rx)', Nbps, 'pam')';
                bit_coded_rx_hard = demapping(real(symb_coded_rx)', Nbps, 'pam')';
            end
            
            bits_uncoded_down_scaled = bit_uncoded_rx(1:N_bits);
            bits_coded_hard_down_scaled = bit_coded_rx_hard(1:N_bits/CodeRate);

            BER_uncoded_temp(1,n) = ErrorCalculator(bits_uncoded_down_scaled,bits);
            
            


            %% Hard Decoder
            bit_decoded_rx_hard_1 = LDPC_hard_decoder(bits_coded_hard_down_scaled,H,1);
            BER_HardDecoded_1_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard_1,bits);

            bit_decoded_rx_hard_2 = LDPC_hard_decoder(bits_coded_hard_down_scaled,H,2);
            BER_HardDecoded_2_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard_2,bits);

            bit_decoded_rx_hard_4 = LDPC_hard_decoder(bits_coded_hard_down_scaled,H,4);
            BER_HardDecoded_4_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard_4,bits);

            bit_decoded_rx_hard_8 = LDPC_hard_decoder(bits_coded_hard_down_scaled,H,8);
            BER_HardDecoded_8_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard_8,bits);


        end

        BER_uncoded(1,k) = mean(BER_uncoded_temp);
        BER_HardDecoded_1(1,k) = mean(BER_HardDecoded_1_temp);
        BER_HardDecoded_2(1,k) = mean(BER_HardDecoded_2_temp);
        BER_HardDecoded_4(1,k) = mean(BER_HardDecoded_4_temp);
        BER_HardDecoded_8(1,k) = mean(BER_HardDecoded_8_temp);
        
    end

    figure(m);
    semilogy(EbNoArray,BER_uncoded,'-*');hold on;
    semilogy(EbNoArray,BER_HardDecoded_1,'-*');hold on;
    semilogy(EbNoArray,BER_HardDecoded_2,'-*');hold on;
    semilogy(EbNoArray,BER_HardDecoded_4,'-*');hold on;
    semilogy(EbNoArray,BER_HardDecoded_8,'-*');hold on;
    grid on; title("BER Nbps = " + int2str(m));
    legend('Uncoded','max Iter = 1','max Iter = 2','max Iter = 4','max Iter = 8');
    xlabel("E_b/N_0 [dB]");
    ylabel("BER");

end

