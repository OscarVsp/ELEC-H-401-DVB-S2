clc;clear;close all;
addpath('../1 Optimal communication chain over the ideal channel');

%% Parameters

N_packet = 8;
N_bit_per_pack = 128;
CodeRate = 1/2;
N_bits=N_bit_per_pack*N_packet;
Nbps = 4;               %Nombre of bits per symbol (1 = BPSK, 2 = 4QAM, 4 = 16QAM, 6 = 64QAM)
M = 4;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
fs = 2*M*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
T_symb = 1/(2*(fs/10));
EbNoArray = 1:3:20; %Energy of one by over the PSD of the noise ratio (in dB)
N_taps = 201; %number of taps of the filter
Average = 10;
max_iter_array = [1 2 4 8];
BER_uncoded = zeros(1,length(EbNoArray));
BER_HardDecoded_1 = zeros(1,length(EbNoArray));
BER_HardDecoded_2 = zeros(1,length(EbNoArray));
BER_HardDecoded_4 = zeros(1,length(EbNoArray));
BER_HardDecoded_8 = zeros(1,length(EbNoArray));
H0 = makeLdpc(N_bit_per_pack, N_bit_per_pack/CodeRate,0,1,3); % Create initial parity check matrix of size 128 x 256

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

        bits_tx = randi(2,1,N_bits)-1;



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

        upsampled_symb_coded_tx = UpSampling(symb_coded_tx,N_bits*2,Nbps,M);
        upsampled_symb_uncoded_tx = UpSampling(symb_uncoded_tx,N_bits,Nbps,M);

        %% Transmitter Filter

        filter = HalfrootNyquistFilter(fs,T_symb,N_taps); 
        signal_coded_tx = conv(upsampled_symb_coded_tx,filter);
        signal_uncoded_tx = conv(upsampled_symb_uncoded_tx,filter);
%         
%         signal_coded_tx = upsampled_symb_coded_tx;
%         signal_uncoded_tx = upsampled_symb_uncoded_tx;

        
        %% Transmission Channel
        [signal_coded_rx,No] = NoiseAddition(signal_coded_tx,EbNo,fs,N_bits);    
        signal_uncoded_rx = NoiseAddition(signal_uncoded_tx,EbNo,fs,N_bits); 

        %% Receiver Filter
        matched_filter = flip(filter); 
        upsampled_symb_coded_rx = conv(signal_coded_rx,matched_filter,'valid'); 
        upsampled_symb_uncoded_rx = conv(signal_uncoded_rx,matched_filter,'valid'); 
%         
%         upsampled_symb_coded_rx = signal_coded_rx;
%         upsampled_symb_uncoded_rx = signal_uncoded_rx;

       
        %% Downsampling

        symb_coded_rx = DownSampling(upsampled_symb_coded_rx,N_bits*2,Nbps,M);
        symb_uncoded_rx = DownSampling(upsampled_symb_uncoded_rx,N_bits,Nbps,M);

        %% Demapping
        %Only for the uncoded and hard decoder symb
        
        if (Nbps > 1)
            bit_uncoded_rx = demapping(symb_uncoded_rx', Nbps, 'qam')';
            bit_coded_rx_hard = demapping(symb_coded_rx', Nbps, 'qam')';
        else
            bit_uncoded_rx = demapping(real(symb_uncoded_rx)', Nbps, 'pam')';
            bit_coded_rx_hard = demapping(real(symb_coded_rx)', Nbps, 'pam')';
        end
        
        BER_uncoded_temp(1,n) = ErrorCalculator(bit_uncoded_rx,bits_tx);
        
        
        %% Hard Decoder
        bit_decoded_rx_hard_1 = LDPC_hard_decoder(bit_coded_rx_hard,H,1);
        BER_HardDecoded_1_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard_1,bits_tx);
        
        bit_decoded_rx_hard_2 = LDPC_hard_decoder(bit_coded_rx_hard,H,2);
        BER_HardDecoded_2_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard_2,bits_tx);
        
        bit_decoded_rx_hard_4 = LDPC_hard_decoder(bit_coded_rx_hard,H,4);
        BER_HardDecoded_4_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard_4,bits_tx);
        
        bit_decoded_rx_hard_8 = LDPC_hard_decoder(bit_coded_rx_hard,H,8);
        BER_HardDecoded_8_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard_8,bits_tx);
        

    end
    
    BER_uncoded(1,k) = mean(BER_uncoded_temp);
    BER_HardDecoded_1(1,k) = mean(BER_HardDecoded_1_temp);
    BER_HardDecoded_2(1,k) = mean(BER_HardDecoded_2_temp);
    BER_HardDecoded_4(1,k) = mean(BER_HardDecoded_4_temp);
    BER_HardDecoded_8(1,k) = mean(BER_HardDecoded_8_temp);
end

figure;
semilogy(EbNoArray,BER_uncoded);hold on;
semilogy(EbNoArray,BER_HardDecoded_1);hold on;
semilogy(EbNoArray,BER_HardDecoded_2);hold on;
semilogy(EbNoArray,BER_HardDecoded_4);hold on;
semilogy(EbNoArray,BER_HardDecoded_8);hold on;
grid on; title('BER 64QAM');
legend('Uncoded','max Iter = 1','max Iter = 2','max Iter = 4','max Iter = 8');
xlabel("E_b/N_0 [dB]");
ylabel("BER");


