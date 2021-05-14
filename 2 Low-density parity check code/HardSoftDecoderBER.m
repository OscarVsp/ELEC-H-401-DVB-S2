clc;clear;close all;
addpath('../1 Optimal communication chain over the ideal channel');

%% Parameters

N_packet = 384;
N_bit_per_pack = 128;
CodeRate = 1/2;
N_bits=N_bit_per_pack*N_packet;
Nbps = 1;    %Nombre of bits per symbol (1 = BPSK, 2 = 4QAM, 4 = 16QAM, 6 = 64QAM)
f_cut = 1e6; %Hz Cutoff frequency
fsymb = 2*f_cut; 
T_symb = 1/fsymb;
fsamp = 16*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
M = 24; %Upsampling factor (link to fsamp/fsymb)
EbNo = 8; %Energy of one by over the PSD of the noise ratio (in dB)
N_taps = 101; %number of taps of the filter
beta = 0.3; %Makes the window smoother as beta increases // roll-off factor given in the specifications


EbNoArray = -4:1:16; %Energy of one by over the PSD of the noise ratio (in dB)
Average = 100;
BER_uncoded = zeros(1,length(EbNoArray));
BER_SoftDecoded = zeros(1,length(EbNoArray));
BER_HardDecoded = zeros(1,length(EbNoArray));
H0 = makeLdpc(N_bit_per_pack, N_bit_per_pack/CodeRate,0,1,3); % Create initial parity check matrix of size 128 x 256

for k=1:length(EbNoArray)
    100*(k-1)/length(EbNoArray)
    EbNo=EbNoArray(k);
    
    BER_uncoded_temp = zeros(1,length(Average));
    BER_SoftDecoded_temp = zeros(1,length(Average));
    BER_HardDecoded_temp = zeros(1,length(Average));
    
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

        symb_coded_tx = mapping(codedbits_tx', Nbps, 'pam')';
        symb_uncoded_tx = (mapping(bits_tx',Nbps,'pam'))';

        %% Upsampling

        upsampled_symb_coded_tx = UpSampling(symb_coded_tx,N_bits/CodeRate,Nbps,M);
        upsampled_symb_uncoded_tx = UpSampling(symb_uncoded_tx,N_bits,Nbps,M);

        %% Transmitter Filter

        filter = HalfrootNyquistFilter(fsamp,T_symb,N_taps); 
        signal_coded_tx = conv(upsampled_symb_coded_tx,filter);
        signal_uncoded_tx = conv(upsampled_symb_uncoded_tx,filter);

        %% Transmission Channel
        [signal_coded_rx,No] = NoiseAddition(signal_coded_tx,EbNo,fsamp,N_bits/CodeRate);    
        signal_uncoded_rx = NoiseAddition(signal_uncoded_tx,EbNo,fsamp,N_bits); 

        %% Receiver Filter
        matched_filter = flip(filter); 
        upsampled_symb_coded_rx = conv(signal_coded_rx,matched_filter,'valid'); 
        upsampled_symb_uncoded_rx = conv(signal_uncoded_rx,matched_filter,'valid'); 


        %% Downsampling

        symb_coded_rx = DownSampling(upsampled_symb_coded_rx,N_bits/CodeRate,Nbps,M);
        symb_uncoded_rx = DownSampling(upsampled_symb_uncoded_rx,N_bits,Nbps,M);

        %% Demapping
        %Only for the uncoded and hard decoder symb
        bit_uncoded_rx = demapping(real(symb_uncoded_rx)', Nbps, 'pam')';
        bit_uncoded_down_scaled = bit_uncoded_rx(1:N_bits);
        bit_coded_rx_hard = demapping(real(symb_coded_rx)', Nbps, 'pam')';
        bits_coded_hard_down_scaled = bit_coded_rx_hard(1:N_bits/CodeRate);

        %% Hard Decoder
        bit_decoded_rx_hard = LDPC_hard_decoder(bits_coded_hard_down_scaled,H,10);
        
        %% SoftDecoder
        %Only for the coded symb
            
        bit_decoded_rx_soft = LDPC_soft_decoder(symb_coded_rx,H,No/2,10);
        bit_decoded_rx_soft_down_scaled = bit_decoded_rx_soft(1:N_bits);
        
        BER_uncoded_temp(1,n) = ErrorCalculator(bit_uncoded_down_scaled,bits_tx);
        BER_SoftDecoded_temp(1,n) = ErrorCalculator(bit_decoded_rx_soft_down_scaled,bits_tx);
        BER_HardDecoded_temp(1,n) = ErrorCalculator(bit_decoded_rx_hard,bits_tx);
    end
    
    BER_uncoded(1,k) = mean(BER_uncoded_temp);
    BER_SoftDecoded(1,k) = mean(BER_SoftDecoded_temp);
    BER_HardDecoded(1,k) = mean(BER_HardDecoded_temp);
end

figure;
semilogy(EbNoArray,BER_uncoded,'-*');hold on;
semilogy(EbNoArray,BER_SoftDecoded,'-*');hold on;
semilogy(EbNoArray,BER_HardDecoded,'-*');hold on;
grid on; title('BER BPSK');
legend('Uncoded','Soft Decoded','Hard Decoded');
xlabel("E_b/N_0 [dB]");
ylabel("BER");


