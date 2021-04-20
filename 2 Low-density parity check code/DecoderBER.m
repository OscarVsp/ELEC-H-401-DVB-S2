clc;clear;close all;
addpath('../1 Optimal communication chain over the ideal channel');

%% Parameters

% image_tx = 'cp.png'; %Only 3x8bits images
% [bit_tx,nbl,nbc,nbd] = ImageToBits(image_tx);

N_packet = 12;
N_bit_per_pack = 128;
CodeRate = 1/2;
N_bits=N_bit_per_pack*N_packet;
Nbps = 1;               %Nombre of bits per symbol (1 = BPSK, 2 = 4QAM, 4 = 16QAM, 6 = 64QAM)
M = 4;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
fs = 2*M*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
T_symb = 1/(2*(fs/10));
EbNoArray = 1:20; %Energy of one by over the PSD of the noise ratio (in dB)
N_taps = 201; %number of taps of the filter
Average = 10;
BER_uncoded = zeros(1,length(EbNoArray));
BER_SoftDecoded = zeros(1,length(EbNoArray));
BER_HardDecoded = zeros(1,length(EbNoArray));


for k=1:length(EbNoArray)
    k
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


        %% Encoder

        H0 = makeLdpc(N_bit_per_pack, N_bit_per_pack/CodeRate,0,1,3); % Create initial parity check matrix of size 128 x 256
        blocks=reshape(bits_tx,N_packet,N_bit_per_pack);
        [checkbits, H] = makeParityChk(blocks', H0, 0);
        checkbits = (checkbits)';
        codedbits = (horzcat(checkbits,blocks));
        codedbits_tx=(reshape(codedbits,[],1))';


        %% Mapping

        symb_coded_tx = mapping(codedbits_tx', Nbps, 'pam')';
        symb_uncoded_tx = (mapping(bits_tx',Nbps,'pam'))';

        %% Upsampling

        upsampled_symb_coded_tx = UpSampling(symb_coded_tx,N_bits*2,Nbps,M);
        upsampled_symb_uncoded_tx = UpSampling(symb_uncoded_tx,N_bits,Nbps,M);

        %% Transmitter Filter

        filter = HalfrootNyquistFilter(fs,T_symb,N_taps); 
        signal_coded_tx = conv(upsampled_symb_coded_tx,filter);
        signal_uncoded_tx = conv(upsampled_symb_uncoded_tx,filter);

        %% Transmission Channel
        signal_coded_rx = NoiseAddition(signal_coded_tx,EbNo,fs,N_bits);    
        [signal_uncoded_rx,No] = NoiseAddition(signal_uncoded_tx,EbNo,fs,N_bits); 

        %% Receiver Filter
        matched_filter = flip(filter); 
        upsampled_symb_coded_rx = conv(signal_coded_rx,matched_filter,'valid'); 
        upsampled_symb_uncoded_rx = conv(signal_uncoded_rx,matched_filter,'valid'); 


        %% Downsampling

        symb_coded_rx = DownSampling(upsampled_symb_coded_rx,N_bits*2,Nbps,M);
        symb_uncoded_rx = DownSampling(upsampled_symb_uncoded_rx,N_bits,Nbps,M);

        %% Demapping
        %Only for the uncoded symb
        bit_uncoded_rx = demapping(real(symb_uncoded_rx)', Nbps, 'pam')';


        %% SoftDecoder
        %Only for the coded symb
        codedblocks = reshape(symb_coded_rx,N_packet,N_bit_per_pack/CodeRate);
        bit_decoded_rx = zeros(N_packet,N_bit_per_pack);
        for i = 1:N_packet
            block = codedblocks(i,:);
            bit_coded_rx = LDPC_soft_decoder(H,block,No/2,10);
            bit_decoded_rx(i,:) = bit_coded_rx(1,N_bit_per_pack+1:end);
        end

        bit_decoded_rx = reshape(bit_decoded_rx,1,[]);

        BER_uncoded_temp(1,n) = ErrorCalculator(bit_uncoded_rx,bits_tx);
        BER_SoftDecoded_temp(1,n) = ErrorCalculator(bit_decoded_rx,bits_tx);
        BER_HardDecoded_temp(1,n) = 0;
    end
    
    BER_uncoded(1,k) = mean(BER_uncoded_temp);
    BER_SoftDecoded(1,k) = mean(BER_SoftDecoded_temp);
    BER_HardDecoded(1,k) = mean(BER_HardDecoded_temp);
end

figure;
semilogy(EbNoArray,BER_uncoded);hold on;
semilogy(EbNoArray,BER_SoftDecoded);hold on;
%semilogy(EbNoArray,BER_HardDecoded);hold on;
grid on; title('BER BPSK M=4');
legend('Uncoded','Soft Decoded');
%legend('Uncoded','Soft Decoded','Hard Decoded');
xlabel("E_b/N_0 [dB]");
ylabel("BER");

