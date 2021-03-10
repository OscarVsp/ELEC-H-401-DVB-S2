clear all; clc; close all;
%% Parameters

Ns = 100; %Number of sample to send for each test

Nbit = 12000;
Nbps = 6;               %Nombre of bits per symbol
M = 16;       %Upsampling factor
M_array = [2 4 8 16 32 48 64];
f_cut = 1e6; %Hz Cutoff frequency
T_symb = 1/(2*f_cut);
fs = M*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
EbNo = 10; %Energy of one by over the PSD of the noise ratio (in dB)

test_value_array = (1:1:20); %the values to test for the parameters
BER_array = zeros(4,length(test_value_array));

fig1 = figure(1);
for j = 1:7
    
    M = M_array(j);

    for i = 1:length(test_value_array)
        compteur = i*100/length(test_value_array)/7 +(j-1)*100/7 %Just a percent compteur for us
        EbNo = test_value_array(i); %HERE REMPLACE THE PARAMETER
        BER_samples = zeros(1,Ns);
        for n = 1:Ns
            %% Bit Generator

            bit_tx = randi(2,1,Nbit)-1; 
            Nbit = length(bit_tx);

            %% Mapping

            if (Nbps > 1)
                symb_tx = mapping(bit_tx', Nbps, 'qam')';
            else
                symb_tx = mapping(bit_tx', Nbps, 'pam')';
            end

            %figure(12); grid on; title("Mapped signal in absolute value"); stem(abs(symb_tx));
            %% Upsampling

            upsampled_symb_tx = UpSampling(symb_tx,Nbit,Nbps,M);


            %% Transmitter Filter

            filter = HalfrootNyquistFilter(fs,T_symb); %So we compute de filter only 1 time as this is the same
            N_filter = length(filter);

            signal_tx = conv(upsampled_symb_tx,filter);%Convolution of the signal with the filter



            %figure(7); stem(abs(signal_tx)); title("transmited signal")
            %signal_tx = upsampled_symb_tx;


            %% Transmission Channel

            signal_rx = NoiseAddition(signal_tx,EbNo,fs,Nbit);      %With Noise
            %signal_rx = signal_tx;                              %Without Nosie

            %fig_signal_tx = figure('Name','signal_tx','NumberTitle','off');plot(signal_rx,'b.');grid on;hold on;plot(signal_tx,'ro');



            %% Receiver Filter
            matched_filter = flip(filter); %Get the time reversal of the filter g(-t) 
            upsampled_symb_rx = conv(signal_rx,matched_filter,'valid'); %Matched filter convolved with the signal


            %figure(8);stem(abs(upsampled_symb_rx)); title("upsampled received signal"); grid on;
            %upsampled_symb_rx = upsampled_symb_rx( (N_filter):(length(upsampled_symb_rx)-(N_filter-1)) ); %Removing unecessary parts due to convolution (conv length = N + M -1 and we need to stay at M) --> to be sure to start at t=0
            %figure(9);stem(abs(upsampled_symb_rx)); title("Received signal truncated to start at first tap"); grid on;
            %upsampled_symb_rx = signal_rx;

            %% Downsampling

            symb_rx = DownSampling(upsampled_symb_rx,Nbit,Nbps,M);

            %figure(10); grid on; stem(abs(symb_rx)); title("Received downsampled signals")
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
    semilogy(test_value_array,BER_array(j,:));
    hold on;
end

title('BER fonction of Eb/No.');
xlabel('Eb/No (dB)');
ylabel('BER');
legend('M = 2','M = 4','M = 8','M = 16','M = 32','M = 48','M = 64');

