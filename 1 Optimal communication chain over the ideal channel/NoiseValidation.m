clear all; clc; close all;


%% Parameters

Nbit = 600000;
Nbps = 4;               %Nombre of bits per symbol (1 = BPSK, 2 = 4QAM, 4 = 16QAM, 6 = 64QAM)
M = 4;       %Upsampling factor
f_cut = 1e6; %Hz Cutoff frequency
fs = 2*M*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut) Its the freq on which the conv of the filter and the signal will be done --> has to be the same !!!
T_symb = 1/(2*(fs/10));
N_taps = 201; %number of taps of the filter

EbNo_array = [2 4 8 12 16 24];
ErrorRatio_array = zeros(1,length(EbNo_array));
legend_array = ["test"];

Nbps_array = [1 2 4 6];

for j=1:length(Nbps_array)
    Nbps = Nbps_array(j)
    figure(j);grid on;
    for i=1:length(EbNo_array)
        EbNo = EbNo_array(i);

        %% Bit Generator

        bits = randi(2,1,Nbit)-1;
        check_mult = mod(Nbit,Nbps);
        if check_mult == 0
            bit_tx = bits;
        else
            bit_tx = [bits  zeros(1,Nbps - check_mult)];
        end
        Nbit_tx = length(bit_tx);

        %% Mapping
        %Maps the bits into desired symbols in the complex plane 
        if (Nbps > 1)
            symb_tx = mapping(bit_tx', Nbps, 'qam')';
        else
            symb_tx = mapping(bit_tx', Nbps, 'pam')';
        end


        %% Upsampling

        upsampled_symb_tx = UpSampling(symb_tx,Nbit_tx,Nbps,M);


        %% Transmitter Filter

        signal_tx = upsampled_symb_tx;


        %% Transmission Channel

        signal_rx = NoiseAddition(signal_tx,EbNo,fs,Nbit);


        %% Receiver Filter


        upsampled_symb_rx = signal_rx; 


        %% Downsampling

        symb_rx = DownSampling(upsampled_symb_rx,Nbit_tx,Nbps,M);


        %% Demapping

        if (Nbps > 1)
            bit_rx = demapping(symb_rx', Nbps, 'qam')';
        else
            bit_rx = demapping(real(symb_rx)', Nbps, 'pam')';
        end

        bit_down_scaled = bit_rx(1:Nbit);

        ErrorRatio = ErrorCalculator(bit_down_scaled,bits);
        ErrorRatio_array(i) = ErrorRatio;
        legend_array(i) = int2str(EbNo)+" ("+int2str(ErrorRatio*100)+")";
        plot(symb_rx,'.');hold on;
    end
    title("Symbols with noise, Nbps = "+int2str(Nbps)+", Nbit = "+int2str(Nbit)+", M = "+int2str(M));
    legend(legend_array);
end