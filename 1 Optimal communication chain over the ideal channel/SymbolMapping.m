clear all; clc; close all;


%% Parameters

Nbit = 2000;             %Nombre of bits
Nbps = 4;               %Nombre of bits per symbol


%% Bit Generator

bit_tx = randi(2,Nbit,1)-1;


%% Mapping

if (Nbps > 1)
    symb_tx = mapping(bit_tx, Nbps, 'qam');
else
    symb_tx = mapping(bit_tx, Nbps, 'pam');
end


%% Transmitter Filter

signal_tx = HalfrootNyquistFilter(symb_tx);


%% Transmission Channel

%signal_rx = NoiseAddition(signal_tx,fs,Nbit);      %With Noise
signal_rx = signal_tx;                              %Without Nosie

fig_signal_tx = figure('Name','signal_tx','NumberTitle','off');plot(signal_rx,'b.');grid on;hold on;plot(signal_tx,'ro');



%% Receiver Filter

symb_rx = HalfrootNyquistFilter(signal_rx);


%% Demapping

if (Nbps > 1)
    bit_rx = demapping(symb_rx, Nbps, 'qam').';
else
    bit_rx = demapping(real(symb_rx), Nbps, 'pam').';
end


%% Bits check

%check = norm(bit_tx - bit_rx)