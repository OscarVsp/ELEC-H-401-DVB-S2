clear all; clc; close all;


%% Parameters

Nbit = 200;             %Nombre of bits
Nbps = 4;               %Nombre of bits per symbol
mod = 'qam';            %Type of modulation ('qam' or 'pam')


%% Bit Generator

bit_tx = randi(2,Nbit,1)-1;

%fig_bit_tx = figure('Name','bits','NumberTitle','off');
%plot((0:size(bit_tx)-1),bit_tx);
%hold;


%% Mapping

symb_tx = mapping(bit_tx, Nbps, mod);

fig_symb_tx = figure('Name','Symbols','NumberTitle','off');
plot(symb_tx,'o'); %mettre les valeurs de références
hold;
plot(symb_tx,'.');
hold;

%% Transmitter Filter

signal_tx = HalfrootNyquistFilter(symb_tx);


%% Transmission Channel

signal_rx = signal_tx;


%% Receiver Filter

symb_rx = HalfrootNyquistFilter(signal_rx);


%% Demapping

bit_rx = demapping(symb_rx, Nbps, mod);


%% Bits check

check = norm(bit_tx - bit_rx)