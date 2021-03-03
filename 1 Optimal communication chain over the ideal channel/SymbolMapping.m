%%%----Parameters----%%%

Nbit = 200;             %Nombre of bits
Nbps = 4;               %Nombre of bits per symbol
mod = 'pam';            %Type of modulation ('qam' or 'pam')


%%%----Bit Generator----%%%

bit_tx = randi(2,Nbit,1)-1;


%%%----Mapping----%%%

symb_tx = mapping(bit_tx, Nbps, mod);

%%%----Transmitter Filter----%%%

signal_tx = HalfrootNyquistFilter(symb_tx);


%%%----Transmission Channel----%%%

signal_rx = signal_tx;


%%%----Receiver Filter----%%%

symb_rx = HalfrootNyquistFilter(signal_rx);


%%%----Demapping----%%%

bit_rx = demapping(symb_rx, Nbps, mod);


%%%----Bits check----%%%

check = norm(bit_tx - bit_rx)