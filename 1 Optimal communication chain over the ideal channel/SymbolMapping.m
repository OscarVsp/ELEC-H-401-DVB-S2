%%%%Parameters%%%%

Nbit = 200;         %Nombre of bit
Nbps = 4;           %Nombre of bits per symbol

%%%%Bit Generator%%%%

bit_tx = randi(2,Nbit,1)-1;

%%%%Mapping%%%%

symb_tx = mapping(bit_tx, Nbps, 'qam');

%%%Transmission%%%

symb_rx = symb_tx;

%%%Demapping%%%%

bit_rx = demapping(symb_rx, Nbps, 'qam');

check = norm(bit_tx - bit_rx)