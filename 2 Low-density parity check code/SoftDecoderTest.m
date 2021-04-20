clc;clear;close all;
addpath('../1 Optimal communication chain over the ideal channel');

Nbps = 1
N_packet = 2;
N_error = 2;
N_bit_per_pack = 16;
CodeRate = 1/2;
N_bits=N_bit_per_pack*N_packet;
H0 = makeLdpc(N_bit_per_pack, N_bit_per_pack/CodeRate,0,1,3); % Create initial parity check matrix of size 128 x 256

bits_tx = randi(2,1,N_bits)-1
blocks=reshape(bits_tx,N_bit_per_pack,N_packet)
[checkbits, H] = makeParityChk(blocks, H0, 0)

codedbits = horzcat(checkbits.',blocks.')


for i=length(codedbits)
    block = codedbits(1,:);
    error = zeros(1,N_bit_per_pack*2);
    for j=1:N_error
        error(randi(N_bit_per_pack*2,1,1)) = 1;
    end
    block_received = mod(block+error,2);
    block_received_symb = mapping(block_received.',Nbps,'pam').';
    decoded_bits = SoftDecoder(H,block_received_symb,(1e-10)/2,30);
    error_norm = norm(error)
    received_error = norm(block_received-block)
    decoded_error = norm(decoded_bits-block)
end





