function [u] = LDPC_encoder(bits,parity_matrix) 
% Encodes a bit stream with an LDPC algorithm
%
% bits : bit stream to be encoded needing to be multiple of the parity
% matrix length
% 
% pariry_matrix : parity check matrix from a given parity check matrix
% 
% u : encoded bit stream
L = length(bits);
m= length(parity_matrix); %number of rows of the parity matrix
if fix(L/m) ~= L/m %check if padding zeroes needed
    new_L = ceil(L/m)*m;
    bits = [bits zeros(new_L-L,1)]; %padding zeroes
    L = new_L;
end

n = 2*m;
G = [parity_matrix eye(m)];

u = zeros(1,L*n/m);

for i=1:L/m
    u(( 1+(i-1)*n):((i-1)*n+n) )= mod(bits((1+(i-1)*m):((i-1)*m+m))*G,2);
end
