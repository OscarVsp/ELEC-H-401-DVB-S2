function [bits_tx,nbl,nbc,nbd] = ImageToBits(file)
% INPUTS :
% - file : name of the jpeg file
% OUTPUTS :
% - bits_tx : bits

pres = 1;

image_name = file ; % file.jpg
image_tx = uint8(imread(image_name)); %8 bits/pixel/layer
image_tx_int = image_tx(1:pres:end,1:pres:end,1:end);
[nbl,nbc,nbd] = size(image_tx_int);
image_tx_bin = de2bi(image_tx_int,8);
bits_tx = reshape(double(image_tx_bin),1,nbl*nbc*nbd*8);
end