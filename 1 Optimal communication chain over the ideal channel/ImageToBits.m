function [bits_tx,nbl,nbc] = ImageToBits(file)
% INPUTS :
% - file : name of the jpeg file
% OUTPUTS :
% - bits_tx : bits
image_tx = file ; % file.jpg
image_tx_int = uint8(imread(image_tx,'jpg')) ; % gray, 8 bits/pixel
image_tx_int = image_tx_int(1:end,1:end,1:1);
[nbl,nbc] = size(image_tx_int) ;
image_tx_bin = de2bi(image_tx_int,8);
bits_tx = reshape(image_tx_bin,1,nbl*nbc*8)>0;