function [image_rx] = BitsToImage(bits_rx,nbl,nbc,nbd)
% INPUTS :
% - bits_tx : bits
% OUTPUTS :
% - NA

image_rx_bin = reshape(bits_rx,nbl*nbc*nbd,8);
image_rx_int = bi2de(image_rx_bin);
image_rx = uint8(reshape(image_rx_int,nbl,nbc,nbd));
end