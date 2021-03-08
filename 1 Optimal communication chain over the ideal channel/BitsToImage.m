function [image_rx] = BitsToImage(bits_rx,nbl,nbc)
% INPUTS :
% - bits_tx : bits
% OUTPUTS :
% - NA

image_rx_bin = reshape(bits_rx,nbl*nbc,8);
image_rx_int = bi2de(image_rx_bin);
image_rx = reshape(image_rx_int,nbl,nbc);
figure(1) ; colormap gray ; image(image_rx) ; % plot result