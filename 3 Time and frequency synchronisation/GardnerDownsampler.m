function [symb_rx, estimated_t_shift] = GardnerDownsampler(upsampled_symb_rx,Nbit,Nbps,M,time_shift, kappa)
%This function Downsamples a given data sequence and adds a time shift and
%correct the time shift induced with the gardner algorithm. 
%
% INPUTS
% upsampled_symb_rx : upsampled signal to downsample
% Nbit : number of sended bits
% M : oversampling factor
% time_shift : time shift in indices; The perturbation
% kappa : slope correction factor
%
%  OUTPUTS
% symb_rx : symbols, downsampled
% estimated_t_shift : estimated time shift by the gardner algorithm


eps_hat = zeros(1,Nnit/Nbps);%Initialized at zero
for i = 1:Nbit-1 %for each symbol
    eps_hat(i+1) = eps_hat(i) + 2*kappa*real( upsampled_symb_rx(1+time_shift+round(M*(i-1)/2) )*( conj(upsampled_symb_rx(1+time_shift+M*(i-1))- conj(upsampled_symb_rx(1+time_shift+M*(i-2))) )) );          
end
end

