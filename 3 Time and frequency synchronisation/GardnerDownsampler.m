function [symb_rx, estimated_t_shift] = GardnerDownsampler(upsampled_symb_rx,Nbit,Nbps,M,time_shift, kappa)
%This function Downsamples a given data sequence and adds a time shift and
%correct the time shift induced with the gardner algorithm. 
%
% INPUTS
% upsampled_symb_rx : upsampled signal to downsample
% Nbit : number of sended bits
% M : oversampling factor
% time_shift : time shift in indices
% kappa : slope correction factor
%
%  OUTPUTS
% symb_rx : symbols, downsampled
% estimated_t_shift : estimated time shift by the gardner algorithm

for i = 1:Nbit/Nbps %for each symbol
end
end

