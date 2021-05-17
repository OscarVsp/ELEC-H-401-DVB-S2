function [symb_rx, eps_hat] = Gardner(upsampled_symb_rx,M,kappa)
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


%Gardner algo
Nbit =length(upsampled_symb_rx);
symb_rx = zeros(1,Nbit/M);
eps_hat = zeros(1,Nbit/M); symb_rx(1)= upsampled_symb_rx(1);%Initialized at zero
for i = 1:Nbit/M -1 %for each symbol
    time_window = 1+(i-1)*M:M*i+M/2; %Taking the time window of 1 symbol
    symbol = upsampled_symb_rx(time_window);
    n = M*i+1 -eps_hat(i);
    n_half = M*i+1-M/2 - eps_hat(i);
    symb_rx(i+1) = interp1(time_window,symbol,n,'linear');
    symb_half = interp1(time_window,symbol,n_half,'linear');
    
    eps_hat(i+1) = eps_hat(i) + 2*kappa*real( symb_half*( conj(symb_rx(i+1)- conj(symb_rx(i)) )) );
    
    %upsampled_symb_rx = circshift(upsampled_symb_rx,round(eps_hat(i+1))*M/Nbps);
end

end

