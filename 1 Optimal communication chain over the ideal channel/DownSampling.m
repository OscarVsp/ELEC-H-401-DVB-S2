function [symb_rx] = DownSampling(upsampled_symb_rx,Nbit,Nbps,M,varargin)

% INPUTS:
% - symb_tx : vector of input symbol to upsample
% - Nbit : nomber of bit
% - Nbps : nomber of symbol per bits
% - M : umsampling factor
% - varargin : only one value : the time shift
% OUTPUTS:
% - upsampled_symb_tx : vector of ouput symbol upsampled
time_shift=0;
if length(varargin) ==1
    time_shift=varargin{1};
end

symb_rx = zeros(1,Nbit/Nbps);
for i = 1:Nbit/Nbps
    symb_rx(i)=upsampled_symb_rx(1+time_shift+M*(i-1));
end


end
