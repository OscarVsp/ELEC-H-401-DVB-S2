function [symb_rx] = DownSampling(upsampled_symb_rx,Nbit,Nbps,M)

% INPUTS:
% - symb_tx : vector of input symbol to upsample
% - Nbit : nomber of bit
% - Nbps : nomber of symbol per bits
% - M : umsampling factor
%
% OUTPUTS:
% - upsampled_symb_tx : vector of ouput symbol upsampled


symb_rx = zeros(1,Nbit/Nbps);
for i = 1:Nbit/Nbps
    symb_rx(i)=upsampled_symb_rx(1+M*(i-1));
end


end



