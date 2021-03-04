function [upsampled_symb_tx] = UpSampling(symb_tx,Nbit,Nbps,M)

% INPUTS:
% - symb_tx : vector of input symbol to upsample
% - Nbit : nomber of bit
% - Nbps : nomber of symbol per bits
% - M : umsampling factor
%
% OUTPUTS:
% - upsampled_symb_tx : vector of ouput symbol upsampled


upsampled_symb_tx = zeros(1,Nbit/Nbps*M);
for i = 1:M:Nbit/Nbps*M
    upsampled_symb_tx(i)=symb_tx((i-1)/M +1);
end

end

