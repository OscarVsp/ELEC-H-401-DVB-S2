function [ErrorRatio] = ErrorCalculator(bit_rx,bit_tx)
% INPUTS :
% - bit_rx : bits received
% - bit_tx : bits sended
% OUTPUTS :
% - error : error ratio

%Errors calculation 
errors = 0;
sub = bit_tx - bit_rx;

for i = 1:length(bit_tx)
    if sub(i) ~= 0
        errors = errors + 1;
    end
end
ErrorRatio = errors / length(bit_tx);
end