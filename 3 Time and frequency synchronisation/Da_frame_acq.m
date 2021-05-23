function [symb_out,n_hat, df_hat,Dk] = Da_frame_acq(y,pilot,data_symbol_length,K,T_symb,Nbps,Nbit)

%This function gets the ToA, the CFO and the time shift of the received signal.
%The algorithm uses the differential cross correlation and removes the pilots.

%INPUTS
% y : Input sequence
% pilot : Known pilot mapped into symbols to correlate with
% data_length : length of the data packet preceded by the pilots
% K : cross correlation window length

%OUTPUTS
% symb_out : symbols with ToA, CFO and phase shift corrected

pilot_length = length(pilot);
codeblock_length =  pilot_length + data_symbol_length;
seq_length = length(y);
window_ratio = 1; %ratio of the window with respect to the codeblock_length (typically 1.5)

y = [y zeros(1,2*codeblock_length)]; %padding zeros to avoid indexing out
n_hat = zeros(1, floor(seq_length/codeblock_length)); %estimated indices of the pilots
df_hat = zeros(1,floor(seq_length/codeblock_length)); %estimated phase shift of the sequences
symb_out = zeros(1,Nbit/Nbps ); %number of symbols at the output
for i= 1:codeblock_length:seq_length-codeblock_length %-window_ratio*codeblock_length %for each packet with a pilot
	indices = i:i+window_ratio*codeblock_length-1;	
	window_y = [y(indices) zeros(1,pilot_length)];
	Dk = zeros(K,length(indices));
	for n = 1:length(indices)
		for k=1:K
			for l=k+1:pilot_length
				 Dk(k,n) = Dk(k,n) + (conj(y(n+l))*pilot(l) )*conj( conj(window_y(n+l-k))*pilot(l-k) );
			end
			Dk(k,n) = Dk(k,n)/(pilot_length-k);
		end
	end
	[~ , n_arg_max] = max(sum(Dk)); %sum of all the columns gives n_hat
	n_hat(1+(i-1)/codeblock_length) = indices(n_arg_max); % Get the good index in the window 
        df_hat(1+(i-1)/codeblock_length) = -sum(unwrap(angle(Dk(:,n_arg_max)))./([1:K]') )/K; %Get the phase value
        
        n_max_y =  indices(n_arg_max);
        data = y(n_max_y+pilot_length+1:n_max_y + pilot_length +data_symbol_length);  %phase interpolation
        n=i+pilot_length:codeblock_length + i-1; % taking the good indices where the CFO was
        
        if i ==1
            phase_interp = 0;
            phase_diff =0;
        else
            phase_diff = df_hat(1+(i-1)/codeblock_length)-df_hat( (i-1+codeblock_length)/codeblock_length -1);
            phase_interp= df_hat( (i-1+codeblock_length)/codeblock_length -1) :phase_diff/length(data): df_hat(1+(i-1)/codeblock_length)-phase_diff/length(data);
        end
        
        data = data.*exp(-1j*df_hat(1)*n-1j*phase_interp); %get rid of the CFO by multiplying back the complex exponential
        
        symb_out(1+(i-1)*data_symbol_length/codeblock_length:(i+codeblock_length-1)*data_symbol_length/codeblock_length) = data;
end

	
end		
	