function [bits] = LDPC_hard_decoder(r,H)
% Hard decoder for the LDPC algorithm
% r : noisy received bits 
% H : parity check matrix
% 
% bits : decoded bits

[m, n] = size(H);
L = length(r);
bits = zeros(1,L*m/n);

for b = 1:L/n
    dr = r( 1+(b-1)*n:(b-1)*n+n);%taking a n size block
    new_dr = dr;
    vote_result = zeros(1,n);
    stop_count=0;
    while (norm(mod(new_dr*H',2)) ~=0) && (stop_count<10) %checking if syndrome equals 0 vector
    stop_count =stop_count +1;
    syndrome = mod(new_r*H',2);
    for i=1:n %for each variable node
        check_node = H(:,i); %taking the link between check and variable nodes
        vote = new_r(i);
        count =1;%Because the received vector is taken for the majority vote
        for j=1:m %for each check node for this variable node
            if check_node(j) == 1 %decision of the check node needed
                count =count+ 1;
                if syndrome(j) == 1 %changing corresponding bits if check node equals 1 from syndrome
                    vote = vote+ mod(r(i) + 1,2);  %inversing the value because assumed to be an error
                else
                    vote = vote + r(i);
                end
            end
        end
        vote_result(i) = round(vote/count); %Vote
    end
    new_dr = vote_result;
    end
    bits(1+(b-1)*m:(b-1)*m+m) = new_dr(m+1:end); %the last bits are the transmitted ones
end