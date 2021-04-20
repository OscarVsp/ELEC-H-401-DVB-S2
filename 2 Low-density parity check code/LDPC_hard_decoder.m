function [bits] = LDPC_hard_decoder(r,H,max_iter)
% Hard decoder for the LDPC algorithm
% r : noisy received bits 
% H : parity check matrix
% max_iter : maximal number of iterations
%
% bits : decoded bits

[m, n] = size(H);
L = length(r);
bits = zeros(1,L*m/n);

for block = 1:L/n
    dr = r( 1+(block-1)*n:(block-1)*n+n);%taking a n size block
    v_nodes = dr;    %Initial value of v nodes

    L_q = zeros(m,n);           %Initialize the value sent from v nodes to c nodes
    u = ones(1,n);              %Initialize the output value

    %% Step 0
    %We take the probability at the v nodes and we send it to the c nodes.
    for l=1:n       %For each v nodes
        nodes_index = find(H(:,l));
        L_q(nodes_index,l)=v_nodes(l);
    end
    n_iter = 0;
    while (n_iter < max_iter && norm(mod(u*H',2))~=0)
        syndrome = mod(u*H',2)
        n_iter + 1 %compteur


        %% Step 1
        %We compute the probability to send back to each v nodes from the
        %previous probability
        L_r = zeros(m,n);   %Reset L_r
        for l=1:m          %For each c nodes
            nodes_index = find(H(l,:));
            for index=nodes_index     %For each v nodes connected to this c nodes
                temp_index = nodes_index;
                temp_index(temp_index == index) = [];%Make a temp index to avoid taking into account the probability sent by this v nodes 
                L_r(l,index)=mod(sum(L_q(l,temp_index)),2);     
            end
        end

        %% Step 2
        %We update the v nodes with the last value and the received message
        %from c nodes and send it back to the c nodes
        L_q = zeros(m,n);
        for l=1:n       %For each v node
            %Mazjority vote
            nodes_index = find(H(:,l));
            u(l) = round( (v_nodes(l) + sum(L_r(nodes_index,l)) )/(length(nodes_index)+1) );        

            %Send value to each c nodes
            for index=nodes_index
                temp_index = nodes_index;
                temp_index(temp_index == index) = [];
                L_q(index,l) = vote - L_r(index,l);         %Don't take into account the last received prob from one c node in the new value for this node
            end
        end

        n_iter = n_iter + 1;

    end

%     while (norm(mod(new_dr*H',2)) ~=0) && (stop_count<10) %checking if syndrome equals 0 vector
%     stop_count =stop_count +1;
%     syndrome = mod(new_r*H',2);
%     for i=1:n %for each variable node
%         check_node = H(:,i); %taking the link between check and variable nodes
%         vote = new_r(i);
%         count =1;%Because the received vector is taken for the majority vote
%         for j=1:m %for each check node for this variable node
%             if check_node(j) == 1 %decision of the check node needed
%                 count =count+ 1;
%                 if syndrome(j) == 1 %changing corresponding bits if check node equals 1 from syndrome
%                     vote = vote+ mod(r(i) + 1,2);  %inversing the value because assumed to be an error
%                 else
%                     vote = vote + r(i);
%                 end
%             end
%         end
%         vote_result(i) = round(vote/count); %Vote
%     end
%     new_dr = vote_result;
%     end
    bits(1+(block-1)*m:(block-1)*m+m) = u(m+1:end); %the last bits are the transmitted ones
end