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
                L_q(nodes_index,l) = round(v_nodes(l)+sum(L_r(temp_index,l)))/(length(temp_index)+1);        %Don't take into account the last received prob from one c node in the new value for this node
            end
        end

        n_iter = n_iter + 1;

    end
    bits(1+(block-1)*m:(block-1)*m+m) = u(m+1:end); %the last bits are the transmitted ones
end