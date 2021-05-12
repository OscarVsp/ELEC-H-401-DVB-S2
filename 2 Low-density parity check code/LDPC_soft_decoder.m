function [bits] = LDPC_soft_decoder(r,H,variance,max_iter)


%% small size LDPC encoder and hard decoder

 
[m n] = size(H);
L = length(r);

bits = zeros(1,L*m/n);

for block = 1:L/n
    dr = r( 1+(block-1)*n:(block-1)*n+n);%taking a n size block
    v_nodes = -2*dr/variance;    %Initial value of v nodes
    
    
    u = ones(1,n);              %Initialize the output value

    for i=1:n
        if v_nodes(i)<0
            u(i)=1;
        else
            u(i)=0;
        end
    end
    syndrome = norm(mod(u*H',2));
    
    if syndrome ~=0
        L_q = zeros(m,n);           %Initialize the value sent from v nodes to c nodes
        

        %% Step 0
        %We take the probability at the v nodes and we send it to the c nodes.
        for l=1:n       %For each v nodes
            nodes_index = find(H(:,l));
            L_q(nodes_index,l)=v_nodes(l);
        end

        n_iter = 0;
        syndrome = norm(mod(u*H',2));
        while (n_iter < max_iter && syndrome~=0)

            n_iter + 1; %compteur


            %% Step 1
            %We compute the probability to send back to each v nodes from the
            %previous probability
            L_r = zeros(m,n);   %Reset L_r
            for l=1:m          %For each c nodes
                nodes_index = find(H(l,:));     
                for j=1:length(nodes_index)
                    temp_index = nodes_index;
                    temp_index(j) = [];%Make a temp index to avoid taking into account the probability sent by this v nodes
                    if length(L_q(1,temp_index)) ~= 0 %%%Needed ?? %%%
                        K = prod(sign(L_q(l,temp_index)));  %Compute the khi factor
                        A = min(abs(L_q(l,temp_index)));  %Compute the approx alpha factor
                        L_r(l,nodes_index(j))=K*A;     %[Check is index works here]
                    end
                end
            end

            %% Step 2
            %We update the v nodes with the last value and the received message
            %from c nodes and send it back to the c nodes
            L_q = zeros(m,n);
            for l=1:n       %For each v node
                %Mazjority vote
                nodes_index = find(H(:,l));
                vote = v_nodes(l) + sum(L_r(nodes_index,l)); 
                u(l) = vote < 0;
                %Send probability to each c nodes
                for j=1:length(nodes_index)
                    L_q(nodes_index(j),l) = vote - L_r(nodes_index(j),l);         %Don't take into account the last received prob from one c node in the new prob for this node
                end
            end
            syndrome = norm(mod(u*H',2));
            n_iter = n_iter + 1;

        end
    end
    bits(1+(block-1)*m:(block-1)*m+m) = u(m+1:end); %the last bits are the transmitted ones
end

