function [u] = SoftDecoder(H,r,variance,max_iter)


%% small size LDPC encoder and hard decoder

 
[m n] = size(H);
v_nodes = -2*r/variance    %Initial value of v nodes

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
    syndrome = norm(mod(u*H',2))
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
            K = prod(sign(L_q(l,temp_index)));  %Compute the khi factor
            A = min(abs(L_q(l,temp_index)));  %Compute the approx alpha factor
            L_r(l,index)=K*A;     %[Check is index works here]
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
        if vote<0
            u(l)=1;
        else
            u(l)=0;
        end
        %Send probability to each c nodes
        for index=nodes_index
            L_q(index,l) = vote - L_r(index,l);         %Don't take into account the last received prob from one c node in the new prob for this node
        end
    end
    
    n_iter = n_iter + 1;

end
n_iter

