%% small size LDPC encoder and hard decoder

%Parity check sparse matrix given in the project
H= [1 1 0 1 1 0 0 1 0 0;
    0 1 1 0 1 1 1 0 0 0;
    0 0 0 1 0 0 0 1 1 1;
    1 1 0 0 0 1 1 0 1 0;
    0 0 1 0 0 1 0 1 0 1]; 
%We need to have an identity matrix on the left side to deduce the parity
%matrix

I = eye(5);
Pt = [0 1 1 1 0;
      1 0 1 0 0;
      1 0 1 0 1;
      0 0 1 1 1;
      1 1 0 0 1]; %hand computed modulo 2 with row exchange allowed
P = Pt';

H_true = [I Pt]; %true parity check matrix

G = [P I]; %Generator matrix


%use mod( "operation",2) to make modulo 2 operations

%Using an arbitrary vector
d= [ 1 0 0 0 1]

u0 = mod(d*G,2)

r = mod(u0+[0 0 0 0 0 0 0 0 0 1],2) %received signal with 1 bit error

%% small size LDPC hard decoder

max_iter = 10;
[m n] = size(H);
v_nodes = r;    %Initial value of v nodes

L_q = zeros(m,n);           %Initialize the value sent from v nodes to c nodes
u = v_nodes;              %Initialize the output value

%% Step 0
%We take the value at the v nodes and we send it to the c nodes.
for l=1:n       %For each v nodes
    nodes_index = find(H(:,l));
    L_q(nodes_index,l)=v_nodes(l);
end
L_q
n_iter = 0;
syndrome = norm(mod(u*H',2));
while (n_iter < max_iter && syndrome~=0)
    %% Step 1
    %We compute the probability to send back to each v nodes from the
    %previous probability
    L_r = zeros(m,n);   %Reset L_r
    for l=1:m          %For each c nodes
        nodes_index = find(H(l,:));
        for j=1:length(nodes_index)
            index = nodes_index(j);
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
        vote_u = (v_nodes(l) + sum(L_r(nodes_index,l)) )/(length(nodes_index)+1); 
        if vote_u>0.5
            u(l) = 1;
        else 
            u(l) = 0;
        end
        
        %Send value to each c nodes
        for j=1:length(nodes_index)
            index = nodes_index(j);
            temp_index = nodes_index;
            temp_index(temp_index == index) = [];
            sum_L_r = (v_nodes(l)+sum(L_r(temp_index,l)));
            vote_L_q = (v_nodes(l)+sum(L_r(temp_index,l)))/(length(temp_index)+1);
            if vote_L_q > 0.5
                L_q(index,l) = 1;
            else
                L_q(index,l) = 0;
            end
        
            
            %Don't take into account the last received prob from one c node in the new value for this node  
        end
        %v_nodes(l) = u(l);  %Not Sure
    end
    syndrome = norm(mod(u*H',2));
    n_iter = n_iter + 1;
end
n_iter
u0
u