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

H_true = [I Pt] %true parity check matrix

G = [P I] %Generator matrix

%use mod( "operation",2) to make modulo 2 operations

%Using an arbitrary vector
d= [ 1 0 0 0 1]

u0 = mod(d*G,2)

r = mod(u0+[0 0 0 0 0 0 0 0 0 0],2) %received signal with 1 bit error

s1 = mod(u0*H_true',2); %syndrome 0 if no error 
s2= mod(r*H_true',2); % different from zero if error


%% small size LDPC hard decoder

max_iter = 10;
[m n] = size(H);
v_nodes = r;    %Initial value of v nodes

L_q = zeros(m,n);           %Initialize the value sent from v nodes to c nodes
u = ones(1,n);              %Initialize the output value

%% Step 0
%We take the value at the v nodes and we send it to the c nodes.
for l=1:n       %For each v nodes
    nodes_index = find(H(:,l));
    L_q(nodes_index,l)=v_nodes(l);
end

n_iter = 0;
while (n_iter < max_iter && norm(mod(u*H',2))~=0)
    syndrome = mod(u*H',2);
    n_iter + 1; %compteur
    
    
    %% Step 1
    %We compute the value to send back to each v nodes from the
    %previous value without taking into account the previous value from
    %that specific v-node
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
    %from c nodes and compute the next value for reach c-node, again without the previous message froml that node
    L_q = zeros(m,n);
    for l=1:n       %For each v node
        %Mazjority vote
        nodes_index = find(H(:,l));
        u(l) = round( (v_nodes(l) + sum(L_r(nodes_index,l)) )/(length(nodes_index)+1) );        

        %Send value to each c nodes
        for index=nodes_index'
            temp_index = nodes_index;
            temp_index(temp_index == index) = [];
            L_q(nodes_index,l) = round( (v_nodes(l)+sum(L_r(temp_index,l)))/(length(temp_index)+1) );         %Don't take into account the last received prob from one c node in the new value for this node
        end
        %v_nodes(l) = u(l);
    end
    u
    n_iter = n_iter + 1;
    
end
n_iter
u0
u