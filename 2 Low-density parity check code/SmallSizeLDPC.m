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
%% Step 0 
r = mod(u0+[0 0 0 1 0 0 0 0 0 0],2) %received signal with 1 bit error

s1 = mod(u0*H_true',2); %syndrome 0 if no error 
s2= mod(r*H_true',2); % different from zero if error

%% Step 1 

%Receiving messages from check nodes and vote the result 
% [m n] = size(H_true);
% vote_result = zeros(1,n);
% stop_count=0;
% new_r = r;
% while (norm(mod(new_r*H_true',2)) ~=0) && (stop_count<10)
%     stop_count =stop_count +1
%     syndrome = mod(new_r*H_true',2)
%     for i=1:n %for each variable node
%         check_node = H_true(:,i); %taking the link between check and variable nodes
%         vote = r(i);
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
%     new_r = vote_result
% end
% u
%% small size LDPC encoder and hard decoder

 max_iter = 10;
[m n] = size(H);
v_nodes = r;    %Initial value of v nodes

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
            L_q(index,l) = u(index)-         %Don't take into account the last received prob from one c node in the new value for this node
        end
    end
    
    n_iter = n_iter + 1;
    
end
n_iter
u0
u