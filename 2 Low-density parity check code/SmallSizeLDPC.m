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
d= [ 1 1 1 0 1]

u = mod(d*G,2)
%% Step 0 
r = mod(u+[0 0 0 0 0 0 0 1 0 0],2) %received signal with 1 bit error

s1 = mod(u*H_true',2); %syndrome 0 if no error 
s2= mod(r*H_true',2); % different from zero if error

%% Step 1 

%Receiving messages from check nodes and vote the result 
[m n] = size(H_true);
vote_result = zeros(1,n);
stop_count=0;
new_r = r;
while (norm(mod(new_r*H_true',2)) ~=0) && (stop_count<10)
    stop_count =stop_count +1
    syndrome = mod(new_r*H_true',2)
    for i=1:n %for each variable node
        check_node = H_true(:,i); %taking the link between check and variable nodes
        vote = r(i);
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
    new_r = vote_result
end
u
