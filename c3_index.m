function c3_index(file1,file2,delimeter)


tic; % start timer
fileID = fopen(file1);
A= dlmread(file1,delimeter);
A2= dlmread(file2,delimeter);

a=0.5; % damping factor
n=max ( max(A(:)), max(A2(:)) ); % number of nodes
n2=size(A,1); % number of links
n3=size(A2,1);
w2=size(A,2); % Width of A
Kout=zeros(1,n); 
Kout2=zeros(1,n);
Kout3=Kout2;
Kout4=Kout2;
Pa= spalloc(double(n),double(n),10000); % P(i,j)= 0 if i dangling, 1/Kout(i) else
Pa2= spalloc(double(n),double(n),10000);
Pa3= Pa2;
Pa4= Pa2;
e= ones(1,n);

fprintf('number of nodes: %d, number of links: %d %d \n',n,n2,n3); 



end