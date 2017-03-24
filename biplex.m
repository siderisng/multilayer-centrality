tic; % start timer
fileID = fopen('lesmis2.txt');
C = textscan(fileID,'%f32 %f32');
fclose(fileID);
A = cell2mat(C);
a=0.85; %damping factor
n=max(A(:)); %number of nodes
n2=size(A,1); % number of links
Kout=zeros(1,n); 
Pa= zeros(n,n); % P(i,j)= 0 if i dangling, 1/Kout(i) else
e= ones(1,n);

v = (1/n) * ones(n, 1); % personalization vector
for i = 1:n2
    x=A(i,1);
    Kout(x)= Kout(x)+1; % how many outgoing links
end
for i = 1:n2
    x=A(i,1);
    y=A(i,2);
    Pa(x,y)=1; % there is a link x->y
end
for i = 1:n
    for j=1:n
        if Pa(i,j)~=0
            Pa(i,j)=Pa(i,j)/Kout(i);
        end
    end
end 


d=zeros(n,1);
u=zeros(n,1);
dangling=0;

for i=1:n
    if Kout(i)==0
        d(i,1)=1;
        dangling= dangling+1;
    end
end

for i=1:n
    if d(i,1)==1
        u(i,1)=1/dangling;
    end
end

if dangling~=0
    Pa= Pa + d*u';
end


e=ones(n,1); % e vector
tol=1e-9;   % tolerance
pu=(1/(2*n)) * ones(n, 1); % starting values for pagerank
pd=(1/(2*n)) * ones(n, 1);
i=0;
dp=1;
while (dp >= tol)   % computation of pagerank, stops when pi converges
    prevpi=pu+pd;
    prevpu=pu;
    prevpd=pd;
    pu = [prevpu'*a*Pa + prevpd'*a]';
    pd = [(1-a)*prevpu' + (1-a)*prevpd'*e*v']';
    pi=pu+pd;
    dp=norm(pi-prevpi,1); % diffence of this iteration compared to the previous
    i=i+1;
end
fprintf('Total number of nodes: %d\n',n);
fprintf('Total number of links: %d\n',n2);
fprintf('Total number of iterations: %d\n',i);
fprintf('Final pagerank:\n')
disp(pi')
fprintf('Sum of pagerank: %f \n', sum(pi(:,1)));
fprintf('Final Ranking:\n \tValue: \t  Index:\n')
top= min(n,10);
c=flipud(unique(sort(pi)));    
result=c(1:top);         %top ten
ind=find(pi>=c(top));    %their indices
result=flipud(sortrows([pi(ind) ind],1));
disp(result)


toc; % end timer
beep; % sound when finished