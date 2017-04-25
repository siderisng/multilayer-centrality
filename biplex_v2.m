function biplex_v2(file1,file2)

%% First initializations
tic; % start timer
fileID = fopen(file1);
if nargin==1
    C = textscan(fileID,'1 %f32 %f32');
    C2 = textscan(fileID,'2 %f32 %f32');
else 
    fileID2 = fopen(file2);
    C = textscan(fileID,'%f32 %f32');
    C2 = textscan(fileID2,'%f32 %f32');
    fclose(fileID2);
end
fclose(fileID);

A = cell2mat(C);
A2 = cell2mat(C2);
a=0.85; %damping factor
n=max ( max(A(:)), max(A2(:)) ); %number of nodes
n2=size(A,1); % number of links
n3=size(A2,1);
Kout=zeros(1,n); 
Kout2=zeros(1,n);
Pa= spalloc(double(n),double(n),10000); % P(i,j)= 0 if i dangling, 1/Kout(i) else
Pa2= spalloc(double(n),double(n),10000);
e= ones(1,n);

fprintf('number of nodes: %d, number of links: %d %d \n',n,n2,n3); 

%% Initializations of Pa

maxreal=max( size(unique(A),1), size(unique(A2),1));
v = zeros(n, 1) ;   % personalization vector
v2 = zeros(n, 1);
for i = 1:n2
    x=A(i,1);
    Kout(x)= Kout(x)+1; % how many outgoing links
end
for i = 1:n3
    x=A2(i,1);
    Kout2(x)= Kout2(x)+1; % how many outgoing links
end
for i = 1:n2
    x=A(i,1);
    y=A(i,2);
    Pa(x,y)=1/Kout(x); % there is a link x->y
end
for i = 1:n3
    x=A2(i,1);
    y=A2(i,2);
    Pa2(x,y)=1/Kout2(x); % there is a link x->y
end

fprintf('maxReal is %d \n', maxreal);

%% third

uni=unique(A);
uni2=unique(A2);
RealCounter= size(uni,1);
RealCounter2= size(uni2,1);


fprintf('Realcounters: %d %d \n',RealCounter,RealCounter2);


%% fourth
NoRealCounter=n-RealCounter;
NoRealCounter2=n-RealCounter2;
fprintf('NoRealcounters: %d %d \n',NoRealCounter,NoRealCounter2);

%% extension to work with dangling nodes
dangling=0;
dangling2=0;
d=zeros(double(n),1);
for i=1:n
    if Kout(i)==0
        d(i,1)=1;
        dangling= dangling+1;
    end
end
u=0;
for i=1:n
    if d(i,1)==1
        u=1/dangling;
    end
    break;
end
clearvars A A2 C C2 file1 fileId x y
whos 

% du=d*u';
% clearvars d u
d2=zeros(double(n),1);
for i=1:n
    if Kout2(i)==0
        d2(i,1)=1;
        dangling2= dangling2+1;
    end
end
u2=0;
for i=1:n
    if d2(i,1)==1
        u2=1/dangling2;
        break
    end
end
 
% du2=d2*u2';
% clearvars d2 u2

% if dangling~=0
%     Pa= Pa + du;
% end
% 
% if dangling2~=0
%     Pa2= Pa2 + du2;
% end
whos

%% initialization of personalization vectors
v=zeros(n,1);
v2=v;


for i=1:RealCounter
    if n==RealCounter
        v(i)= 1/maxreal;
    else
        x=uni(i);
        v(x)=a/maxreal;
    end
end
for i=1:RealCounter2
    if n==RealCounter2
        v2(i)= 1/maxreal;
    else
        x=uni2(i);
        v2(x)=a/maxreal;
    end
end


summ= 1 - sum(v(:,1));
summ2= 1 - sum(v2(:,1));
for i=1:n
    if v(i)==0
        v(i)=summ/NoRealCounter;
    end 
    if v2(i)==0
        v2(i)=summ2/NoRealCounter2;
    end 
end
v;
v2;
sum(v(:,1));
sum(v2(:,1));



%% preparing for iterative


e=ones(n,1); % e vector
tol=1e-10;   % tolerance
pu=(1/(2*n)) * ones(n, 1); % starting values for pagerank
pu= double(pu);
pu2=pu; pd=pu; pd2=pu;
pi=(pu+pu2+pd+pd2)/2;
sum(pi(:,1));
i=0;
dp=1;
%v=(1/n)*ones(n,1);
%v2=v;

%% iterative procedure

while (dp >= tol)   % computation of pagerank, stops when pi converges
    prevpi=(pu+pu2+pd+pd2)/2;
    prevpu=pu;
    prevpu2=pu2;
    prevpd2=pd2;
    prevpd=pd;
    test=prevpu'*a;
    test2=prevpu2'*a;
    C=zeros(1,n);
    D=C;
	for i = 1:1
        for j = 1:n
            P=full(Pa(:,j));
            P2=full(Pa2(:,j));
            if (d(j)~=1)
                for k = 1:n
                    C(i,j) = C(i,j) + test(i,k)*P(k);
                end
            else
                for k= 1:n
                    if (d(k)==1)
                        C(i,j) = C(i,j) + test(i,k)*u;
                    else
                        C(i,j) = C(i,j) + test(i,k)*P(k);
                    end
                end
            end
            if (d2(j)~=1)
                for k = 1:n
                    D(i,j) = D(i,j) + test2(i,k)*P2(k);
                end
            else
                for k = 1:n
                    if (d2(k)==1)
                        D(i,j) = D(i,j) + test2(i,k)*u2;
                    else
                        D(i,j) = D(i,j) + test2(i,k)*P2(k);
                    end
                end
            end
        end
	end

    pu = [C + prevpu2' + 2*a*prevpd2']'/2;
    pu2 = [prevpu' + D + 2*a*prevpd2']'/2;
    pd = [(1-a)*(prevpu' + prevpd'*e*v' + prevpd2'*e*v')]'/2;
    pd2 = [(1-a)*(prevpu2' + prevpd'*e*v2' + prevpd2'*e*v2')]'/2;
    pi=(pu+pu2+pd+pd2)/2;
    dp=norm(pi-prevpi,1) % difference of this iteration compared to the previous
    i=i+1;
end

%% results

fprintf('Total number of nodes: %d\n',n);
fprintf('Total number of links: %d\n',n2+n3);
fprintf('Total number of iterations: %d\n',i);
fprintf('Final pagerank:\n')
% disp(pi')
fprintf('Sum of elements:%f\n',sum(pi));
fprintf('Final Ranking:\n \tValue: \t  Node:\n')

c=flipud((sort(pi))); 
top= min(n,10);
result=c(1:top);         %top ten
ind=find(pi>=c(top));    %their indices
result=flipud(sortrows([pi(ind) ind],1));
for i=1:top
    fprintf('\t%f32\t%d \n',result(i,1),result(i,2));
end

toc; % end timer
beep; % sound when finished