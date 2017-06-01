function biplex_data(file1,file2,delimeter)

%% First initializations
tic; % start timer
fileID = fopen(file1);
A= dlmread(file1,delimeter);
A2= dlmread(file2,delimeter);

a=0.85; % damping factor
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

%% Initializations of Pa

maxreal=max( size(unique(nonzeros(A)),1), size(unique(nonzeros(A2)),1))
maxreal=maxreal-1; % remove zeros
v = zeros(n,1) ;   % personalization vector
v2 = zeros(n,1);

for i=1:n2
    i/n2*100
    for j=2:w2
        x=A(i,j);
        if x==0
           break;
        else
            for k=2:w2
                y=A(i,k);
                if y==0
                    break;
                elseif x==y
                    continue;
                elseif Pa(x,y)==1
                    continue;
                else
                    Pa(x,y)= 1;
                    Kout(x)= Kout(x) + 1;
                end
            end
            y= A(i,1);
            Pa3(x,y)=1;
            Pa4(y,x)=1;
            Kout3(x)= Kout3(x) + 1;
            Kout4(y)= Kout4(y) + 1;
        end
    end
end
for i = 1:n3
    x=A2(i,1);
    Kout2(x)= Kout2(x)+1; % how many outgoing links
end
for i = 1:n3
    i/n3*100
    x=A2(i,1);
    y=A2(i,2);
    Pa2(x,y)=1/Kout2(x); % there is a link x->y
    for j = 1:n2
        if A(j,1)==x
            for k = 1:n3
                if A(k,1)==y
                    % j is the index of first paper , k of second
                    % we create links from every author of j to every
                    % author of k
                    for m = 2:w2
                        x2= A(j,m);
                        if x2==0
                            break;
                        else
                            for n = 2:w2
                                y2= A(k,n);
                                if y2==0
                                    break;
                                elseif x2==y2
                                    continue;
                                elseif Pa(x2,y2)==1
                                    continue;
                                else
                                    Pa(x2,y2)=1;
                                    Kout(x2)=Kout(x2) + 1;
                                end
                            end
                        end          
                    end
                end
            end
        end
    end
end

for i = 1:n
    for j = 1:n
        if Pa(i,j)==1
            Pa(i,j)=Pa(i,j)/Kout(i);
        end
        if Pa3(i,j)==1
            Pa3(i,j)=Pa3(i,j)/Kout3(i);
        end
        if Pa4(i,j)==1
            Pa4(i,j)=Pa4(i,j)/Kout4(i);
        end

    end
end

Pa
Pa2
Pa3

% Kout4 = fliplr(Kout3);


Kout
Kout2
Kout3
Kout4
% file=strrep(strrep(file1,'tests/matlab/',''),'.txt','');
% file
% dlmwrite(strcat('outputs/Pa_',file,'_'),Pa);
% dlmwrite(strcat('Pa2_',file,'_'),Pa2,'	');
% dlmwrite(strcat('Pa3_',file,'_'),Pa3,'	');
% dlmwrite(strcat('Pa4_',file,'_'),Pa4,'	');
% 
% dlmwrite(strcat('Kout_',file,'_'),Kout,'	');
% dlmwrite(strcat('Kout2_',file,'_'),Kout2,'	');
% dlmwrite(strcat('Kout3_',file,'_'),Kout3,'	');
% dlmwrite(strcat('Kout4_',file,'_'),Kout4,'	');
% 
% fprintf('maxReal is %d \n', maxreal);

%% third

uni=unique(nonzeros(A));
uni2=unique(nonzeros(A2));
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
dangling3=0;
dangling4=0;
d=zeros(double(n),1);
for i=1:n
    if Kout(i)==0
        d(i,1)=1;
        dangling= dangling+1;
    end
end
u=zeros(n,1);
for i=1:n
    if d(i,1)==1
        u(i)=1/dangling;
    end
end
clearvars A A2 C C2 fileId x y
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
u2=zeros(n,1);
for i=1:n
    if d2(i,1)==1
        u2(i)=1/dangling2;
        
    end
end

d3=zeros(double(n),1);
for i=1:n
    if Kout3(i)==0
        d3(i,1)=1;
        dangling3= dangling3+1;
    end
end
u3=zeros(n,1);
for i=1:n
    if d3(i,1)==1
        u3(i)=1/dangling3;
        
    end
end


d4=zeros(double(n),1);
for i=1:n
    if Kout4(i)==0
        d4(i,1)=1;
        dangling4= dangling4+1;
    end
end
u4=zeros(n,1);
for i=1:n
    if d4(i,1)==1
        u4(i)=1/dangling4;
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
% whos

% %% initialization of personalization vectors
% v=zeros(n,1);
% v2=v;
% 
% 
% for i=1:RealCounter
%     if n==RealCounter
%         v(i)= 1/maxreal;
%     else
%         x=uni(i);
%         v(x)=a/maxreal;
%     end
% end
% for i=1:RealCounter2
%     if n==RealCounter2
%         v2(i)= 1/maxreal;
%     else
%         x=uni2(i);
%         v2(x)=a/maxreal;
%     end
% end
% 
% 
% summ= 1 - sum(v(:,1));
% summ2= 1 - sum(v2(:,1));
% for i=1:n
%     if v(i)==0
%         v(i)=summ/NoRealCounter;
%     end 
%     if v2(i)==0
%         v2(i)=summ2/NoRealCounter2;
%     end 
% end
% v;
% v2;
% sum(v(:,1));
% sum(v2(:,1));



%% preparing for iterative


e=ones(n,1); % e vector
tol=1e-7;   % tolerance
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
P=sparse(1,double(n));
P2=sparse(1,double(n));

while (dp >= tol)   % computation of pagerank, stops when pi converges
    prevpi=pi;
    prevpu=pu;
    prevpu2=pu2;
    prevpd2=pd2;
    prevpd=pd;
    test=prevpu'*a;
    test2=prevpu2'*a;
    test3=prevpd';
    test4=prevpd2';

    parfor j = 1:n
%             fprintf('start %d\n',j);
%             tic;
%         P=(Pa(:,j));
%         P2=(Pa2(:,j));
% %             toc;
% %             tic;
%         P = P + d(j)*u;
%         P2 = P2 + d2(j)*u2;
%             toc;
%             tic;

        C(j) = test*(Pa(:,j) + d(j)*u);
        D(j) = test2*(Pa2(:,j)+ d2(j)*u2);
        E(j) = test3*(Pa3(:,j) + d3(j)*u3);
        F(j) = test4*(Pa3(:,j) + d3(j)*u3);
        G(j) = test3*(Pa4(:,j) + d4(j)*u4);
        H(j) = test4*(Pa4(:,j) + d4(j)*u4);
%             toc;
    end
    

    pu = [C + prevpu2' + 2*a*prevpd2']'/2;
    pu2 = [prevpu' + D + 2*a*prevpd2']'/2;
    pd = [(1-a)*(prevpu' + E + F)]'/2;
    pd2 = [(1-a)*(prevpu2' + G + H)]'/2;
    pi=(pu+pu2+pd+pd2)/2;
    dp=norm(pi-prevpi,1) % difference of this iteration compared to the previous one
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

date=strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
file=strrep(strrep(file1,'tests/matlab/',''),'.txt','');
dlmwrite(strcat('outputs/biplex_data_',file,'_',date),result,'	')

toc; % end timer
beep; % sound when finished