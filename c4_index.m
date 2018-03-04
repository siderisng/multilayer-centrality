function c4_index(coauthorship,citations,bi_results,delimiter)
    

tic; % start timer
A= dlmread(coauthorship,delimiter);
A2= dlmread(citations,delimiter);

%% parse biplex results
n=max ( max(A(:)), max(A2(:)) ); % number of nodes
B = dlmread(bi_results,'	');
B = sortrows(B,1);
pi = zeros(n,1);
n2 = size(B,1)
for i=1:n2
    if B(i,1)~=0
        pi(B(i,1))=B(i,2);
    end
end
mean1 = mean(pi);
dev = std(pi);
pi = (pi-mean1)./dev;


%% C3 INDEX



    
a=0; % credit distribution among author, 0 for uniform, >0 for more credit for author with bigger c3
theta= 0.5; % damping factor, suggested 0.5

n2=size(A,1); % number of links
n3=size(A2,1);
w2=size(A,2); % Width of A
Kout=zeros(1,n); 
Kout2=zeros(1,n);
Kout3=Kout2;
Pa= spalloc(double(n),double(n),10000^2); % P(i,j)= 0 if i dangling, 1/Kout(i) else
Pa2= spalloc(double(n),double(n),10000^2);
Pa3= Pa2;

fprintf('number of nodes: %d, number of links: %d %d \n',n,n2,n3); 

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
                else
                    Pa2(x,y)= Pa2(x,y) + 1;
                    Kout2(x)= Kout2(x) + 1;
                end
            end
        end
    end
end

for i = 1:n3
    i/n3*100
    x=A2(i,1);
    y=A2(i,2);
    Pa3(x,y)=1; % there is a link x->y
    Kout3(x)=Kout3(x)+1;
    for j = 1:n2
        if A(j,1)==x
            for k = 1:n2
                if A(k,1)==y
                    % j is the index of first paper , k of second
                    % we create links from every author of j to every
                    % author of k
                    for m = 2:w2
                        x2= A(j,m);
                        if x2==0
                            break;
                        else
                            for l = 2:w2
                                y2= A(k,l);
                                if y2==0
                                    break;
                                elseif x2==y2
                                    continue;
                                else
                                    Pa(x2,y2)=Pa(x2,y2)+1;
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

PQI=zeros(1,n);


prevPQI=PQI;
for j = 1:n 
    sum = 0;
    Pa3k = full(Pa3(:,j));
    for k= 1:n      
        if Pa3k(k)~=0 % find a paper that cites j'th paper
            sum = sum + prevPQI(k)/Kout3(k);
        end
    end
    PQI(j) = (1-theta) + theta* sum;
end



PQI;
PCI= zeros(1,n);
C4= zeros(1,n);
AAI= zeros(1,n);
ACI= zeros(1,n);
dp=10000;
tol=1e-7;   % tolerance
iterations=0;
while (dp >= tol)
   prevPQI=PQI;
   prevC4=C4;
   prevAAI=AAI;
   prevACI=ACI;
   for j= 1:n
       Pa3k = full(Pa3(:,j));
       sum = 0;
       for k= 1:n
           if Pa3k(k)~=0 % find a paper that cites i'th paper
               sum = sum + prevPQI(k)/Kout3(k);
           end
       end
       PQI(j) = (1-theta) + theta* sum;
       sum = 0;
       for i = 1:n2
          for k=2:w2
              if A(i,k)==j % paper written by j                  
                  sumsum=0;
                  for l=2:w2
                      author=A(i,l);
                      if author~=0
                        sumsum= sumsum+ prevC4(author)^a;
                      end
                  end
                  sum = sum + prevPQI(i)/sumsum;
              end
          end
       end
       PCI(j)=prevC4(j)^a * sum;
       
       sum=0;
       Pa2k= full(Pa2(:,j));
       for k=1:n
           if Pa2k(k)>0
               sum = sum + prevAAI(k)/Kout2(k);
           end
       end
       AAI(j)=sum;
       
       sum=0;
       Pak= full(Pa(:,j));
       for k=1:n
           if Pak(k)>0
               sum = sum + prevACI(k)/Kout(k);
           end
       end
       ACI(j)= (1-theta) + theta*sum;
       C4(j) = (1-theta) + theta*(ACI(j)+AAI(j)+PCI(j)+pi(j));
   end
   
   
   dp=norm(C4-prevC4,1) % difference of this iteration compared to the previous one
   iterations= iterations+1;
end

fprintf('Total number of iterations: %d\n',iterations);
fprintf('Final C^4 index:\n')
disp(C4)

C4=C4';

dlmwrite('c4_result_500',C4,' ');
c=flipud((sort(C4))); 
%top= min(n,10);
top = n;
result=c(1:top);         %top ten
ind=find(C4>=c(top));    %their indices
result=flipud(sortrows([ind C4(ind)],2));
for i=1:top
    %fprintf('\t%f32\t%d \n',result(i,1),result(i,2));
end

date=strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
file=strrep(strrep(coauthorship,'tests/matlab/',''),'.txt','');
dlmwrite(strcat('outputs/c4_index_',file,'_',date),result,'	')

toc;
beep;



end