function h_index(file1,file2)

tic; % start timer
fileID = fopen(file1,'r');

n=10;
papers=zeros(1,n);
A = dlmread(file1,' ');

nR=max(A(:));
n2=size(A,1);
n3=size(A,2);
NumOfPapers=zeros(1,nR);

for i=1:n2
   for j=2:n3
       if A(i,j)==0
           continue;
       else
           NumOfPapers(A(i,j))= NumOfPapers(A(i,j)) + 1 ;
       end
   end
end

NumOfPapers

B= dlmread(file2,' ')
nP=max(B(:));
nC=size(B,1);
NumOfCitations = zeros(1,nP);

for i=1:nC
    NumOfCitations(B(i,2))=NumOfCitations(B(i,2)) +1;
end

NumOfCitations
n5=max(NumOfPapers(:));
CPR=zeros(nR,n5); %% citations for every paper of an author (authors are on rows) 
for i=1:n2
    for j=2:n3
      if A(i,j)==0
        continue
      else
          for k=1:n5
              if CPR(A(i,j),k)==0
                CPR(A(i,j),k)= NumOfCitations(A(i,1));
                break
              end
          end
      end
    end
end


CPR
CPR= sort(CPR,2,'descend')

hIndex=zeros(1,nR);
for i=1:nR
    for j=1:n5
        if j==n5
            hIndex(1,i)=j
            break
        end
        if CPR(i,j)>=j
            continue
        else
            hIndex(1,i)=max(j-1,0);
            break
        end
    end
end

res=hIndex';
fprintf('Final Ranking:\n \th-Index: \tAuthor:\n')
c=flipud((sort(res))); 
top= min(nR,10);
result=c(1:top);         %top ten
ind=find(res>=c(top));    %their indices
result=flipud(sortrows([res(ind) ind],1));
for i=1:top
    fprintf('\t%f\t%d \n',result(i,1),result(i,2));
end