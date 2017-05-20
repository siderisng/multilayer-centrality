function h_index(file1,file2,delimeter)

tic; % start timer
fileID = fopen(file1,'r');

n=10;
papers=zeros(1,n);
A= dlmread(file1,delimeter);

nR=max(A(:));
n2=size(A,1);
n3=size(A,2);
NumOfPapers=zeros(1,nR);
fprintf('parsed file 1\n');
for i=1:n2
   for j=2:n3
       if A(i,j)==0
           continue;
       else
           NumOfPapers(A(i,j))= NumOfPapers(A(i,j)) + 1 ; % initialize number or papers for every author
       end
   end
end

NumOfPapers;
fprintf('initialized NumOfPapers\n');
B= dlmread(file2,delimeter);
fprintf('parsed file 2\n');
nP=max(B(:));
fprintf('number of papers: %d , number of researchers: %d \n', nP,nR);
nC=size(B,1);
NumOfCitations = zeros(1,max(nR,nP));

for i=1:nC
    NumOfCitations(1,B(i,2))=NumOfCitations(1,B(i,2)) +1; % citations for every paper
end
fprintf('initialized NumOfCitations\n');
NumOfCitations;
n5=max(NumOfPapers(:));
CPR=sparse(double(nR),double(n5)); %% citations for every paper of an author (authors are on rows) 
for i=1:n2
    i/n2*100
    for j=2:n3
      if A(i,j)==0
        continue
      else
          for k=1:n5
              k;
              if CPR(A(i,j),k)==0
                CPR(A(i,j),k)= NumOfCitations(1,A(i,1));
                break
              end
          end
      end
    end
end

fprintf('initialized CPR\n');
CPR;
CPR= sort(CPR,2,'descend');
fprintf('sorted CPR\n');
hIndex=zeros(1,nR);
for i=1:nR  % calculate h-index
    i/nR*100
    for j=1:n5
        if j==n5
            hIndex(1,i)=j;
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
    fprintf('\t%d\t\t\t%d \n',result(i,1),result(i,2));
end
date=strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
dlmwrite(strcat('outputs/h_index_',date),result,'	')
toc; % end timer
beep;