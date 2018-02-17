function filter500(file1,file2,file3,delimeter)


tic; % start timer

top500=dlmread(file3,delimeter);

A= dlmread(file1,delimeter);


A2= dlmread(file2,delimeter);

unA2 = unique(A2);



n=size(A,1);
w=size(A,2);

toDelete=[];

for i=1:n
    i/n*100
    flag=1;
    for j=2:w
        if (A(i,j)==0)
            continue;
        end
        if (flag==0)
            break;
        end
        if ((ismember(A(i,j), top500(:)))==0)
            flag=0;
            toDelete(end+1)=i;
        end
    end 
end

A(toDelete,:)=[];
unA = unique(A);

n=size(A2,1);
w=size(A2,2);

toDelete=[];
for i=1:n
    i/n*100
    flag=1;
    for j=2:w
        if (A2(i,j)==0)
            continue;
        end
        if (flag==0)
            break;
        end
        if (~ismember(A2(i,j), unA(:)))
            flag=0;
            toDelete(end+1)=i;
        end
    end 
end

A2(toDelete,:)=[];

dlmwrite('paper_authors_filtered',A,' ')
dlmwrite('citations_filtered',A2,' ')

toc;

end