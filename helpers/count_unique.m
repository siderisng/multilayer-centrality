function count_unique(file1,file2,delimeter)


tic; % start timer
A= dlmread(file1,delimeter);

A = unique(A);
A2= dlmread(file2,delimeter);

A2 = unique(A2);
A= unique(vertcat(A,A2));
size(A)

n=size(A,1) % number of links


for i=1:n
    B(i,1) = i;
    B(i,2) = A(i);
end

A= dlmread(file1,delimeter);
A2= dlmread(file2,delimeter);

n=size(A,1); % number of links
n2=size(A2,1);
w2=size(A2,2); % Width of A
w=size(A,2);
n3=size(B,1);



for i=1:n
    i/n*100
    for j=1:w
        if (A(i,j)~=0)
            for k=1:n3
                if (A(i,j)==B(k,2))
                    A(i,j)=B(k,1);
                end
            end
        end
    end
end

for i=1:n2
    i/n2*100
    for j=1:w2
        if (A2(i,j)~=0)
            for k=1:n3
                if (A2(i,j)==B(k,2))
                    A2(i,j)=B(k,1);
                end
            end
        end
    end
end

dlmwrite('paper_authors_filtered_mapped',A,' ')
dlmwrite('citations_filtered_mapped',A2,' ')
dlmwrite('mapping',B,' ')

toc;

end