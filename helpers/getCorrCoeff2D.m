function getCorrCoeff2D(file1,file2,file3,file4,delimeter)
    A= dlmread(file1,delimeter); %biplex
    A2= dlmread(file2,delimeter); %h-index
    A3= dlmread(file3,delimeter); %c3-index
    A4= dlmread(file4,delimeter); %c4-index
        
    n=size(A,1); % length of A 
    n2=size(A2,1); % length of A2
    n3=size(A3,1); % length of A3
    n4=size(A4,2); % length of A4
    for i=1:n
        if (A(i,1)~=0)
            v(A(i,1))=A(i,2);
        end
    end
    
    for i=1:n2
        if (A2(i,1)~=0)
            v2(A2(i,1))=A2(i,2);
        end
    end
    
    for i=1:n3
        if (A3(i,1)~=0)
            v3(A3(i,1))=A3(i,2);
        end
    end
    
    for i=1:n4
        v4(i) = A4(1,i);
    end
    
    
    fprintf('Correlation between Biplex PR and C3-index is: %f\n\n',corr2(v,v4))
    fprintf('Correlation between H-Index and C3-index is: %f\n\n',corr2(v2,v4))
    fprintf('Correlation between C3-Index and C4-index is: %f\n\n',corr2(v3,v4))
end