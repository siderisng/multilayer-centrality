function create_mapping_arrays(file1,file2,delimeter)


tic; % start timer
A= dlmread(file1,delimeter);
A2= dlmread(file2,delimeter);

n=size(A,1); % number of links
n2=size(A2,1);
w2=size(A,2); % Width of A
w= 2;

for i=1:n
    for j=1:w
        
    end
end

for i=1:n2
    for j=1:w2
    end
end

end