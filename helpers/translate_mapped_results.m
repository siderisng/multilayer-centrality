function translate_mapped_results(file,mapping,names,delimeter)

tic;
A= dlmread(file,delimeter);
A2= dlmread(mapping,' ');
T = readtable(names,'Delimiter','	','Format','%d	%s	%s	%s');
T = table2cell(T);
A = num2cell(A);
A2 = num2cell(A2);

n=size(A,1); % length of A 
n2=size(A2,1); % length of A2
n3=size(names,1); % length of A2
w=size(A,2); % width of A
w2=2; % width of A2;

fid = fopen('outputs/c4_index_final','w')

for i=1:n
    i/n*100
    q = A{i,1};
    for j=1:n2
        if A2{j,1}==q
            A{i,1} = A2{j,2};
        end
    end
    fprintf(fid,'%d	%2.10f\n',A{i,1},A{i,2});
end

fclose(fid)
toc; % end timer
beep;

end