function keepOnlyAuthors(file,citations)

tic; % start timer
A= dlmread(file,'	');
citations = dlmread(citations,' ');

n = size(A,1);
fid = fopen('outputs/biplex_only_authors_filtered_mapped','w');
for i=1:n
    if ismember(A(i,1),citations)==0
        fprintf('%d : IS AUTHOR!\n',A(i,1));
        fprintf(fid,'%d	%2.10f\n',A(i,1),A(i,2));
    end
end
fclose(fid);
toc;
beep;

end