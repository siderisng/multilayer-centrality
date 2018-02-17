function translate_mapped_results(file,mapping,names,delimeter)

tic;
A= dlmread(file,delimeter);
A2= dlmread(mapping,' ');
%names = dlmread(names,'	');
%%fileID = fopen(names);
%names = textscan(fileID,'%d %s %s %s','HeaderLines', 1,'delimiter','\t');
%fclose(fileID);
T = readtable(names,'Delimiter','	','Format','%d	%s	%s	%s');
T = table2cell(T);
A = num2cell(A);
A2 = num2cell(A2);

n=size(A,1); % length of A 
n2=size(A2,1); % length of A2
n3=size(names,1); % length of A2
w=size(A,2); % width of A
w2=2; % width of A2;

fid = fopen('outputs/biplex_data_final','w')

for i=1:n
    i/n*100
    q = A{i,1};
    for j=1:n2
        if A2{j,1}==q
            A{i,1} = A2{j,2};
        end
    end
    flag = 0;
    if (A{1,2}~=0)
        for k=1:91751
            if T{k,1}==A{i,1}
                firstname = T{k,3};
                lastname = T{k,2};
                A{i,3} = strcat(firstname,{' '},lastname);
                flag = 1;

            end
        end
    end
    if flag==0
        A{i,3} = '';
    end
    fprintf(fid,'%d	%2.10f	%s\n',A{i,1},A{i,2},char(A{i,3}));
end

fclose(fid)
%date=strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
%dlmwrite(strcat('outputs/h_index_final_',date),A,'delimiter','	','precision','%d')
toc; % end timer
beep;

end