function plotmatrices(file1,file2,delimeter)
    A= dlmread(file1,delimeter);
    A2= dlmread(file2,delimeter);
    
    maximum = max(max(max(A),max(A2)));
    x = A(:,1);
    y = A(:,2);
    x2 = A2(:,1);
    y2 = A2(:,2);
    myx = linspace(0,maximum,size(A,1));
    plot(x,48*10^6*y,'b.',x,y2,'r.')
    xlabel('id'); ylabel('value');title('Comparison');
    legend('48*10^6 * Biplex PR','H-Index')
end