function plotmatrices(file1,file2,file3,file4,delimeter)
    A= dlmread(file1,delimeter);
    A2= dlmread(file2,delimeter);
    A3= dlmread(file3,delimeter);
    A4= dlmread(file3,delimeter);
    
    maximum = max(max(max(A),max(max(A2),max(A3))));
    x = A(:,1);
    y = A(:,2);
    x2 = A2(:,1);
    y2 = A2(:,2);
    y3 = A3(:,2);
    y4 = A4(:,2);
    mean1 = mean(y);
    mean2 = mean(y2);
    mean3 = mean(y3);
    mean4 = mean(y4);
    dev = std(y);
    dev2 = std(y2);
    dev3 = std(y3);
    dev4 = std(y4);
    f = (y-mean1)./dev;
    y = arrayfun(@(x) (x-mean1)/dev,y);
    y2 = arrayfun(@(x) (x-mean2)/dev2,y2);
    y3 = arrayfun(@(x) (x-mean3)/dev3,y3);
    y4 = arrayfun(@(x) (x-mean4)/dev4,y3);
    plot(x,y3,'b.',x,y4,'r.')
    xlabel('id'); ylabel('value');title('Standard score');
    legend('C3-Index','C4-Index')
    set(gcf,'renderer','painters')
    print('c3-c4','-depsc','-r300');
end