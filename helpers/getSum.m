function getSum(file,delimeter)
    A= dlmread(file,delimeter);
    tot = sum(A)
    fprintf('Sum of values is: %f\n\n',tot(2));
end