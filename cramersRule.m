function [x, determinants] = cramersRule(textFile)
%% CRAMERS RULE FOR SOLVING SYSTEMS OF LINEAR EQUATIONS %%
%takes a text file for input, textFile
%textFile contains (in order), the dimension of the matrix A, n,
%the cell values of A(reading off all cells from row 1, then row 2,...,row n),
%and the values of the vector, b, for Ax = b (reading off b-values the same as A)

%open file
fileID = fopen(textFile);
%set specification to floating point numbers
fileSpec = '%f';
%read in all contents of file into row vector
A = fscanf(fileID, fileSpec, [1 inf]);
%extract n, remove from A
n = A(1);
A(1) = [];
%convert A to matrix from vector form, will have b as bottom row
A = vec2mat(A, n);
%extract b, force into column vector, remove b from A
b = (A(n+1,:))';
A(n+1,:) = [];

%get determinant of A, add to return list of determinants
determinants = zeros(n+1,1);
determinants(1) = det(A);
%% MAIN LOOP(S) %%
%compute determinate for every Aj
for j = 1:n
    Aj = A;
    Aj(:,j) = b;
    %add to (j+1)th spot of determinant list for order
    determinants(j+1) = det(Aj);
end
%use determinants to compute all xi 
x = zeros(n,1);
for i = 2:n+1
    x(i-1) = determinants(i)/determinants(1);
end
%display for user
disp("Determinant A = " + determinants(1));
for k = 2:n+1
    disp("Determinant A"+(k-1) + " = " + determinants(k));
end
for c = 1:n
    disp("x"+c + " = " + x(c));
end
end

