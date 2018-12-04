function [p_x, p_prime_x] = hornersMethod(textFile)
%% HORNERS METHOD FOR POLYNOMIAL EVALUATION %%
%takes a text file for input, textFile
%textFile contains the degree  of our p(x), n, a degree ascending 
%list of the coefficients of p(x), and a real number, x0,
%to evaluate p(x0) and p'(x0) efficiently in O(n) time

%open file 
fileID = fopen(textFile);
%set format specification to floating point numbers
formatSpec = '%f';
%get all contents of file into row vector
coeff = fscanf(fileID, formatSpec, [1 inf]);
%last element has x0, extract and delete
x0 = coeff(end);
coeff(end) = [];
%first element is degree n, extract and delete
n = coeff(1);
coeff(1) = []; 
%set p_x and p_prime_x to first element, highest degree coefficient
p_x = coeff(1);
p_prime_x = coeff(1);
%% MAIN LOOP %%
for i = 2:n
    p_x = p_x*x0 + coeff(i);
    p_prime_x = p_prime_x*x0 + p_x;
end
%final computation for p(x0)
p_x = p_x*x0 + coeff(end);
%display results
disp("P(x0) = " + p_x)
disp("P'(x0) = " + p_prime_x)
end

