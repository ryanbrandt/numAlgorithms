function x1 = newtonHorner(textFile)
%% NEWTON-RAPHSON METHOD WITH HORNER %%
%takes a text file for input, textFile
%textFile contains (in order) degree of p(x), n, degree-descending
%list of coefficients of p(x), initial guess, x0, error tolerance,                              
%epsilon and iteration tolerance, N

%open file
fileID = fopen(textFile);
%set specification to floating point numbers
fileSpec = '%f';
%get all contents of file into row vector
coeff = fscanf(fileID, fileSpec);
%extract degree n, delete from vector
n = coeff(1);
coeff(1) = [];
%extract x0, delete from vector
x0 = coeff(end-2);
coeff(end-2) = [];
%extract epsilon, delete from vector
eps = coeff(end-1);
coeff(end-1) = [];
%extract N, delete from vector
N = coeff(end);
coeff(end) = [];
%set initial error and number of iterations
error = inf;
count = 0;
%% MAIN LOOP %%
while(error > eps && count <= N)
    count = count + 1;
    %set p_x and p_prime_x, which correspond to p(x0) and p'(x0) for Horner
    %to first element in coeff, which is the highest degree coefficient
    p_x = coeff(1);
    p_prime_x = coeff(1);
    
    %Horners method to evaluate
    for i = 2:n
        p_x = p_x*x0 + coeff(i);
        p_prime_x = p_prime_x*x0 + p_x;
    end
    %final computation for p(x)
    p_x = p_x*x0 + coeff(end);
    
    %get newest approximation to root, x1, by newtons formula
    x1 = x0 - (p_x/p_prime_x);
    %update error
    error = abs(x1-x0);
    %set x0 to x1 for next iteration
    x0 = x1;
end
%make sure iteration tolerance not exceeded
if(count > N)
    error("Iteration tolerance exceeded");
end
%else, display answer
disp("Root Approximation: " + x1);
end




    