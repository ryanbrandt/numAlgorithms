function [coeffs, plt] = chebyshevLeastSquares(textFile)
%% LEAST SQUARES FIT WITH CHEBYSHEV POLYNOMIALS %%
%Takes on argument, a textfile, textFile
%textFile contains (in order) the number of points we are working with, N,
%the degree of the desired fit, n, and the points where x precedes y (e.g.
%for (1,0) 1 then 0)
%returns the coefficients a0,a1,..an+1 obtained for the fit (e.g. degree 3
%means we have 4 polynomials with coefficients ai) and a plot of the real
%data against our fit

%open file
fileID = fopen(textFile);
%set specification to floating point numbers
fileSpec = '%f';
%read entire file in array
all = fscanf(fileID, fileSpec, [1 inf]);
%extract N and n, delete from array
N = all(1);
n = all(2);
all(1) = [ ];
all(1) = [ ];
%get x and y into seperate column vectors
x = zeros(N, 1);
y = zeros(N, 1);
xIndexer = 0;
for i = 1:2*N
    %is y value if even index
    if rem(i,2) == 0
        y(i - (i/2)) = all(i);
    else
        x(i-xIndexer) = all(i);
        xIndexer = xIndexer+1;
    end
end
%known first two chebyshev polynomials
T0 = ones(N,1);
T1 = x/N-1;
chebyMat = vec2mat([T0 T1],2);
%recursively generate chebyshev polynomials up to n+1 (e.g. degree 3 means
%4 polynomials) and evaluate them at all x, add to matrix of polynomials
for i = 3:n+1
    T_cur = 2.*(x/N-1).*T1 - T0;
    chebyMat = vec2mat([chebyMat T_cur], i);
    T0 = T1;
    T1 = T_cur;
end
%make A matrix to solve for coefficients
A = zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
        A(i,j) = sum(chebyMat(:,i).*chebyMat(:,j));
    end
end
%get b-vector
b = zeros(n+1,1);
for i = 1:n+1
    b(i) = sum(y.*chebyMat(:,i));
end
%solve for coefficients a0,a1...,an+1
coeffs = A\b;
%get actual Pn(x) evaluation
P_n_x = 0;
for i = 1:n+1
    P_n_x = P_n_x + coeffs(i)*chebyMat(:,i);
end
%make plot
plt = figure('Name','Chebyshev Least Squares'); scatter(x,y); hold on; plot(x,P_n_x); legend("Actual", "Fit");

end

