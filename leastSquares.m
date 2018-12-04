function [coeffs,plt] = leastSquares(textFile)
%% LEAST SQUARES METHOD FOR POLYNOMIAL FITTING %%
%Takes a text file for input, textFile
%textFile contains (in order) the number of points we are basing the fit
%off, N, the degree of the desired polynomial fit, k, and the points, where
%x precedes y (e.g. for (1,0) we have 1 then 0, then the next point in the
%same fashion, and so on)
%returns the coefficients (ascending order) for the desired fit and a plot

%open file
fileID = fopen(textFile);
%set specification to floating point numbers
fileSpec = '%f';
%put all file contents into row vector
all = fscanf(fileID, fileSpec, [1 inf]);
%extract N and k, delete from all
N = all(1);
k = all(2);
all(1) = [ ];
all(1) = [ ];
%get x and y into seperate vectors
x = zeros(1, N);
y = zeros(1, N);
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

%make matrix and b-vector to solve, size (k+1,k+1) and (1,k+1), respectively
A = zeros(k+1,k+1);
b = zeros(k+1,1);
%populate b-vector
for j = 1:k+1
    b(j) = sum((x.^(j-1)).*y);
end
%populate matrix, columnwise (e.g. do column 1, then 2,...,then k+1)
for n = 1:k+1
    power = n-1;
    for m = 1:k+1
        if n == 1 && m == 1
            A(m,n) = N;
            power = power+1;
        else
            A(m,n) = sum(x.^((power)));
            power = power+1;
        end
    end
end

coeffs = zeros(1, k+1);
A_det = det(A);
%% MAIN LOOP %%
for c = 1:k+1
    %use cramers rule to solve system for all a0,a1,...,ak
    Ac = A;
    Ac(:,c) = b;
    %add to coeffs, ascending order [a0,a1,...,ak]
    coeffs(c) = (det(Ac)/A_det);
end

%make plot to return, real values versus fit
pred = polyval(fliplr(coeffs),x);
plt = figure('Name','Least Squares Fit'); scatter(x,y); hold on; plot(x,pred); legend('Actual', 'Fit'); 

end

