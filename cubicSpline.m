function [coeffs, plt] = cubicSpline(textFile)
%% CUBIC SPLINE METHOD FOR POLYNOMIAL INTERPOLATION %%
%Takes one argument, a text file, textFile
%textFile contains (in order) the number of points (x,y), n, and the points
%where x precedes y (e.g. for (1,0) 1 then 0, and so on)
%returns a coefficient matrix for all s(x) and a plot to show our
%inteprolation against our actual points

%open file
fileID = fopen(textFile);
%set speciifcation to floating point numbers
fileSpec = '%f';
%read file into row vector
all = fscanf(fileID, fileSpec, [1 inf]);
%get n, delete
n = all(1);
all(1) = [ ];
%get x and y into seperate row vectors
x = zeros(1, n);
y = zeros(1, n);
xIndexer = 0;
for k = 1:2*n
    if rem(k,2) == 0
        y(k-(k/2)) = all(k);
    else
        x(k-xIndexer) = all(k);
        xIndexer = xIndexer + 1;
    end
end
%initialize space
h = zeros(n-1,1); b = zeros(n-1,1); d = zeros(n-1,1); a = y;
%get hi
for i = 1:n-1
        h(i) = x(i+1)-x(i);
end
%make tridiagonal matrix to solve for ci
A = zeros(n,n);
A(1,1) = 1; A(n,n) = 1; A(2,1) = h(1); A(n-1,n) = h(end);
j = 2;
%main diagonal
for i = 2:n-1
    A(i,j) = 2*(h(i-1)+h(i));
    j = j+1;
end
%upper and lower diagonal
j = 3;
for i = 2:n-2
    A(i,j) = h(i-1);
    j = j+1;
end
j = 2;
for i = 3:n-1
    A(i,j) = h(j-1);
    j = j+1;
end
%get solution vector (e.g. Ax = solution)
sol = zeros(n,1);
sol(1) = 0; sol(n) = 0;
for i = 2:n-1
    sol(i) = (3/h(i))*(a(i+1)-a(i)) - (3/h(i-1))*(a(i)-a(i-1));
end
%solve for ci
c = A\sol;
%get bi by substituting ci
for i = 1:n-1
    b(i) = ((a(i+1)-a(i))/h(i)) - (((2*c(i)+c(i+1))*h(i))/3);
end
%get di by substituting ci
for i = 1:n-1
    d(i) = (1/3*h(i))*(c(i+1)-c(i));
end
%get s(x) matrix for output
coeffs = zeros(n-1,4);
for i = 1:n-1
    coeffs(i,1) = a(i);
    coeffs(i,2) = b(i);
    coeffs(i,3) = c(i);
    coeffs(i,4) = d(i);
end
%get a plot for output, need to plot each s(x) in its subdomain
%need to have a temp placeholderfor x since fplot evaluates x
vals = x;
plt = figure('Name', 'Cubic Spline Fit'); scatter(x,y); hold on;
for i = 1:n-1
    subDomain = [x(i) x(i+1)];
    %add si(x) in its subdomain to plot
    plt  = fplot(@(x) coeffs(i,1) + coeffs(i,2).*(x-vals(i)) + coeffs(i,3).*(x-vals(i)).^2 + coeffs(i,4).*(x-vals(i)).^3, subDomain);
    hold on;
end

