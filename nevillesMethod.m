function p_x_0 = nevillesMethod(textFile)
%% NEVILLE'S METHOD FOR EVALUATING INTERPOLATED POLYNOMIALS %%
%Takes a text file for input, textFile
%textFile contains (in order), the degree n of Pn(x) (n = # points - 1),
%the points, where x precedes y (e.g. for (1,0) input = 1 then 0, and so on
%for all points) and x0, for which we are evaluating Pn(x0)
%returns the value of Pn(x0)

%open file
fileID = fopen(textFile);
%set specification to floating point numbers
fileSpec = '%f';
%read in all contents into row vector
all = fscanf(fileID, fileSpec, [1 inf]);
%get n, x0, delete from all
n = all(1);
x0 = all(end);
all(1) = [ ];
all(end) = [ ];

%get x, get y into seperate vectors
y = zeros(1,n+1);
x = zeros(1,n+1);
xIndexer = 0;
for i = 1:length(all)
    if rem(i,2) == 0
        y(i - (i/2)) = all(i);
    else
        x(i-xIndexer) = all(i);
        xIndexer = xIndexer+1;
    end
end
%Make Neville's table/matrix
table = diag(y);
%% MAIN LOOP %%
for d = 2:n+1
    for i = 1:(n+1)-(d-1)
        %get column index, do nevilles formula 
        j = i + (d - 1);
        table(i,j) = (((x0 - x(i))*table(i+1,j)) - ((x0-x(j))*table(i,j-1)))/(x(j)-x(i));
    end
end
p_x_0 = table(1,n+1);
disp("P" + n + "(" + x0 + ") = "  + p_x_0);
end

