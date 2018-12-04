function approx_intgrl = simpsonsRule(textFile)
%% SIMPSONS RULE FOR NUMERICAL INTEGRATION %%
%Takes one argument, a txt file, textFile
%textFile contains (in order) the function f(x) to be integrated, the
%bounds of integration, a and b (a precedes b) and the number of
%subintervals we are applying to the composite rule, n

%open file
fileID = fopen(textFile);
%specify as string, since we have f(x)
all = textscan(fileID, '%s', 'Delimiter', ' ');
f_x = all{1}{1};
%get variable working with, only one variable
var = string(symvar(all{1}{1}));
%make into feval evaluatable function (e.g. add @(<var>) to front)
var = strcat('@(', var);
var = strcat(var, ')');
f_x = strcat(var, f_x);
f_x = str2func(f_x);
%to num
a = str2num(all{1}{2});
b = str2num(all{1}{3});
n = str2num(all{1}{4});

%get non iterative comptations, h, f(x0) and f(xn) (add these to integral)
h = (b-a)/n;
approx_intgrl = feval(f_x, a) + feval(f_x,b);
%to check if odd/even sum index
parity = 1;
%% MAIN LOOP %%
for i = a+h:h:b-h
    %even sum indexes by 2
    if(rem(parity,2) == 0)
       approx_intgrl = approx_intgrl + 2*feval(f_x,i);
    %odd sum indexes by 4    
    else
        approx_intgrl = approx_intgrl + 4*feval(f_x,i);
    end
    parity = parity+1;
end
%lastly, divide by h/3
approx_intgrl = (h/3)*(approx_intgrl);

end

