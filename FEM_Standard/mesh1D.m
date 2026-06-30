function [x,h] = mesh1D(n)
h = 1/(n+1);
x = linspace(0,1,n+2);
end
