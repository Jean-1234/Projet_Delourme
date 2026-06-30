function [A,F] = assembleFEM(x,h,epsi,f,aeps)
n = length(x) - 2;
A = sparse(n,n);
I0 = integral(aeps, x(1), x(2));
I1 = integral(aeps, x(2), x(3));
A(1,1) = (I0 + I1)/h^2;
A(1,2) = -I1/h^2;
for i = 2:n-1
    xim1 = x(i);
    xi   = x(i+1);
    xip1 = x(i+2);
    I1 = integral(aeps, xim1, xi);
    I2 = integral(aeps, xi, xip1);
    A(i,i)   = (I1 + I2)/h^2;
    A(i,i-1) = -I1/h^2;
    A(i,i+1) = -I2/h^2;
    F(i)     = h * f(xi);
end
In_1 = integral(aeps, x(n), x(n+1));
In   = integral(aeps, x(n+1), x(n+2));
A(n,n)   = (In_1 + In)/h^2;
A(n,n-1) = -In_1/h^2;
#second membre
F = h * f(x(2:end-1))' ;
end
