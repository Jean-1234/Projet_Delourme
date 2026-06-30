function [A,F] = assembleFEMv2(x,h,epsi,f,aeps,choice)

n = length(x) - 2;

switch choice

  case 1

    I = [1:n , 2:n , 1:n-1];
    J = [1:n , 1:n-1 , 2:n];

    diag_main = 2*ones(n,1);
    sub_diag  = -ones(n-1,1);

    V = [diag_main ;sub_diag ;sub_diag];

    A = (1/h^2) * sparse(I,J,V,n,n);

  case 2

    % primitive exacte
    Fint = @(x) (epsi/(pi*sqrt(3))) * ...
        atan( sqrt(3) * tan( pi*x/epsi ) );

    val_element = zeros(n+1,1);

    for i = 1:n+1
        val_element(i) = (Fint(x(i+1)) - Fint(x(i))) / (h*h);
    end

    diag_main = val_element(1:n) + val_element(2:n+1);
     sub_diag  = -val_element(2:n);
    I = [1:n , 2:n , 1:n-1];
    J = [1:n , 1:n-1 , 2:n];
    V = [diag_main ; sub_diag ; sub_diag];

    A = sparse(I,J,V,n,n);

  otherwise
    error('choix invalide');

end

% second membre
F = zeros(n,1);
for i = 1:n
    F(i)=(h/2)*(f(x(i))+f(x(i+1))) ;
end

end
