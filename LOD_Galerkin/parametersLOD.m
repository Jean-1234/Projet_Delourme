function [aeps, f, epsi] = parametersLOD(choice)

switch choice

    case 1
        epsi = 0.05;
        aeps = @(x) ones(size(x));
        f    = @(x) pi^2 .* sin(pi*x);

    case 2
        epsi = [2^-2, 2^-4, 2^-6, 2^-8];
        aeps = @(x, ep) 1 ./ (2 - cos(2*pi*x ./ ep));
        f    = @(x) ones(size(x));

    case 3
    epsi = [2^-2, 2^-4, 2^-6, 2^-8];

    % coefficient discontinu
    a    = @(y) 1 .*(y < 0.5) + 2 .*(y >= 0.5);
    aeps = @(x, ep) a(mod(x./ep, 1));

    % u_ex = cos(2*pi*x) - 1,  CL : u(0)=u(1)=0
    % f = -(a_eps * u_ex')' = 4*pi^2 * a_eps * cos(2*pi*x)  p.p.
    f = @(x,ep) 4*pi^2 * a(mod(x./ep,1)) .* cos(2*pi*x);

    otherwise
        error('Choix invalide : 1 ou 2');

end
end
