function [aeps,f,epsi] = parameters(choice)

switch choice

    case 1

        epsi = 0.05;
        aeps = @(x) ones(size(x));
        f = @(x) pi^2 * sin(pi*x);

    case 2

        epsi = [2^-2, 2^-4, 2^-6, 2^-8];
        aeps = @(x,ep) 1 ./ (2 - cos(2*pi*x/ep));
        f = @(x) ones(size(x));

    case 3

        epsi = [2^-2, 2^-4, 2^-6, 2^-8];

        % fonction de base sur [0,1]
        a = @(y) 1 .* (y <= 0.5) + 2 .* (y > 0.5);

        % version oscillante
        aeps = @(x,ep) a(mod(x./ep,1));

        f = @(x) ones(size(x));

    otherwise
        error('Choix invalide');

end

end
