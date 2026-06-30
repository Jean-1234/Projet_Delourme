function Uex = exactSolution(x, epsi, choice)

switch choice

    case 1
        Uex = sin(pi*x);

    case 2
        Uex = x - x.^2 ...
             - epsi * ( ...
             (1/(4*pi)) * sin(2*pi*x/epsi) ...
             - (x/(2*pi)) .* sin(2*pi*x/epsi) ...
             - (epsi/(4*pi^2)) * cos(2*pi*x/epsi) ...
             + (epsi/(4*pi^2)) );



    otherwise
        error('Choix invalide');
end

end
