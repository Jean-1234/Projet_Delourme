function [phi_j, varphi_j] = phiLOD(j, x, x_fin, aeps)
% phiLOD : fonction de base multiscale Lambda_j = (1-C) lambda_j
%          et fonction chapeau classique varphi_j = lambda_j.
%
% En 1D, la formule exacte (eq. 4.20, Verfurth p.73) donne :
%
%   Lambda_j(x) = int_{x_{j-1}}^x  a^{-1}  /  int_{x_{j-1}}^{x_j} a^{-1}
%                                                   si x dans T_{j-1}
%
%   Lambda_j(x) = 1 - int_{x_j}^x a^{-1} / int_{x_j}^{x_{j+1}} a^{-1}
%                                                   si x dans T_j
%
% Entrees :
%   j      : indice du noeud interieur (1 <= j <= n, noeud = x(j+1))
%   x      : noeuds grossiers (n+2 x 1), inclut 0 et 1
%   x_fin  : points d evaluation fins
%   aeps   : handle, coefficient a_epsilon
% Sorties :
%   phi_j    : Lambda_j evaluee sur x_fin  (base multiscale)
%   varphi_j : lambda_j evaluee sur x_fin  (base chapeau classique)

H              = x(2) - x(1);
noeud_central  = x(j+1);
borne_gauche   = x(j);
borne_droite   = x(j+2);

ainv = @(t) 1 ./ aeps(t);

% ---- Base chapeau classique varphi_j ----

varphi_j = zeros(size(x_fin));

idx_G = (x_fin >= borne_gauche) & (x_fin <= noeud_central);
idx_D = (x_fin >= noeud_central) & (x_fin <= borne_droite);

varphi_j(idx_G) = (x_fin(idx_G) - borne_gauche) / H;
varphi_j(idx_D) = (borne_droite - x_fin(idx_D)) / H;

% ---- Base multiscale Lambda_j (formule exacte 1D) ----

phi_j = zeros(size(x_fin));

Int_G = integral(ainv, borne_gauche, noeud_central);
Int_D = integral(ainv, noeud_central, borne_droite);

for k = 1:length(x_fin)

    t = x_fin(k);

    if t >= borne_gauche && t <= noeud_central

        phi_j(k) = integral(ainv, borne_gauche, t) / Int_G;

    elseif t >= noeud_central && t <= borne_droite

        phi_j(k) = 1 - integral(ainv, noeud_central, t) / Int_D;

    end

end

end
