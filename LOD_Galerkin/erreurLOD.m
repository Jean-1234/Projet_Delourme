function [err_l2, err_H1] = erreurLOD(U_lod, Uex, M, K)
% erreurLOD : calcul des erreurs L2 et H1 via matrices fines.
%
% Normes matricielles :
%   ||e||_L2^2  = e' * M * e
%   ||e||_H1^2  = e' * (K + M) * e     (norme H1 complete)
%
% Entrees :
%   U_lod : solution LOD sur grille fine (N_fin+2 x 1)
%   Uex   : solution exacte sur grille fine (N_fin+2 x 1)
%   M     : matrice de masse fine (N_fin+2 x N_fin+2)
%   K     : matrice de rigidite fine avec aeps (N_fin+2 x N_fin+2)
% Sorties :
%   err_l2 : erreur L2
%   err_H1 : erreur H1

e      = U_lod(:) - Uex(:);
err_l2 = sqrt(e' * M * e);
err_H1 = sqrt(e' * (K + M) * e);

end
