function U_fin = reconstructLOD(U_lod, x, x_fin, aeps)
% reconstructLOD : reconstruit la solution LOD sur grille fine.
%
% u^ms = sum_j alpha_j * Lambda_j
%
% U_lod contient les valeurs nodales grossieres alpha_j = U_lod(j+1)
% (indices 2 a n+1, les bords 0 et 1 sont nuls).
%
% On reconstruit :
%   U_fin(k) = sum_{j=1}^{n} alpha_j * Lambda_j(x_fin(k))
%
% Entrees :
%   U_lod : solution LOD grossiere (n+2 x 1), avec CL
%   x     : noeuds grossiers (n+2 x 1)
%   x_fin : grille fine (N_fin x 1)
%   aeps  : handle, coefficient a_epsilon
% Sortie :
%   U_fin : solution reconstruite sur grille fine (N_fin x 1)

n     = length(x) - 2;
alpha = U_lod(2:end-1);      % coefficients interieurs

U_fin = zeros(length(x_fin), 1);

for j = 1:n
    [phi_j, ~] = phiLOD(j, x, x_fin, aeps);
    U_fin = U_fin + alpha(j) * phi_j(:);
end

end
