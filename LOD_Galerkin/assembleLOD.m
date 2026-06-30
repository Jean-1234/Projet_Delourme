function [A_lod, F_lod] = assembleLOD(x, h, x_fin, aeps, f)
% assembleLOD : matrice de rigidite et second membre LOD Galerkin.
%
% LOD GALERKIN : les deux termes utilisent les bases multiscales Lambda_j.
%
%   A_lod(i,j) = b(Lambda_i, Lambda_j) = int_Omega a * Lambda_i' * Lambda_j' dx
%
%   F_lod(j)   = int_Omega f * Lambda_j dx
%                                       <- CRUCIAL : Lambda_j, pas lambda_j
%
% En 1D, la formule exacte de Lambda_j (eq. 4.20, Verfurth p.73) donne
% pour la matrice (p.74, Verfurth) :
%
%   A_lod(i,i)   = 1/Int_i   +  1/Int_{i+1}
%   A_lod(i,i-1) = -1/Int_i
%   A_lod(i,i+1) = -1/Int_{i+1}
%
% avec Int_j = int_{T_j} a_eps^{-1} dx  (moyenne harmonique sur T_j).
%
% Pour le second membre, on integre f*Lambda_j numeriquement sur x_fin.
%
% Entrees :
%   x     : noeuds grossiers (n+2 x 1)
%   h     : pas grossier
%   x_fin : grille fine d evaluation (N_fin x 1)
%   aeps  : handle, coefficient a_epsilon
%   f     : handle, second membre
% Sorties :
%   A_lod : matrice de rigidite LOD (n x n), sparse
%   F_lod : second membre LOD (n x 1)

n    = length(x) - 2;
ainv = @(t) 1 ./ aeps(t);

% ---- Integrales harmoniques par element T_j = [x(j), x(j+1)] ----

Int = zeros(n+1, 1);
for j = 1:n+1
    Int(j) = integral(ainv, x(j), x(j+1));
end

% ---- Matrice de rigidite LOD ----

A_lod = sparse(n, n);

A_lod(1, 1) = 1/Int(1) + 1/Int(2);
A_lod(1, 2) = -1/Int(2);

for i = 2:n-1
    A_lod(i, i-1) = -1/Int(i);
    A_lod(i, i)   =  1/Int(i) + 1/Int(i+1);
    A_lod(i, i+1) = -1/Int(i+1);
end

A_lod(n, n-1) = -1/Int(n);
A_lod(n, n)   =  1/Int(n) + 1/Int(n+1);

% ---- Second membre LOD : F(j) = int f * Lambda_j dx ----
% On evalue Lambda_j sur x_fin et on integre par la regle des trapezes.

h_fin = x_fin(2) - x_fin(1);
f_fin = f(x_fin);

F_lod = zeros(n, 1);

for j = 1:n

    [phi_j, ~] = phiLOD(j, x, x_fin, aeps);

    F_lod(j) = h_fin * ( 0.5*f_fin(1)*phi_j(1) ...
                       + sum(f_fin(2:end-1) .* phi_j(2:end-1)) ...
                       + 0.5*f_fin(end)*phi_j(end) );
end

end
