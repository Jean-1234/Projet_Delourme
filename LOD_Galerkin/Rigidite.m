function K = Rigidite(x, h, aeps)
% Rigidite : matrice de rigidite P1 sur la grille x,
%            avec coefficient variable aeps.
%
% K(i,i)   = (1/h^2) * [ int_{T_{i-1}} aeps  +  int_{T_i} aeps ]
% K(i,i+1) = -(1/h^2) * int_{T_i} aeps
%
% Entrees :
%   x    : noeuds (n+2 x 1), inclut les bords 0 et 1
%   h    : pas du maillage
%   aeps : handle de fonction, coefficient a_epsilon
% Sortie :
%   K : matrice de rigidite (n+2 x n+2), sparse

n = length(x) - 2;
K = sparse(n+2, n+2);

for k = 1:n+1
    Ik = integral(aeps, x(k), x(k+1));
    K(k:k+1, k:k+1) = K(k:k+1, k:k+1) + (Ik/h^2) * [ 1  -1 ; ...
                                                         -1   1];
end

end
