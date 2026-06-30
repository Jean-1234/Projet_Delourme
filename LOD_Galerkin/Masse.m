function M = Masse(x, h)
% Masse : matrice de masse P1 sur la grille x.
% Entrees :
%   x : noeuds (n+2 x 1), inclut les bords 0 et 1
%   h : pas du maillage
% Sortie :
%   M : matrice de masse (n+2 x n+2), sparse

n = length(x) - 2;
M = sparse(n+2, n+2);

for k = 1:n+1
    M(k:k+1, k:k+1) = M(k:k+1, k:k+1) + h * [1/3  1/6 ; ...
                                                 1/6  1/3];
end

end
