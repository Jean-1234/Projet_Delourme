function [Xh, XH] = construire_maillages(N_fine, N_coarse)
    % Génère les coordonnées des nœuds
    Xh = linspace(0, 1, N_fine + 1);
    XH = linspace(0, 1, N_coarse + 1);
end
