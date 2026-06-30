function [patch_nodes, kernel_patch_indices] = construire_patch(j, k, Xh, XH, IH)

    x_left  = XH(max(1, j-k));
    x_right = XH(min(length(XH), j+k));

    % patch sur la grille fine
    patch_nodes = find(Xh >= x_left & Xh <= x_right);

    % nœuds grossiers
    coarse_nodes = find(ismember(Xh, XH));

    % espace Wh(patch)
    kernel_patch_indices = setdiff(patch_nodes, coarse_nodes);

end
