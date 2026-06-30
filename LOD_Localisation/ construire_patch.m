function [patch_nodes, kernel_patch_indices] = construire_patch(j, k, Xh, XH, IH)
    % ÉTAPE 4 & 5 : Identifie le domaine local (k-patch) autour du nœud grossier j
    NH = length(XH) - 1;
    Nf = length(Xh) - 1;
    ratio = Nf / NH;

    % Déterminer les bornes des éléments grossiers inclus dans le k-patch
    elem_gauche = max(1, j - k);
    elem_droite = min(NH, j + k - 1);

    % Convertir en indices de nœuds fins
    node_fine_debut = (elem_gauche - 1) * ratio + 1;
    node_fine_fin = elem_droite * ratio + 1;

    patch_nodes = node_fine_debut : node_fine_fin;

    % Trouver Wh(omega) = Nœuds fins du patch qui s'annulent aux nœuds grossiers
    % En 1D, ce sont les nœuds du patch sauf ceux qui sont multiples du ratio
    coarse_nodes_in_fine = 1 : ratio : (Nf + 1);
    kernel_patch_indices = setdiff(patch_nodes, coarse_nodes_in_fine);
end
