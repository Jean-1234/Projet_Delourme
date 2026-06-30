function [Phi_LOD] = phase_offline(Xh, XH, Kh, k)

NH = length(XH) - 1;
Nh = length(Xh) - 1;

Phi_LOD = zeros(Nh+1, NH+1);

for j = 1:NH+1


    [~, kernel_patch] = construire_patch(j, k, Xh, XH);


    K_omega = Kh(kernel_patch, kernel_patch);


    Phi_j_std = zeros(Nh+1, 1);


    x_center = XH(j);

    if j > 1
        x_left = XH(j-1);

        idx_gauche = find(Xh >= x_left & Xh <= x_center);
        Phi_j_std(idx_gauche) = (Xh(idx_gauche) - x_left) / (x_center - x_left);
    end

    if j < NH+1
        x_right = XH(j+1);

        idx_droite = find(Xh >= x_center & Xh <= x_right);
        Phi_j_std(idx_droite) = (x_right - Xh(idx_droite)) / (x_right - x_center);
    end


    b_global = -Kh * Phi_j_std;
    bj = b_global(kernel_patch);


    if ~isempty(kernel_patch)
        cj_local = K_omega \ bj;


        Qj = zeros(Nh+1, 1);
        Qj(kernel_patch) = cj_local;
    else
        Qj = zeros(Nh+1, 1);
    end


    Phi_LOD(:,j) = Phi_j_std + Qj;

end

end
