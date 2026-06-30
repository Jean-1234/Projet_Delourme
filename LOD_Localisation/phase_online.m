function [u_H, u_LOD] = phase_online(Xh, XH, Phi_LOD, Kh, Mh, f)

NH = length(XH) - 1;
Nh = length(Xh) - 1;

A_LOD = Phi_LOD' * Kh * Phi_LOD;

f_fine = f(Xh)';

% Correction de l'intégration : on utilise la matrice de masse fine Mh
F_LOD = Phi_LOD' * (Mh * f_fine);

idx = 2:NH;

u_H = zeros(NH+1,1);
u_H(idx) = A_LOD(idx,idx) \ F_LOD(idx);

u_LOD = Phi_LOD * u_H;

end
