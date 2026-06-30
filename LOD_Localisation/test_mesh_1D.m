clear; close all; clc;

NH = 4;
Nh = 10;

fprintf('Maillage grossier : %d noeuds (NH+2 = %d)\n', NH+2, NH+2);
[xH, h_grossier] = mesh1D(NH)

fprintf('Maillage fin : %d noeuds (Nh+2 = %d)\n', Nh+2, Nh+2);
[xh, h_fin] = mesh1D(Nh)
