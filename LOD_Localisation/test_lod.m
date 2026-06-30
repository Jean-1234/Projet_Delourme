clear; close all; clc ;
fprintf('\n--- TEST phiLOD ---\n')

NH  = 4;
Nh=199;
[xH, h_grossier]   = mesh1D(NH);
[xh, h_fin]  = mesh1D(Nh);
aeps_1   = @(t) ones(size(t));


j = 2;
[phi, varphi] = phiLOD(j, xH, xh, aeps_1);

j = 2;
[phi1, var1] = phiLOD(j, xH, xh, aeps_1);

plot(xh, var1, 'r-d', 'LineWidth', 1.5);
 hold on
plot(xh, phi1, 'b-v',  'LineWidth', 1.5)
xlabel('x');
grid on



