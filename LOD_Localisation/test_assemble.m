clear; close all; clc;
fprintf('\n--- TEST assembleLOD ---\n');

% maillages
NH = 6;
Nh = 200;

[xH, ~] = mesh1D(NH);   % maillage grossier LOD
[xh, ~] = mesh1D(Nh);   % maillage fin

% coefficient homogène
aeps = @(t) ones(size(t));

% second membre simple (test standard)
f = @(x) sin(pi*x);

% assemblage
[A_lod, M_lod, K_lod, F_lod] = assembleLOD(xH, xh, aeps, f);

% affichage structure
disp('A_lod ='); full(A_lod)
disp('K_lod ='); full(K_lod)
disp('M_lod ='); full(M_lod)
disp('F_lod ='); F_lod
