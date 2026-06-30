function plotComparaison(ep, choice)
% plotComparaison : compare FEM standard (fin et grossier), LOD et solution exacte.
%
% Reproduit l effet de la Figure 1.2 et 1.4 de Verfurth pour H = 2^{-3}.
%
% Entrees :
%   ep     : epsilon (defaut 2^{-5})
%   choice : cas test (defaut 2)
clear; close all;clc;
if nargin < 1, ep     = 2^-5; end
if nargin < 2, choice = 2;    end

% ---- Parametres ----

N_H   = 2^3 - 1;        % noeuds interieurs grossiers : H = 2^{-3}
N_fin = 1999;            % grille fine

aeps = @(t) 1 ./ (2 - cos(2*pi*t ./ ep));
f    = @(t) ones(size(t));

% ---- Grilles ----

[x_H, H]          = mesh1D(N_H);
[x_fin_bords, h_fin] = mesh1D(N_fin);
x_fin = x_fin_bords(:);

% ---- Solution exacte ----

Uex = exactSolution(x_fin_bords, ep, choice)';
Uex = Uex(:);

% ---- FEM grossier (sur grille H) ----

[A_fem_H, F_fem_H] = assembleFEM(x_H, H, ep, f, aeps);
U_fem_H            = solveFEM(A_fem_H, F_fem_H);

% interpolation sur grille fine pour affichage
U_fem_H_fin = interp1(x_H, U_fem_H, x_fin, 'linear');

% ---- FEM fin (sur grille fine h_fin) ----

[A_fem_h, F_fem_h] = assembleFEM(x_fin_bords, h_fin, ep, f, aeps);
U_fem_h            = solveFEM(A_fem_h, F_fem_h);

% ---- LOD Galerkin ----

[A_lod, F_lod] = assembleLOD(x_H, H, x_fin, aeps, f);
U_lod          = solveLOD(A_lod, F_lod);
U_lod_fin      = reconstructLOD(U_lod, x_H, x_fin, aeps);

% ---- Figure vue globale ----

figure('Name', 'Comparaison FEM / LOD / Exacte', ...
       'Position', [100 100 1000 420])

subplot(1, 2, 1)
plot(x_fin, Uex,          'k-',  'LineWidth', 2.0); hold on
plot(x_fin, U_fem_h(:),   'b--', 'LineWidth', 1.5)
plot(x_fin, U_fem_H_fin,  'r-.',  'LineWidth', 1.5)
plot(x_fin, U_lod_fin,    'g-',  'LineWidth', 1.5)
legend('u exacte', ...
       sprintf('FEM fin  h=2^{-%d}', round(-log2(h_fin))), ...
       sprintf('FEM grossier  H=2^{-3}'), ...
       'LOD Galerkin  H=2^{-3}', ...
       'Location', 'northeast')
xlabel('x'); grid on; box on
title(sprintf('Vue globale   ', round(log2(ep))))

% ---- Figure zoom ----

subplot(1, 2, 2)
plot(x_fin, Uex,          'k-',  'LineWidth', 2.0); hold on
plot(x_fin, U_fem_h(:),   'b--', 'LineWidth', 1.5)
plot(x_fin, U_fem_H_fin,  'r-.',  'LineWidth', 1.5)
plot(x_fin, U_lod_fin,    'g-',  'LineWidth', 1.5)
xlim([0.47 0.52])
xlabel('x'); grid on; box on
title('Zoom')

% ---- Affichage erreurs ----

M  = Masse(x_fin_bords, h_fin);
K  = Rigidite(x_fin_bords, h_fin, aeps);

e_fem_h   = U_fem_h(:)   - Uex;
e_fem_H   = U_fem_H_fin(:) - Uex;
e_lod     = U_lod_fin(:)  - Uex;

fprintf('\n--- H = 2^{-3},  epsilon = %.2e ---\n', ep)
fprintf('FEM fin   : L2 = %.3e   H1 = %.3e\n', ...
        sqrt(e_fem_h' * M * e_fem_h), ...
        sqrt(e_fem_h' * (K+M) * e_fem_h))
fprintf('FEM gros. : L2 = %.3e   H1 = %.3e\n', ...
        sqrt(e_fem_H' * M * e_fem_H), ...
        sqrt(e_fem_H' * (K+M) * e_fem_H))
fprintf('LOD       : L2 = %.3e   H1 = %.3e\n', ...
        sqrt(e_lod'   * M * e_lod), ...
        sqrt(e_lod'   * (K+M) * e_lod))

end
