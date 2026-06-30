
%#  u_eps (solution exacte), u_0 , I_H u_eps
#   eps = 2^{-5},  H = 2^{-4}

clear; clc; close all


# PARAMETRES


ep  = 2^-5;
H   = 2^-4;
N_H = floor(1/H);

% maillage fin pour tracer u_eps et u_0 en continu
N_fin = 2000;
x_fin = linspace(0, 1, N_fin)';


# SOLUTION EXACTE

u_eps = exactSolution(x_fin, ep, 2)';
u_eps = u_eps(:);


# PARTIE MACROSCOPIQUE u_0
#A_0 = 1/2,  f = 1  =>  -A_0 u_0'' = 1  avec u_0(0)=u_0(1)=0
# u_0(x) = x - x^2   marche (EDO )

u_0 = x_fin - x_fin.^2;



[x_H, h_H] = mesh1D(N_H);         # noeuds grossiers x_H

# valeurs de u_eps aux noeuds grossiers
u_eps_nodes = exactSolution(x_H, ep, 2)';
u_eps_nodes = u_eps_nodes(:);

# interpoler lineairement I_H u_eps sur x_fin
IH_ueps = interp1(x_H, u_eps_nodes, x_fin, 'linear');


figure()

#gauche
subplot(1, 2, 1)
plot(x_fin, u_eps,   'b-',  'LineWidth', 1.5); hold on
plot(x_fin, u_0,     'r-',  'LineWidth', 1.5)
plot(x_fin, IH_ueps, 'g-',  'LineWidth', 1.5)
xlim([0 1]); ylim([-0.01 0.26])
xlabel(''); ylabel('')
legend('u_\epsilon', 'u_0', 'I_H u_\epsilon', ...
       'Location', 'northeast', 'FontSize', 10)
grid on; box on

# zoom autour de x = 0.5
subplot(1, 2, 2)
plot(x_fin, u_eps,   'b-',  'LineWidth', 1.5); hold on
plot(x_fin, u_0,     'r-',  'LineWidth', 1.5)
plot(x_fin, IH_ueps, 'g-',  'LineWidth', 1.5)
xlim([0.45 0.55])
ylim([0.2400 0.2605])
xlabel(''); ylabel('')
grid on; box on

title(sprintf('u_\\epsilon,  u_0  et  I_H u_\\epsilon     \\epsilon = 2^{-5},  H = 2^{-4}'), ...
        'FontSize', 11)
