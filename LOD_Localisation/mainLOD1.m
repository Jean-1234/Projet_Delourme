
tic()
clear all; close all; clc;

disp('Choisir un cas :')
disp('1 : a_eps = 1  (cas constant, pas d''oscillations)')
disp('2 : a_eps = 1/(2 - cos(2*pi*x/eps))  (cas oscillant périodique)')
choice = input('Votre choix = ');


Nh = 1024;
k_patch = 0;

% Tailles de grilles grossières pour l'étude de convergence
NH_vec = [4, 8, 16, 32, 64];
H_val = 1 ./ NH_vec;

if choice == 1
    epsilon_vec = 0.02; % non significatif ici
    a_eps_base = @(x, ep) ones(size(x));
    func = @(x) ones(size(x));

elseif choice == 2
    epsilon_vec = [2^-2, 2^-4, 2^-6,2^-8];
    a_eps_base = @(x, ep) 1 ./ (2 - cos(2*pi*x ./ ep));
    func = @(x) ones(size(x));

else
    error('Choix invalide');
end

err_L2 = zeros(length(epsilon_vec), length(NH_vec));
err_H1 = zeros(length(epsilon_vec), length(NH_vec));

ordreL2 = zeros(length(epsilon_vec),1);
ordreH1 = zeros(length(epsilon_vec),1);

for k = 1:length(epsilon_vec)

    ep = epsilon_vec(k);
    a_eps = @(x) a_eps_base(x, ep);

    fprintf('\n--- Calcul de la référence fine pour epsilon = %.4f ---\n', ep);

    % Solution de référence fine
    [Xh, ~] = construire_maillages(Nh, 4);

    [Kh, Mh] = assembler_systeme_fin(Xh, a_eps);

    F_fine = Mh * func(Xh)';

    u_ref = zeros(length(Xh),1);
    u_ref(2:end-1) = ...
        Kh(2:end-1,2:end-1) \ F_fine(2:end-1);

    for i = 1:length(NH_vec)

        NH = NH_vec(i);

        fprintf('  LOD Grossière N = %d...\n', NH);

        [~, XH] = construire_maillages(Nh, NH);

        Phi_ms = phase_offline(Xh, XH, Kh, k_patch);

        [u_H, u_LOD] = ...
            phase_online(Xh, XH, Phi_ms, Kh, Mh, func);

        error_vector = u_ref - u_LOD;

        % Erreur L2
        err_L2(k,i) = ...
            sqrt(error_vector' * Mh * error_vector);

        % Erreur H1
        err_H1(k,i) = ...
            sqrt(error_vector' * (Kh + Mh) * error_vector);

    end

    % Ordres de convergence
    p = polyfit(log(H_val), log(err_L2(k,:)),1);
    q = polyfit(log(H_val), log(err_H1(k,:)),1);

    ordreL2(k) = p(1);
    ordreH1(k) = q(1);

    fprintf('  => Ordre L2 = %.3f\n', ordreL2(k));
    fprintf('  => Ordre H1 = %.3f\n', ordreH1(k));

end


%% Création des légendes

legL2 = cell(1,length(epsilon_vec));
legH1 = cell(1,length(epsilon_vec));

for k = 1:length(epsilon_vec)

    legL2{k} = sprintf('\\epsilon = %.3f (p = %.2f)', ...
                        epsilon_vec(k), ordreL2(k));

    legH1{k} = sprintf('\\epsilon = %.3f (p = %.2f)', ...
                        epsilon_vec(k), ordreH1(k));

end


%% Affichage des graphes


figure('Name','Pentes de convergence LOD','Position',[100 100 1200 400]);
colors = {'r-o','b-s','g-^','m-d'};


%% Graphe L2

subplot(1,2,1)
hold on
grid on

for k = 1:length(epsilon_vec)

    loglog(H_val,...
           err_L2(k,:)*10,...
           colors{k},...
           'LineWidth',2,...
           'MarkerSize',7);

end

% Courbe théorique O(H²)
loglog(H_val,...
       H_val*10,...
       'c-d',...
       'LineWidth',2);
loglog(H_val,...
       H_val.^2*10,...
       'k-d',...
       'LineWidth',2);

xlabel('H ')
ylabel('Erreur relative L2')

title('Convergence Norme L2')

legend([legL2,{'O(H)'}, {'O(H^2)'}],...
       'Location','southeast', ...
       'FontSize',6)



%% Graphe H1

subplot(1,2,2)
hold on
grid on

for k = 1:length(epsilon_vec)

    loglog(H_val,...
           err_H1(k,:)*10,...
           colors{k},...
           'LineWidth',2,...
           'MarkerSize',7);

end

% Courbe théorique O(H)
loglog(H_val,...
       H_val,...
       'k-v',...
       'LineWidth',2);
loglog(H_val,...
       H_val.^2*10,...
       'c-d',...
       'LineWidth',2);
xlabel('H')
ylabel('Erreur relative H1')

title('Convergence Norme H1')

legend([legH1, {'O(H)'},{'O(H^2)'}], ...
       'Location','southeast', ...
       'FontSize',6);


temps=toc()
