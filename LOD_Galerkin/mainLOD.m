clear; clc; close all
tic()
disp('Choisir un cas :')
disp('1 : A_eps = 1  ')
disp('2 : A_eps = 1/(2-cos(2*pi*x/eps))  ')
disp('2 : A_eps discontinu')
choice = input('Votre choix = ');

[aeps, f, epsi] = parametersLOD(choice);

% Grille grossiere : N_val noeuds interieurs
% Grille fine fixe pour evaluation et erreurs

N_val  = [7, 15, 31, 63, 127];
N_fin  = 200;                   % noeuds interieurs de la grille fine
[x_fin_bords, h_fin] = mesh1D(N_fin);
x_fin = x_fin_bords(:);         % inclut 0 et 1


% CAS 1

if choice == 1

    err_l2 = zeros(size(N_val));
    err_H1 = zeros(size(N_val));
    H_val  = zeros(size(N_val));

    for i = 1:length(N_val)

        n       = N_val(i);
        [x, h]  = mesh1D(n);
        H_val(i) = h;

        [A_lod, F_lod] = assembleLOD(x, h, x_fin, aeps, f);
        U_lod          = solveLOD(A_lod, F_lod);

        % Reconstruction sur grille fine
        U_fin = reconstructLOD(U_lod, x, x_fin, aeps);

        % Matrices fines pour les erreurs
        M = Masse(x_fin_bords, h_fin);
        K = Rigidite(x_fin_bords, h_fin, aeps);
        Uex = exactSolution(x_fin_bords, epsi, 1)';

        [err_l2(i), err_H1(i)] = erreurLOD(U_fin, Uex(:), M, K);

    end

    % Solutions
    figure('Name', 'Solutions LOD Galerkin - cas 1')
    n     = N_val(end);
    [x,h] = mesh1D(n);
    [A_lod, F_lod] = assembleLOD(x, h, x_fin, aeps, f);
    U_lod = solveLOD(A_lod, F_lod);
    U_fin = reconstructLOD(U_lod, x, x_fin, aeps);
    Uex   = exactSolution(x_fin_bords, epsi, 1)';
    plot(x_fin, U_fin, 'b-', 'LineWidth', 1.5); hold on
    plot(x_fin, Uex,   'k-', 'LineWidth', 1.5)
    legend('LOD Galerkin', 'Exacte'); grid on; xlabel('x')
    title('LOD Galerkin vs Exacte - cas 1')

    % Convergence
    p = polyfit(log(H_val), log(err_l2), 1);
    q = polyfit(log(H_val), log(err_H1), 1);
    fprintf('Ordre L2 = %.3f\n', p(1))
    fprintf('Ordre H1 = %.3f\n', q(1))

    figure('Name', 'Convergence LOD Galerkin - cas 1')

    subplot(1, 2, 1)
    loglog(H_val, err_l2,       'b-o', 'LineWidth', 2); hold on
    loglog(H_val, H_val.^2 * 10, 'g-d', 'LineWidth', 1.5)
    loglog(H_val, H_val    * 10, 'k-v', 'LineWidth', 1.5)
    legend(sprintf('LOD L2 (ordre=%.2f)', p(1)), 'O(H^2)', 'O(H)')
    xlabel('H'); title('Erreur L2'); grid on

    subplot(1, 2, 2)
    loglog(H_val, err_H1,       'r-o', 'LineWidth', 2); hold on
    loglog(H_val, H_val    * 10, 'g-v', 'LineWidth', 1.5)
    loglog(H_val, H_val.^2 * 10, 'k-d', 'LineWidth', 1.5)
    legend(sprintf('LOD H1 (ordre=%.2f)', q(1)), 'O(H)', 'O(H^2)')
    xlabel('H'); title('Erreur H1'); grid on


% CAS 2

   elseif choice==2

    err_l2 = zeros(length(epsi), length(N_val));
    err_H1 = zeros(length(epsi), length(N_val));
    H_val  = zeros(size(N_val));

    for k = 1:length(epsi)
        ep = epsi(k);
        fprintf('\n--- epsilon = %.4f ---\n', ep)

        for i = 1:length(N_val)

            n        = N_val(i);
            [x, h]   = mesh1D(n);
            H_val(i) = h;

            aeps_ep = @(t) aeps(t, ep);

            [A_lod, F_lod] = assembleLOD(x, h, x_fin, aeps_ep, f);
            U_lod          = solveLOD(A_lod, F_lod);

            % Reconstruction sur grille fine
            U_fin = reconstructLOD(U_lod, x, x_fin, aeps_ep);

            % Matrices fines pour les erreurs
            M   = Masse(x_fin_bords, h_fin);
            K   = Rigidite(x_fin_bords, h_fin, aeps_ep);
            Uex = exactSolution(x_fin_bords, ep, 2)';

            [err_l2(k,i), err_H1(k,i)] = erreurLOD(U_fin, Uex(:), M, K);

        end
    end

    colors = {'r-o','b-o','g-o','k-o'};
    labels = {'\epsilon=2^{-2}','\epsilon=2^{-4}','\epsilon=2^{-6}','\epsilon=2^{-8}'};

    % Solutions
    n      = N_val(end);
    [x, h] = mesh1D(n);

    figure('Name', 'Solutions LOD Galerkin - cas 2')
    for k = 1:length(epsi)
        ep      = epsi(k);
        aeps_ep = @(t) aeps(t, ep);
        [A_lod, F_lod] = assembleLOD(x, h, x_fin, aeps_ep, f);
        U_lod = solveLOD(A_lod, F_lod);
        U_fin = reconstructLOD(U_lod, x, x_fin, aeps_ep);
        Uex   = exactSolution(x_fin_bords, ep, 2)';
        subplot(2, 2, k)
        plot(x_fin, U_fin, 'b-', 'LineWidth', 1.5); hold on
        plot(x_fin, Uex,   'k-', 'LineWidth', 1.5)
        legend('LOD Galerkin', 'Exacte', 'Location', 'northeast')
        title(['\epsilon = ', num2str(ep)])
        xlabel('x'); grid on
    end

    % Convergence L2
    figure('Name', 'Convergence LOD Galerkin - cas 2')

    subplot(1, 2, 1)
    hold on; grid on
    for k = 1:length(epsi)
        p = polyfit(log(H_val), log(err_l2(k,:)), 1);
        fprintf('eps=%.2e  L2 ordre=%.2f\n', epsi(k), p(1))
        loglog(H_val, err_l2(k,:), colors{k}, 'LineWidth', 2)
    end
    loglog(H_val, H_val.^2 * 10, 'k-d', 'LineWidth', 1.5)
    loglog(H_val, H_val    * 10, 'g-v', 'LineWidth', 1.5)
    title('Erreur L2 - LOD Galerkin'); xlabel('H')
    legend([labels, {'O(H^2)','O(H)'}], 'Location', 'northwest', 'FontSize', 7)

    % Convergence H1
    subplot(1, 2, 2)
    hold on; grid on
    for k = 1:length(epsi)
        q = polyfit(log(H_val), log(err_H1(k,:)), 1);
        fprintf('eps=%.2e  H1 ordre=%.2f\n', epsi(k), q(1))
        loglog(H_val, err_H1(k,:), colors{k}, 'LineWidth', 2)
    end
    loglog(H_val, H_val    * 10, 'm-d', 'LineWidth', 1.5)
    loglog(H_val, H_val.^2 * 10, 'c-v', 'LineWidth', 1.5)
    title('Erreur H1 - LOD Galerkin'); xlabel('H')
    legend([labels, {'O(H)','O(H^2)'}], 'Location', 'northwest', 'FontSize', 7)

    % CAS 2

     else

    err_l2 = zeros(length(epsi), length(N_val));
    err_H1 = zeros(length(epsi), length(N_val));
    H_val  = zeros(size(N_val));

    for k = 1:length(epsi)
        ep = epsi(k);
        fprintf('\n--- epsilon = %.4f ---\n', ep)

        for i = 1:length(N_val)

            n        = N_val(i);
            [x, h]   = mesh1D(n);
            H_val(i) = h;

            aeps_ep = @(t) aeps(t, ep);
            f_ep = @(t) f(t, ep);

            [A_lod, F_lod] = assembleLOD(x, h, x_fin, aeps_ep, f_ep);

            U_lod          = solveLOD(A_lod, F_lod);

            % Reconstruction sur grille fine
            U_fin = reconstructLOD(U_lod, x, x_fin, aeps_ep);

            % Matrices fines pour les erreurs
            M   = Masse(x_fin_bords, h_fin);
            K   = Rigidite(x_fin_bords, h_fin, aeps_ep);
            Uex = exactSolution(x_fin_bords, ep, 3)';

            [err_l2(k,i), err_H1(k,i)] = erreurLOD(U_fin, Uex(:), M, K);

        end
    end

    colors = {'r-o','b-o','g-o','k-o'};
    labels = {'\epsilon=2^{-2}','\epsilon=2^{-4}','\epsilon=2^{-6}','\epsilon=2^{-8}'};

    % Solutions
    n      = N_val(end);
    [x, h] = mesh1D(n);

    figure('Name', 'Solutions LOD Galerkin - cas 2')
    for k = 1:length(epsi)
        ep      = epsi(k);
        aeps_ep = @(t) aeps(t, ep);
        [A_lod, F_lod] = assembleLOD(x, h, x_fin, aeps_ep, f);
        U_lod = solveLOD(A_lod, F_lod);
        U_fin = reconstructLOD(U_lod, x, x_fin, aeps_ep);
        Uex   = exactSolution(x_fin_bords, ep, 3)';
        subplot(2, 2, k)
        plot(x_fin, U_fin, 'b-', 'LineWidth', 1.5); hold on
        plot(x_fin, Uex,   'k-', 'LineWidth', 1.5)
        legend('LOD Galerkin', 'Exacte', 'Location', 'northeast')
        title(['\epsilon = ', num2str(ep)])
        xlabel('x'); grid on
    end

    % Convergence L2
    figure('Name', 'Convergence LOD Galerkin - cas 2')

    subplot(1, 2, 1)
    hold on; grid on
    for k = 1:length(epsi)
        p = polyfit(log(H_val), log(err_l2(k,:)), 1);
        fprintf('eps=%.2e  L2 ordre=%.2f\n', epsi(k), p(1))
        loglog(H_val, err_l2(k,:), colors{k}, 'LineWidth', 2)
    end
    loglog(H_val, H_val.^2 * 10, 'k-d', 'LineWidth', 1.5)
    loglog(H_val, H_val    * 10, 'g-v', 'LineWidth', 1.5)
    title('Erreur L2 - LOD Galerkin'); xlabel('H')
    legend([labels, {'O(H^2)','O(H)'}], 'Location', 'northwest', 'FontSize', 7)

    % Convergence H1
    subplot(1, 2, 2)
    hold on; grid on
    for k = 1:length(epsi)
        q = polyfit(log(H_val), log(err_H1(k,:)), 1);
        fprintf('eps=%.2e  H1 ordre=%.2f\n', epsi(k), q(1))
        loglog(H_val, err_H1(k,:), colors{k}, 'LineWidth', 2)
    end
    loglog(H_val, H_val    * 10, 'm-d', 'LineWidth', 1.5)
    loglog(H_val, H_val.^2 * 10, 'c-v', 'LineWidth', 1.5)
    title('Erreur H1 - LOD Galerkin'); xlabel('H')
    legend([labels, {'O(H)','O(H^2)'}], 'Location', 'northwest', 'FontSize', 7)

end



temps=toc()
