clear
clc
close all

% CHOIX DU TEST
disp('Choisir un cas :')
disp('1 : Cas constant')
disp('2 : Cas avec a_epsilon variable')

choice = input('Votre choix = ');

% PARAMÈTRES
[aeps, f, epsi] = parameters(choice);

N_val = 20:20:600;


% CAS 1


if choice == 1

    errl2 = zeros(size(N_val));
    H_val = zeros(size(N_val));

    for i = 1:length(N_val)

        n = N_val(i);

        [x, h] = mesh1D(n);
        H_val(i) = h;

        [A, F] = assembleFEMv2(x, h, epsi, f, aeps, 1);

        U_int = A \ F;
        U = [0; U_int; 0];

        M = Masse(x, h);

        Uex = exactSolution(x, epsi, choice)';

        U = U(:);
        Uex = Uex(:);

        errl2(i) = sqrt((U - Uex)' * M * (U - Uex));

    end

    plotSolutions(x, U, Uex)
    title('Comparaison U_h vs U_{ex} pour a_{\epsilon} constant');

    fprintf('Erreur L2 = %.6e\n', errl2(end));

    p = polyfit(log(H_val), log(errl2), 1);

    fprintf('\nOrdre L2 = %.3f\n', p(1));

    figure;
    loglog(H_val, errl2, 'r-o', 'LineWidth', 2);
    hold on;

    loglog(H_val, H_val.^2*10, 'm-d', 'LineWidth', 2);
    loglog(H_val, H_val*10, 'k-v', 'LineWidth', 2);

    legend(sprintf('Erreur L2 (ordre = %.3f)', p(1)), ...
           'O(h^2)', ...
           'O(h)', ...
           'Location', 'northwest');

    xlabel('h');
    ylabel('Erreur L2');
    grid on;

% CAS 2


else

    err_all = zeros(length(epsi), length(N_val));
    H_val = zeros(size(N_val));

    % ERREUR
    for k = 1:length(epsi)

        ep = epsi(k);

        for i = 1:length(N_val)

            n = N_val(i);

            [x, h] = mesh1D(n);
            H_val(i) = h;

            [A, F] = assembleFEMv2(x, h, ep, f, aeps, 2);

            U_int = A \ F;
            U = [0; U_int; 0];

            M = Masse(x, h);

            Uex = exactSolution(x, ep, choice)';

            U = U(:);
            Uex = Uex(:);

            err_all(k,i) = sqrt((U - Uex)' * M * (U - Uex));

        end
    end


    % AFFICHAGE SOLUTIONS


    figure;
    n = N_val(end);
    [x, h] = mesh1D(n);

    for k = 1:length(epsi)

        ep = epsi(k);

        [A, F] = assembleFEMv2(x, h, ep, f, aeps, 2);

        U_int = A \ F;
        U = [0; U_int; 0];

        Uex = exactSolution(x, ep, choice)';

        subplot(2,2,k)
        plot(x, U, 'r--', 'LineWidth', 1);
        hold on;
        plot(x, Uex, 'k-', 'LineWidth', 1);

        grid on;
        title(['\epsilon = ', num2str(ep)]);

        legend('U_h','U_{ex}');
        xlabel('x');
        ylabel('u(x)');

    end


    % CONVERGENCE


    figure;
    hold on;
    grid on;

    colors = {'r-o','b-o','g-o','c-o'};
    labels = {'2^{-2}','2^{-4}','2^{-6}','2^{-8}'};

    for k = 1:length(epsi)

        erreur = err_all(k,:);

        q = polyfit(log(H_val), log(erreur), 1);

        fprintf('\nPour epsilon= %.3e, Ordre L2 = %.3f\n', epsi(k), q(1));

        loglog(H_val, erreur, colors{k}, 'LineWidth', 2);

    end

    loglog(H_val, H_val.^2/10, 'm-v', 'LineWidth', 2);
    loglog(H_val, H_val/10, 'k-d', 'LineWidth', 2);

    legend([labels, {'O(h^2)','O(h)'}]);

    xlabel('h');
    ylabel('Erreur L2');
    title('Convergence pour différents \epsilon');

end
