tic()

clear
clc
close all


% CHOIX DU TEST

disp('Choisir un cas :')
disp('1 : Cas  constant')
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

        [A, F] = assembleFEM(x, h, epsi, f, aeps);

        U_int = A \ F;
        U = [0; U_int; 0];

        M = Masse(x, h);
        K=Rigidite(x,h);
        Uex = exactSolution(x, epsi, choice)';

        U = U(:);
        Uex = Uex(:);


        errH1(i) = sqrt( (U-Uex)'*(K)*(U-Uex) );
        errl2(i) = sqrt((U - Uex)' * M * (U - Uex));

    end

    plotSolutions(x, U, Uex)
    title('Comparaison U_h vs U_{ex} pour a_{\epsilon} constant');
    fprintf('Erreur L2 = %.6e\n', errl2(end));

    p = polyfit(log(H_val), log(errl2), 1);
    q= polyfit(log(H_val), log(errH1), 1);
    fprintf('\nOrdre L2 = %.3f\n', p(1));
    fprintf('\nOrdre H1 = %.3f\n', q(1));
    figure;
    loglog(H_val, errl2, 'r-o', 'LineWidth', 2);
    hold on;

    loglog(H_val, H_val.^2*10, 'm-d', 'LineWidth', 2);
    loglog(H_val, H_val*10, 'k-v', 'LineWidth', 2);
    loglog(H_val, errH1, 'b-o', 'LineWidth', 2);
    legend(sprintf('Erreur L2 (ordre = %.3f)', p(1)), ...
           'O(h^2)', ...
           'O(h)', ...
           sprintf('Erreur H1 (ordre = %.3f)', q(1)),...
           'Location', 'northwest');

    xlabel('h');
    ylabel('Erreur L2');
    grid on;

% CAS 2

else

    err_all = zeros(length(epsi), length(N_val));
    H_val = zeros(size(N_val));
    err_all_H1=zeros(length(epsi), length(N_val));
    % ERREUR

    for k = 1:length(epsi)

        ep = epsi(k);

        for i = 1:length(N_val)

            n = N_val(i);

            [x, h] = mesh1D(n);
            H_val(i) = h;

            [A, F] = assembleFEM(x, h, ep, f, @(x) aeps(x, ep));

            U_int = A \ F;
            U = [0; U_int; 0];

            M = Masse(x, h);

        K=Rigidite(x,h);
            Uex = exactSolution(x, ep, choice)';

            U = U(:);
            Uex = Uex(:);

            err_all(k, i) = sqrt((U - Uex)' * M * (U - Uex));
            err_all_H1(k, i) = sqrt((U - Uex)' * (M+K) * (U - Uex));

        end

    end


    figure;
    title('Comparaison U_h vs U_{ex} pour différents \epsilon');

    n = N_val(end);
    [x, h] = mesh1D(n);

    for k = 1:length(epsi)

        ep = epsi(k);

        [A, F] = assembleFEM(x, h, ep, f, @(x) aeps(x, ep));

        U_int = A \ F;
        U = [0; U_int; 0];

        Uex = exactSolution(x, ep, choice)';

        subplot(2,2,k)

        plot(x, U, 'r--', 'LineWidth', 1);
        hold on;
        plot(x, Uex, 'k--', 'LineWidth', 1);

        grid on;

        title(['\epsilon = ', num2str(ep)]);

        legend('U_h','U_{ex}');
        xlabel('x');
        ylabel('u(x)');

    end


    figure;

%L2
subplot(1,2,1)
hold on;
grid on;

colors = {'r-o','b-o','g-o','c-o'};
labels = {'2^{-2}','2^{-4}','2^{-6}','2^{-8}'};

for k = 1:length(epsi)

    erreurl2 = err_all(k,:);
    q1 = polyfit(log(H_val), log(erreurl2), 1);

    fprintf('\nPour epsilon= %.3e, Ordre L2 = %.3f\n', epsi(k), q1(1));

    loglog(H_val, erreurl2, colors{k}, 'LineWidth', 2);

end

loglog(H_val, H_val.^2*10, 'm-v', 'LineWidth', 2);
loglog(H_val, H_val*10, 'k-d', 'LineWidth', 2);

title('Convergence L2');
xlabel('h');
ylabel('Erreur L2');
legend([labels, {'O(h^2)','O(h)'}], 'Location','southeast');


%  H1
subplot(1,2,2)
hold on;
grid on;

for k = 1:length(epsi)

    erreurH1 = err_all_H1(k,:);
    q2 = polyfit(log(H_val), log(erreurH1), 1);

    fprintf('\nPour epsilon= %.3e, Ordre H1 = %.3f\n', epsi(k), q2(1));

    loglog(H_val, erreurH1, colors{k}, 'LineWidth', 2);

end

loglog(H_val, H_val*10, 'k-v', 'LineWidth', 2);
loglog(H_val, H_val.^2*10, 'm-d', 'LineWidth', 2);

title('Convergence H1');
xlabel('h');
ylabel('Erreur H1');
legend([labels, {'O(h)','O(h^2)'}], 'Location','northeast');
end


temps=toc()
