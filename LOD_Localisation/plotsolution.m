% plotSolution.m


%% Parametres
NH_plot   = 16;
n_eps     = length(epsilon_vec);
fsize     = 14;  % taille police generale
fsize_leg = 12;  % taille police legende
lw        = 2;   % epaisseur des courbes

%% Figure
figure('Name','Comparaison LOD vs Reference', ...
       'Position',[100 100 900 700]);

for k = 1:n_eps

    ep    = epsilon_vec(k);
    a_eps = @(x) a_eps_base(x, ep);

    [Xh, ~]        = construire_maillages(Nh, 4);
    [Kh, Mh]       = assembler_systeme_fin(Xh, a_eps);
    F_fine         = Mh * func(Xh)';
    u_ref          = zeros(length(Xh), 1);
    u_ref(2:end-1) = Kh(2:end-1,2:end-1) \ F_fine(2:end-1);

    [~, XH]        = construire_maillages(Nh, NH_plot);
    Phi_ms         = phase_offline(Xh, XH, Kh, k_patch);
    [~, u_LOD]     = phase_online(Xh, XH, Phi_ms, Kh, Mh, func);

    subplot(ceil(n_eps/2), 2, k)
    hold on; grid on;
    plot(Xh, u_ref,  'k-', 'LineWidth', lw);
    plot(Xh, u_LOD,  'b-', 'LineWidth', lw);
    xlabel('x',  'FontSize', fsize);
    ylabel('u',  'FontSize', fsize);
    ylim([-0.05, 0.35]);
    title(sprintf('\\epsilon = %.7g', ep), ...
          'FontWeight', 'bold', 'FontSize', fsize);
    legend({'Exacte', 'LOD Galerkin'}, ...
           'Location', 'north', 'FontSize', fsize_leg);
    set(gca, 'FontSize', fsize);

end

title(sprintf('N_H = %d,  k_{patch} = %d', NH_plot, k_patch), ...
        'FontSize', fsize + 2, 'FontWeight', 'bold');
