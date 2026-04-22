%  ISAC System: Dinkelbach-SCA Algorithm for Multi-Objective EE Optimization
%  BOUNDED uncertainty model  –  Figures 3 to 8
%  Benchmarks per figure:
%    (1) Robust          : r_k = p.r_k_vec,  i_b = p.i_b,  i_k = p.i_k
%    (2) Perfect CSI     : r_k = 0,           i_b = p.i_b,  i_k = p.i_k
%    (3) Perfect Hardware: r_k = p.r_k_vec,  i_b = IMP_EPS, i_k = IMP_EPS
%    (4) Perfect Both    : r_k = 0,           i_b = IMP_EPS, i_k = IMP_EPS
%  All four curves share ONE set of normalization bounds per figure.
clc; clear; close all;

%% =========================================================================
%  LOCK RANDOM SEED
rng(42);

%% =========================================================================
%  LOAD SYSTEM PARAMETERS
p = parameters();

%  Small but nonzero impairment floor to keep solver stable when i_b=i_k~0
IMP_EPS = 0;

%% =========================================================================
%  GENERATE BASELINE CHANNEL  H_hat  (N x K)
H_hat = generate_channels(p.N, p.K, 1);
save('saved_channels.mat', 'H_hat', 'p');
fprintf('Channel realisation saved to saved_channels.mat\n\n');

%% =========================================================================
%  COMPUTE UNCERTAINTY RADII  r_k  FOR THE BASELINE CHANNEL
for k = 1:p.K
    p.xi_k_vec(k) = p.delta * norm(H_hat(:, k));
    p.r_k_vec(k)  = sqrt( (p.xi_k_vec(k)^2 / 2) * ...
                     chi2inv(1 - p.delta_outage, 2*p.N) );
end
fprintf('Baseline uncertainty radii r_k computed.\n\n');

%  Shorthand benchmark parameter sets
r_k_rob  = p.r_k_vec;
r_k_perf = zeros(1, p.K);

v_rk  = {r_k_rob,  r_k_perf, r_k_rob,  r_k_perf};
v_ib  = {p.i_b,    p.i_b,    IMP_EPS,  IMP_EPS };
v_ik  = {p.i_k,    p.i_k,    IMP_EPS,  IMP_EPS };
v_lbl = {'Robust','Perfect CSI','Perfect HW','Perfect Both'};

%  Shared plot styles
styles    = {'r-s','k--^','b-o','m--d'};
styles_ec = {'r-s','k-^','b-o','m-d'};
styles_es = {'r-s','k-^','b-o','m-d'};
mfc       = {'r',  'k',   'b',  'm'  };

ec_lbl = cellfun(@(x) ['EE_c ' x], v_lbl, 'UniformOutput', false);
es_lbl = cellfun(@(x) ['EE_s ' x], v_lbl, 'UniformOutput', false);

%% =========================================================================
%  GENERATE H_bank FOR FIGURE 7
N_vec_f7 = 6:1:9;
H_bank   = generate_channels(max(N_vec_f7), p.K, p.N_MC);
fprintf('H_bank generated: %d x %d x %d.\n\n', max(N_vec_f7), p.K, p.N_MC);

%% =========================================================================
%  PHASE 1 – GLOBAL NORMALISATION CONSTANTS  (baseline parameters)
fprintf('=== Phase 1: Global normalisation constants ===\n');

[EEcmax, EEcmin, EEsmax, EEsmin] = get_norm_constants( ...
    H_hat, p.theta_targets, p.N, p.K, p.M, ...
    p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
    p.i_b, p.i_k, p.P_static, p.r_k_vec, p.T_max, p.epsilon, p.N_rand);

fprintf('\nEE_c,max = %.4f  |  EE_c,min = %.4f\n', EEcmax, EEcmin);
fprintf('EE_s,max = %.4f  |  EE_s,min = %.4f\n\n', EEsmax, EEsmin);
%% =========================================================================
%  FIGURE 2 – IMPACT OF WEIGHTING COEFFICIENT omega
% --------------------------------------------------------------------------
fprintf('=== Figure 2: Impact of omega ===\n');

omega_vec   = 0:0.05:1;
n_om        = numel(omega_vec);

% Weighted normalised contributions  (what the objective actually uses)
wEEc_om = zeros(n_om, 1);   %  omega       * (EEc - EEcmin)/(EEcmax - EEcmin)
wEEs_om = zeros(n_om, 1);   %  (1-omega)   * (EEs - EEsmin)/(EEsmax - EEsmin)

for idx = 1:n_om
    om = omega_vec(idx);
    fprintf('  omega=%.2f\n', om);
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
        p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
        p.i_b, p.i_k, p.P_static, p.r_k_vec, om, ...
        p.T_max, p.epsilon, EEcmax, EEsmax, p.N_rand, EEcmin, EEsmin);

    % Normalise using the global bounds computed in Phase 1
    ec_norm = (ec - EEcmin) / (EEcmax - EEcmin);
    es_norm = (es - EEsmin) / (EEsmax - EEsmin);

    % Apply the omega weighting  ← this is the fix
    wEEc_om(idx) = om       * ec_norm;
    wEEs_om(idx) = (1-om)   * es_norm;
end

% ------------------------------------------------------------------
%  Find crossover:  omega * EEc_norm  =  (1-omega) * EEs_norm
%  i.e. the zero-crossing of  (wEEc_om - wEEs_om)
% ------------------------------------------------------------------
diff_vec  = wEEc_om - wEEs_om;
cross_idx = find(diff(sign(diff_vec)) ~= 0, 1);   % last sign change

if ~isempty(cross_idx)
    % Linear interpolation between cross_idx and cross_idx+1
    om1 = omega_vec(cross_idx);   d1 = diff_vec(cross_idx);
    om2 = omega_vec(cross_idx+1); d2 = diff_vec(cross_idx+1);
    omega_cross = om1 - d1*(om2-om1)/(d2-d1);
    wEEc_cross  = interp1(omega_vec, wEEc_om, omega_cross, 'linear');
    fprintf('  Crossover omega* = %.4f  (weighted EEc_norm = EEs_norm = %.4f)\n', ...
            omega_cross, wEEc_cross);
else
    omega_cross = NaN;
    fprintf('  No crossover found in [0,1].\n');
end
p.omega = omega_cross; 
% ------------------------------------------------------------------
%  Plot
% ------------------------------------------------------------------
figure(2); clf;
set(gcf, 'Position', [100 100 750 480]);
 
yyaxis left;
plot(omega_vec, wEEc_om, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, ...
    'MarkerFaceColor', 'b', 'DisplayName', '\omega \cdot EE_c^{norm}');
ylabel('\omega \cdot EE_c^{norm}', 'FontSize', 13, 'Color', 'b');
set(gca, 'YColor', 'b');
 
yyaxis right;
plot(omega_vec, wEEs_om, 'r-s', 'LineWidth', 2, 'MarkerSize', 5, ...
    'MarkerFaceColor', 'r', 'DisplayName', '(1-\omega) \cdot EE_s^{norm}');
ylabel('(1-\omega) \cdot EE_s^{norm}', 'FontSize', 13, 'Color', 'r');
set(gca, 'YColor', 'r');
 
% Mark crossover point
if ~isempty(cross_idx)
    yyaxis left;
    hold on;
    plot(omega_cross, wEEc_cross, 'kp', 'MarkerSize', 12, ...
        'MarkerFaceColor', 'k', ...
        'DisplayName', sprintf('Crossover \\omega^* = %.3f', omega_cross));
    hold off;
end
 
grid on;
xlabel('Weighting Coefficient  \omega', 'FontSize', 13);
set(gca, 'XTick', omega_vec);
title('Fig. 2: Weighted EE_c^{norm} and EE_s^{norm} vs. \omega', 'FontSize', 13);
legend('show', 'Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 12);
drawnow;
fprintf('    Figure 2 rendered.\n\n');

%% =========================================================================
%  FIGURE 1 – CONVERGENCE OF DINKELBACH-SCA
%  Now plots: total WEE, omega*EEc_norm, and (1-omega)*EEs_norm per iteration
% --------------------------------------------------------------------------
fprintf('=== Figure 1: Convergence ===\n');

[Wc_f1, Ws_f1, EEc_f1, EEs_f1, obj_hist,EEc_hist, EEs_hist] = ...
    solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
               p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.r_k_vec, ...
               p.omega, p.T_max, p.epsilon, EEcmax, EEsmax, p.N_rand, EEcmin, EEsmin);


% ---- plot ----------------------------------------------------------------
n_iter = numel(obj_hist);
iter_ax = 1:n_iter;
EEc_hist_norm = p.omega     * (EEc_hist - EEcmin) ./ (EEcmax - EEcmin);
EEs_hist_norm = (1-p.omega) * (EEs_hist - EEsmin) ./ (EEsmax - EEsmin);

figure(1);
plot(iter_ax, obj_hist,       'b-o',  'LineWidth',2, 'MarkerSize',6); hold on;
plot(iter_ax, EEc_hist_norm,  'r--s', 'LineWidth',1.8,'MarkerSize',5);hold on; 
plot(iter_ax, EEs_hist_norm,  'g--^', 'LineWidth',1.8,'MarkerSize',5);hold on;
hold off; grid on;
xlabel('Number of Iterations','FontSize',13);
ylabel('Normalised Value','FontSize',13);
title(sprintf('Fig. 1: Convergence of Dinkelbach-SCA  (\\omega=%.2f)', p.omega),'FontSize',13);
%str = sprintf('EE_c = %.4f bit/J/Hz\nEE_s = %.4f', EEc_f1, EEs_f1);
%text(0.97, 0.15, str, 'Units','normalized','HorizontalAlignment','right', ...
%    'VerticalAlignment','bottom','FontSize',11,'BackgroundColor','white', ...
%    'EdgeColor','black','Margin',5);
legend(sprintf('WEE (\\omega=%.2f)', p.omega), ...
       sprintf('\\omega \\cdot EE_c^{norm}'), ...
       sprintf('(1-\\omega) \\cdot EE_s^{norm}'), ...
       'Location','best','FontSize',8);
set(gca,'FontSize',12); drawnow;
fprintf('    Figure 1 rendered.\n\n');


%% =========================================================================
%  FIGURE 3 – IMPACT OF SENSING CONSTRAINT Gamma_min  (4 benchmarks)
%  Layout: 3 subplots – (a) WEE, (b) EEc, (c) EEs
% --------------------------------------------------------------------------
fprintf('=== Figure 3: Gamma_min sweep (4 benchmarks) ===\n');

Gamma_dBm_vec = 16:2:24;
Gamma_W_vec   = db2pow(Gamma_dBm_vec) * 1e-3;
n_gam         = numel(Gamma_dBm_vec);

% --- Pass 1: shared bounds ---
ec_max_all = zeros(n_gam,4);  ec_min_all = zeros(n_gam,4);
es_max_all = zeros(n_gam,4);  es_min_all = zeros(n_gam,4);

for v = 1:4
    for idx = 1:n_gam
        [ec_max_all(idx,v), ec_min_all(idx,v), ...
         es_max_all(idx,v), es_min_all(idx,v)] = ...
            get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
                p.P_max, p.sigma2, Gamma_W_vec(idx), p.gamma_min, ...
                v_ib{v}, v_ik{v}, p.P_static, v_rk{v}, p.T_max, p.epsilon, p.N_rand);
    end
end

[EEc_max_f3, best_idx_f3] = max(ec_max_all(:,1));
EEc_min_f3 = ec_min_all(best_idx_f3, 1);
EEs_max_f3 = es_max_all(best_idx_f3, 1);
EEs_min_f3 = es_min_all(best_idx_f3, 1);
fprintf('  Fig3 shared EEc bounds: [%.4f, %.4f]\n', EEc_min_f3, EEc_max_f3);
fprintf('  Fig3 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f3, EEs_max_f3);

% --- Pass 2: solve all 4 variants ---
WEE_f3 = zeros(n_gam,4);
EEc_f3 = zeros(n_gam,4);
EEs_f3 = zeros(n_gam,4);

for v = 1:4
    for idx = 1:n_gam
        fprintf('  [Solve] %-16s  Gamma=%d dBm\n', v_lbl{v}, Gamma_dBm_vec(idx));
        [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, Gamma_W_vec(idx), p.gamma_min, ...
            v_ib{v}, v_ik{v}, p.P_static, v_rk{v}, p.omega, ...
            p.T_max, p.epsilon, EEc_max_f3, EEs_max_f3, p.N_rand, EEc_min_f3, EEs_min_f3);
        EEc_f3(idx,v) = ec;  EEs_f3(idx,v) = es;
        WEE_f3(idx,v) = p.omega     * ((ec - EEc_min_f3)/(EEc_max_f3 - EEc_min_f3)) + ...
                        (1-p.omega) * ((es - EEs_min_f3)/(EEs_max_f3 - EEs_min_f3));
    end
end

% --- Plot: 3 subplots ---
figure(3); set(gcf,'Position',[100 100 1600 450]);

subplot(1,3,1);
for v = 1:4
    plot(Gamma_dBm_vec, WEE_f3(:,v), styles{v}, 'LineWidth',2, ...
        'MarkerSize',7,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Sensing Threshold  \Gamma_{min}  (dBm)','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('(a) WEE  (\\omega=%.2f)', p.omega),'FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

subplot(1,3,2);
for v = 1:4
    plot(Gamma_dBm_vec, EEc_f3(:,v), styles_ec{v}, 'LineWidth',2, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Sensing Threshold  \Gamma_{min}  (dBm)','FontSize',13);
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13);
title('(b) EE_c','FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

subplot(1,3,3);
for v = 1:4
    plot(Gamma_dBm_vec, EEs_f3(:,v), styles_es{v}, 'LineWidth',2, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Sensing Threshold  \Gamma_{min}  (dBm)','FontSize',13);
ylabel('Sensing EE  (EE_s)','FontSize',13);
title('(c) EE_s','FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

sgtitle(sprintf('Fig. 3: Impact of \\Gamma_{min}  (\\omega=%.2f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 3 rendered.\n\n');

%% =========================================================================
%  FIGURE 4 – IMPACT OF gamma_min (SINR)  (4 benchmarks)
%  Layout: 3 subplots – (a) WEE, (b) EEc, (c) EEs
% --------------------------------------------------------------------------
fprintf('=== Figure 4: gamma_min sweep (4 benchmarks) ===\n');

gamma_dB_vec  = 3:1:8;
gamma_lin_vec = db2pow(gamma_dB_vec);
n_sinr        = numel(gamma_dB_vec);

% --- Pass 1: shared bounds ---
ec_max_all = zeros(n_sinr,4);  ec_min_all = zeros(n_sinr,4);
es_max_all = zeros(n_sinr,4);  es_min_all = zeros(n_sinr,4);

for v = 1:4
    for idx = 1:n_sinr
        [ec_max_all(idx,v), ec_min_all(idx,v), ...
         es_max_all(idx,v), es_min_all(idx,v)] = ...
            get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
                p.P_max, p.sigma2, p.Gamma_min, gamma_lin_vec(idx), ...
                v_ib{v}, v_ik{v}, p.P_static, v_rk{v}, p.T_max, p.epsilon, p.N_rand);
    end
end

[EEc_max_f4, best_idx_f4] = max(ec_max_all(:,1));
EEc_min_f4 = ec_min_all(best_idx_f4, 1);
EEs_max_f4 = es_max_all(best_idx_f4, 1);
EEs_min_f4 = es_min_all(best_idx_f4, 1);
fprintf('  Fig4 shared EEc bounds: [%.4f, %.4f]\n', EEc_min_f4, EEc_max_f4);
fprintf('  Fig4 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f4, EEs_max_f4);

% --- Pass 2: solve ---
WEE_f4 = zeros(n_sinr,4);
EEc_f4 = zeros(n_sinr,4);
EEs_f4 = zeros(n_sinr,4);

for v = 1:4
    for idx = 1:n_sinr
        fprintf('  [Solve] %-16s  gamma=%d dB\n', v_lbl{v}, gamma_dB_vec(idx));
        [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, gamma_lin_vec(idx), ...
            v_ib{v}, v_ik{v}, p.P_static, v_rk{v}, p.omega, ...
            p.T_max, p.epsilon, EEc_max_f4, EEs_max_f4, p.N_rand, EEc_min_f4, EEs_min_f4);
        EEc_f4(idx,v) = ec;  EEs_f4(idx,v) = es;
        WEE_f4(idx,v) = p.omega     * ((ec - EEc_min_f4)/(EEc_max_f4 - EEc_min_f4)) + ...
                        (1-p.omega) * ((es - EEs_min_f4)/(EEs_max_f4 - EEs_min_f4));
    end
end

% --- Plot: 3 subplots ---
figure(4); set(gcf,'Position',[100 100 1600 450]);

subplot(1,3,1);
for v = 1:4
    plot(gamma_dB_vec, WEE_f4(:,v), styles{v}, 'LineWidth',2, ...
        'MarkerSize',7,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Minimum SINR  \gamma_{min}  (dB)','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('(a) WEE  (\\omega=%.2f)', p.omega),'FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

subplot(1,3,2);
for v = 1:4
    plot(gamma_dB_vec, EEc_f4(:,v), styles_ec{v}, 'LineWidth',2, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Minimum SINR  \gamma_{min}  (dB)','FontSize',13);
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13);
title('(b) EE_c','FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

subplot(1,3,3);
for v = 1:4
    plot(gamma_dB_vec, EEs_f4(:,v), styles_es{v}, 'LineWidth',2, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Minimum SINR  \gamma_{min}  (dB)','FontSize',13);
ylabel('Sensing EE  (EE_s)','FontSize',13);
title('(c) EE_s','FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

sgtitle(sprintf('Fig. 4: Impact of \\gamma_{min}  (\\omega=%.2f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 4 rendered.\n\n');

%% =========================================================================
%  FIGURE 5 – IMPACT OF HARDWARE IMPAIRMENT COEFFICIENTS (i_b, i_k)
%  (unchanged – single benchmark, bar chart layout)
% --------------------------------------------------------------------------
fprintf('=== Figure 5: Hardware impairments sweep ===\n');

i_b_vec    = [0.01, 0.02, 0.03];
i_k_vec    = [0.02, 0.04, 0.06];
n_imp      = numel(i_b_vec);
imp_labels = arrayfun(@(a,b) sprintf('(%.2f,%.2f)',a,b), ...
                      i_b_vec, i_k_vec, 'UniformOutput', false);

% --- Pass 1: shared bounds ---
ec_max_imp = zeros(n_imp,1);  ec_min_imp = zeros(n_imp,1);
es_max_imp = zeros(n_imp,1);  es_min_imp = zeros(n_imp,1);

for idx = 1:n_imp
    fprintf('  [Norm] i_b=%.3f  i_k=%.3f\n', i_b_vec(idx), i_k_vec(idx));
    [ec_max_imp(idx), ec_min_imp(idx), es_max_imp(idx), es_min_imp(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
            i_b_vec(idx), i_k_vec(idx), p.P_static, ...
            p.r_k_vec, p.T_max, p.epsilon, p.N_rand);
end

EEc_max_f5 = max(ec_max_imp);  EEc_min_f5 = min(ec_min_imp);
EEs_max_f5 = max(es_max_imp);  EEs_min_f5 = min(es_min_imp);
fprintf('  Fig5 EEc bounds: [%.4f, %.4f]\n', EEc_min_f5, EEc_max_f5);
fprintf('  Fig5 EEs bounds: [%.4f, %.4f]\n', EEs_min_f5, EEs_max_f5);

% --- Pass 2: solve ---
WEE_imp = zeros(n_imp,1);
EEc_imp = zeros(n_imp,1);
EEs_imp = zeros(n_imp,1);

for idx = 1:n_imp
    fprintf('  [Solve] i_b=%.3f  i_k=%.3f\n', i_b_vec(idx), i_k_vec(idx));
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
        p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
        i_b_vec(idx), i_k_vec(idx), p.P_static, p.r_k_vec, p.omega, ...
        p.T_max, p.epsilon, EEc_max_f5, EEs_max_f5, p.N_rand, EEc_min_f5, EEs_min_f5);
    EEc_imp(idx) = ec;  EEs_imp(idx) = es;
    WEE_imp(idx) = p.omega     * ((ec - EEc_min_f5)/(EEc_max_f5 - EEc_min_f5)) + ...
                   (1-p.omega) * ((es - EEs_min_f5)/(EEs_max_f5 - EEs_min_f5));
    fprintf('    EEc=%.4f  EEs=%.4f  WEE=%.4f\n', ec, es, WEE_imp(idx));
end

x_ticks = 1:n_imp;
figure(5); set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
bar(x_ticks, WEE_imp, 0.5, 'FaceColor',[0.2 0.5 0.8],'EdgeColor','k');
grid on;
set(gca,'XTick',x_ticks,'XTickLabel',imp_labels,'FontSize',12);
xlabel('(i_b,  i_k)  pairs','FontSize',13);
ylabel('Normalised WEE','FontSize',13);
title(sprintf('(a) WEE  (\\omega=%.2f)', p.omega),'FontSize',13);
xtickangle(15);

subplot(1,2,2);
yyaxis left;
plot(x_ticks, EEc_imp, 'b-o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b');
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13,'Color','b');
set(gca,'YColor','b','XTick',x_ticks,'XTickLabel',imp_labels,'FontSize',12);
yyaxis right;
plot(x_ticks, EEs_imp, 'r-s','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r');
ylabel('Sensing EE  (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r'); grid on;
xlabel('(i_b,  i_k)  pairs','FontSize',13);
title('(b) Raw EE_c and EE_s','FontSize',13);
legend('EE_c','EE_s','Location','best','FontSize',11);
xtickangle(15);

sgtitle(sprintf('Fig. 5: Impact of Hardware Impairments  (\\omega=%.2f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 5 rendered.\n\n');

%% =========================================================================
%  FIGURE 6 – WEIGHTED EE vs P_max  (4 benchmarks)
%  Layout: 3 subplots – (a) WEE, (b) EEc, (c) EEs
% --------------------------------------------------------------------------
fprintf('=== Figure 6: P_max sweep (4 benchmarks) ===\n');

Pmax_dBm_vec = [20 25 30 35];
Pmax_W_vec   = db2pow(Pmax_dBm_vec) * 1e-3;
n_pmax       = numel(Pmax_dBm_vec);

% --- Pass 1: shared bounds ---
ec_max_all = zeros(n_pmax,4);  ec_min_all = zeros(n_pmax,4);
es_max_all = zeros(n_pmax,4);  es_min_all = zeros(n_pmax,4);

for v = 1:4
    for idx = 1:n_pmax
        [ec_max_all(idx,v), ec_min_all(idx,v), ...
         es_max_all(idx,v), es_min_all(idx,v)] = ...
            get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
                Pmax_W_vec(idx), p.sigma2, p.Gamma_min, p.gamma_min, ...
                v_ib{v}, v_ik{v}, p.P_static, v_rk{v}, p.T_max, p.epsilon, p.N_rand);
    end
end

[EEc_max_f6, best_idx_f6] = max(ec_max_all(:,1));
EEc_min_f6 = ec_min_all(best_idx_f6, 1);
EEs_max_f6 = es_max_all(best_idx_f6, 1);
EEs_min_f6 = es_min_all(best_idx_f6, 1);
fprintf('  Fig6 shared EEc bounds: [%.4f, %.4f]\n', EEc_min_f6, EEc_max_f6);
fprintf('  Fig6 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f6, EEs_max_f6);

% --- Pass 2: solve ---
WEE_f6 = zeros(n_pmax,4);
EEc_f6 = zeros(n_pmax,4);
EEs_f6 = zeros(n_pmax,4);

for v = 1:4
    for idx = 1:n_pmax
        fprintf('  [Solve] %-16s  P_max=%d dBm\n', v_lbl{v}, Pmax_dBm_vec(idx));
        [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            Pmax_W_vec(idx), p.sigma2, p.Gamma_min, p.gamma_min, ...
            v_ib{v}, v_ik{v}, p.P_static, v_rk{v}, p.omega, ...
            p.T_max, p.epsilon, EEc_max_f6, EEs_max_f6, p.N_rand, EEc_min_f6, EEs_min_f6);
        EEc_f6(idx,v) = ec;  EEs_f6(idx,v) = es;
        WEE_f6(idx,v) = p.omega     * ((ec - EEc_min_f6)/(EEc_max_f6 - EEc_min_f6)) + ...
                        (1-p.omega) * ((es - EEs_min_f6)/(EEs_max_f6 - EEs_min_f6));
    end
end

% --- Plot: 3 subplots ---
figure(6); set(gcf,'Position',[100 100 1600 450]);

subplot(1,3,1);
for v = 1:4
    plot(Pmax_dBm_vec, WEE_f6(:,v), styles{v}, 'LineWidth',2, ...
        'MarkerSize',7,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Maximum Transmit Power  P_{max}  (dBm)','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('(a) WEE  (\\omega=%.2f)', p.omega),'FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

subplot(1,3,2);
for v = 1:4
    plot(Pmax_dBm_vec, EEc_f6(:,v), styles_ec{v}, 'LineWidth',2, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Maximum Transmit Power  P_{max}  (dBm)','FontSize',13);
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13);
title('(b) EE_c','FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

subplot(1,3,3);
for v = 1:4
    plot(Pmax_dBm_vec, EEs_f6(:,v), styles_es{v}, 'LineWidth',2, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Maximum Transmit Power  P_{max}  (dBm)','FontSize',13);
ylabel('Sensing EE  (EE_s)','FontSize',13);
title('(c) EE_s','FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);

sgtitle(sprintf('Fig. 6: Weighted EE vs. P_{max}  (\\omega=%.2f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 6 rendered.\n\n');

%% =========================================================================
%  FIGURE 7 – WEIGHTED EE vs NUMBER OF ANTENNAS N  (4 benchmarks)
%  Only subplot (a): WEE vs N  –  EEc/EEs subplot removed
% --------------------------------------------------------------------------
fprintf('=== Figure 7: N sweep (4 benchmarks) ===\n');

n_N = numel(N_vec_f7);

% --- Pass 1: shared bounds ---
ec_max_all = zeros(n_N,4);  ec_min_all = zeros(n_N,4);
es_max_all = zeros(n_N,4);  es_min_all = zeros(n_N,4);

for v = 1:4
    for idx = 1:n_N
        N_cur  = N_vec_f7(idx);
        H_norm = H_bank(1:N_cur, :, 1);

        if ismember(v, [1 3])
            r_k_cur = zeros(1, p.K);
            for k = 1:p.K
                xi_k       = p.delta * norm(H_norm(:,k));
                r_k_cur(k) = sqrt((xi_k^2/2) * chi2inv(1-p.delta_outage, 2*N_cur));
            end
        else
            r_k_cur = zeros(1, p.K);
        end

        [ec_max_all(idx,v), ec_min_all(idx,v), ...
         es_max_all(idx,v), es_min_all(idx,v)] = ...
            get_norm_constants(H_norm, p.theta_targets, N_cur, p.K, p.M, ...
                p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
                v_ib{v}, v_ik{v}, p.P_static, r_k_cur, p.T_max, p.epsilon, p.N_rand);
    end
end

[EEc_max_f7, best_idx_f7] = max(ec_max_all(:,1));
EEc_min_f7 = ec_min_all(best_idx_f7, 1);
EEs_max_f7 = es_max_all(best_idx_f7, 1);
EEs_min_f7 = es_min_all(best_idx_f7, 1);
fprintf('  Fig7 shared EEc bounds: [%.4f, %.4f]\n', EEc_min_f7, EEc_max_f7);
fprintf('  Fig7 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f7, EEs_max_f7);

% --- Pass 2: solve ---
WEE_f7 = zeros(n_N,4);
EEc_f7 = zeros(n_N,4);   %#ok<NASGU>  kept for possible future use
EEs_f7 = zeros(n_N,4);   %#ok<NASGU>

for v = 1:4
    for idx = 1:n_N
        N_cur = N_vec_f7(idx);
        H_cur = H_bank(1:N_cur, :, 1);

        if ismember(v, [1 3])
            r_k_cur = zeros(1, p.K);
            for k = 1:p.K
                xi_k       = p.delta * norm(H_cur(:,k));
                r_k_cur(k) = sqrt((xi_k^2/2) * chi2inv(1-p.delta_outage, 2*N_cur));
            end
        else
            r_k_cur = zeros(1, p.K);
        end

        fprintf('  [Solve] %-16s  N=%d\n', v_lbl{v}, N_cur);
        [~,~,ec,es,~] = solve_ISAC(H_cur, p.theta_targets, N_cur, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
            v_ib{v}, v_ik{v}, p.P_static, r_k_cur, p.omega, ...
            p.T_max, p.epsilon, EEc_max_f7, EEs_max_f7, p.N_rand, EEc_min_f7, EEs_min_f7);
        EEc_f7(idx,v) = ec;
        EEs_f7(idx,v) = es;
        WEE_f7(idx,v) = p.omega     * ((ec - EEc_min_f7)/(EEc_max_f7 - EEc_min_f7)) + ...
                        (1-p.omega) * ((es - EEs_min_f7)/(EEs_max_f7 - EEs_min_f7));
        fprintf('    WEE=%.4f\n', WEE_f7(idx,v));
    end
end

% --- Plot: single subplot (WEE only) ---
figure(7); set(gcf,'Position',[100 100 650 480]);

for v = 1:4
    plot(N_vec_f7, WEE_f7(:,v), styles{v}, 'LineWidth',2, ...
        'MarkerSize',7,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Number of Antennas  N','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('Fig. 7: Weighted EE vs. N  (\\omega=%.2f)', p.omega), ...
      'FontSize',14,'FontWeight','bold');
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);
drawnow;
fprintf('    Figure 7 rendered.\n\n');

%% =========================================================================
%  FIGURE 8 – IMPACT OF CHANNEL UNCERTAINTY LEVEL delta  (unchanged)
% --------------------------------------------------------------------------
fprintf('=== Figure 8: Channel uncertainty level delta ===\n');

delta_vec = [0, 0.005, 0.01, 0.02, 0.03];
n_delta   = numel(delta_vec);

% --- Pass 1: shared bounds ---
ec_max_all = zeros(n_delta,1);  ec_min_all = zeros(n_delta,1);
es_max_all = zeros(n_delta,1);  es_min_all = zeros(n_delta,1);

for idx = 1:n_delta
    delta_cur = delta_vec(idx);
    r_k_cur   = zeros(1, p.K);
    for k = 1:p.K
        xi_k       = delta_cur * norm(H_hat(:,k));
        r_k_cur(k) = sqrt((xi_k^2/2) * chi2inv(1-p.delta_outage, 2*p.N));
    end
    fprintf('  [Norm] delta=%.4f\n', delta_cur);
    [ec_max_all(idx), ec_min_all(idx), es_max_all(idx), es_min_all(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
            p.i_b, p.i_k, p.P_static, r_k_cur, p.T_max, p.epsilon, p.N_rand);
end

EEc_max_f8 = max(ec_max_all);  EEc_min_f8 = min(ec_min_all);
EEs_max_f8 = max(es_max_all);  EEs_min_f8 = min(es_min_all);
fprintf('  Fig8 shared EEc bounds: [%.4f, %.4f]\n', EEc_min_f8, EEc_max_f8);
fprintf('  Fig8 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f8, EEs_max_f8);

% --- Pass 2: solve ---
WEE_delta = zeros(n_delta,1);
EEc_delta = zeros(n_delta,1);
EEs_delta = zeros(n_delta,1);

for idx = 1:n_delta
    delta_cur = delta_vec(idx);
    r_k_cur   = zeros(1, p.K);
    for k = 1:p.K
        xi_k       = delta_cur * norm(H_hat(:,k));
        r_k_cur(k) = sqrt((xi_k^2/2) * chi2inv(1-p.delta_outage, 2*p.N));
    end
    fprintf('  [Solve] delta=%.4f\n', delta_cur);
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
        p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, r_k_cur, p.omega, ...
        p.T_max, p.epsilon, EEc_max_f8, EEs_max_f8, p.N_rand, EEc_min_f8, EEs_min_f8);
    EEc_delta(idx) = ec;
    EEs_delta(idx) = es;
    WEE_delta(idx) = p.omega     * ((ec - EEc_min_f8)/(EEc_max_f8 - EEc_min_f8)) + ...
                     (1-p.omega) * ((es - EEs_min_f8)/(EEs_max_f8 - EEs_min_f8));
    fprintf('    EEc=%.4f  EEs=%.4f  WEE=%.4f\n', ec, es, WEE_delta(idx));
end

perf_idx = 1;
rob_idx  = 2:n_delta;

WEE_perf_val = WEE_delta(perf_idx);
EEc_perf_val = EEc_delta(perf_idx);
EEs_perf_val = EEs_delta(perf_idx);

figure(8); set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
plot(delta_vec(rob_idx), WEE_delta(rob_idx), ...
     'b-o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b');
hold on;
yline(WEE_perf_val, 'k--', 'LineWidth',1.5);
hold off; grid on;
xlabel('Channel Uncertainty Level  \delta','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('(a) WEE  (\\omega=%.2f)', p.omega),'FontSize',13);
legend('Robust (Bounded)','Perfect Channel (\delta=0)', ...
       'Location','south','FontSize',11);
set(gca,'FontSize',12);

subplot(1,2,2);
yyaxis left;
plot(delta_vec(rob_idx), EEc_delta(rob_idx), ...
     'b-o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','b');
hold on;
yline(EEc_perf_val, 'b--', 'LineWidth',1.5);
hold off;
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13,'Color','b');
set(gca,'YColor','b');
yyaxis right;
plot(delta_vec(rob_idx), EEs_delta(rob_idx), ...
     'r-s','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','r');
hold on;
yline(EEs_perf_val, 'r--', 'LineWidth',1.5);
hold off;
ylabel('Sensing EE  (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r'); grid on;
xlabel('Channel Uncertainty Level  \delta','FontSize',13);
title('(b) Raw EE_c and EE_s','FontSize',13);
legend('EE_c Robust','EE_c Perfect (\delta=0)', ...
       'EE_s Robust','EE_s Perfect (\delta=0)', ...
       'Location','best','FontSize',11);
set(gca,'FontSize',12);

sgtitle(sprintf('Fig. 8: Impact of Channel Uncertainty Level \\delta  (\\omega=%.2f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 8 rendered.\n\n');
fprintf('=== All figures complete. ===\n');
