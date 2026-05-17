%  ISAC System: Dinkelbach-SCA Algorithm for Multi-Objective EE Optimization
%  STATISTICAL uncertainty model  –  Figures 1 to 8
%
%  Channel uncertainty model:
%    Delta_h_k = xi_k * p_k,  p_k ~ CN(0, I_N)
%    xi_k = delta * ||h_hat_k||
%
%  Benchmarks per figure:
%    (1) Robust          : xi_k = p.xi_k_vec,  i_b = p.i_b,   i_k = p.i_k
%    (2) Perfect CSI     : xi_k = 0,           i_b = p.i_b,   i_k = p.i_k
%    (3) Perfect Hardware: xi_k = p.xi_k_vec,  i_b = IMP_EPS, i_k = IMP_EPS
%    (4) Perfect Both    : xi_k = 0,           i_b = IMP_EPS, i_k = IMP_EPS
%  All four curves share ONE set of normalization bounds per figure.
clc; clear; close all;

%% =========================================================================
%  LOCK RANDOM SEED
%rng(42);

%% =========================================================================
%  LOAD SYSTEM PARAMETERS
p = parameters();

%  Small but nonzero impairment floor to keep solver stable when i_b=i_k~0
IMP_EPS = 0;

%% =========================================================================
%  GENERATE BASELINE CHANNEL  H_hat  (N x K)
%{
H_hat = generate_channels(p.N, p.K, 1);
save('saved_channels.mat', 'H_hat', 'p');
fprintf('Channel realisation saved to saved_channels.mat\n\n');
%}
load('saved_channels.mat');
%% =========================================================================
%  COMPUTE PER-USER UNCERTAINTY STANDARD DEVIATION
%  xi_k = delta * ||h_hat_k||   (statistical uncertainty radius)
% --------------------------------------------------------------------------
for k = 1:p.K
    p.xi_k_vec(k) = p.delta * norm(H_hat(:, k));
end
fprintf('=== ISAC Dinkelbach-SCA  (Statistical CSI) ===\n');
fprintf('N=%d  K=%d  M=%d  P_max=30 dBm\n', p.N, p.K, p.M);
fprintf('delta=%.3f  delta_outage=%.3f\n\n', p.delta, p.delta_outage);

%  Shorthand benchmark parameter sets
xi_k_rob  = p.xi_k_vec;
xi_k_perf = zeros(1, p.K);

%  4 variant definitions reused across figures
v_xik = {xi_k_rob,  xi_k_perf, xi_k_rob,  xi_k_perf};
v_ib  = {p.i_b,     p.i_b,     IMP_EPS,   IMP_EPS  };
v_ik  = {p.i_k,     p.i_k,     IMP_EPS,   IMP_EPS  };
v_lbl = {'Robust','Perfect CSI','Perfect HW','Perfect Both'};

%  Shared plot styles
styles    = {'r-s','k--^','b-o','m--d'};
styles_ec = {'r-s','k-^','b-o','m-d'};
styles_es = {'r--s','k--^','b--o','m--d'};
mfc       = {'r',  'k',   'b',  'm'  };

%  Legend label helpers
ec_lbl = cellfun(@(x) ['EE_c ' x], v_lbl, 'UniformOutput', false);
es_lbl = cellfun(@(x) ['EE_s ' x], v_lbl, 'UniformOutput', false);

%% =========================================================================
%  GENERATE H_bank FOR FIGURE 7
N_vec_f7 = 6:1:9;
H_bank   = generate_channels(max(N_vec_f7), p.K, p.N_MC);
fprintf('H_bank generated: %d x %d x %d.\n\n', max(N_vec_f7), p.K, p.N_MC);

%% =========================================================================
%  PHASE 1 – GLOBAL NORMALISATION CONSTANTS  (baseline parameters)
% --------------------------------------------------------------------------
fprintf('=== Phase 1: Global normalisation constants ===\n');

[EEcmax, EEcmin, EEsmax, EEsmin] = get_norm_constants( ...
    H_hat, p.theta_targets, p.N, p.K, p.M, ...
    p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
    p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, ...
    p.T_max, p.epsilon, p.N_rand);

fprintf('\nEE_c,max = %.4f  |  EE_c,min = %.4f\n', EEcmax, EEcmin);
fprintf('EE_s,max = %.4f  |  EE_s,min = %.4f\n\n', EEsmax, EEsmin);

%% =========================================================================
%  FIGURE 2 – IMPACT OF WEIGHTING COEFFICIENT omega  (run FIRST)
%  Finds crossover omega* where omega*EEc_norm = (1-omega)*EEs_norm
% --------------------------------------------------------------------------
fprintf('=== Figure 2: Impact of omega ===\n');

omega_vec = 0:0.05:1;
n_om      = numel(omega_vec);

wEEc_om = zeros(n_om, 1);   %  omega     * (EEc - EEcmin)/(EEcmax - EEcmin)
wEEs_om = zeros(n_om, 1);   %  (1-omega) * (EEs - EEsmin)/(EEsmax - EEsmin)

for idx = 1:n_om
    om = omega_vec(idx);
    fprintf('  omega=%.3f\n', om);
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
        p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
        p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, om, ...
        p.T_max, p.epsilon, EEcmax, EEsmax, p.N_rand, EEcmin, EEsmin);

    ec_norm = (ec - EEcmin) / (EEcmax - EEcmin);
    es_norm = (es - EEsmin) / (EEsmax - EEsmin);

    wEEc_om(idx) = om       * ec_norm;
    wEEs_om(idx) = (1-om)   * es_norm;
     % Break if either weighted component leaves [0,1]
    if wEEc_om(idx) > 1 || wEEc_om(idx) < 0 || ...
       wEEs_om(idx) > 1 || wEEs_om(idx) < 0
        fprintf('  WARNING: omega=%.2f caused out-of-range value. Stopping sweep.\n', om);
        wEEc_om(idx+1:end) = NaN;
        wEEs_om(idx+1:end) = NaN;
        break;
    end
end

% Find crossover: omega*EEc_norm = (1-omega)*EEs_norm
diff_vec  = wEEc_om - wEEs_om;
cross_idx = find(diff(sign(diff_vec)) ~= 0, 1);

if ~isempty(cross_idx)
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

% Update omega to crossover value for all subsequent figures
p.omega = omega_cross;

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
%  Plots: total WEE, omega*EEc_norm, (1-omega)*EEs_norm per iteration
%  Uses omega* found in Figure 2
% --------------------------------------------------------------------------
fprintf('=== Figure 1: Convergence ===\n');

[Wc_f1, Ws_f1, EEc_f1, EEs_f1, obj_hist, EEc_hist, EEs_hist] = ...
    solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
               p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, ...
               p.delta_outage, p.xi_k_vec, ...
               p.omega, p.T_max, p.epsilon, EEcmax, EEsmax, p.N_rand, EEcmin, EEsmin);

figure(1);
plot( obj_hist,   'b-o',  'LineWidth',2, 'MarkerSize',6); hold on;
plot(EEc_hist,   'r--s', 'LineWidth',1.8,'MarkerSize',5); hold on;
plot(EEs_hist,   'g--^', 'LineWidth',1.8,'MarkerSize',5);
hold off; grid on;
xlabel('Number of Iterations','FontSize',13);
ylabel('Normalised Value','FontSize',13);
title(sprintf('Fig. 1: Convergence of Dinkelbach-SCA  (\\omega=%.3f)', p.omega),'FontSize',13);
legend('WEE', '\omega\cdotEE_c^{norm}', '(1-\omega)\cdotEE_s^{norm}', ...
       'Location','best','FontSize',11);
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

% --- Pass 1: shared bounds from ALL 4 variants x ALL sweep points ---
ec_max_all = zeros(n_gam,4);  ec_min_all = zeros(n_gam,4);
es_max_all = zeros(n_gam,4);  es_min_all = zeros(n_gam,4);

for v = 1:4
    for idx = 1:n_gam
        [ec_max_all(idx,v), ec_min_all(idx,v), ...
         es_max_all(idx,v), es_min_all(idx,v)] = ...
            get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
                p.P_max, p.sigma2, Gamma_W_vec(idx), p.gamma_min, ...
                v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, v_xik{v}, ...
                p.T_max, p.epsilon, p.N_rand);
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
            v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, v_xik{v}, p.omega, ...
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
title(sprintf('(a) WEE  (\\omega=%.3f)', p.omega),'FontSize',13);
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

sgtitle(sprintf('Fig. 3: Impact of \\Gamma_{min}  (\\omega=%.3f)', p.omega), ...
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
                v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, v_xik{v}, ...
                p.T_max, p.epsilon, p.N_rand);
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
            v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, v_xik{v}, p.omega, ...
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
title(sprintf('(a) WEE  (\\omega=%.3f)', p.omega),'FontSize',13);
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

sgtitle(sprintf('Fig. 4: Impact of \\gamma_{min}  (\\omega=%.3f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 4 rendered.\n\n');

%% =========================================================================
%  FIGURE 5 – IMPACT OF HARDWARE IMPAIRMENT COEFFICIENTS (i_b, i_k)
%  4 benchmarks (Robust, Perfect CSI, Perfect HW, Perfect Both)
%  For "Perfect HW" variants: i_b/i_k are fixed at IMP_EPS (flat line).
%  For "Robust" / "Perfect CSI" variants: i_b/i_k are swept.
% --------------------------------------------------------------------------
fprintf('=== Figure 5: Hardware impairments sweep (4 benchmarks) ===\n');

i_b_vec    = [0.01, 0.02, 0.03, 0.04];
i_k_vec    = [0.02, 0.04, 0.06, 0.08];
n_imp      = numel(i_b_vec);
imp_labels = arrayfun(@(a,b) sprintf('(%.2f,%.2f)',a,b), ...
                      i_b_vec, i_k_vec, 'UniformOutput', false);

% For each variant, build effective i_b / i_k used at each sweep point:
%   v=1 Robust       : sweep i_b_vec(idx), i_k_vec(idx);  xi_k = xi_k_rob
%   v=2 Perfect CSI  : sweep i_b_vec(idx), i_k_vec(idx);  xi_k = 0
%   v=3 Perfect HW   : fixed IMP_EPS, IMP_EPS;             xi_k = xi_k_rob
%   v=4 Perfect Both : fixed IMP_EPS, IMP_EPS;             xi_k = 0
ib_use = @(v,idx) deal( ...
    (v<=2)*i_b_vec(idx) + (v>2)*IMP_EPS, ...
    (v<=2)*i_k_vec(idx) + (v>2)*IMP_EPS );

% --- Pass 1: shared bounds across all 4 variants x all sweep points ---
ec_max_all = zeros(n_imp,4);  ec_min_all = zeros(n_imp,4);
es_max_all = zeros(n_imp,4);  es_min_all = zeros(n_imp,4);

for v = 1:4
    for idx = 1:n_imp
        [ib_cur, ik_cur] = ib_use(v, idx);
        fprintf('  [Norm] %-16s  i_b=%.3f  i_k=%.3f\n', v_lbl{v}, ib_cur, ik_cur);
        [ec_max_all(idx,v), ec_min_all(idx,v), ...
         es_max_all(idx,v), es_min_all(idx,v)] = ...
            get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
                p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
                ib_cur, ik_cur, p.P_static, p.delta_outage, v_xik{v}, ...
                p.T_max, p.epsilon, p.N_rand);
    end
end

[EEc_max_f5, best_idx_f5] = max(ec_max_all(:,1));
EEc_min_f5 = ec_min_all(best_idx_f5, 1);
EEs_max_f5 = es_max_all(best_idx_f5, 1);
EEs_min_f5 = es_min_all(best_idx_f5, 1);
fprintf('  Fig5 shared EEc bounds: [%.4f, %.4f]\n', EEc_min_f5, EEc_max_f5);
fprintf('  Fig5 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f5, EEs_max_f5);

% --- Pass 2: solve all 4 variants ---
WEE_imp = zeros(n_imp, 4);
EEc_imp = zeros(n_imp, 4);
EEs_imp = zeros(n_imp, 4);

for v = 1:4
    for idx = 1:n_imp
        [ib_cur, ik_cur] = ib_use(v, idx);
        fprintf('  [Solve] %-16s  i_b=%.3f  i_k=%.3f\n', v_lbl{v}, ib_cur, ik_cur);
        [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
            ib_cur, ik_cur, p.P_static, p.delta_outage, v_xik{v}, p.omega, ...
            p.T_max, p.epsilon, EEc_max_f5, EEs_max_f5, p.N_rand, EEc_min_f5, EEs_min_f5);
        EEc_imp(idx,v) = ec;  EEs_imp(idx,v) = es;
        WEE_imp(idx,v) = p.omega     * ((ec - EEc_min_f5)/(EEc_max_f5 - EEc_min_f5)) + ...
                         (1-p.omega) * ((es - EEs_min_f5)/(EEs_max_f5 - EEs_min_f5));
        fprintf('    EEc=%.4f  EEs=%.4f  WEE=%.4f\n', ec, es, WEE_imp(idx,v));
    end
end

x_ticks = 1:n_imp;
figure(5); set(gcf,'Position',[100 100 700 480]);

for v = 1:4
    plot(x_ticks, WEE_imp(:,v), styles{v}, 'LineWidth',2, ...
        'MarkerSize',7, 'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
set(gca,'XTick',x_ticks,'XTickLabel',imp_labels,'FontSize',12);
xlabel('(i_b,  i_k)  pairs','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('Fig. 5: Impact of Hardware Impairments  (\\omega=%.3f)', p.omega), ...
      'FontSize',13,'FontWeight','bold');
legend(v_lbl,'Location','best','FontSize',10);
xtickangle(15);
drawnow;
fprintf('    Figure 5 rendered.\n\n');

%% =========================================================================
%  FIGURE 6 – WEIGHTED EE vs P_max  (4 benchmarks)
%  Layout: 3 subplots – (a) WEE, (b) EEc, (c) EEs
% --------------------------------------------------------------------------
fprintf('=== Figure 6: P_max sweep (4 benchmarks) ===\n');

Pmax_dBm_vec = 5:5:40;
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
                v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, v_xik{v}, ...
                p.T_max, p.epsilon, p.N_rand);
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
            v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, v_xik{v}, p.omega, ...
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
title(sprintf('(a) WEE  (\\omega=%.3f)', p.omega),'FontSize',13);
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

sgtitle(sprintf('Fig. 6: Weighted EE vs. P_{max}  (\\omega=%.3f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 6 rendered.\n\n');

%% =========================================================================
%  FIGURE 7 – WEIGHTED EE vs NUMBER OF ANTENNAS N  (4 benchmarks)
%  Single subplot (WEE only).
%  Note: xi_k is N-dependent and recomputed locally per N.
%        Variants 1 & 3 use statistical xi_k; variants 2 & 4 use xi_k = 0.
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
            xi_k_cur = zeros(1, p.K);
            for k = 1:p.K
                xi_k_cur(k) = p.delta * norm(H_norm(:,k));
            end
        else
            xi_k_cur = zeros(1, p.K);
        end

        [ec_max_all(idx,v), ec_min_all(idx,v), ...
         es_max_all(idx,v), es_min_all(idx,v)] = ...
            get_norm_constants(H_norm, p.theta_targets, N_cur, p.K, p.M, ...
                p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
                v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, xi_k_cur, ...
                p.T_max, p.epsilon, p.N_rand);
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
EEc_f7 = zeros(n_N,4);  
EEs_f7 = zeros(n_N,4);  

for v = 1:4
    for idx = 1:n_N
        N_cur = N_vec_f7(idx);
        H_cur = H_bank(1:N_cur, :, 1);

        if ismember(v, [1 3])
            xi_k_cur = zeros(1, p.K);
            for k = 1:p.K
                xi_k_cur(k) = p.delta * norm(H_cur(:,k));
            end
        else
            xi_k_cur = zeros(1, p.K);
        end

        fprintf('  [Solve] %-16s  N=%d\n', v_lbl{v}, N_cur);
        [~,~,ec,es,~] = solve_ISAC(H_cur, p.theta_targets, N_cur, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
            v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, xi_k_cur, p.omega, ...
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
title(sprintf('Fig. 7: Weighted EE vs. N  (\\omega=%.3f)', p.omega), ...
      'FontSize',14,'FontWeight','bold');
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);
drawnow;
fprintf('    Figure 7 rendered.\n\n');

%% =========================================================================
%  FIGURE 8 – IMPACT OF CHANNEL UNCERTAINTY LEVEL delta  (4 benchmarks)
%  Statistical CSI: xi_k is swept via delta; "Perfect CSI" uses xi_k=0.
%  v=1 Robust       : xi_k swept via delta;  i_b=p.i_b,   i_k=p.i_k
%  v=2 Perfect CSI  : xi_k = 0 (flat);       i_b=p.i_b,   i_k=p.i_k
%  v=3 Perfect HW   : xi_k swept via delta;  i_b=IMP_EPS, i_k=IMP_EPS
%  v=4 Perfect Both : xi_k = 0 (flat);       i_b=IMP_EPS, i_k=IMP_EPS
% --------------------------------------------------------------------------
fprintf('=== Figure 8: Channel uncertainty level delta (4 benchmarks) ===\n');

delta_vec = [0, 0.005, 0.01, 0.02, 0.03];
n_delta   = numel(delta_vec);

% Helper: compute xi_k_vec for a given delta
xik_from_delta = @(d) arrayfun(@(k) d * norm(H_hat(:,k)), 1:p.K);

% For variants 2 & 4, xi_k = 0 regardless of delta
xik_for_v = @(v, d) (ismember(v,[1 3])) * xik_from_delta(d);

% --- Pass 1: shared bounds ---
ec_max_all = zeros(n_delta,4);  ec_min_all = zeros(n_delta,4);
es_max_all = zeros(n_delta,4);  es_min_all = zeros(n_delta,4);

for v = 1:4
    for idx = 1:n_delta
        xi_k_cur = xik_for_v(v, delta_vec(idx));
        fprintf('  [Norm] %-16s  delta=%.4f\n', v_lbl{v}, delta_vec(idx));
        [ec_max_all(idx,v), ec_min_all(idx,v), ...
         es_max_all(idx,v), es_min_all(idx,v)] = ...
            get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
                p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
                v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, xi_k_cur, ...
                p.T_max, p.epsilon, p.N_rand);
    end
end

[EEc_max_f8, best_idx_f8] = max(ec_max_all(:,1));
EEc_min_f8 = ec_min_all(best_idx_f8, 1);
EEs_max_f8 = es_max_all(best_idx_f8, 1);
EEs_min_f8 = es_min_all(best_idx_f8, 1);
fprintf('  Fig8 shared EEc bounds: [%.4f, %.4f]\n', EEc_min_f8, EEc_max_f8);
fprintf('  Fig8 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f8, EEs_max_f8);

% --- Pass 2: solve ---
WEE_delta = zeros(n_delta, 4);
EEc_delta = zeros(n_delta, 4);
EEs_delta = zeros(n_delta, 4);

for v = 1:4
    for idx = 1:n_delta
        xi_k_cur = xik_for_v(v, delta_vec(idx));
        fprintf('  [Solve] %-16s  delta=%.4f\n', v_lbl{v}, delta_vec(idx));
        [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
            p.Gamma_min, p.gamma_min, v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, xi_k_cur, ...
            p.omega, p.T_max, p.epsilon, EEc_max_f8, EEs_max_f8, p.N_rand, EEc_min_f8, EEs_min_f8);
        EEc_delta(idx,v) = ec;
        EEs_delta(idx,v) = es;
        WEE_delta(idx,v) = p.omega     * ((ec - EEc_min_f8)/(EEc_max_f8 - EEc_min_f8)) + ...
                           (1-p.omega) * ((es - EEs_min_f8)/(EEs_max_f8 - EEs_min_f8));
        fprintf('    EEc=%.4f  EEs=%.4f  WEE=%.4f\n', ec, es, WEE_delta(idx,v));
    end
end

figure(8); set(gcf,'Position',[100 100 700 480]);

for v = 1:4
    plot(delta_vec, WEE_delta(:,v), styles{v}, 'LineWidth',2, ...
        'MarkerSize',7, 'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Channel Uncertainty Level  \delta','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('Fig. 8: Impact of \\delta  (\\omega=%.3f)', p.omega), ...
      'FontSize',13,'FontWeight','bold');
legend(v_lbl,'Location','best','FontSize',10);
set(gca,'FontSize',12);
drawnow;
fprintf('    Figure 8 rendered.\n\n');
%% =========================================================================
%  FIGURE 9 – TRANSMIT BEAMPATTERN GAIN
% --------------------------------------------------------------------------
fprintf('=== Figure 9: Beampattern ===\n');

theta_scan = -90 : 0.5 : 90;
n_theta    = numel(theta_scan);
sl_bp      = @(th) exp(1j*pi*((0:p.N-1)'-(p.N-1)/2)*sind(th));

% --- Solve once per variant at baseline parameters ----------------------
Wc_bp = cell(1,4);
Ws_bp = cell(1,4);

for v = 1:4
    fprintf('  [Beampattern] variant: %s\n', v_lbl{v});
    [Wc_bp{v}, Ws_bp{v}, ~, ~, ~] = ...
        solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
                   p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
                   v_ib{v}, v_ik{v}, p.P_static, p.delta_outage, v_xik{v}, ...  
                   p.omega, p.T_max, p.epsilon, ...
                   EEcmax, EEsmax, p.N_rand, EEcmin, EEsmin);

end

% --- Compute beampattern gain -------------------------------------------
BP_dB = zeros(n_theta, 4);

for v = 1:4
    W_tot = zeros(p.N, p.N);
    for j = 1:p.K; W_tot = W_tot + Wc_bp{v}{j}; end
    for l = 1:p.M; W_tot = W_tot + Ws_bp{v}{l};  end
    tr_W  = real(trace(W_tot));

    BP_lin = zeros(n_theta, 1);
    for idx = 1:n_theta
        vm = sl_bp(theta_scan(idx));
        BP_lin(idx) = real(vm' * W_tot * vm) + v_ib{v} * tr_W;
    end

    BP_lin     = max(BP_lin, 1e-30);
    BP_dB(:,v) = 10*log10(BP_lin / max(BP_lin));
end

% --- Plot ---------------------------------------------------------------
figure(9); set(gcf,'Position',[100 100 750 500]);

for v = 1:4
    plot(theta_scan, BP_dB(:,v), styles{v}, ...
         'LineWidth', 2, 'MarkerSize', 4, ...
         'MarkerIndices', 1:20:n_theta, ...
         'MarkerFaceColor', mfc{v});
    hold on;
end

% Target angle markers
for m = 1:p.M
    xline(p.theta_targets(m), 'k:', 'LineWidth', 1.2, 'HandleVisibility','off');
end
hold off; grid on;

xlim([-90 90]);
ylim([max(min(BP_dB(:)), -40)  max(BP_dB(:))+3]);
set(gca, 'XTick', -90:15:90, 'FontSize', 12);
xlabel('Angle  \theta  (degrees)',        'FontSize', 13);
ylabel('Normalised Beampattern Gain (dB)','FontSize', 13);
title(sprintf('Fig. 9: Transmit Beampattern  (\\omega=%.3f)', p.omega), ...
      'FontSize', 13, 'FontWeight', 'bold');
legend(v_lbl, 'Location','best', 'FontSize', 10);
drawnow;
fprintf('    Figure 9 rendered.\n\n');
%% =========================================================================
    %  SAVE SIMULATION DATA
    % --------------------------------------------------------------------------
    fprintf('=== Saving Results ===\n');

    % 1. Mirror the directory logic from save_figure.m
    timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
    [parentDir, currentFolder, ~] = fileparts(pwd);
    baseDir = fullfile(parentDir, 'Results', currentFolder);
    saveDir = fullfile(baseDir, timestamp);

    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    % 2. Define the targeted save path
    matFilePath = fullfile(saveDir, 'isac_results.mat');
    fprintf('Saving data to: %s\n', matFilePath);

    % 3. Targeted save (avoids crashing from parallel pool objects)
% 3. Targeted save (avoids crashing from parallel pool objects)
    save(matFilePath, 'p', 'H_hat', 'Gamma_dBm_vec', 'gamma_dB_vec', ...
        'i_b_vec', 'i_k_vec', 'n_imp', 'imp_labels', 'Pmax_dBm_vec', 'N_vec_f7', ...
        'delta_vec', 'omega_vec', 'wEEc_om', 'wEEs_om', 'omega_cross', 'wEEc_cross', ... 
        'obj_hist', 'EEc_hist', 'EEs_hist', 'Wc_f1', 'Ws_f1', ...
        '-regexp', '^WEE_f', '^EEc_f', '^EEs_f', '^WEE_delta', '^EEc_delta', '^EEs_delta');

    fprintf('Data successfully saved.\n');
