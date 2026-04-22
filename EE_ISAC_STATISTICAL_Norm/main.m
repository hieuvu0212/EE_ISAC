%  ISAC System: Dinkelbach-SCA Algorithm for Multi-Objective EE Optimization
%  STATISTICAL uncertainty model  –  Figures 1 to 8
%
%  Channel uncertainty model:
%    Delta_h_k = xi_k * p_k,  p_k ~ CN(0, I_N)
%    xi_k = delta * ||h_hat_k||
%
%  Hardware uncertainty model:
%    i_b = transmitter distortion coefficient
%    i_k = receiver  distortion coefficient
%
%  Benchmarks per figure:
%    (1) Robust          : xi_k = p.xi_k_vec,  i_b = p.i_b,   i_k = p.i_k
%    (2) Perfect CSI     : xi_k = 0,            i_b = p.i_b,   i_k = p.i_k
%    (3) Perfect Hardware: xi_k = p.xi_k_vec,  i_b = IMP_EPS, i_k = IMP_EPS
%    (4) Perfect Both    : xi_k = 0,            i_b = IMP_EPS, i_k = IMP_EPS
%  All four curves share ONE set of normalization bounds anchored to Robust.
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
xi_k_rob  = p.xi_k_vec;     % statistical channel uncertainty
xi_k_perf = zeros(1, p.K);  % perfect channel knowledge: xi_k = 0

%  4 variant definitions reused across Fig 3,4,6,7:
%    v_xik{v} : channel uncertainty vector
%    v_ib{v}  : transmitter hardware distortion
%    v_ik{v}  : receiver  hardware distortion
v_xik = {xi_k_rob,  xi_k_perf, xi_k_rob,  xi_k_perf};
v_ib  = {p.i_b,     p.i_b,     IMP_EPS,   IMP_EPS  };
v_ik  = {p.i_k,     p.i_k,     IMP_EPS,   IMP_EPS  };
v_lbl = {'Robust','Perfect CSI','Perfect HW','Perfect Both'};

%  Shared plot styles (1=Robust, 2=PerfCSI, 3=PerfHW, 4=PerfBoth)
styles = {'r-s','k--^','b-o','m--d'};
mfc    = {'r',  'k',   'b',  'm'  };

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
    p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, p.T_max, p.epsilon, p.N_rand);

fprintf('\nEE_c,max = %.4f  |  EE_c,min = %.4f\n', EEcmax, EEcmin);
fprintf('EE_s,max = %.4f  |  EE_s,min = %.4f\n\n', EEsmax, EEsmin);

%% =========================================================================
%  FIGURE 1 – CONVERGENCE OF DINKELBACH-SCA
% --------------------------------------------------------------------------
fprintf('=== Figure 1: Convergence ===\n');

[Wc_f1, Ws_f1, EEc_f1, EEs_f1, obj_hist] = ...
    solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
               p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, ...
               p.omega, p.T_max, p.epsilon, EEcmax, EEsmax, p.N_rand, EEcmin, EEsmin);

figure(1);
plot(1:numel(obj_hist), obj_hist, 'b-o','LineWidth',2,'MarkerSize',6);
grid on;
xlabel('Number of Iterations','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('Fig. 1: Convergence of Dinkelbach-SCA  (\\omega=%.2f)', p.omega),'FontSize',13);
str = sprintf('EE_c = %.4f bit/J/Hz\nEE_s = %.4f', EEc_f1, EEs_f1);
<<<<<<< HEAD
text(0.97, 0.15, str, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 11, ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'Margin', 5);
legend(sprintf('\\omega = %.2f', p.omega), 'Location', 'southeast');
set(gca,'FontSize',12);
drawnow;
=======
text(0.97, 0.15, str, 'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontSize',11,'BackgroundColor','white', ...
    'EdgeColor','black','Margin',5);
legend(sprintf('\\omega = %.2f', p.omega), 'Location','southeast');
set(gca,'FontSize',12); drawnow;
>>>>>>> Newworld
fprintf('    Figure 1 rendered.\n\n');

%% =========================================================================
%  FIGURE 2 – IMPACT OF WEIGHTING COEFFICIENT omega
% --------------------------------------------------------------------------
fprintf('=== Figure 2: Impact of omega ===\n');

omega_vec = 0:0.05:1;
n_om      = numel(omega_vec);
EEc_om    = zeros(n_om,1);
EEs_om    = zeros(n_om,1);

for idx = 1:n_om
    om = omega_vec(idx);
    fprintf('  omega=%.2f\n', om);
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
        p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, om, ...
        p.T_max, p.epsilon, EEcmax, EEsmax, p.N_rand, EEcmin, EEsmin);
    EEc_om(idx) = ec;
    EEs_om(idx) = es;
end

figure(2);
yyaxis left;
plot(omega_vec, EEc_om, 'b-o','LineWidth',2,'MarkerSize',5,'MarkerFaceColor','b');
ylabel('Communication EE (EE_c)','FontSize',13,'Color','b');
set(gca,'YColor','b');
yyaxis right;
plot(omega_vec, EEs_om, 'r-s','LineWidth',2,'MarkerSize',5,'MarkerFaceColor','r');
ylabel('Sensing EE (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r');
grid on;
xlabel('Weighting Coefficient  \omega','FontSize',13);
<<<<<<< HEAD
set(gca, 'XTick', omega_vec, ...
         'XTickLabelRotation', 45, ...
         'FontSize', 10);
=======
set(gca,'XTick',omega_vec,'XTickLabelRotation',45,'FontSize',10);
>>>>>>> Newworld
title('Fig. 2: Trade-off between EE_c and EE_s','FontSize',13);
set(gca,'FontSize',12); drawnow;
fprintf('    Figure 2 rendered.\n\n');

%% =========================================================================
%  FIGURE 3 – IMPACT OF SENSING CONSTRAINT Gamma_min  (4 benchmarks)
% --------------------------------------------------------------------------
fprintf('=== Figure 3: Gamma_min sweep (4 benchmarks) ===\n');

Gamma_dBm_vec = 16:2:24;
Gamma_W_vec   = db2pow(Gamma_dBm_vec) * 1e-3;
n_gam         = numel(Gamma_dBm_vec);

<<<<<<< HEAD
% Khởi tạo biến cho Robust (Statistical)
ec_max_rob = zeros(n_gam,1); ec_min_rob = zeros(n_gam,1);
es_max_rob = zeros(n_gam,1); es_min_rob = zeros(n_gam,1);
% Khởi tạo biến cho Perfect
ec_max_perf = zeros(n_gam,1); ec_min_perf = zeros(n_gam,1);
es_max_perf = zeros(n_gam,1); es_min_perf = zeros(n_gam,1);

xi_k_perf = zeros(1, p.K); % Vector 0 ép thành Perfect CSI

for idx = 1:n_gam
    Gam_cur = Gamma_W_vec(idx);
    fprintf('  [Norm pass] Gamma_min=%d dBm\n', Gamma_dBm_vec(idx));
    
    % Chuẩn hóa Robust CSI
    [ec_max_rob(idx), ec_min_rob(idx), es_max_rob(idx), es_min_rob(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, Gam_cur, p.gamma_min, ...
            p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, p.T_max, p.epsilon, p.N_rand);
    
    % Chuẩn hóa Perfect CSI
    [ec_max_perf(idx), ec_min_perf(idx), es_max_perf(idx), es_min_perf(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, Gam_cur, p.gamma_min, ...
            p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_perf, p.T_max, p.epsilon, p.N_rand);
end

[EEc_max_f3_rob, best_idx_ec3_rob] = max(ec_max_rob); EEc_min_f3_rob = ec_min_rob(best_idx_ec3_rob);
EEs_max_f3_rob = es_max_rob(best_idx_ec3_rob);        EEs_min_f3_rob = es_min_rob(best_idx_ec3_rob);

[EEc_max_f3_perf, best_idx_ec3_perf] = max(ec_max_perf); EEc_min_f3_perf = ec_min_perf(best_idx_ec3_perf);
EEs_max_f3_perf = es_max_perf(best_idx_ec3_perf);        EEs_min_f3_perf = es_min_perf(best_idx_ec3_perf);

% Giải bài toán 
WEE_rob = zeros(n_gam,1); EEc_rob = zeros(n_gam,1); EEs_rob = zeros(n_gam,1);
WEE_perf = zeros(n_gam,1); EEc_perf = zeros(n_gam,1); EEs_perf = zeros(n_gam,1);

for idx = 1:n_gam
    Gam_cur = Gamma_W_vec(idx);
    fprintf('  [Solve pass] Gamma_min=%d dBm\n', Gamma_dBm_vec(idx));

    % Robust
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
        Gam_cur, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, p.omega, ...
        p.T_max, p.epsilon, EEc_max_f3_rob, EEs_max_f3_rob, p.N_rand, EEc_min_f3_rob, EEs_min_f3_rob);
    EEc_rob(idx) = ec; EEs_rob(idx) = es;
    WEE_rob(idx) = p.omega*((ec - EEc_min_f3_rob)/(EEc_max_f3_rob - EEc_min_f3_rob)) + ...
                   (1-p.omega)*((es - EEs_min_f3_rob)/(EEs_max_f3_rob - EEs_min_f3_rob));

    % Perfect
    [~,~,ec_p,es_p,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
        Gam_cur, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_perf, p.omega, ...
        p.T_max, p.epsilon, EEc_max_f3_perf, EEs_max_f3_perf, p.N_rand, EEc_min_f3_perf, EEs_min_f3_perf);
    EEc_perf(idx) = ec_p; EEs_perf(idx) = es_p;
    WEE_perf(idx) = p.omega*((ec_p - EEc_min_f3_perf)/(EEc_max_f3_perf - EEc_min_f3_perf)) + ...
                    (1-p.omega)*((es_p - EEs_min_f3_perf)/(EEs_max_f3_perf - EEs_min_f3_perf));
end

figure(3); set(gcf,'Position',[100 100 1200 500]);
subplot(1,2,1);
plot(Gamma_dBm_vec, WEE_perf, 'k--^','LineWidth',2,'MarkerSize',8); hold on;
plot(Gamma_dBm_vec, WEE_rob, 'r-s','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r');
grid on; xlabel('Sensing Threshold \Gamma_{min} (dBm)'); ylabel('Normalised WEE');
title(sprintf('(a) WEE (\\omega=%.2f)', p.omega)); legend('Perfect CSI', 'Statistical CSI');

subplot(1,2,2);
yyaxis left;
plot(Gamma_dBm_vec, EEc_perf, 'b--^','LineWidth',2,'MarkerSize',7); hold on;
plot(Gamma_dBm_vec, EEc_rob, 'b-o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','b');
ylabel('Communication EE (EE_c)'); set(gca,'YColor','b');
=======
% --- Pass 1: collect bounds from ALL 4 variants x ALL sweep points ---
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

% Bounds anchored to best Robust scenario (column 1)
[EEc_max_f3, best_idx_f3] = max(ec_max_all(:,1));
EEc_min_f3 = ec_min_all(best_idx_f3, 1);
EEs_max_f3 = es_max_all(best_idx_f3, 1);
EEs_min_f3 = es_min_all(best_idx_f3, 1);
fprintf('  Fig3 shared EEc bounds: [%.4f, %.4f]  (best Gamma=%d dBm)\n', ...
        EEc_min_f3, EEc_max_f3, Gamma_dBm_vec(best_idx_f3));
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

% --- Plot ---
figure(3); set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
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

subplot(1,2,2);
yyaxis left;
for v = 1:4
    plot(Gamma_dBm_vec, EEc_f3(:,v), styles{v}, 'LineWidth',1.8, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off;
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13,'Color','b');
set(gca,'YColor','b');
yyaxis right;
for v = 1:4
    plot(Gamma_dBm_vec, EEs_f3(:,v), styles{v}, 'LineWidth',1.8, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off;
ylabel('Sensing EE  (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r'); grid on;
xlabel('Sensing Threshold  \Gamma_{min}  (dBm)','FontSize',13);
title('(b) Raw EE_c and EE_s','FontSize',13);
legend([ec_lbl, es_lbl],'Location','best','FontSize',9);
set(gca,'FontSize',12);
>>>>>>> Newworld

yyaxis right;
plot(Gamma_dBm_vec, EEs_perf, 'r--^','LineWidth',2,'MarkerSize',7); hold on;
plot(Gamma_dBm_vec, EEs_rob, 'r-s','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','r');
ylabel('Sensing EE (EE_s)'); set(gca,'YColor','r');
grid on; xlabel('Sensing Threshold \Gamma_{min} (dBm)'); title('(b) Raw EE_c & EE_s');
legend('EE_c (Perf)','EE_c (Stat)','EE_s (Perf)','EE_s (Stat)','Location','best');
sgtitle(sprintf('Fig. 3: Impact of \\Gamma_{min} (\\omega=%.2f)', p.omega), 'FontWeight','bold');
drawnow;

%% =========================================================================
<<<<<<< HEAD
%%  FIGURE 4  –  Impact of gamma_min (SINR)
fprintf('=== Figure 4: Computing global norms across all gamma_min values ===\n');
gamma_dB_vec  = 3 : 1 : 8;
gamma_lin_vec = db2pow(gamma_dB_vec);
n_sinr        = numel(gamma_dB_vec);

ec_max_rob = zeros(n_sinr,1); ec_min_rob = zeros(n_sinr,1);
es_max_rob = zeros(n_sinr,1); es_min_rob = zeros(n_sinr,1);
ec_max_perf = zeros(n_sinr,1); ec_min_perf = zeros(n_sinr,1);
es_max_perf = zeros(n_sinr,1); es_min_perf = zeros(n_sinr,1);

for idx = 1:n_sinr
    gam_cur = gamma_lin_vec(idx);
    fprintf('  [Norm pass] gamma_min=%d dB\n', gamma_dB_vec(idx));
    
    [ec_max_rob(idx), ec_min_rob(idx), es_max_rob(idx), es_min_rob(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, p.Gamma_min, gam_cur, ...
        p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, p.T_max, p.epsilon, p.N_rand);
        
    [ec_max_perf(idx), ec_min_perf(idx), es_max_perf(idx), es_min_perf(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, p.Gamma_min, gam_cur, ...
        p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_perf, p.T_max, p.epsilon, p.N_rand);
end

[EEc_max_f4_rob, best_ec4_rob] = max(ec_max_rob); EEc_min_f4_rob = ec_min_rob(best_ec4_rob);
EEs_max_f4_rob = es_max_rob(best_ec4_rob);        EEs_min_f4_rob = es_min_rob(best_ec4_rob);

[EEc_max_f4_perf, best_ec4_perf] = max(ec_max_perf); EEc_min_f4_perf = ec_min_perf(best_ec4_perf);
EEs_max_f4_perf = es_max_perf(best_ec4_perf);        EEs_min_f4_perf = es_min_perf(best_ec4_perf);

WEE_rob = zeros(n_sinr,1); EEc_rob = zeros(n_sinr,1); EEs_rob = zeros(n_sinr,1);
WEE_perf = zeros(n_sinr,1); EEc_perf = zeros(n_sinr,1); EEs_perf = zeros(n_sinr,1);

for idx = 1:n_sinr
    gam_cur = gamma_lin_vec(idx);
    fprintf('  [Solve pass] gamma_min=%d dB\n', gamma_dB_vec(idx));

    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, p.Gamma_min, gam_cur, ...
        p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, p.omega, p.T_max, p.epsilon, EEc_max_f4_rob, EEs_max_f4_rob, p.N_rand, EEc_min_f4_rob, EEs_min_f4_rob);
    EEc_rob(idx) = ec; EEs_rob(idx) = es;
    WEE_rob(idx) = p.omega*((ec - EEc_min_f4_rob)/(EEc_max_f4_rob - EEc_min_f4_rob)) + (1-p.omega)*((es - EEs_min_f4_rob)/(EEs_max_f4_rob - EEs_min_f4_rob));

    [~,~,ec_p,es_p,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, p.Gamma_min, gam_cur, ...
        p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_perf, p.omega, p.T_max, p.epsilon, EEc_max_f4_perf, EEs_max_f4_perf, p.N_rand, EEc_min_f4_perf, EEs_min_f4_perf);
    EEc_perf(idx) = ec_p; EEs_perf(idx) = es_p;
    WEE_perf(idx) = p.omega*((ec_p - EEc_min_f4_perf)/(EEc_max_f4_perf - EEc_min_f4_perf)) + (1-p.omega)*((es_p - EEs_min_f4_perf)/(EEs_max_f4_perf - EEs_min_f4_perf));
end

figure(4); set(gcf,'Position',[100 100 1200 500]);
subplot(1,2,1);
plot(gamma_dB_vec, WEE_perf, 'k--^','LineWidth',2,'MarkerSize',8); hold on;
plot(gamma_dB_vec, WEE_rob, 'g-^','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','g');
grid on; xlabel('Minimum SINR \gamma_{min} (dB)'); ylabel('Normalised WEE');
title(sprintf('(a) WEE (\\omega=%.2f)', p.omega)); legend('Perfect CSI', 'Statistical CSI');

subplot(1,2,2);
yyaxis left;
plot(gamma_dB_vec, EEc_perf, 'b--^','LineWidth',2,'MarkerSize',7); hold on;
plot(gamma_dB_vec, EEc_rob, 'b-o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','b');
ylabel('Communication EE (EE_c)'); set(gca,'YColor','b');
=======
%  FIGURE 4 – IMPACT OF gamma_min (SINR)  (4 benchmarks)
% --------------------------------------------------------------------------
fprintf('=== Figure 4: gamma_min sweep (4 benchmarks) ===\n');

gamma_dB_vec  = 3:1:8;
gamma_lin_vec = db2pow(gamma_dB_vec);
n_sinr        = numel(gamma_dB_vec);

% --- Pass 1 ---
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
fprintf('  Fig4 shared EEc bounds: [%.4f, %.4f]  (best gamma=%d dB)\n', ...
        EEc_min_f4, EEc_max_f4, gamma_dB_vec(best_idx_f4));
fprintf('  Fig4 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f4, EEs_max_f4);

% --- Pass 2 ---
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

% --- Plot ---
figure(4); set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
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

subplot(1,2,2);
yyaxis left;
for v = 1:4
    plot(gamma_dB_vec, EEc_f4(:,v), styles{v}, 'LineWidth',1.8, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off;
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13,'Color','b');
set(gca,'YColor','b');
yyaxis right;
for v = 1:4
    plot(gamma_dB_vec, EEs_f4(:,v), styles{v}, 'LineWidth',1.8, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off;
ylabel('Sensing EE  (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r'); grid on;
xlabel('Minimum SINR  \gamma_{min}  (dB)','FontSize',13);
title('(b) Raw EE_c and EE_s','FontSize',13);
legend([ec_lbl, es_lbl],'Location','best','FontSize',9);
set(gca,'FontSize',12);
>>>>>>> Newworld

yyaxis right;
plot(gamma_dB_vec, EEs_perf, 'g--^','LineWidth',2,'MarkerSize',7); hold on;
plot(gamma_dB_vec, EEs_rob, 'g-^','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','g');
ylabel('Sensing EE (EE_s)'); set(gca,'YColor','g');
grid on; xlabel('Minimum SINR \gamma_{min} (dB)'); title('(b) Raw EE_c & EE_s');
legend('EE_c (Perf)','EE_c (Stat)','EE_s (Perf)','EE_s (Stat)','Location','best');
sgtitle(sprintf('Fig. 4: Impact of \\gamma_{min} (\\omega=%.2f)', p.omega), 'FontWeight','bold');
drawnow;

%% =========================================================================
<<<<<<< HEAD
%%  FIGURE 5  –  Impact of Hardware Impairment Coefficients (i_b, i_k)
fprintf('=== Figure 5: Impact of hardware impairments (i_b, i_k) ===\n');

i_b_vec = [ 0.01, 0.02,0.03];   % Transmitter distortion
i_k_vec = [0.02, 0.04,0.06];   % Receiver  distortion
n_imp   = numel(i_b_vec);

% --- Labels for x-axis ticks ---
=======
%  FIGURE 5 – IMPACT OF HARDWARE IMPAIRMENT COEFFICIENTS (i_b, i_k)
%  Single Robust benchmark. Hardware IS the swept variable here.
% --------------------------------------------------------------------------
fprintf('=== Figure 5: Hardware impairments sweep ===\n');

i_b_vec    = [0.01, 0.02, 0.03];
i_k_vec    = [0.02, 0.04, 0.06];
n_imp      = numel(i_b_vec);
>>>>>>> Newworld
imp_labels = cell(1, n_imp);
for idx = 1:n_imp
    imp_labels{idx} = sprintf('(%.2f, %.2f)', i_b_vec(idx), i_k_vec(idx));
end

<<<<<<< HEAD
% --- Pass 1: per-point normalisation bounds ---
=======
% --- Pass 1: bounds anchored to best impairment point ---
>>>>>>> Newworld
ec_max_imp = zeros(n_imp,1);  ec_min_imp = zeros(n_imp,1);
es_max_imp = zeros(n_imp,1);  es_min_imp = zeros(n_imp,1);

for idx = 1:n_imp
    fprintf('  [Norm pass] i_b=%.3f  i_k=%.3f\n', i_b_vec(idx), i_k_vec(idx));
    [ec_max_imp(idx), ec_min_imp(idx), es_max_imp(idx), es_min_imp(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
            i_b_vec(idx), i_k_vec(idx), p.P_static, ...
<<<<<<< HEAD
            p.delta_outage, p.xi_k_vec, p.T_max, p.epsilon, p.N_rand);
end

% Global bounds across all impairment levels
=======
            p.delta_outage, xi_k_rob, p.T_max, p.epsilon, p.N_rand);
end

>>>>>>> Newworld
[EEc_max_f5, best_idx5] = max(ec_max_imp);
EEc_min_f5 = ec_min_imp(best_idx5);
EEs_max_f5 = es_max_imp(best_idx5);
EEs_min_f5 = es_min_imp(best_idx5);
<<<<<<< HEAD

fprintf('  Fig5 EEc bounds: [%.4f, %.4f]\n', EEc_min_f5, EEc_max_f5);
fprintf('  Fig5 EEs bounds: [%.4f, %.4f]\n\n', EEs_min_f5, EEs_max_f5);

% --- Pass 2: solve for each (i_b, i_k) pair ---
=======
fprintf('  Fig5 EEc bounds: [%.4f, %.4f]\n', EEc_min_f5, EEc_max_f5);
fprintf('  Fig5 EEs bounds: [%.4f, %.4f]\n\n', EEs_min_f5, EEs_max_f5);

% --- Pass 2: solve ---
>>>>>>> Newworld
WEE_imp = zeros(n_imp,1);
EEc_imp = zeros(n_imp,1);
EEs_imp = zeros(n_imp,1);

for idx = 1:n_imp
    fprintf('  [Solve pass] i_b=%.3f  i_k=%.3f\n', i_b_vec(idx), i_k_vec(idx));
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, ...
        p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
<<<<<<< HEAD
        i_b_vec(idx), i_k_vec(idx), p.P_static, ...
        p.delta_outage, p.xi_k_vec, p.omega, ...
        p.T_max, p.epsilon, EEc_max_f5, EEs_max_f5, p.N_rand, EEc_min_f5, EEs_min_f5);

    EEc_imp(idx) = ec;
    EEs_imp(idx) = es;
=======
        i_b_vec(idx), i_k_vec(idx), p.P_static, p.delta_outage, xi_k_rob, p.omega, ...
        p.T_max, p.epsilon, EEc_max_f5, EEs_max_f5, p.N_rand, EEc_min_f5, EEs_min_f5);
    EEc_imp(idx) = ec;  EEs_imp(idx) = es;
>>>>>>> Newworld
    WEE_imp(idx) = p.omega     * ((ec - EEc_min_f5)/(EEc_max_f5 - EEc_min_f5)) + ...
                   (1-p.omega) * ((es - EEs_min_f5)/(EEs_max_f5 - EEs_min_f5));
    fprintf('    EEc=%.4f  EEs=%.4f  WEE=%.4f\n', ec, es, WEE_imp(idx));
end

<<<<<<< HEAD
% --- Plot ---
x_ticks = 1:n_imp;

figure(5); set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
bar(x_ticks, WEE_imp, 0.5, 'FaceColor',[0.2 0.5 0.8], 'EdgeColor','k');
grid on;
set(gca, 'XTick', x_ticks, 'XTickLabel', imp_labels, 'FontSize', 12);
=======
x_ticks = 1:n_imp;
figure(5); set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
bar(x_ticks, WEE_imp, 0.5, 'FaceColor',[0.2 0.5 0.8],'EdgeColor','k');
grid on;
set(gca,'XTick',x_ticks,'XTickLabel',imp_labels,'FontSize',12);
>>>>>>> Newworld
xlabel('(i_b,  i_k)  pairs','FontSize',13);
ylabel('Normalised WEE','FontSize',13);
title(sprintf('(a) WEE  (\\omega=%.2f)', p.omega),'FontSize',13);
xtickangle(15);

subplot(1,2,2);
yyaxis left;
plot(x_ticks, EEc_imp, 'b-o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b');
<<<<<<< HEAD
ylabel('Communication EE (EE_c)  [bps/Hz/W]','FontSize',13,'Color','b');
set(gca,'YColor','b', 'XTick', x_ticks, 'XTickLabel', imp_labels, 'FontSize', 12);

yyaxis right;
plot(x_ticks, EEs_imp, 'r-s','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r');
ylabel('Sensing EE (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r');

grid on;
=======
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13,'Color','b');
set(gca,'YColor','b','XTick',x_ticks,'XTickLabel',imp_labels,'FontSize',12);
yyaxis right;
plot(x_ticks, EEs_imp, 'r-s','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r');
ylabel('Sensing EE  (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r'); grid on;
>>>>>>> Newworld
xlabel('(i_b,  i_k)  pairs','FontSize',13);
title('(b) Raw EE_c and EE_s','FontSize',13);
legend('EE_c','EE_s','Location','best','FontSize',11);
xtickangle(15);

sgtitle(sprintf('Fig. 5: Impact of Hardware Impairments  (\\omega=%.2f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 5 rendered.\n\n');

%% =========================================================================
<<<<<<< HEAD
%%  FIGURE 6  –  Weighted EE vs P_max (Bỏ comment để chạy)
=======
%  FIGURE 6 – WEIGHTED EE vs P_max  (4 benchmarks)
% --------------------------------------------------------------------------
fprintf('=== Figure 6: P_max sweep (4 benchmarks) ===\n');
>>>>>>> Newworld

Pmax_dBm_vec = [20 25 30 35];
Pmax_W_vec   = db2pow(Pmax_dBm_vec) * 1e-3;
n_pmax       = numel(Pmax_dBm_vec);

<<<<<<< HEAD
ec_max_rob = zeros(n_pmax,1); ec_min_rob = zeros(n_pmax,1);
es_max_rob = zeros(n_pmax,1); es_min_rob = zeros(n_pmax,1);
ec_max_perf = zeros(n_pmax,1); ec_min_perf = zeros(n_pmax,1);
es_max_perf = zeros(n_pmax,1); es_min_perf = zeros(n_pmax,1);

for idx = 1:n_pmax
    P_cur = db2pow(Pmax_dBm_vec(idx)) * 1e-3;
    fprintf('  [Norm pass] P_max=%d dBm\n', Pmax_dBm_vec(idx));
    
    [ec_max_rob(idx), ec_min_rob(idx), es_max_rob(idx), es_min_rob(idx)] = get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, P_cur, p.sigma2, p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, p.T_max, p.epsilon, p.N_rand);
    [ec_max_perf(idx), ec_min_perf(idx), es_max_perf(idx), es_min_perf(idx)] = get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, P_cur, p.sigma2, p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_perf, p.T_max, p.epsilon, p.N_rand);
end

[EEc_max_f6_rob, best_ec6_rob] = max(ec_max_rob); EEc_min_f6_rob = ec_min_rob(best_ec6_rob);
EEs_max_f6_rob = es_max_rob(best_ec6_rob);        EEs_min_f6_rob = es_min_rob(best_ec6_rob);

[EEc_max_f6_perf, best_ec6_perf] = max(ec_max_perf); EEc_min_f6_perf = ec_min_perf(best_ec6_perf);
EEs_max_f6_perf = es_max_perf(best_ec6_perf);        EEs_min_f6_perf = es_min_perf(best_ec6_perf);

WEE_rob = zeros(n_pmax,1); EEc_rob = zeros(n_pmax,1); EEs_rob = zeros(n_pmax,1);
WEE_perf = zeros(n_pmax,1); EEc_perf = zeros(n_pmax,1); EEs_perf = zeros(n_pmax,1);

for idx = 1:n_pmax
    P_cur = db2pow(Pmax_dBm_vec(idx)) * 1e-3;
    fprintf('  [Solve pass] P_max=%d dBm\n', Pmax_dBm_vec(idx));

    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, P_cur, p.sigma2, p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, p.xi_k_vec, p.omega, p.T_max, p.epsilon, EEc_max_f6_rob, EEs_max_f6_rob, p.N_rand, EEc_min_f6_rob, EEs_min_f6_rob);
    EEc_rob(idx) = ec; EEs_rob(idx) = es;
    WEE_rob(idx) = p.omega*((ec - EEc_min_f6_rob)/(EEc_max_f6_rob - EEc_min_f6_rob)) + (1-p.omega)*((es - EEs_min_f6_rob)/(EEs_max_f6_rob - EEs_min_f6_rob));

    [~,~,ec_p,es_p,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, P_cur, p.sigma2, p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_perf, p.omega, p.T_max, p.epsilon, EEc_max_f6_perf, EEs_max_f6_perf, p.N_rand, EEc_min_f6_perf, EEs_min_f6_perf);
    EEc_perf(idx) = ec_p; EEs_perf(idx) = es_p;
    WEE_perf(idx) = p.omega*((ec_p - EEc_min_f6_perf)/(EEc_max_f6_perf - EEc_min_f6_perf)) + (1-p.omega)*((es_p - EEs_min_f6_perf)/(EEs_max_f6_perf - EEs_min_f6_perf));
end

figure(6); set(gcf,'Position',[100 100 1200 500]);
subplot(1,2,1);
plot(Pmax_dBm_vec, WEE_perf, 'k--^','LineWidth',2,'MarkerSize',8); hold on;
plot(Pmax_dBm_vec, WEE_rob, 'm-d','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','m');
grid on; xlabel('Maximum Transmit Power P_{max} (dBm)'); ylabel('Normalised WEE');
title(sprintf('(a) WEE (\\omega=%.2f)', p.omega)); legend('Perfect CSI', 'Statistical CSI');

subplot(1,2,2);
yyaxis left;
plot(Pmax_dBm_vec, EEc_perf, 'b--^','LineWidth',2,'MarkerSize',7); hold on;
plot(Pmax_dBm_vec, EEc_rob, 'b-o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','b');
ylabel('Communication EE (EE_c)'); set(gca,'YColor','b');
=======
% --- Pass 1 ---
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
fprintf('  Fig6 shared EEc bounds: [%.4f, %.4f]  (best P_max=%d dBm)\n', ...
        EEc_min_f6, EEc_max_f6, Pmax_dBm_vec(best_idx_f6));
fprintf('  Fig6 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f6, EEs_max_f6);

% --- Pass 2 ---
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

% --- Plot ---
figure(6); set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
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

subplot(1,2,2);
yyaxis left;
for v = 1:4
    plot(Pmax_dBm_vec, EEc_f6(:,v), styles{v}, 'LineWidth',1.8, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off;
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13,'Color','b');
set(gca,'YColor','b');
yyaxis right;
for v = 1:4
    plot(Pmax_dBm_vec, EEs_f6(:,v), styles{v}, 'LineWidth',1.8, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off;
ylabel('Sensing EE  (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r'); grid on;
xlabel('Maximum Transmit Power  P_{max}  (dBm)','FontSize',13);
title('(b) Raw EE_c and EE_s','FontSize',13);
legend([ec_lbl, es_lbl],'Location','best','FontSize',9);
set(gca,'FontSize',12);
>>>>>>> Newworld

yyaxis right;
plot(Pmax_dBm_vec, EEs_perf, 'm--^','LineWidth',2,'MarkerSize',7); hold on;
plot(Pmax_dBm_vec, EEs_rob, 'm-d','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','m');
ylabel('Sensing EE (EE_s)'); set(gca,'YColor','m');
grid on; xlabel('Maximum Transmit Power P_{max} (dBm)'); title('(b) Raw EE_c & EE_s');
legend('EE_c (Perf)','EE_c (Stat)','EE_s (Perf)','EE_s (Stat)','Location','best');
sgtitle(sprintf('Fig. 6: Weighted EE vs. P_{max} (\\omega=%.2f)', p.omega), 'FontWeight','bold');
drawnow;

%% =========================================================================
<<<<<<< HEAD
%%  FIGURE 7  –  Weighted EE vs N
fprintf('=== Figure 7: Computing global norms across all N values ===\n');
N_vec = 6:1:9 ;
n_N   = numel(N_vec);

ec_max_rob = zeros(n_N,1); ec_min_rob = zeros(n_N,1);
es_max_rob = zeros(n_N,1); es_min_rob = zeros(n_N,1);
ec_max_perf = zeros(n_N,1); ec_min_perf = zeros(n_N,1);
es_max_perf = zeros(n_N,1); es_min_perf = zeros(n_N,1);

for idx = 1:n_N
    N_cur   = N_vec(idx);
    H_norm  = H_bank(1:N_cur, :, 1);
    
    xi_k_cur = zeros(1,p.K);
    for k = 1:p.K; xi_k_cur(k) = p.delta * norm(H_norm(:,k)); end
    
    fprintf('  [Norm pass] N=%d\n', N_cur);
    [ec_max_rob(idx), ec_min_rob(idx), es_max_rob(idx), es_min_rob(idx)] = get_norm_constants(H_norm, p.theta_targets, N_cur, p.K, p.M, p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_cur, p.T_max, p.epsilon, p.N_rand);
    [ec_max_perf(idx), ec_min_perf(idx), es_max_perf(idx), es_min_perf(idx)] = get_norm_constants(H_norm, p.theta_targets, N_cur, p.K, p.M, p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_perf, p.T_max, p.epsilon, p.N_rand);
end

[EEc_max_f7_rob, best_ec7_rob] = max(ec_max_rob); EEc_min_f7_rob = ec_min_rob(best_ec7_rob);
EEs_max_f7_rob = es_max_rob(best_ec7_rob);        EEs_min_f7_rob = es_min_rob(best_ec7_rob);

[EEc_max_f7_perf, best_ec7_perf] = max(ec_max_perf); EEc_min_f7_perf = ec_min_perf(best_ec7_perf);
EEs_max_f7_perf = es_max_perf(best_ec7_perf);        EEs_min_f7_perf = es_min_perf(best_ec7_perf);
fprintf('  Fig7 Robust  : EEc=[%.4f, %.4f]  EEs=[%.4f, %.4f]  (best N=%d)\n', ...
        EEc_min_f7_rob,  EEc_max_f7_rob,  EEs_min_f7_rob,  EEs_max_f7_rob,  N_vec(best_ec7_rob));
fprintf('  Fig7 Perfect : EEc=[%.4f, %.4f]  EEs=[%.4f, %.4f]  (best N=%d)\n', ...
        EEc_min_f7_perf, EEc_max_f7_perf, EEs_min_f7_perf, EEs_max_f7_perf, N_vec(best_ec7_perf));
WEE_rob = zeros(n_N,1); EEc_rob = zeros(n_N,1); EEs_rob = zeros(n_N,1);
WEE_perf = zeros(n_N,1); EEc_perf = zeros(n_N,1); EEs_perf = zeros(n_N,1);

for idx = 1:n_N
    N_cur = N_vec(idx);
    H_cur = H_bank(1:N_cur, :, 1);   
    
    xi_k_cur = zeros(1, p.K);
    for k = 1:p.K; xi_k_cur(k) = p.delta * norm(H_cur(:,k)); end

    fprintf('  [Solve pass] N=%d\n', N_cur);

    [~,~,ec,es,~] = solve_ISAC(H_cur, p.theta_targets, N_cur, p.K, p.M, p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_cur, p.omega, p.T_max, p.epsilon, EEc_max_f7_rob, EEs_max_f7_rob, p.N_rand, EEc_min_f7_rob, EEs_min_f7_rob);
    EEc_rob(idx) = ec; EEs_rob(idx) = es;
    WEE_rob(idx) = p.omega*((ec - EEc_min_f7_rob)/(EEc_max_f7_rob - EEc_min_f7_rob)) + (1-p.omega)*((es - EEs_min_f7_rob)/(EEs_max_f7_rob - EEs_min_f7_rob));

    [~,~,ec_p,es_p,~] = solve_ISAC(H_cur, p.theta_targets, N_cur, p.K, p.M, p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_perf, p.omega, p.T_max, p.epsilon, EEc_max_f7_perf, EEs_max_f7_perf, p.N_rand, EEc_min_f7_perf, EEs_min_f7_perf);
    EEc_perf(idx) = ec_p; EEs_perf(idx) = es_p;
    WEE_perf(idx) = p.omega*((ec_p - EEc_min_f7_perf)/(EEc_max_f7_perf - EEc_min_f7_perf)) + (1-p.omega)*((es_p - EEs_min_f7_perf)/(EEs_max_f7_perf - EEs_min_f7_perf));
end

figure(7); set(gcf,'Position',[100 100 1200 500]);
subplot(1,2,1);
plot(N_vec, WEE_perf, 'k--^','LineWidth',2,'MarkerSize',8); hold on;
plot(N_vec, WEE_rob,'k-p','LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[0.35 0.35 0.35]);
grid on; xlabel('Number of Antennas N'); ylabel('Normalised WEE');
title(sprintf('(a) WEE (\\omega=%.2f)', p.omega)); legend('Perfect CSI', 'Statistical CSI');

subplot(1,2,2);
yyaxis left;
plot(N_vec, EEc_perf, 'b--^','LineWidth',2,'MarkerSize',7); hold on;
plot(N_vec, EEc_rob, 'b-o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','b');
ylabel('Communication EE (EE_c)'); set(gca,'YColor','b');

yyaxis right;
plot(N_vec, EEs_perf, 'k--^','LineWidth',2,'MarkerSize',7); hold on;
plot(N_vec, EEs_rob, 'k-p','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0.35 0.35 0.35]);
ylabel('Sensing EE (EE_s)'); set(gca,'YColor','k');
grid on; xlabel('Number of Antennas N'); title('(b) Raw EE_c & EE_s');
legend('EE_c (Perf)','EE_c (Stat)','EE_s (Perf)','EE_s (Stat)','Location','best');
sgtitle(sprintf('Fig. 7: Weighted EE vs. N (\\omega=%.2f)', p.omega), 'FontWeight','bold');
drawnow;

%% =========================================================================
%%  FIGURE 8 – IMPACT OF CSI UNCERTAINTY LEVEL delta
%
%  delta = 0  --> perfect CSI benchmark (xi_k = 0, no robustness penalty)
%  delta > 0  --> robust design with increasing uncertainty
% --------------------------------------------------------------------------
fprintf('=== Figure 8: Impact of CSI uncertainty level delta ===\n');
delta_vec = [0, 0.005, 0.01, 0.02, 0.04, 0.05];
n_delta   = numel(delta_vec);

% --- Pass 1: per-point bounds ---
ec_max_all = zeros(n_delta,1);  ec_min_all = zeros(n_delta,1);
es_max_all = zeros(n_delta,1);  es_min_all = zeros(n_delta,1);

for idx = 1:n_delta
    delta_cur = delta_vec(idx);
    
    % Tính vector phương sai xi_k_cur cho mức delta hiện tại
    xi_k_cur = zeros(1, p.K);
    for k = 1:p.K
        xi_k_cur(k) = delta_cur * norm(H_hat(:, k));
    end
    
    fprintf('  [Norm pass] delta=%.4f\n', delta_cur);
    [ec_max_all(idx), ec_min_all(idx), es_max_all(idx), es_min_all(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
            p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_cur, p.T_max, p.epsilon, p.N_rand);
end

[EEc_max_f8, best_idx8] = max(ec_max_all); EEc_min_f8 = ec_min_all(best_idx8);
EEs_max_f8 = es_max_all(best_idx8);        EEs_min_f8 = es_min_all(best_idx8);

fprintf('  Fig8 EEc bounds: [%.4f, %.4f]  (best delta=%.4f)\n', ...
        EEc_min_f8, EEc_max_f8, delta_vec(best_idx8));
fprintf('  Fig8 EEs bounds: [%.4f, %.4f]  (best delta=%.4f)\n', ...
        EEs_min_f8, EEs_max_f8, delta_vec(best_idx8));

% --- Pass 2: solve for all delta values (delta=0 IS the perfect-CSI point) ---
WEE_delta = zeros(n_delta,1);
EEc_delta = zeros(n_delta,1);
EEs_delta = zeros(n_delta,1);

for idx = 1:n_delta
    delta_cur = delta_vec(idx);
    
    % Tính lại vector phương sai xi_k_cur
    xi_k_cur = zeros(1, p.K);
    for k = 1:p.K
        xi_k_cur(k) = delta_cur * norm(H_hat(:, k));
    end
    
    fprintf('  [Solve pass] delta=%.4f\n', delta_cur);
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
        p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_cur, p.omega, ...
        p.T_max, p.epsilon, EEc_max_f8, EEs_max_f8, p.N_rand, EEc_min_f8, EEs_min_f8);
        
    EEc_delta(idx) = ec;
    EEs_delta(idx) = es;
    WEE_delta(idx) = p.omega     * ((ec - EEc_min_f8)/(EEc_max_f8 - EEc_min_f8)) + ...
                     (1-p.omega) * ((es - EEs_min_f8)/(EEs_max_f8 - EEs_min_f8));
    fprintf('    EEc=%.4f  EEs=%.4f  WEE=%.4f\n', ec, es, WEE_delta(idx));
end

% delta=0 is the perfect-CSI reference point; delta>0 are the robust points
perfect_idx   = 1;
imperfect_idx = 2:n_delta;

WEE_perf_val = WEE_delta(perfect_idx);
EEc_perf_val = EEc_delta(perfect_idx);
EEs_perf_val = EEs_delta(perfect_idx);

figure(8);
set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
plot(delta_vec(imperfect_idx), WEE_delta(imperfect_idx), ...
     'b-o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b');
hold on;
yline(WEE_perf_val, 'k--', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('CSI Uncertainty Level  \delta','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('(a) WEE  (\\omega=%.2f)', p.omega),'FontSize',13);
legend('Statistical CSI (Robust)', 'Perfect CSI (\delta=0)', ...
       'Location','northeast','FontSize',11);
=======
%  FIGURE 7 – WEIGHTED EE vs NUMBER OF ANTENNAS N  (4 benchmarks)
%  Note: xi_k is N-dependent and recomputed locally per N.
%        Variants 1 & 3 use statistical xi_k; variants 2 & 4 use xi_k = 0.
% --------------------------------------------------------------------------
fprintf('=== Figure 7: N sweep (4 benchmarks) ===\n');

N_vec = N_vec_f7;
n_N   = numel(N_vec);

% --- Pass 1 ---
ec_max_all = zeros(n_N,4);  ec_min_all = zeros(n_N,4);
es_max_all = zeros(n_N,4);  es_min_all = zeros(n_N,4);

for v = 1:4
    for idx = 1:n_N
        N_cur  = N_vec(idx);
        H_norm = H_bank(1:N_cur, :, 1);

        % xi_k scales with N; variants 1 & 3 = statistical, 2 & 4 = perfect CSI
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
fprintf('  Fig7 shared EEc bounds: [%.4f, %.4f]  (best N=%d)\n', ...
        EEc_min_f7, EEc_max_f7, N_vec(best_idx_f7));
fprintf('  Fig7 shared EEs bounds: [%.4f, %.4f]\n', EEs_min_f7, EEs_max_f7);

% --- Pass 2 ---
WEE_f7 = zeros(n_N,4);
EEc_f7 = zeros(n_N,4);
EEs_f7 = zeros(n_N,4);

for v = 1:4
    for idx = 1:n_N
        N_cur = N_vec(idx);
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
        EEc_f7(idx,v) = ec;  EEs_f7(idx,v) = es;
        WEE_f7(idx,v) = p.omega     * ((ec - EEc_min_f7)/(EEc_max_f7 - EEc_min_f7)) + ...
                        (1-p.omega) * ((es - EEs_min_f7)/(EEs_max_f7 - EEs_min_f7));
        fprintf('    WEE=%.4f\n', WEE_f7(idx,v));
    end
end

% --- Plot ---
figure(7); set(gcf,'Position',[100 100 1200 500]);

subplot(1,2,1);
for v = 1:4
    plot(N_vec, WEE_f7(:,v), styles{v}, 'LineWidth',2, ...
        'MarkerSize',7,'MarkerFaceColor',mfc{v}); hold on;
end
hold off; grid on;
xlabel('Number of Antennas  N','FontSize',13);
ylabel('Normalised Weighted EE','FontSize',13);
title(sprintf('(a) WEE  (\\omega=%.2f)', p.omega),'FontSize',13);
legend(v_lbl,'Location','best','FontSize',10);
>>>>>>> Newworld
set(gca,'FontSize',12);

subplot(1,2,2);
yyaxis left;
<<<<<<< HEAD
plot(delta_vec(imperfect_idx), EEc_delta(imperfect_idx), ...
     'b-o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','b');
hold on;
yline(EEc_perf_val, 'b--', 'LineWidth', 1.5);
=======
for v = 1:4
    plot(N_vec, EEc_f7(:,v), styles{v}, 'LineWidth',1.8, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
>>>>>>> Newworld
hold off;
ylabel('Communication EE  (EE_c)  [bps/Hz/W]','FontSize',13,'Color','b');
set(gca,'YColor','b');

yyaxis right;
<<<<<<< HEAD
plot(delta_vec(imperfect_idx), EEs_delta(imperfect_idx), ...
     'r-s','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','r');
hold on;
yline(EEs_perf_val, 'r--', 'LineWidth', 1.5);
hold off;
ylabel('Sensing EE  (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r');
grid on;
xlabel('CSI Uncertainty Level  \delta','FontSize',13);
title('(b) Raw EE_c and EE_s','FontSize',13);
legend('EE_c (Statistical)', 'EE_c (Perfect CSI)', ...
       'EE_s (Statistical)', 'EE_s (Perfect CSI)', ...
       'Location','best','FontSize',11);
=======
for v = 1:4
    plot(N_vec, EEs_f7(:,v), styles{v}, 'LineWidth',1.8, ...
        'MarkerSize',6,'MarkerFaceColor',mfc{v}); hold on;
end
hold off;
ylabel('Sensing EE  (EE_s)','FontSize',13,'Color','r');
set(gca,'YColor','r'); grid on;
xlabel('Number of Antennas  N','FontSize',13);
title('(b) Raw EE_c and EE_s','FontSize',13);
legend([ec_lbl, es_lbl],'Location','best','FontSize',9);
>>>>>>> Newworld
set(gca,'FontSize',12);

sgtitle(sprintf('Fig. 8: Impact of CSI Uncertainty (Statistical)  (\\omega=%.2f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
<<<<<<< HEAD
fprintf('    Figure 8 rendered.\n\n');
=======
fprintf('    Figure 7 rendered.\n\n');

%% =========================================================================
%  FIGURE 8 – IMPACT OF CHANNEL UNCERTAINTY LEVEL delta
%
%  delta = 0  --> perfect channel knowledge (xi_k=0, no robustness cost)
%  delta > 0  --> statistical robust design: xi_k = delta * ||h_hat_k||
%  Hardware impairments fixed at p.i_b / p.i_k throughout.
% --------------------------------------------------------------------------
fprintf('=== Figure 8: Channel uncertainty level delta ===\n');

delta_vec = [0, 0.005, 0.01, 0.02, 0.03];
n_delta   = numel(delta_vec);

% --- Pass 1: shared bounds across all delta values ---
ec_max_all = zeros(n_delta,1);  ec_min_all = zeros(n_delta,1);
es_max_all = zeros(n_delta,1);  es_min_all = zeros(n_delta,1);

for idx = 1:n_delta
    delta_cur = delta_vec(idx);
    xi_k_cur  = zeros(1, p.K);
    for k = 1:p.K
        xi_k_cur(k) = delta_cur * norm(H_hat(:,k));
    end
    fprintf('  [Norm pass] delta=%.4f\n', delta_cur);
    [ec_max_all(idx), ec_min_all(idx), es_max_all(idx), es_min_all(idx)] = ...
        get_norm_constants(H_hat, p.theta_targets, p.N, p.K, p.M, ...
            p.P_max, p.sigma2, p.Gamma_min, p.gamma_min, ...
            p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_cur, p.T_max, p.epsilon, p.N_rand);
end

[EEc_max_f8, best_idx8] = max(ec_max_all);
EEc_min_f8 = ec_min_all(best_idx8);
EEs_max_f8 = es_max_all(best_idx8);
EEs_min_f8 = es_min_all(best_idx8);
fprintf('  Fig8 EEc bounds: [%.4f, %.4f]  (best delta=%.4f)\n', ...
        EEc_min_f8, EEc_max_f8, delta_vec(best_idx8));
fprintf('  Fig8 EEs bounds: [%.4f, %.4f]  (best delta=%.4f)\n', ...
        EEs_min_f8, EEs_max_f8, delta_vec(best_idx8));

% --- Pass 2: solve for all delta values ---
WEE_delta = zeros(n_delta,1);
EEc_delta = zeros(n_delta,1);
EEs_delta = zeros(n_delta,1);

for idx = 1:n_delta
    delta_cur = delta_vec(idx);
    xi_k_cur  = zeros(1, p.K);
    for k = 1:p.K
        xi_k_cur(k) = delta_cur * norm(H_hat(:,k));
    end
    fprintf('  [Solve pass] delta=%.4f\n', delta_cur);
    [~,~,ec,es,~] = solve_ISAC(H_hat, p.theta_targets, p.N, p.K, p.M, p.P_max, p.sigma2, ...
        p.Gamma_min, p.gamma_min, p.i_b, p.i_k, p.P_static, p.delta_outage, xi_k_cur, p.omega, ...
        p.T_max, p.epsilon, EEc_max_f8, EEs_max_f8, p.N_rand, EEc_min_f8, EEs_min_f8);
    EEc_delta(idx) = ec;
    EEs_delta(idx) = es;
    WEE_delta(idx) = p.omega     * ((ec - EEc_min_f8)/(EEc_max_f8 - EEc_min_f8)) + ...
                     (1-p.omega) * ((es - EEs_min_f8)/(EEs_max_f8 - EEs_min_f8));
    fprintf('    EEc=%.4f  EEs=%.4f  WEE=%.4f\n', ec, es, WEE_delta(idx));
end

% delta=0 is the perfect channel knowledge reference; delta>0 are robust points
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
legend('Statistical CSI (Robust)','Perfect Channel (\delta=0)', ...
       'Location','northeast','FontSize',11);
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
legend('EE_c (Statistical)','EE_c (Perfect, \delta=0)', ...
       'EE_s (Statistical)','EE_s (Perfect, \delta=0)', ...
       'Location','best','FontSize',11);
set(gca,'FontSize',12);

sgtitle(sprintf('Fig. 8: Impact of Channel Uncertainty Level \\delta  (\\omega=%.2f)', p.omega), ...
        'FontSize',14,'FontWeight','bold');
drawnow;
fprintf('    Figure 8 rendered.\n\n');
fprintf('=== All figures complete. ===\n');
>>>>>>> Newworld
