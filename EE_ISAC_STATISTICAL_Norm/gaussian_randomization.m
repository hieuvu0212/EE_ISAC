% -------------------------------------------------------------------------
%  gaussian_randomization : Algorithm 2  (Statistical CSI version)
%
%  SIGNATURE CHANGE vs original:
%    xi_k_vec  added as 16th argument (was missing — needed for Bernstein).
%    delta_outage remains the 15th argument (unchanged).
%
%  FIX: Randomized rank-1 candidates are now checked against ALL three
%  Bernstein robust SINR constraints before being accepted:
%
%    Build Qk for each user k from the rank-1 matrices, then compute
%    the tightest feasible Bernstein slack pair (u_k*, v_k*):
%
%    (C3) v_k* = max(0,  -xi_k^2 * lambda_min(Qk))
%         Ensures  v_k*I + xi_k^2*Qk >= 0.
%
%    (C2) u_k* = || [xi_k^2 * vec(Qk) ;  sqrt(2)*xi_k * Qk * h_hat_k] ||
%         Minimum u_k satisfying the SOC constraint.
%
%    (C1) Check:
%         xi_k^2*Tr(Qk) - sqrt(2*ln(1/delta_outage))*u_k*
%             + ln(delta_outage)*v_k*  + c_k >= 0
%         where c_k = real(h_hat_k' * Qk * h_hat_k) - (1+i_k)*sigma2.
%
%  A candidate is accepted only when ALL three hold for ALL users.
%  The sensing threshold is checked first (cheap gate).
%
%  CALLERS: solve_ISAC must pass xi_k_vec as the 16th argument.
% -------------------------------------------------------------------------
function [Wc_r1, Ws_r1] = gaussian_randomization( ...
        Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
        N, K, M, P_max, sigma2, Gamma_min, gamma_min, ...
        i_b, i_k, P_static, delta_outage, xi_k_vec, ...
        omega, EEcmax, EEsmax, N_rand, EEcmin, EEsmin)

    sl = @(th) exp(1j*pi*((0:N-1)'-(N-1)/2)*sind(th));

    % --- Pre-compute Bernstein scalar constants --------------------------
    % ln(delta_outage) < 0,  sqrt(2*ln(1/delta)) > 0
    ln_delta    = log(delta_outage);
    sqrt_factor = sqrt(-2 * ln_delta);   % = sqrt(2*ln(1/delta_outage))

    % --- Eigen-decompose each SDR matrix ---------------------------------
    Uc = cell(1,K);  Sc_half = cell(1,K);
    Us = cell(1,M);  Ss_half = cell(1,M);
    for j = 1:K
        [U,S] = eig((Wc_SDR{j}+Wc_SDR{j}')/2);
        s = max(real(diag(S)), 0);
        Uc{j} = U;  Sc_half{j} = diag(sqrt(s));
    end
    for l = 1:M
        [U,S] = eig((Ws_SDR{l}+Ws_SDR{l}')/2);
        s = max(real(diag(S)), 0);
        Us{l} = U;  Ss_half{l} = diag(sqrt(s));
    end

    % --- Initialise best with the SDR solution itself -------------------
    nc = max(EEcmax - EEcmin, 1e-9);
    ns = max(EEsmax - EEsmin, 1e-9);
    [EEc_sdr0, EEs_sdr0] = compute_EE(Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
                        sigma2, i_b, i_k, P_static, N, K, M);
    best_obj = omega*((EEc_sdr0 - EEcmin)/nc) + (1-omega)*((EEs_sdr0 - EEsmin)/ns);
    Wc_r1    = Wc_SDR;
    Ws_r1    = Ws_SDR;

    % =====================================================================
    for n_trial = 1:N_rand

        % --- Generate candidate beamforming vectors ----------------------
        wc = cell(1,K);  ws = cell(1,M);
        for j = 1:K
            z = (randn(N,1)+1j*randn(N,1))/sqrt(2);
            wc{j} = Uc{j} * Sc_half{j} * z;
        end
        for l = 1:M
            z = (randn(N,1)+1j*randn(N,1))/sqrt(2);
            ws{l} = Us{l} * Ss_half{l} * z;
        end

        % Scale to satisfy total power budget
        P_cand = 0;
        for j = 1:K; P_cand = P_cand + norm(wc{j})^2; end
        for l = 1:M; P_cand = P_cand + norm(ws{l})^2; end
        if P_cand > P_max
            sc = sqrt(P_max / P_cand);
            for j = 1:K; wc{j} = wc{j}*sc; end
            for l = 1:M; ws{l} = ws{l}*sc; end
        end

        % Rebuild rank-1 covariance matrices
        Wc_n = cell(1,K);  Ws_n = cell(1,M);
        for j = 1:K; Wc_n{j} = wc{j}*wc{j}'; end
        for l = 1:M; Ws_n{l} = ws{l}*ws{l}'; end

        % Aggregate total beamforming matrix
        W_tot_n = zeros(N,N);
        for j = 1:K; W_tot_n = W_tot_n + Wc_n{j}; end
        for l = 1:M; W_tot_n = W_tot_n + Ws_n{l};  end
        tr_tot  = real(trace(W_tot_n));

        % =================================================================
        % GATE 1: Sensing threshold  (cheap — check first)
        % =================================================================
        ok = true;
        for m = 1:M
            vm = sl(theta_targets(m));
            if real(vm'*W_tot_n*vm) + i_b*tr_tot < Gamma_min
                ok = false; break;
            end
        end
        if ~ok; continue; end

        % =================================================================
        % GATE 2: Bernstein robust SINR constraints  (C1), (C2), (C3)
        %
        %  For each user k, build Qk from the candidate rank-1 matrices,
        %  derive the minimum-feasible slack pair (u_k*, v_k*), and verify
        %  the scalar Bernstein inequality (C1).
        % =================================================================
        for k = 1:K
            hk   = H_hat(:,k);
            xi_k = xi_k_vec(k);

            % Build Qk  (same definition as in solve_ISAC CVX block)
            Qk = (1/gamma_min) * Wc_n{k};
            for j = 1:K
                if j ~= k; Qk = Qk - Wc_n{j}; end
            end
            for l = 1:M
                Qk = Qk - Ws_n{l};
            end
            Qk = Qk - i_k * W_tot_n ...
                     - (1+i_k)*i_b * diag(diag(W_tot_n));

            % c_k  (deterministic SINR residual at nominal channel)
            c_k = real(hk' * Qk * hk) - (1+i_k)*sigma2;

            % --- (C3) Minimum v_k* satisfying  v_k*I + xi_k^2*Qk >= 0 --
            %     Requires  v_k >= -xi_k^2 * lambda_min(Qk)
            %     (if lambda_min(Qk) >= 0 then v_k* = 0 suffices)
            lambda_min_Qk = min(real(eig((Qk+Qk')/2)));
            v_k_star = max(0, -xi_k^2 * lambda_min_Qk);

            % --- (C2) Minimum u_k* satisfying the SOC --------------------
            %     u_k >= || [xi_k^2 * vec(Qk) ;  sqrt(2)*xi_k * Qk * hk] ||
            vec_Qk   = Qk(:);   % N^2 x 1 complex vector
            lhs_vec  = [xi_k^2 * vec_Qk ; sqrt(2)*xi_k * (Qk * hk)];
            u_k_star = norm(lhs_vec);

            % --- (C1) Bernstein scalar inequality check ------------------
            %     xi_k^2*Tr(Qk) - sqrt(2*ln(1/delta))*u_k*
            %         + ln(delta)*v_k*  + c_k >= 0
            bernstein_val = xi_k^2 * real(trace(Qk)) ...
                          - sqrt_factor * u_k_star ...
                          + ln_delta    * v_k_star ...
                          + c_k;

            if bernstein_val < 0
                ok = false; break;
            end
        end
        if ~ok; continue; end

        % =================================================================
        % Candidate passed all checks — evaluate objective
        % =================================================================
        [EEc_n, EEs_n] = compute_EE(Wc_n, Ws_n, H_hat, theta_targets, ...
                            sigma2, i_b, i_k, P_static, N, K, M);
        f_n = omega*((EEc_n - EEcmin)/nc) + (1-omega)*((EEs_n - EEsmin)/ns);

        if f_n > best_obj
            best_obj = f_n;
            Wc_r1 = Wc_n;
            Ws_r1 = Ws_n;
        end
    end
end
