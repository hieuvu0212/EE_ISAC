% -------------------------------------------------------------------------
%  solve_ISAC : Dinkelbach-SCA solver + Gaussian Randomization
%  Statistical CSI version: Bernstein-type probabilistic SINR constraints.
%
%  Channel error model:
%    h_k = h_hat_k + Delta_h_k,
%    Delta_h_k ~ CN(0, xi_k^2 * I_N),
%    xi_k^2 = delta^2 * || vec(h_hat_k) ||^2
%
%  Robust SINR constraint (Bernstein inequality form):
%    xi_k^2 * Tr(Q_k) - sqrt(2*ln(1/delta_outage)) * u_k + ln(delta_outage) * v_k
%        + c_k >= 0,                                            
%    || [xi_k^2 * vec(Q_k) ; sqrt(2)*xi_k * Q_k * h_hat_k] || <= u_k,  
%    v_k * I_N + xi_k^2 * Q_k >= 0,   v_k >= 0,                
%
%  where  Q_k = (1/gamma_min)*W_{c,k} - sum_{j~=k} W_{c,j}
%              - sum_l W_{s,l} - i_k * W_total
%              - (1+i_k)*i_b * diag(diag(W_total))
%  and    c_k = h_hat_k' * Q_k * h_hat_k - (1+i_k)*sigma2.
%
%  Returns rank-1 physical matrices and EE values from those matrices.
% -------------------------------------------------------------------------
function [Wc_r1, Ws_r1, EEc, EEs, obj_hist] = ...
    solve_ISAC(H_hat, theta_targets, N, K, M, P_max, sigma2, ...
               Gamma_min, gamma_min, i_b, i_k, P_static, ...
               delta_outage, xi_k_vec, ...
               omega, T_max, epsilon, EEcmax, EEsmax, N_rand,EEcmin, EEsmin)
    if nargin < 21; EEcmin = 0; end
    if nargin < 22; EEsmin = 0; end 
    sl    = @(th) exp(1j*pi*((0:N-1)'-(N-1)/2)*sind(th));
    log_2 = log(2);

    % --- Initialise beamformers -------------------------------------------
    p0 = P_max / (K + M);
    Wc = cell(1,K);  Ws = cell(1,M);
    for j = 1:K
        w = (randn(N,1)+1j*randn(N,1))/sqrt(2);
        w = w/norm(w)*sqrt(p0);
        Wc{j} = w*w';
    end
    for l = 1:M
        w = sl(theta_targets(l))*sqrt(p0);
        Ws{l} = w*w';
    end

    nc = max(EEcmax - EEcmin, 1e-9);
    ns = max(EEsmax - EEsmin, 1e-9);
    q  = 0;
    obj_hist = zeros(T_max,1);
    Wc_SDR = Wc;  Ws_SDR = Ws;

    % --- Dinkelbach-SCA main loop ----------------------------------------
    for t = 1:T_max

        % Compute linearisation point U_k^(t) at current iterate
        W_total = zeros(N,N);
        for j = 1:K; W_total = W_total + Wc{j}; end
        for l = 1:M; W_total = W_total + Ws{l};  end

        Uk_t = zeros(K,1);
        for k = 1:K
            hk  = H_hat(:,k);  hh = hk*hk';
            inter = 0;
            for j = 1:K; if j~=k; inter=inter+real(hk'*Wc{j}*hk); end; end
            for l = 1:M;          inter=inter+real(hk'*Ws{l}*hk);       end
            Uk_t(k) = inter ...
                    + i_k*real(hk'*W_total*hk) ...
                    + (1+i_k)*i_b*real(hk'*diag(diag(W_total))*hk) ...
                    + (1+i_k)*sigma2;
            Uk_t(k) = max(Uk_t(k), 1e-30);
        end

        % -----------------------------------------------------------------
        % CVX SDP
        % -----------------------------------------------------------------
        cvx_begin quiet
            cvx_solver mosek 
            variable Wc_v(N,N,K) complex semidefinite
            variable Ws_v(N,N,M) complex semidefinite
            variable u_v(K) nonnegative          % Bernstein slack u_k >= 0
            variable v_v(K) nonnegative          % Bernstein slack v_k >= 0

            W_sum     = sum(Wc_v,3) + sum(Ws_v,3);
            tot_power = real(trace(W_sum));

            % --- SCA lower bound of sum-rate ---
            Rk_lb = 0;
            for k = 1:K
                hk = H_hat(:,k);  hh = hk*hk';
                Tk_cvx = (1+i_k)*real(trace(W_sum*hh)) ...
                       + (1+i_k)*i_b*real(trace(diag(diag(W_sum))*hh)) ...
                       + (1+i_k)*sigma2;

                inter_k = 0;
                for j = 1:K
                    if j~=k
                        inter_k = inter_k + real(trace(Wc_v(:,:,j)*hh));
                    end
                end
                for l = 1:M
                    inter_k = inter_k + real(trace(Ws_v(:,:,l)*hh));
                end
                Uk_cvx = inter_k ...
                       + i_k*real(trace(W_sum*hh)) ...
                       + (1+i_k)*i_b*real(trace(diag(diag(W_sum))*hh)) ...
                       + (1+i_k)*sigma2;

                Rk_lb = Rk_lb + ( log(Tk_cvx)/log_2 ...
                    - log(Uk_t(k))/log_2 ...
                    - (1/(Uk_t(k)*log_2)) * (Uk_cvx - Uk_t(k)) );
            end

            % --- Sensing beampattern gain ---
            sense_gain = 0;
            for m = 1:M
                vm = sl(theta_targets(m));
                sense_gain = sense_gain + real(trace((vm*vm')*W_sum)) + i_b*tot_power;
            end

            % --- Dinkelbach objective ---
            P_tot_expr = P_static + (1+i_b)*tot_power;
            obj_expr = (omega/nc)*(Rk_lb - EEcmin*P_tot_expr) + ((1-omega)/ns)*(sense_gain - EEsmin*P_tot_expr) - q * P_tot_expr;

            maximize( obj_expr )

            subject to
                % Power budget
                tot_power <= P_max;

                % Sensing threshold  forall m
                for m = 1:M
                    vm = sl(theta_targets(m));
                    real(trace((vm*vm')*W_sum)) + i_b*tot_power >= Gamma_min;
                end

                % -------------------------------------------------------
                % Statistical CSI robust SINR constraints (Bernstein form)
                % -------------------------------------------------------
                for k = 1:K
                    hk   = H_hat(:,k);
                    xi_k = xi_k_vec(k);   % per-user uncertainty std

                    % Build Q_k (effective SINR numerator matrix)
                    Qk = (1/gamma_min)*Wc_v(:,:,k);
                    for j = 1:K; if j~=k; Qk = Qk - Wc_v(:,:,j); end; end
                    for l = 1:M;          Qk = Qk - Ws_v(:,:,l);       end
                    Qk = Qk - i_k*W_sum ...
                             - (1+i_k)*i_b*diag(diag(W_sum));

                    % c_k  (deterministic part of SINR at nominal channel)
                    c_k = real(hk'*Qk*hk) - (1+i_k)*sigma2;

                    % (C1) Bernstein scalar constraint
                    xi_k^2 * real(trace(Qk)) ...
                        - sqrt(2*log(1/delta_outage)) * u_v(k) ...
                        + log(delta_outage) * v_v(k) ...
                        + c_k >= 0;

                    % (C2) Second-order cone constraint
                    % || [xi_k^2 * vec(Qk) ; sqrt(2)*xi_k * Qk * hk] || <= u_k
                    % vec(Qk) for an N x N matrix -> N^2 x 1 real+imag stack
                    vec_Qk = Qk(:);   % complex N^2-vector
                    lhs_vec = [xi_k^2 * vec_Qk ; sqrt(2)*xi_k * Qk * hk];
                    norm(lhs_vec) <= u_v(k);

                    % (C3) Matrix inequality: v_k*I + xi_k^2 * Q_k >= 0
                    v_v(k)*eye(N) + xi_k^2 * Qk == hermitian_semidefinite(N);
                end

        cvx_end

        if strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Failed')
            warning('CVX: %s at t=%d. Using last feasible point.', ...
                    cvx_status, t);
            obj_hist = obj_hist(1:max(t-1,1));
            break;
        end

        % Store SDR solution (symmetrised)
        for j = 1:K
            Wc_SDR{j} = (Wc_v(:,:,j)+Wc_v(:,:,j)')/2;
        end
        for l = 1:M
            Ws_SDR{l} = (Ws_v(:,:,l)+Ws_v(:,:,l)')/2;
        end
        Wc = Wc_SDR;  Ws = Ws_SDR;

        % Update Dinkelbach parameter using SDR solution
        [EEc_sdr, EEs_sdr, R_sum_sdr, p_sense_sdr] = ...
            compute_EE(Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
                       sigma2, i_b, i_k, P_static, N, K, M);
        % Total power consumed at this iterate
        W_tmp = zeros(N,N);
        for j=1:K; W_tmp=W_tmp+Wc_SDR{j}; end
        for l=1:M; W_tmp=W_tmp+Ws_SDR{l}; end
        P_cons_sdr = P_static + (1+i_b)*real(trace(W_tmp));

        % Raw numerator matching exactly what CVX maximises:
        %   omega/nc * R_sum  +  (1-omega)/ns * p_sense
        num_raw = (omega/nc)*(R_sum_sdr - EEcmin*P_cons_sdr) + ((1-omega)/ns)*(p_sense_sdr - EEsmin*P_cons_sdr);

        % Dinkelbach parameter update must match the CVX objective form
        q_new = num_raw / max(P_cons_sdr, 1e-30);

        % Normalised weighted EE for convergence history, clamped to [0,1]
        obj_hist(t) = omega*((EEc_sdr-EEcmin)/nc) + (1-omega)*((EEs_sdr-EEsmin)/ns);

        % ------------------------------------------------------------------
        % Iteration log
        % ------------------------------------------------------------------
        fprintf(' Iter %2d: Rate = %6.2f bps/Hz | SensGain = %6.4f | P_tot = %6.4f W | q = %6.4f | EEc = %6.4f bit/J/Hz | EEs = %6.4f \n', ...
        t, R_sum_sdr, p_sense_sdr, P_cons_sdr, q_new, EEc_sdr, EEs_sdr);
        if abs(q_new - q) <= epsilon
            obj_hist = obj_hist(1:t);
            break;
        end
        q = q_new;
    end

    % Trim trailing zeros if loop hit T_max without triggering break
    nz = find(obj_hist ~= 0, 1, 'last');
    if ~isempty(nz); obj_hist = obj_hist(1:nz); end

    [Wc_r1, Ws_r1] = gaussian_randomization( ...
        Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
        N, K, M, P_max, sigma2, Gamma_min, gamma_min, ...
        i_b, i_k, P_static, delta_outage, omega, EEcmax, EEsmax, N_rand, EEcmin, EEsmin);

    [EEc, EEs] = compute_EE(Wc_r1, Ws_r1, H_hat, theta_targets, ...
                            sigma2, i_b, i_k, P_static, N, K, M);
end
