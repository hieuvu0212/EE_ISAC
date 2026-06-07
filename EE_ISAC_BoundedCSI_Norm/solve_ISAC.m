% -------------------------------------------------------------------------
%  solve_ISAC : Dinkelbach-SCA solver + Gaussian Randomization
%  Returns rank-1 physical matrices and EE values from those matrices.
%
%  Normalization uses min-max scaling:
%  nc = EEcmax - EEcmin,  ns = EEsmax - EEsmin
%  EEc and EEs are on the same [0,1] scale in the weighted objective.
% -------------------------------------------------------------------------
function [Wc_r1, Ws_r1, EEc, EEs, obj_hist, EEc_hist, EEs_hist] = ...
    solve_ISAC(H_hat, theta_targets, N, K, M, P_max, sigma2, ...
               Gamma_min, gamma_min, i_b, i_k, P_static, r_k_vec, ...
               omega, T_max, epsilon, EEcmax, EEsmax, N_rand, ...
               EEcmin, EEsmin)

    % ------------------------------------------------------------------
    % Default EEcmin / EEsmin to 0 when not supplied (Phase-1 calls)
    % ------------------------------------------------------------------
    if nargin < 21; EEcmin = 0; end
    if nargin < 22; EEsmin = 0; end

    sl    = @(th) exp(1j*pi*((0:N-1)'-(N-1)/2)*sind(th));
    log_2 = log(2);

    % Min-max normalization ranges (guard against zero range)
    nc = EEcmax - EEcmin;
    ns = EEsmax - EEsmin;

% --- Initialise beamformers via worst-case-margin feasibility SDP ------
cvx_begin quiet
    cvx_solver mosek
    variable Wc0(N,N,K) complex semidefinite
    variable Ws0(N,N,M) complex semidefinite
    variable lambda0(K) nonnegative
    variable tslack

    W_sum0   = sum(Wc0,3) + sum(Ws0,3);
    tot_pow0 = real(trace(W_sum0));

    maximize( tslack )
    subject to
        tot_pow0 <= P_max;

        % sensing, slack scaled by Gamma_min  ->  >= (1+t) Gamma_min
        for m = 1:M
            vm = sl(theta_targets(m));
            real(trace((vm*vm')*W_sum0)) + i_b*tot_pow0 >= (1+tslack)*Gamma_min;
        end

        % robust SINR (S-procedure LMI), slack scaled by the noise term
        for k = 1:K
            hk  = H_hat(:,k);
            r_k = r_k_vec(k);
            Qk = (1/gamma_min)*Wc0(:,:,k);
            for j = 1:K; if j~=k; Qk = Qk - Wc0(:,:,j); end; end
            for l = 1:M;          Qk = Qk - Ws0(:,:,l);       end
            Qk = Qk - i_k*W_sum0 - (1+i_k)*i_b*diag(diag(W_sum0));
            if r_k > 0
                sc_val = real(hk'*Qk*hk) - (1+tslack)*(1+i_k)*sigma2 - lambda0(k)*r_k^2;
                [Qk+lambda0(k)*eye(N), Qk*hk;
                 (Qk*hk)',             sc_val] == hermitian_semidefinite(N+1);
            else
                real(hk'*Qk*hk) - (1+tslack)*(1+i_k)*sigma2 >= 0;
            end
        end
cvx_end

if strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Failed') || tslack < 0
    error('solve_ISAC:infeasible', ...
          'Feasibility SDP gave t = %.3g (status: %s); infeasible at these thresholds.', ...
          tslack, cvx_status);
end

Wc = cell(1,K);  Ws = cell(1,M);
for j = 1:K; Wc{j} = (Wc0(:,:,j)+Wc0(:,:,j)')/2; end
for l = 1:M; Ws{l} = (Ws0(:,:,l)+Ws0(:,:,l)')/2; end


q        = 0;
obj_hist = zeros(T_max, 1);
EEc_hist = zeros(T_max, 1);
EEs_hist = zeros(T_max, 1);
Wc_SDR   = Wc;  Ws_SDR = Ws;

    % --- Dinkelbach-SCA main loop ----------------------------------------
    for t = 1:T_max

        % Compute U_k^(t) at current iterate
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

        % CVX SDP
        cvx_begin quiet
            cvx_solver mosek
            variable Wc_v(N,N,K) complex semidefinite
            variable Ws_v(N,N,M) complex semidefinite
            variable lambda_v(K) nonnegative
            
            

            W_sum     = sum(Wc_v,3) + sum(Ws_v,3);
            tot_power = real(trace(W_sum));

            % SCA lower bound of sum-rate
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

            % Sensing beampattern gain
            sense_gain = 0;
            for m = 1:M
                vm = sl(theta_targets(m));
                sense_gain = sense_gain + real(trace((vm*vm')*W_sum)) + i_b*tot_power;
            end

            % Power consumption expression
            P_tot_expr = P_static + (1+i_b)*tot_power;

            % -----------------------------------------------------------------
            % Dinkelbach linearised objective (CVX-compliant):
            %   maximise  omega/nc * Rk_lb + (1-omega)/ns * sense_gain
            %             - q_aug * P_tot_expr
            % -----------------------------------------------------------------
            obj_expr = (omega/nc)*(Rk_lb - EEcmin*P_tot_expr) + ...
                       ((1-omega)/ns)*(sense_gain - EEsmin*P_tot_expr) - ...
                       q * P_tot_expr;

            maximize( obj_expr )

            subject to
                % Power budget
                tot_power <= P_max;

                % Sensing threshold  forall m
                for m = 1:M
                    vm = sl(theta_targets(m));
                    real(trace((vm*vm')*W_sum)) + i_b*tot_power >= Gamma_min;
                end

                % Robust SINR LMI  forall k
                for k = 1:K
                    hk = H_hat(:,k);
                    r_k = r_k_vec(k);
                    Qk = (1/gamma_min)*Wc_v(:,:,k);
                    for j = 1:K; if j~=k; Qk=Qk-Wc_v(:,:,j); end; end
                    for l = 1:M;          Qk=Qk-Ws_v(:,:,l);       end
                    Qk = Qk - i_k*W_sum ...
                             - (1+i_k)*i_b*diag(diag(W_sum));
                    %{
                    sc_val = real(hk'*Qk*hk) ...
                             - (1+i_k)*sigma2 - lambda_v(k)*r_k^2;
                    [Qk+lambda_v(k)*eye(N),  Qk*hk ; ...
                    (Qk*hk)',               sc_val ] == hermitian_semidefinite(N+1);
                    %}
                    if r_k > 0
                     % Full S-procedure LMI with lambda
                     sc_val = real(hk'*Qk*hk) - (1+i_k)*sigma2 - lambda_v(k)*r_k^2;
                    [Qk+lambda_v(k)*eye(N), Qk*hk;
                    (Qk*hk)',              sc_val] == hermitian_semidefinite(N+1);
                    else
                    % Perfect CSI: nominal SINR constraint only
                    real(hk'*Qk*hk) - (1+i_k)*sigma2 >= 0;
                    end
                end
        cvx_end

        if strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Failed')
            warning('CVX: %s at t=%d. Using last feasible point.', ...
                    cvx_status, t);
            last = max(t-1, 1);
            obj_hist = obj_hist(1:last);
            EEc_hist = EEc_hist(1:last);   
            EEs_hist = EEs_hist(1:last);   
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

        % ------------------------------------------------------------------
        % Evaluate SDR solution metrics for logging and Dinkelbach update
        % ------------------------------------------------------------------
        [EEc_sdr, EEs_sdr, R_sum_sdr, p_sense_sdr] = ...
            compute_EE(Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
                       sigma2, i_b, i_k, P_static, N, K, M);

        % Total power consumed at this iterate
        W_tmp = zeros(N,N);
        for j=1:K; W_tmp=W_tmp+Wc_SDR{j}; end
        for l=1:M; W_tmp=W_tmp+Ws_SDR{l}; end
        P_cons_sdr = P_static + (1+i_b)*real(trace(W_tmp));

        % Raw numerator matching exactly what CVX maximises
        num_raw = (omega/nc)*(R_sum_sdr - EEcmin*P_cons_sdr) + ...
                  ((1-omega)/ns)*(p_sense_sdr - EEsmin*P_cons_sdr);

        % Dinkelbach parameter update
        q_new = num_raw / max(P_cons_sdr, 1e-30);

        % Normalised weighted EE for convergence history
        obj_hist(t) = omega*((EEc_sdr-EEcmin)/nc) + (1-omega)*((EEs_sdr-EEsmin)/ns);

        % Per-component history (raw bits/J/Hz and sensing EE)  <-- NEW
        EEc_hist(t) = omega*((EEc_sdr-EEcmin)/nc);
        EEs_hist(t) = (1-omega)*((EEs_sdr-EEsmin)/ns);

        % ------------------------------------------------------------------
        % Iteration log
        % ------------------------------------------------------------------
        fprintf([' Iter %2d: Rate=%6.2f bps/Hz | SensGain=%8.6f |' ...
                 ' P_tot=%6.4f W | q=%8.4f | EEc=%8.4f | EEs=%8.4f\n'], ...
                t, R_sum_sdr, p_sense_sdr, P_cons_sdr, q_new, EEc_sdr, EEs_sdr);

        if abs(q_new - q) <= epsilon
            obj_hist = obj_hist(1:t);
            EEc_hist = EEc_hist(1:t);   
            EEs_hist = EEs_hist(1:t);   
            break;
        end
        q = q_new;
    end

    % Trim trailing zeros if loop hit T_max without triggering break
    nz = find(obj_hist ~= 0, 1, 'last');
    if ~isempty(nz)
        obj_hist = obj_hist(1:nz);
        EEc_hist = EEc_hist(1:nz);   
        EEs_hist = EEs_hist(1:nz);   
    end

    [Wc_r1, Ws_r1] = gaussian_randomization( ...
        Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
        N, K, M, P_max, sigma2, Gamma_min, gamma_min, ...
        i_b, i_k, P_static, r_k_vec, omega, EEcmax, EEsmax, N_rand, EEcmin, EEsmin);

    [EEc, EEs] = compute_EE(Wc_r1, Ws_r1, H_hat, theta_targets, ...
                            sigma2, i_b, i_k, P_static, N, K, M);
end
