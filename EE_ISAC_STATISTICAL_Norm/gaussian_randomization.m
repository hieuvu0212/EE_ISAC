% -------------------------------------------------------------------------
%  gaussian_randomization : Algorithm 2
% -------------------------------------------------------------------------
function [Wc_r1, Ws_r1] = gaussian_randomization( ...
        Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
        N, K, M, P_max, sigma2, Gamma_min, gamma_min, ...
        i_b, i_k, P_static, delta_outage, omega, EEcmax, EEsmax, N_rand, EEcmin, EEsmin)

    sl = @(th) exp(1j*pi*((0:N-1)'-(N-1)/2)*sind(th));

    % Eigen-decompose each SDR matrix
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

    best_obj = -inf;
    Wc_r1    = Wc_SDR;   % fallback = SDR solution
    Ws_r1    = Ws_SDR;
    nc = max(EEcmax - EEcmin, 1e-9);
    ns = max(EEsmax - EEsmin, 1e-9);
    [EEc_sdr0, EEs_sdr0] = compute_EE(Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
                        sigma2, i_b, i_k, P_static, N, K, M);
    best_obj = omega*((EEc_sdr0 - EEcmin)/nc) + (1-omega)*((EEs_sdr0 - EEsmin)/ns);    
    for n_trial = 1:N_rand

        % Generate candidate vectors
        wc = cell(1,K);  ws = cell(1,M);
        for j = 1:K
            z = (randn(N,1)+1j*randn(N,1))/sqrt(2);
            wc{j} = Uc{j} * Sc_half{j} * z;
        end
        for l = 1:M
            z = (randn(N,1)+1j*randn(N,1))/sqrt(2);
            ws{l} = Us{l} * Ss_half{l} * z;
        end

        % Scale to satisfy total power constraint
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

        % Aggregate
        W_tot_n = zeros(N,N);
        for j = 1:K; W_tot_n = W_tot_n + Wc_n{j}; end
        for l = 1:M; W_tot_n = W_tot_n + Ws_n{l};  end
        tr_tot  = real(trace(W_tot_n));

        % Check sensing threshold
        ok = true;
        for m = 1:M
            vm = sl(theta_targets(m));
            if real(vm'*W_tot_n*vm) + i_b*tr_tot < Gamma_min
                ok = false; break;
            end
        end
        if ~ok; continue; end

        % Check SINR at nominal channel (necessary condition from LMI)
        for k = 1:K
            hk = H_hat(:,k);
            Qk = (1/gamma_min)*Wc_n{k};
            for j = 1:K; if j~=k; Qk=Qk-Wc_n{j}; end; end
            for l = 1:M;          Qk=Qk-Ws_n{l};       end
            Qk = Qk - i_k*W_tot_n - (1+i_k)*i_b*diag(diag(W_tot_n));
            if real(hk'*Qk*hk) - (1+i_k)*sigma2 < 0
                ok = false; break;
            end
        end
        if ~ok; continue; end

        % Evaluate weighted objective
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
