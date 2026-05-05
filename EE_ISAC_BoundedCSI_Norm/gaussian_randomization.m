% -------------------------------------------------------------------------
%  gaussian_randomization : Algorithm 2  (Bounded CSI version)
%
%  FIX: Randomized rank-1 candidates are now checked against the full
%  S-procedure robust SINR LMI before being accepted.
%
%  The bounded-CSI robust SINR constraint for user k is the LMI:
%
%    S_k(lambda_k) = [ Qk + lambda_k*I ,   Qk*h_hat_k       ]  >= 0
%                   [ (Qk*h_hat_k)^H   ,   sc_k - lambda_k*r_k^2 ]
%
%  where  Qk  =  (1/gamma_min)*Wck - sum_{j~=k} Wcj - sum_l Wsl
%                  - i_k*W_total - (1+i_k)*i_b*diag(diag(W_total))
%         sc_k = real(h_hat_k' * Qk * h_hat_k) - (1+i_k)*sigma2
%
%  For a fixed Qk (rank-1 candidate), we find the MINIMUM lambda_k >= 0
%  that makes S_k PSD by a scalar bisection, then accept if that
%  minimum lambda_k exists (i.e. the LMI is feasible for some lambda >= 0).
%
%  Analytically, S_k >= 0  iff:
%    (A)  Qk + lambda_k*I >= 0   (i.e. lambda_k >= -lambda_min(Qk))
%    (B)  Schur complement: sc_k - lambda_k*r_k^2
%           - (Qk*h_hat_k)^H * (Qk + lambda_k*I)^{-1} * (Qk*h_hat_k) >= 0
%  We sweep lambda_k to find the smallest feasible value.
%
%  CALLERS: solve_ISAC signature is unchanged (r_k_vec already present).
% -------------------------------------------------------------------------
function [Wc_r1, Ws_r1] = gaussian_randomization( ...
        Wc_SDR, Ws_SDR, H_hat, theta_targets, ...
        N, K, M, P_max, sigma2, Gamma_min, gamma_min, ...
        i_b, i_k, P_static, r_k_vec, omega, EEcmax, EEsmax, N_rand, EEcmin, EEsmin)

    sl = @(th) exp(1j*pi*((0:N-1)'-(N-1)/2)*sind(th));

    % Bisection parameters for lambda search
    N_BISECT   = 40;        % bisection iterations (40 -> ~1e-12 relative error)
    LAM_MAX    = 1e6;       % upper search bound for lambda_k

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
        % GATE 2: S-procedure robust SINR LMI for all users
        %
        %  For each k, check whether there exists lambda_k >= 0 such that
        %  the (N+1) x (N+1) Hermitian block matrix
        %
        %    S_k = [ Qk + lambda_k*I ,  Qk*hk        ]  >= 0
        %          [ hk^H*Qk         ,  sc_k - lam*rk^2 ]
        %
        %  We use a bisection over lambda_k:
        %    - Lower bound:  max(0, -lambda_min(Qk))  [ensures Qk+lam*I >= 0]
        %    - Upper bound:  LAM_MAX
        %    - Schur complement of bottom-right element:
        %        f(lam) = sc_k - lam*rk^2 - hk^H*Qk*(Qk+lam*I)^{-1}*Qk*hk
        %      S_k >= 0  iff  Qk+lam*I >= 0  AND  f(lam) >= 0.
        % =================================================================
        for k = 1:K
            hk = H_hat(:,k);
            r_k = r_k_vec(k);

            % Build Qk from rank-1 candidate
            Qk = (1/gamma_min) * Wc_n{k};
            for j = 1:K
                if j ~= k; Qk = Qk - Wc_n{j}; end
            end
            for l = 1:M
                Qk = Qk - Ws_n{l};
            end
            Qk = Qk - i_k * W_tot_n ...
                     - (1+i_k)*i_b * diag(diag(W_tot_n));
            Qk = (Qk + Qk') / 2;   % symmetrise for numerical safety

            sc_k = real(hk' * Qk * hk) - (1+i_k)*sigma2;

            % Minimum lambda needed to make Qk + lambda*I >= 0
            lam_min_needed = max(0, -min(real(eig(Qk))));

            if r_k < 1e-12
                % Perfect CSI: degenerate to deterministic check
                % sc_k >= 0  and Qk >= 0  (or standard SINR check)
                if sc_k < -1e-9
                    ok = false; break;
                end
                continue;
            end

            % Bisection: find lambda in [lam_min_needed, LAM_MAX]
            % such that the Schur complement f(lambda) >= 0.
            lam_lo = lam_min_needed;
            lam_hi = LAM_MAX;

            % Quick check: if even at lam_hi the LMI fails, reject
            if ~schur_ok(Qk, hk, sc_k, r_k, lam_hi, N)
                ok = false; break;
            end

            % If lam_lo already works, no bisection needed
            if ~schur_ok(Qk, hk, sc_k, r_k, lam_lo, N)
                % Bisect to find smallest feasible lambda
                feasible_found = false;
                for b = 1:N_BISECT
                    lam_mid = (lam_lo + lam_hi) / 2;
                    if schur_ok(Qk, hk, sc_k, r_k, lam_mid, N)
                        lam_hi = lam_mid;
                        feasible_found = true;
                    else
                        lam_lo = lam_mid;
                    end
                end
                if ~feasible_found
                    ok = false; break;
                end
            end
            % If we reach here, a feasible lambda was found for user k
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


% =========================================================================
%  HELPER: check Schur complement of S_k(lambda) >= 0
%  Returns true if S_k(lambda) is PSD (to within numerical tolerance).
% =========================================================================
function result = schur_ok(Qk, hk, sc_k, r_k, lambda, N)
    TOL = -1e-8;   % small negative slack for numerical noise

    A = Qk + lambda * eye(N);
    A = (A + A') / 2;   % force symmetry

    % Check A >= 0 first (necessary for S_k >= 0)
    if min(real(eig(A))) < TOL
        result = false;
        return;
    end

    % Schur complement of bottom-right block:
    %   f = sc_k - lambda*r_k^2 - hk^H * A^{-1} * Qk * hk
    % Use linear solve for numerical stability
    Qk_hk = Qk * hk;
    f = sc_k - lambda*r_k^2 - real(Qk_hk' * (A \ Qk_hk));

    result = (f >= TOL);
end
