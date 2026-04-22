
% -------------------------------------------------------------------------
%  compute_EE : evaluate EE from (possibly rank-1) covariance matrices
% -------------------------------------------------------------------------
function [EEc, EEs, R_sum, p_sense] = compute_EE(Wc, Ws, H_hat, theta_targets, ...
                                  sigma2, i_b, i_k, P_static, N, K, M)
    sl = @(th) exp(1j*pi*((0:N-1)'-(N-1)/2)*sind(th));

    W_total = zeros(N,N);
    for j = 1:K; W_total = W_total + Wc{j}; end
    for l = 1:M; W_total = W_total + Ws{l};  end
    tot_pwr = real(trace(W_total));

    % Communication sum-rate
    R_sum = 0;
    for k = 1:K
        hk = H_hat(:,k);  hh = hk*hk';
        Tk = (1+i_k)*real(trace(W_total*hh)) ...
           + (1+i_k)*i_b*real(trace(diag(diag(W_total))*hh)) ...
           + (1+i_k)*sigma2;
        inter = 0;
        for j = 1:K
            if j ~= k; inter = inter + real(trace(Wc{j}*hh)); end
        end
        for l = 1:M; inter = inter + real(trace(Ws{l}*hh)); end
        Uk = inter ...  
           + i_k*real(trace(W_total*hh)) ...
           + (1+i_k)*i_b*real(trace(diag(diag(W_total))*hh)) ...
           + (1+i_k)*sigma2;
        SINR_k = max((Tk - Uk) / max(Uk, 1e-30), 0);
        R_sum  = R_sum + log2(1 + SINR_k);
    end

    % Sensing gain
    p_sense = 0;
    for m = 1:M
        vm = sl(theta_targets(m));
        p_sense = p_sense + real(vm' * W_total * vm) + i_b*real(trace(W_total)); %\mathbf{v}_m^H \operatorname{diag}(\mathbf{W}_{\text{total}}) \mathbf{v}_m = \operatorname{Tr}(\mathbf{W}_{\text{total}})$$
    end

    P_cons = max(P_static + (1+i_b)*tot_pwr, 1e-30);
    EEc = R_sum   / P_cons;
    EEs = p_sense / P_cons;
end