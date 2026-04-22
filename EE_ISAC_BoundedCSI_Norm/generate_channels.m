% -------------------------------------------------------------------------
%  generate_channels : generate Rayleigh-fading channel matrices.
%
%  Inputs
%    N        – number of transmit antennas
%    K        – number of users
%    n_real   – number of independent realisations
%
%  Output
%    H  – N x K x n_real complex matrix  (n_real=1 returns an N x K matrix)
%
%  NOTE: Call rng(seed) in the main script BEFORE calling this function
%        to guarantee reproducibility.
% -------------------------------------------------------------------------
function H = generate_channels(N, K, n_real)

    if nargin < 3
        n_real = 1;
    end

    % Rayleigh fading: CN(0, I) normalised so E[|h_nk|^2] = 1
    H = (randn(N, K, n_real) + 1j*randn(N, K, n_real)) / sqrt(2);

    % Squeeze to N x K when only one realisation is requested
    if n_real == 1
        H = squeeze(H);
    end

end