% -------------------------------------------------------------------------
%  parameters : returns a struct containing all system parameters.
%  Edit values here; the rest of the code reads from this struct.
% -------------------------------------------------------------------------
function p = parameters()

    % Dimensions
    p.N = 8;     % Number of transmit antennas
    p.K = 6;     % Number of communication users
    p.M = 4;     % Number of sensing targets
    p.omega = 0.5;
    % Power / noise
    p.P_max   = db2pow(30)  * 1e-3;   % Max transmit power       [30 dBm -> W]
    p.sigma2  = db2pow(0) * 1e-3;   % Noise power per user    [0 dBm -> W]
    p.P_static = 1e-3;                % Static circuit power          [W]

    % QoS thresholds
    p.Gamma_min = db2pow(20) * 1e-3;  % Min sensing power threshold [20 dBm -> W]
    p.gamma_min = db2pow(5);          % Min SINR [5 dB -> linear]

    % Hardware impairments
    p.i_b = 0.01;   % Transmitter distortion coefficient
    p.i_k = 0.02;   % Receiver  distortion coefficient

    % Sensing target angles [degrees]
    p.theta_targets = [-30, -10, 10, 30];

    % CSI uncertainty
    p.delta  = 0.01;  % The chanel uncertainty level
    p.delta_outage = 0.05; % Outage probability (Pr{ SINR_k < gamma_min } <= delta_outage) 
    p.xi_k_vec = zeros(1,p.K);
    p.r_k_vec = zeros(1, p.K);


    % Algorithm settings
    p.T_max   = 30;      % Max Dinkelbach/SCA iterations
    p.epsilon = 1e-3;    % Convergence tolerance
    p.N_rand  = 1000;    % Gaussian-randomisation trials
  
    % Monte-Carlo
    p.N_MC = 1;         % Channel realisations per antenna size

end
