function [cdf] = paoi_cdf_gps(step, M, N, mu, nu, b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         function: paoi_cdf_gps                          %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the CDF of the PAoI, over a support with a given granularity   %
% and over M cycle periods, in the GPS system                             %
%                                                                         %
% Inputs:                                                                 %
% -step:    the minimum step size for the CDF support [scalar]            % 
% -N:       the number of clients in each batch [1*B]                     %
% -M:       the number of cycles [scalar]                                 %
% -mu:      the service rate [scalar]                                     %
% -nu:      the period between one batch and the next [1*B]               %
% -b:       the current batch index [scalar]                              %
%                                                                         %
% Outputs:                                                                %
% -cdf:     the CDF of the PAoI [1*(M/step)]                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility variables
B = length(N);
tau = sum(nu);
T = transition_matrix_gps(N, mu, nu);
denom = 0;
ts = ceil(tau / step);

% Find steady state distribution for the batch
steady = steady_state(T);
states = size(T, 1);
idx = 1 : states;
valid = idx(mod(idx - 1, B) + 1 == b);
step_matrix = T ^ B;
xs = (0 : ts - 1) * step;
cdf = zeros(1, ts * M + 1);

% Compute success probability and latency CDF for all states
success = zeros(1, states);
lat_cdf = zeros(states, ts);
for state = valid
    success(state) = success_prob_state_gps(N, state, T);
    lat_cdf(state, :) = latency_cdf_state_gps(xs, N, mu, nu, state, T) * success(state);
end

% Compute transition matrices from 1 to M-1 steps
empty_trans = zeros(M - 1, states, states);
for m = 1 : M - 2
    empty_trans(m, :, :) = failure_transition_gps(N, m, b, T);
end

% Iterate over all states
for state = valid
    mult = steady(state) * success(state);
    denom = denom + mult;
    low_idx = ts + 2;
    high_idx = ts * 2 + 1;
    % First step: two consecutive successes
    for next_state = valid
        sn = state2vec_gps(next_state, N);
        if (sn(b) < N(b))
            prob_succ = (N(b) - sn(b)) / N(b) * step_matrix(state, next_state) / success(state);
            cdf(low_idx : high_idx) = cdf(low_idx : high_idx) + prob_succ * mult * lat_cdf(next_state, :);
        end
    end
    % Next steps: success, then M-2 failures, then another success
    for jumps = 1 : M - 2
        low_idx = (jumps + 1) * ts + 2;
        high_idx = (jumps + 2) * ts + 1;
        for next_state = valid
            prob_succ = mult * empty_trans(jumps, state, next_state);
            cdf(low_idx : high_idx) = cdf(low_idx : high_idx) + prob_succ * lat_cdf(next_state, :);
        end
    end
end

% Sum previous successes to CDF
for jumps = 1 : M - 2
    low_idx = (jumps + 1) * ts + 2;
    high_idx = (jumps + 2) * ts + 1;
    p0 = cdf(low_idx - 1);
    cdf(low_idx : high_idx) = cdf(low_idx : high_idx) + p0;
end

% Normalize over conditional state distribution
cdf = cdf / denom;

end