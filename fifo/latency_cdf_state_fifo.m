function [cdf] = latency_cdf_state_fifo(xs, N, mu, nu, state, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    function: latency_cdf_state_fifo                     %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the CDF of the latency starting in a specific state in a FIFO  %
% system                                                                  %
%                                                                         %
% Inputs:                                                                 %
% -xs:      the values over which to compute the CDF [1*X]                %
% -N:       the number of clients in each batch [1*B]                     %
% -mu:      the service rate of the server [scalar]                       %
% -nu:      the period between one batch and the next [1*B]               %
% -state:   the current state [scalar]                                    %
% -T:       the Markov transition matrix [S*S]                            %
%                                                                         %
% Outputs:                                                                %
% -cdf:     the CDF vector [1*X]                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% utility variables
B = length(N);
cdf = zeros(1, length(xs));

% Find steady state distribution for the batch
sv = state2vec_fifo(state, N);
b = sv(2);
sigma = success_prob_state_fifo(N, state, T);
states = size(T, 1);
idx = 1 : states;
nu_b = [nu(b : end), nu(1 : b - 1)];
idx_b = [b : B, 1 : b];

% Compute indices of batch in the support vector
sv(1) = min(sum(N), sv(1) + N(b));
low_idx = 1;
high_idx = find(xs >= nu(b), 1);
if (isempty(high_idx))
    high_idx = length(xs);
end
done = false;
jumps = 0;
while(~done)
    % Compute transition matrix for next jump and reachable states
    valid_step = idx(mod(idx - 1, B) + 1 == idx_b(jumps + 1));
    step_matrix = T ^ jumps;
    red_xs = xs(low_idx : high_idx) - sum(nu_b(1 : jumps));
    if (jumps == 0)
        % First step
        prob = noarr_latency_cdf_fifo(sv, N, b, mu, red_xs);
        cdf(low_idx : high_idx) = cdf(low_idx : high_idx) + prob;
    else
        % Next steps: iterate over states
        for next_state = valid_step
            sn = state2vec_fifo(next_state, N);
            q = behind_queue_fifo(N, b, sn(2));
            sn(1) = min(sum(N), sn(1) + N(sn(2)));
            % All frames in the batch have already been served
            if (sn(1) <= q)
                cdf(low_idx : high_idx) = cdf(low_idx : high_idx) + step_matrix(state, next_state);
            else
                mult = 1;
                % Some frames in the batch have already been served
                if (sn(1) < q + N(b))
                    cdf(low_idx : high_idx) = cdf(low_idx : high_idx) + step_matrix(state, next_state) * (N(b) - sn(1) + q) / N(b);
                    mult = (sn(1) - q) / N(b);
                end
                prob = step_matrix(state, next_state) * noarr_latency_cdf_fifo(sn, N, b, mu, red_xs) * mult;
                cdf(low_idx : high_idx) = cdf(low_idx : high_idx) + prob;
            end
        end
    end
    % Go to next hop
    jumps = jumps + 1;
    low_idx = high_idx + 1;
    if (low_idx >= length(xs))
        done = true;
    else
        high_idx = find(xs >= sum(nu_b(1 : jumps + 1)), 1);
        if (isempty(high_idx))
            high_idx = length(xs);
        end
    end
end

% Normalize (only successful frames
cdf = cdf / sigma;

end