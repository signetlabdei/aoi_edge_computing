function [T] = expected_latency_state_fifo(N, mu, nu, state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  function: expected_latency_state_fifo                  %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected latency for the GPS system, starting from a given %
% state and considering only successful frames                            %
%                                                                         %
% Inputs:                                                                 %
% -N:       the number of clients in each batch [1*B]                     %
% -mu:      the service rate [scalar]                                     %
% -nu:      the period between one batch and the next [1*B]               %
% -state:   the current state [scalar]                                    %
%                                                                         %
% Outputs:                                                                %
% -T:       the expected latency [scalar]                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility variables
B = length(N);
T = 0;
trans = transition_matrix_fifo(N, mu, nu);
sv = state2vec_fifo(state, N);
b = sv(2);
sigma = success_prob_state_fifo(N, state, trans);
states = size(trans, 1);
idx = 1 : states;
nu_b = [nu(b : end), nu(1 : b - 1)];
idx_b = [b : B, 1 : b];
N_b = [N(b + 1 : end), N(1 : b)];

% Iterate over possible states
for jumps = 0 : B - 1
    T = T + nu_b(jumps + 1);
    valid_step = idx(mod(idx - 1, B) + 1 == idx_b(jumps + 1));
    step_matrix = trans ^ jumps;
    % First interval: the frame has just arrived
    if (jumps == 0)
        sv(1) = min(sv(1) + N(b), sum(N));
        Q = sv(1) - N(b);
        X = sv(1);
        T = T - nu_b(1) / sigma;
        mult = 1 / (mu * sigma);
        % Frames queued ahead
        for k = 0 : Q
            T = T + mult * gammainc(mu * nu_b(1), k + 1);
        end
        % Frames in the same batch
        for k = Q + 1 : X - 1
            mult = (X - k) / mu / X / sigma;
            T = T + mult * gammainc(mu * nu_b(1), k + 1);
        end
    else
        % Next steps: iterate over states
        bn = mod(b + jumps, B);
        if (bn == 0)
            bn = B;
        end
        Qb = behind_queue_fifo(N, b, bn);
        for next_state = valid_step
            sn = state2vec_fifo(next_state, N);
            sn(1) = min(sn(1) + N_b(jumps), sum(N));
            Ts = - nu_b(jumps + 1) / sigma;
            if (sn(1) > Qb)
                X = sn(1) - Qb;
                Q = max(0, X - N(b));
                % Multiply by the remaining frames in the batch (some may
                % already have been served)
                mult = min(1, X / N(b)) / (mu * sigma);
                % Frames queued ahead
                for k = 0 : Q
                    Ts = Ts + mult * gammainc(mu * nu_b(jumps + 1), k + 1);
                end
                % Frames in the same batch
                for k = Q + 1 : X - 1
                    mult = min(1, X / N(b)) * (X - k) / mu / N(b) / sigma;
                    Ts = Ts + mult * gammainc(mu * nu_b(jumps + 1), k + 1);
                end
            end
            T = T + step_matrix(state, next_state) * Ts;
        end
    end
end

end