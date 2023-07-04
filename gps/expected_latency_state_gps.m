function [T] = expected_latency_state_gps(N, mu, nu, state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  function: expected_latency_state_gps                   %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected latency for the FIFO system starting from the     %
% given state, considering only successful frames                         %
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
trans = transition_matrix_gps(N, mu, nu);
sv = state2vec_gps(state, N);
b = sv(B + 1);
p_succ = success_prob_state_gps(N, state, trans);
states = size(trans, 1);
idx = 1 : states;
nu_b = [nu(b : end), nu(1 : b - 1)];
idx_b = [b : B, 1 : b];

% Iterate over possible states
for jumps = 0 : B - 1
    T = T + nu_b(jumps + 1);
    valid_step = idx(mod(idx - 1, B) + 1 == idx_b(jumps + 1));
    step_matrix = trans ^ jumps;
    % First interval: the frame has just arrived
    if (jumps == 0)
        sv(b) = N(b);
        X = sum(sv(1 : B));
        T = T - nu_b(1) / p_succ;
        for k = 0 : X - 1
            mult = (X - k) / mu / X / p_succ;
            T = T + mult * gammainc(mu * nu_b(1), k + 1);
        end
    else
        % Next steps: iterate over states
        for next_state = valid_step
            sn = state2vec_gps(next_state, N);
            sn(idx_b(jumps + 1)) = N(idx_b(jumps + 1));
            X = sum(sn(1 : B));
            Ts = - nu_b(jumps + 1) / p_succ;
            for k = 0 : X - 1
                % Multiplier: s_b/N_b (the others have already been served) 
                mult = (X - k) / mu / X / p_succ;
                Ts = Ts + mult * sn(b) / N(b) * gammainc(mu * nu_b(jumps + 1), k + 1);
            end
            T = T + step_matrix(state, next_state) * Ts;
        end
    end
end

end