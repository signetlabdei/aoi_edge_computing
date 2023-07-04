function [E] = failure_transition_fifo(N, M, b, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    function: failure_transition_gps                     %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the transition matrix over M+1 steps in which all frames fail, %
% except for the first one, which succeeds, in the FIFO system            %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients in each batch [1*B]                         %
% -M:   the number of cycles [scalar]                                     %
% -b:   the current batch index [scalar]                                  %
% -T:   the Markov transition matrix [S*S]                                %
%                                                                         %
% Outputs:                                                                %
% -E:   the failure-conditional transition matrix over M cycles [S*S]     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility variables
B = length(N);
step_matrix = T ^ B;
states = size(T, 1);
idx = 1 : states;
valid = idx(mod(idx - 1, B) + 1 == b);

% First step: success-conditional transition matrix over a cycle
first_step = zeros(states);
for state = valid
    ms = success_prob_state_fifo(N, state, T);
    for next_state = valid
        sn = state2vec_fifo(next_state, N);
        if (sn(1) <= sum(N) - N(b))
            first_step(state, next_state) = step_matrix(state, next_state) / ms;
        else
            first_step(state, next_state) = (sum(N) - sn(1)) / N(b) / ms * step_matrix(state, next_state);
        end
    end
end

% Failure-conditional transition matrix over a cycle
empty_step = zeros(states);
for state = valid
    for next_state = valid
        sn = state2vec_fifo(next_state, N);
        if (sn(1) > sum(N) - N(b))
            empty_step(state, next_state) = (sn(1) - sum(N) + N(b)) / N(b) * step_matrix(state, next_state);
        end
    end
end

% Compute transition matrix: one success, then M failures
E = zeros(states);
empty_step = empty_step ^ M;
for state = valid
    for next_state = valid
        E(state, :) = E(state, :) + first_step(state, next_state) * empty_step(next_state, :);
    end
end

end