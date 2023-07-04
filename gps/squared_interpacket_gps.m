function [Y] = squared_interpacket_gps(N, mu, nu, b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   function: squared_interpacket_gps                     %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected value of the square of the time between two       %
% consecutive successful frames in the GPS system                         %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients in each batch [1*B]                         %
% -mu:  the service rate [scalar]                                         %
% -nu:  the period between one batch and the next [1*B]                   %
% -b:   the current batch index [scalar]                                  %
%                                                                         %
% Outputs:                                                                %
% -Y:   the expected squared interpacket time [scalar]                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary variables
B = length(N);
T_onestep = transition_matrix_gps(N, mu, nu);
Y = 0;

% Find steady state distribution for the batch
steady = steady_state(T_onestep);
states = size(T_onestep, 1);
idx = 1 : states;
valid = idx(mod(idx - 1, B) + 1 == b);

% Compute state weight multiplier
mult = zeros(1, states);
success = zeros(1, states);
for state = valid
    success(state) =  success_prob_state_gps(N, state, T_onestep);
    mult(state) = steady(state) * success(state);
end
mult = mult / sum(mult);

% Compute conditional transition matrix (failures)
T = T_onestep ^ B;
S = zeros(states);
for state = 1 : states
    for next_state = 1 : states
        ms = state2vec_gps(state, N);
        mn = state2vec_gps(next_state, N);
        c = ms(end);
        S(state, next_state) = T(state, next_state) / N(c) * mn(c);
    end
end

% Approximate diagonalization and series solution
[V, D] = eig(S + diag(randn(1, states)) * 1e-6);
for state = 1 : states
    D(state, state) = (1 + D(state, state)) / (1 - D(state, state)) ^ 3;
end
G = V * D / V;

% Iterate over state possibilities
for state = valid
    for next_state = valid
        for med_state = valid
            mx = state2vec_gps(med_state, N);
            zeta = T(state, med_state) * (N(b) - mx(b)) / success(state) / N(b);
            Y = Y + mult(state) * zeta * success(next_state) * G(med_state, next_state) * sum(nu) ^ 2;
        end
    end
end