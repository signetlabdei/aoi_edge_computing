function [T] = expected_latency_gps(N, mu, nu, b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     function: expected_latency_gps                      %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected latency for the GPS system, considering only      %
% successful frames                                                       %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients in each batch [1*B]                         %
% -mu:  the service rate [scalar]                                         %
% -nu:  the period between one batch and the next [1*B]                   %
% -b:   the current batch index [scalar]                                  %
%                                                                         %
% Outputs:                                                                %
% -T:   the expected latency [scalar]                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility variables
B = length(N);
trans = transition_matrix_gps(N, mu, nu);
T = 0;

% Find steady state distribution for the batch
steady = steady_state(trans);
states = size(trans, 1);
idx = 1 : states;
valid = idx(mod(idx - 1, B) + 1 == b);
denom = 0;

% Iterate over next states
for state = valid
    mult = steady(state) * success_prob_state_gps(N, state, trans);
    denom = denom + mult;
    T = T + mult * expected_latency_state_gps(N, mu, nu, state);
end
T = T / denom;

end