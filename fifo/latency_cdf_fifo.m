function [cdf] = latency_cdf_fifo(xs, N, mu, nu, b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       function: latency_cdf_fifo                        %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the CDF of the latency in the FIFO system                      %
%                                                                         %
% Inputs:                                                                 %
% -xs:  the values over which to compute the CDF [1*X]                    %
% -N:   the number of clients in each batch [1*B]                         %
% -mu:  the service rate of the server [scalar]                           %
% -nu:  the period between one batch and the next [1*B]                   %
% -b:   the current batch index [scalar]                                  %
%                                                                         %
% Outputs:                                                                %
% -cdf: the CDF vector [1*X]                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility variables
B = length(N);
T = transition_matrix_fifo(N, mu, nu);
cdf = zeros(1, length(xs));

% Find steady state distribution for the batch
sigma = success_prob_fifo(N, b, T);
steady = steady_state(T);
states = size(T, 1);
idx = 1 : states;
valid = idx(mod(idx - 1, B) + 1 == b);
denom = 0;

% Iterate over states in the batch
for state = valid
    mult = steady(state);
    denom = denom + mult;
    cdf = cdf + mult * latency_cdf_state_fifo(xs, N, mu, nu, state, T) * success_prob_state_fifo(N, state, T);
end

% Normalize (only successful frames)
cdf = cdf / denom;
cdf = cdf / sigma;

end