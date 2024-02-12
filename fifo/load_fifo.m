function [rho] = load_fifo(N, mu, nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           function: load_fifo                           %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected load in the FIFO case                             %
%                                                                         %
% Inputs:                                                                 %
% -N:       the number of clients in each batch [1*B]                     %
% -mu:      the service rate of the server [scalar]                       %
% -nu:      the period between one batch and the next [1*B]               %
%                                                                         %
% Outputs:                                                                %
% -rho: the average system load [scalar]                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary variables
T = transition_matrix_fifo(N, mu, nu);
states = size(T, 1);
rho = 0;

% Compute steady-state distribution
steady = steady_state(T);

for state = 1 : states
    sv = state2vec_fifo(state, N);
    sv(1) = min(sum(N), sv(1) + N(sv(2)));
    rho = rho + steady(state) * nu(sv(2)) * load_empty(sv(1), mu, nu(sv(2)));
end
rho = rho / sum(nu) * length(N);