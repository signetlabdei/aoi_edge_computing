function [rho] = load_gps(N, mu, nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           function: load_gps                            %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected load in the GPS case                              %
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
T = transition_matrix_gps(N, mu, nu);
states = size(T, 1);
rho = 0;
B = length(nu);

% Compute steady-state distribution
steady = steady_state(T);

for state = 1 : states
    sv = state2vec_gps(state, N);
    b = sv(B + 1);
    sv(b) = N(b);
    rho = rho + steady(state) * nu(b) * load_empty(sum(sv(1 : B)), mu, nu(b));
end
rho = rho / sum(nu) * B;