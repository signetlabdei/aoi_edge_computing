function [sigma] = success_prob_state_gps(N, state, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   function: success_prob_state_gps                      %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the success probability for a frame arriving as part of a new  %
% batch in a given state in the GPS system                                %
%                                                                         %
% Inputs:                                                                 %
% -N:       the number of clients in each batch [1*B]                     %
% -state:   the current state [scalar]                                    %
% -T:       the Markov transition matrix [S*S]                            %
%                                                                         %
% Outputs:                                                                %
% -sigma:   the conditional success probability [scalar]                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Utility variables
B = length(N);
sigma = 0;
x = state2vec_gps(state, N);
b = x(B + 1);

% Find set of valid states
states = size(T, 1);
idx = 1 : states;
valid = idx(mod(idx - 1, B) + 1 == b);

% Find the B-step transition matrix
step_matrix = T ^ B;

% Compute success probability for each transition
for next_state = valid
    xn = state2vec_gps(next_state, N);
    sigma = sigma + (N(b) - xn(b)) / N(b) * step_matrix(state, next_state);
end