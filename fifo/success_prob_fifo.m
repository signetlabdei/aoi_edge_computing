function [sigma] = success_prob_fifo(N, b, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       function: success_prob_fifo                       %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the success probability for a frame arriving as part of a new  %
% batch with a given index in the FIFO system                             %
%                                                                         %
% Inputs:                                                                 %
% -N:       the number of clients in each batch [1*B]                     %
% -b:       the current batch index [scalar]                              %
% -T:       the Markov transition matrix [S*S]                            %
%                                                                         %
% Outputs:                                                                %
% -sigma:   the conditional success probability [scalar]                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility variables
B = length(N);
sigma = 0;

% Find steady state distribution for the batch
steady = steady_state(T);
states = size(T, 1);
idx = 1 : states;
valid = idx(mod(idx - 1, B) + 1 == b);
other = idx(mod(idx - 1, B) + 1 ~= b);
steady(other) = 0;
steady = steady / sum(steady);

% Iterate over states
for state = valid
    if (steady(state) > 0)
        sigma = sigma + steady(state) * success_prob_state_fifo(N, state, T);
    end
end