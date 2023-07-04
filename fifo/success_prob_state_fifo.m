function [p_succ] = success_prob_state_fifo(N, state, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    function: success_prob_state_fifo                    %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the success probability for a frame arriving as part of a new  %
% batch in a given state in the FIFO system                               %
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
p_succ = 0;
x = state2vec_fifo(state, N);
b = x(2);

% Find set of valid states
states = size(T, 1);
idx = 1 : states;
valid = idx(mod(idx - 1, B) + 1 == b);

step_matrix = T ^ B;

for next_state = valid
    xn = state2vec_fifo(next_state, N);
    xb = xn(1) - (sum(N) - N(b));
    if (xb <= 0)
        p_succ = p_succ + step_matrix(state, next_state);
    else
        p_succ = p_succ + (N(b) - xb) / N(b) * step_matrix(state, next_state);
    end
end