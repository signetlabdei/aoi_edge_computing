function [matrix] = transition_matrix_fifo(N, mu, nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    function: transition_matrix_fifo                     %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the transition matrix for the FIFO system (one step)           %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients in each batch [1*B]                         %
% -mu:  the service rate of the server [scalar]                           %
% -nu:  the period between one batch and the next [1*B]                   %
%                                                                         %
% Outputs:                                                                %
% -matrix: the Markov transition matrix [S*S]                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = length(N);
states = vec2state_fifo([sum(N), B], N);

matrix = zeros(states);

for state = 1 : states
    x = state2vec_fifo(state, N);
    x(1) = min(sum(N), x(1) + N(x(2)));
    for next_state = 1 : states
        xn = state2vec_fifo(next_state, N);
        % Check batch index
        if ((xn(2) == 1 && x(2) == B) || xn(2) == x(2) + 1)
            % Check if transition is possible
            if (xn(1) <= x(1))
                events = x(1) - xn(1);
                matrix(state, next_state) = events_prob(x(1), mu, nu(x(2)), events);
            end
        end
    end
end

end