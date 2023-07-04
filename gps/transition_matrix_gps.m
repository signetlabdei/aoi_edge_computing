function [matrix] = transition_matrix_gps(N, mu, nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     function: transition_matrix_gps                     %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the transition matrix for the GPS system (one step)            %
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

% Initialize utility variables
B = length(N);
states = vec2state_gps([N, B], N);
matrix = zeros(states);

% Iterate over the states
for state = 1 : states
    x = state2vec_gps(state, N);
    x(x(B + 1)) = N(x(B + 1));
    for next_state = 1 : states
        xn = state2vec_gps(next_state, N);
        % Check batch index
        if ((xn(B + 1) == 1 && x(B + 1) == B) || xn(B + 1) == x(B + 1) + 1)
            % Check if transition is possible
            possible = true;
            for b = 1 : B
                if (xn(b) > x(b))
                    possible = false;
                end
            end
            % Consider event probability for the specific state transition
            if (possible)
                events = sum(x(1 : B) - xn(1 : B));
                prob = events_prob(sum(x(1 : B)), mu, nu(x(B + 1)), events);
                for b = 1 : B
                    prob = prob * nchoosek(x(b), x(b) - xn(b));
                end
                matrix(state, next_state) = prob / nchoosek(sum(x(1 : B)), events);
            end
        end
    end
end

end