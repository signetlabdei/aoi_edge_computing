function [pi] = steady_state(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         function: steady_state                          %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the steady-state probability distribution for a Markov chain   %
% using the eigenvector method                                            %
%                                                                         %
% Inputs:                                                                 %
% -T:   the Markov transition matrix [S*S]                                %
%                                                                         %
% Outputs:                                                                %
% -pi:  the steady-state probability vector [1*S]                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute the left eigenvectors and eigenvalues
N = size(T, 1);
[~, V, D] = eig(T);
idx = 0;
% Find the eigenvector with eigenvalue 1
for i = 1 : N
    if (abs(V(i, i) - 1) < 1e-5)
        idx = i;
    end
end

% Normalize distribution
pi = D(:, idx);
pi = pi / sum(pi);

end