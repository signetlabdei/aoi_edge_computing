function [p_succ] = success_prob_sync(N, mu, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       function: success_prob_sync                       %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the success probability for a frame in the synchronized case   %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients [scalar]                                    %
% -mu:  the service rate [scalar]                                         %
% -tau: the frame generation period [scalar]                              %
%                                                                         %
% Outputs:                                                                %
% -p_succ: the success probability [scalar]                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_succ = 0;

% Iterate over the number of possible Poisson events
for k = 0 : N - 1
    p_succ = p_succ + P_X_sync(N, mu, tau, k) * (N - k) / N;
end

end