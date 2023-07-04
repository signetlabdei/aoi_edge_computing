function [p] = P_X_sync(N, mu, tau, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           function: P_X_sync                            %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the final state probability in the synchronized case (i.e.,    %
% the probability that there are k clients at the end of a frame period)  %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients [scalar]                                    %
% -mu:  the service rate [scalar]                                         %
% -tau: the frame generation period [scalar]                              %
% -k: the number of remaining clients at the end of the frame [scalar]    %
%                                                                         %
% Outputs:                                                                %
% -p: the probability [scalar]                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = mu * tau;

% Consider Poisson events over adjusted time mu*tau
if (k > 0)
    p = (x) ^ (N - k) * exp(-x) / factorial(N - k);
else
    p = gammainc(x, N);
end


end