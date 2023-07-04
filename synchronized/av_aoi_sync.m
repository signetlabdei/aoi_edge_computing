function [aoi] = av_aoi_sync(N, mu, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          function: av_aoi_sync                          %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected Age of Information in the synchronized case       %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients [scalar]                                    %
% -mu:  the service rate [scalar]                                         %
% -tau: the frame generation period [scalar]                              %
%                                                                         %
% Outputs:                                                                %
% -aoi: the average Age of Information [scalar]                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary variables
p_succ = success_prob_sync(N, mu, tau);
x = mu * tau;

% Compute AoI
aoi = tau / 2;
for k = 0 : N - 1
    aoi = aoi + (N - k) * gammainc(x, k + 1) / (N * p_succ * mu);
end