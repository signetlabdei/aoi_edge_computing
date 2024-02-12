function [rho] = load_empty(N, mu, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          function: load_empty                           %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected load in the GPS case (in between states)          %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients [scalar]                                    %
% -mu:  the service rate [scalar]                                         %
% -tau: the frame generation period [scalar]                              %
%                                                                         %
% Outputs:                                                                %
% -rho: the average system load [scalar]                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary variables
x = mu * tau;

% Full system (some frames are still active at the end of the period)
rho = 1 - gammainc(x, N);

% Emptied system (all frames end by time tau)
rho = rho + gammainc(x, N + 1) * N / x;