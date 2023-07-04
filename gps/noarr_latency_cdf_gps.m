function [prob] = noarr_latency_cdf_gps(N, mu, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     function: noarr_latency_cdf_gps                     %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the CDF of the latency for a specific client, conditional on   %
% having no new arrivals in a GPS system                                  %
%                                                                         %
% Inputs:                                                                 %
% -N:       the total number of clients [scalar]                          %
% -mu:      the service rate of the server [scalar]                       %
% -t:       the considered time interval [scalar]                         %
%                                                                         %
% Outputs:                                                                %
% -prob:    the probability of the client being served by time t [scalar] %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Incomplete gamma (all clients are served, the distribution is Erlang)
prob = gammainc(mu .* t, N);

% Loop (N Poisson events, not all clients are served)
for c = 1 : N - 1
    prob = prob + (mu .* t) .^ c .* exp(-mu .* t) / N / factorial(c - 1);
end

end