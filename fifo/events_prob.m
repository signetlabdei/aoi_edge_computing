function [prob] = events_prob(N, mu, t, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          function: events_prob                          %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the probability of having k Poisson events in a given time t,  %
% knowing that there are N active clients and no new arrivals             %
%                                                                         %
% Inputs:                                                                 %
% -N:       the total number of clients [scalar]                          %
% -mu:      the service rate of the server [scalar]                       %
% -t:       the considered time interval [scalar]                         %
% -k:       the number of events [scalar]                                 %
%                                                                         %
% Outputs:                                                                %
% -prob:    the probability of having k events [scalar]                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


prob = 0;

if (k > N)
    return;
end

if (k == N)
    prob = gammainc(mu * t, N);
    return;
end

prob = (mu * t) .^ k .* exp(-mu * t) / factorial(k);

end