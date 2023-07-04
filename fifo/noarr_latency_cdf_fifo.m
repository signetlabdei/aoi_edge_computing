function [prob] = noarr_latency_cdf_dm1(x, N, b, mu, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    function: noarr_latency_cdf_fifo                     %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the CDF of the latency for a specific client, conditional on   %
% having no new arrivals in a FIFO system                                 %
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

% Get the number of packets queued behind the client
q = behind_queue_fifo(N, b, x(2));

% Utility variables
prob = zeros(1, length(t));
remaining = min(N(b), x(1) - q);

% If the client is still active, compute CDF
if (x(1) > q)
    cmin = max(0, x(1) - q - N(b));
    % The batch is partially served
    for c = cmin + 1 : x(1) - q
        final = x(1) - c;
        prob = prob + events_prob(x(1), mu, t, c) * (remaining - (final - q)) / remaining;
    end
    % The whole batch is served
    for c = x(1) - q + 1 : x(1)
        prob = prob + events_prob(x(1), mu, t, c);
    end
end

end