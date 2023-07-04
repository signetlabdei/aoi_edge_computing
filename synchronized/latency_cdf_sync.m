function [cdf] = latency_cdf_sync(xs, N, mu, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       function: latency_cdf_sync                        %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the Cumulative Distribution Function of the latency of a       %
% successful frame over a given set of points in the synchronized case    %
%                                                                         %
% Inputs:                                                                 %
% -xs:  the support over which to compute the CDF [1 by M]                %
% -N:   the number of clients [scalar]                                    %
% -mu:  the service rate [scalar]                                         %
% -tau: the frame generation period [scalar]                              %
%                                                                         %
% Outputs:                                                                %
% -cdf: the conditional latency CDF [1 by M]                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_succ = success_prob_sync(N, mu, tau);

% If there is only one client, latency is exponential
if (N == 1)
    cdf = expcdf(xs, 1 / mu);
else
    % If there is more than one client, latency is more complex
    cdf = gammainc(mu * xs, N);
    for k = 1 : N - 1
        cdf = cdf + (mu .* xs) .^ k .* exp(-mu .* xs) / N / factorial(k - 1);
    end
end

% Condition: successful frames
cdf = cdf / p_succ;

end