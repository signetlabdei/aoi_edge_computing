function [cdf] = paoi_cdf_sync(xs, N, mu, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          function: av_aoi_sync                          %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the Cumulative Distribution Function of the Peak Age of        %
% Information over a given set of points in the synchronized case         %
%                                                                         %
% Inputs:                                                                 %
% -xs:  the support over which to compute the CDF [1 by M]                %
% -N:   the number of clients [scalar]                                    %
% -mu:  the service rate [scalar]                                         %
% -tau: the frame generation period [scalar]                              %
%                                                                         %
% Outputs:                                                                %
% -cdf: the Peak Age of Information CDF [1 by M]                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary variables
p_succ = success_prob_sync(N, mu, tau);
cdf = zeros(1, length(xs));

% Compute geometric part (consecutive failures)
T = floor(xs(end) / tau);
ccdf_step = zeros(1, T);
for t = 1 : T
    ccdf_step(t) = (1 - p_succ) ^ (t - 1);
end

% Compute latency part (successful packet)
for i = 1 : length(xs)
    x = xs(i);
    xf = floor(x / tau);
    if (xf >= 1)
        cdf(i) = 1 - ccdf_step(xf) * (1 - p_succ * latency_cdf_sync(x - xf * tau, N, mu, tau));
    end
end

end