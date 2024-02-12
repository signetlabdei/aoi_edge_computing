%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            script: run_sync                             %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Run a Monte Carlo simulation of the synchronized GPS system and compare %
% the results to the theoretical CDFs                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

% Simulation parameters
N = 5;                          % Number of clients
mu = 4;                         % Service rate
tau = 1.1;                        % Inter-frame period
L = 1E6;                        % Number of Monte Carlo cycles

% CDF support
xs = 0 : 0.0001 : tau;
ds = 0 : 0.0001 : 10 * tau;


% Run Monte Carlo simulation and compute empirical CDFs
[paoi_mc, latency_mc, p_succ_mc, av_aoi_mc, load_mc] = montecarlo_sync(N, mu, tau, L);
latency_mc_cdf = histcounts(latency_mc, xs);
latency_mc_cdf = cumsum(latency_mc_cdf) / sum(latency_mc_cdf);
paoi_mc_cdf = histcounts(paoi_mc, ds);
paoi_mc_cdf = cumsum(paoi_mc_cdf) / length(paoi_mc);

% Compute theoretical CDFs
latency_th_cdf = latency_cdf_sync(xs, N, mu, tau);
paoi_th_cdf = paoi_cdf_sync(ds, N, mu, tau);


% Create latency and PAoI CDF plots
f1 = figure(1);
plot(xs, latency_th_cdf * success_prob_sync(N, mu, tau))
hold on
plot(xs, [0, latency_mc_cdf] * p_succ_mc)
xlabel('Time')
ylabel('Latency CDF')
f2 = figure(2);
plot(ds, paoi_th_cdf)
hold on
plot(ds, [0, paoi_mc_cdf])
xlabel('Time')
ylabel('PAoI CDF')

% Print average AoI (empirical and theoretical)
av_aoi_mc
av_aoi_sync(N, mu, tau)

% Print average load (empirical and theoretical)
load_mc
load_sync(N, mu, tau)