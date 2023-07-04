%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            script: run_fifo                             %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Run a Monte Carlo simulation of the general FIFO system and compare     %
% the results to the theoretical CDFs                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

% Simulation parameters
N = [1,1,2,1,1];            % Number of clients in each batch
mu = 4;                     % Service rate
nu = ones(1, 5) * 1 / 6;    % Inter-batch interval
nu(1) = 1 / 3;
L = 1E6;                    % Number of Monte Carlo cycles
client = 1;                 % Client batch index

% CDF support and utility variables
tau = sum(nu);
xs = 0 : 0.0001 : tau;
ds = 0 : 0.0001 : 10 * tau;

% Compute theoretical CDFs
latency_th_cdf = latency_cdf_fifo(xs, N, mu, nu, client);
paoi_th_cdf = paoi_cdf_fifo(0.0001, 10, N, mu, nu, client);
T_th = transition_matrix_fifo(N, mu, nu);
p_succ_th = success_prob_fifo(N, client, T_th);

% Run Monte Carlo simulation and compute empirical CDFs
[paoi_mc, latency_mc, p_succ_mc, av_aoi_mc] = montecarlo_fifo(N, mu, nu, client, L);
latency_mc_cdf = histcounts(latency_mc, xs);
latency_mc_cdf = cumsum(latency_mc_cdf) / sum(latency_mc_cdf);
paoi_mc_cdf = histcounts(paoi_mc, ds);
paoi_mc_cdf = cumsum(paoi_mc_cdf) / length(paoi_mc);

% Create latency and PAoI CDF plots
f1 = figure(1);
plot(xs, latency_th_cdf * p_succ_th)
hold on
plot(xs, [0, latency_mc_cdf] * p_succ_mc)
xlabel('Time')
ylabel('Latency CDF')
legend('Analytical', 'Monte Carlo')
f2 = figure(2);
plot(ds, paoi_th_cdf)
hold on
plot(ds, [0, paoi_mc_cdf])
xlabel('Time')
ylabel('PAoI CDF')
legend('Analytical', 'Monte Carlo')

% Print average AoI (empirical and theoretical)
av_aoi_mc
Y = expected_aoi_fifo(N, mu, nu, client)