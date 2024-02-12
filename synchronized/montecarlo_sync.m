function [paoi, latency, p_succ, av_aoi, load] = montecarlo_sync(N, mu, tau, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        function: montecarlo_sync                        %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Monte Carlo simulation for synchronized deterministic arrivals (one     %
% batch with N clients) to a GPS server with preemption in service        %
%                                                                         %
% Inputs:                                                                 %
% -N:       the number of clients [scalar]                                %
% -mu:      the service rate [scalar]                                     %
% -tau:     the frame generation period [scalar]                          %
% -L:       the number of Monte Carlo samples [scalar]                    %
%                                                                         %
% Outputs:                                                                %
% -paoi:    the Peak AoI results over time                                %
% -latency: the latency results over time (only successful frames)        %
% -p_succ:  the frame success probability                                 %
% -av_aoi:  the average AoI in the simulation                             %
% -load:    the average system load in the simulation                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility vectors and variable definition
paoi = -ones(1, L);
latency = -ones(1, L);
completed = -ones(1, L);
latest = -1;
av_aoi = 0;
aoi = 0;
load = 0;

% Main Monte Carlo loop
for l = 1 : L
    % Compute service times
    delays = cumsum(exprnd(1 / mu, 1, N));
    completed(l) = length(find(delays < tau));
    % The order of service is random (both in GPS and D^N/M/1)
    index = randi(N);
    % Check if the frame is successful
    if (index <= completed(l))
        latency(l) = delays(index);
        % Update average AoI
        av_aoi = av_aoi + (aoi + latency(l) / 2) * latency(l);
        av_aoi = av_aoi + tau ^ 2 / 2 - latency(l) ^ 2 / 2;
        aoi = tau;
        if (latest > 0)
            % Update Peak AoI
            paoi(l) = tau * (l - latest) + latency(l);
        end
        latest = l;
    else
        % Update average AoI
        av_aoi = av_aoi + (aoi + tau / 2) * tau;
        aoi = aoi + tau;
    end
    % Compute load
    if (delays(N) <= tau)
        load = load + delays(N);
    else
        load = load + tau;
    end
end

% Remove unsuccessful frames and compute success probability
av_aoi  = av_aoi / (tau * L);
paoi = paoi(paoi > 0);
latency = latency(latency > 0);
p_succ = length(latency) / L;
load = load / L / tau;

end