function [paoi, latency, p_succ, av_aoi] = montecarlo_gps(N, mu, nu, client, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        function: montecarlo_gps                         %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Monte Carlo simulation for generalized deterministic arrivals (B        %
% batches with N_b clients) to a GPS server with preemption in service    %
%                                                                         %
% Inputs:                                                                 %
% -N:       the number of clients [1*B]                                   %
% -mu:      the service rate [scalar]                                     %
% -nu:      the period between one batch and the next [1*B]               %
% -client:  the batch index [scalar]                                      %
% -L:       the number of Monte Carlo samples [scalar]                    %
%                                                                         %
% Outputs:                                                                %
% -paoi:    the Peak AoI results over time                                %
% -latency: the latency results over time (only successful frames)        %
% -p_succ:  the frame success probability                                 %
% -av_aoi:  the average AoI in the simulation                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Utility vectors and variable definition
B = length(N);
if (B > 1)
    nu = [nu(client : end), nu(1 : client - 1)];
    N = [N(client : end), N(1 : client - 1)];
end
state = zeros(1, B);
p_succ = 0;
tau = sum(nu);
paoi = -ones(1, L);
latency = -ones(1, L);
completed = -ones(B, L);
latest = -1;
av_aoi = 0;
aoi = 0;
done = 0;
running_latency = 0;

% Main simulation loop
for l = 1 : L
    % Batch index
    for b = 1 : B
        % Check if selected client has been served
        if (b == 1)
            if (l > 1)
                if (done == 0)
                    % Selected client has not been served: update AoI
                    av_aoi = av_aoi + tau * (aoi + tau / 2);
                    aoi = aoi + tau;
                else
                    % Selected client has been served: update AoI and
                    % success probability
                    av_aoi = av_aoi + (aoi + latency(l - 1) / 2) * latency(l - 1);
                    av_aoi = av_aoi + tau ^ 2 / 2 - latency(l - 1) ^ 2 / 2;
                    p_succ = p_succ + 1 / L;
                    aoi = tau;
                end
            end
            done = 0;
            running_latency = 0;
        end
        % New batch: update state and compute frame completion events
        state(b) = N(b);
        X = sum(state);
        delays = cumsum(exprnd(1 / mu, 1, X));
        completed(b, l) = length(find(delays < nu(b)));
        % Assign completion events to batches
        for i = 1 : completed(b, l)
            probs = cumsum(state / sum(state));
            served = find(probs >= rand, 1);
            % Check if the selected batch has been served
            if (served == 1 && done == 0)
                % Check if the selected client is in the batch (first
                % client of the batch)
                if (randi(state(1)) == 1)
                    latency(l) = running_latency + delays(i);
                    done = 1;
                    if (latest > 0)
                        paoi(l) = aoi + latency(l);
                    end
                    latest = l;
                end
            end
            state(served) = state(served) - 1;
        end
        % If client has not been served, update running latency
        if (done == 0)
            running_latency = running_latency + nu(b);
        end
    end
end

% Remove unsuccessful frames and compute success probability
av_aoi  = av_aoi / (tau * L);
latency = latency(latency > 0);
paoi = paoi(paoi > 0);

end