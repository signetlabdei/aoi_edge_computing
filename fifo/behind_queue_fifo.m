function [Q] = behind_queue_fifo(N, b, bn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       function: behind_queue_fifo                       %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the number of frames behind the ones in batch b when batch bn  %
% is reached in the FIFO system                                           %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients in each batch [1*B]                         %
% -b:   the considered batch index [scalar]                               %
% -bn:  the current batch index [scalar]                                  %
%                                                                         %
% Outputs:                                                                %
% -Q:   the number of packets in the queuebehind those from batch b       %
%       [scalar]                                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No packets are queued behind batch b if it has just arrived
Q = 0;
if (b == bn)
    return;
end

% Compute the number of packets in the queue
if (b < bn)
    Q = sum(N(b + 1 : bn));
else
    Q = sum(N(b + 1 : end)) + sum(N(1 : bn));
end

end