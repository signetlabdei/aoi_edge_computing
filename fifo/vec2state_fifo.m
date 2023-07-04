function [state] = vec2state_fifo(x, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        function: vec2state_fifo                         %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Gets the compressed state, used for indexing and vectorial operations,  %
% from the state vector, including the total number of active clients     %
% (element 1) and the current batch index (element 2), in the FIFO system %
%                                                                         %
% Inputs:                                                                 %
% -x:       the state vector [1*(2)]                                      %
% -N:       the number of clients in each batch [1*B]                     %
%                                                                         %
% Outputs:                                                                %
% -state:   the compressed state [scalar]                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility variable
B = length(N);

% Add the batch index to the state
state = x(2) - 1;

% Add the number of clients to the state
state = state + B * x(1);

% Ensure that the state is always at least 1
state = state + 1;

end