function [x] = state2vec_fifo(state, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        function: state2vec_fifo                         %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Gets the state vector, including the total number of active clients     %
% (element 1) and the current batch index (el. 2), from the compressed    %
% state, used for indexing and vectorial operations, for the FIFO system  %
%                                                                         %
% Inputs:                                                                 %
% -state:   the compressed state [scalar]                                 %
% -N:       the number of clients in each batch [1*B]                     %
%                                                                         %
% Outputs:                                                                %
% -x:       the state vector [1*2]                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Utility variables
B = length(N);
x = zeros(1, 2);
state = state - 1;

% Extract state components
x(2) = mod(state, B) + 1;
x(1) = floor(state / B);

end