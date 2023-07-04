function [x] = state2vec_gps(state, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         function: state2vec_gps                         %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Gets the state vector, which includes the number of active clients for  %
% each of the B batches (elements 1 to B) and the current batch index,    %
% from the compressed state, used for indexing and vectorial operations,  %
% for the GPS system                                                      %
%                                                                         %
% Inputs:                                                                 %
% -state:   the compressed state [scalar]                                 %
% -N:       the number of clients in each batch [1*B]                     %
%                                                                         %
% Outputs:                                                                %
% -x:       the state vector [1*(B+1)]                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the utility variables
B = length(N);
x = zeros(1, B + 1);
state = state - 1;

% Get batch index
x(B + 1) = mod(state, B) + 1;
state = floor(state / B);

% Iterate over batches to decompress the state
for i = 1 : B
    x(i) = mod(state, N(i) + 1);
    state = floor(state / (N(i) + 1));
end


end