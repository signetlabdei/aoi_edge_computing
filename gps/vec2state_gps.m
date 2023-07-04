function [state] = vec2state_gps(x, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         function: vec2state_gps                         %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Gets the compressed state, used for indexing and vectorial operations,  %
% from the state vector, which includes the number of active clients for  %
% each of the B batches (elements 1 to B) and the current batch index,    %
% for the GPS system                                                      %
%                                                                         %
% Inputs:                                                                 %
% -x:       the state vector [1*(B+1)]                                    %
% -N:       the number of clients in each batch [1*B]                     %
%                                                                         %
% Outputs:                                                                %
% -state:   the compressed state [scalar]                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize utility variables
state = 0;
B = length(N);
mult = 1;

% Add the batch index (last in the vector)
state = state + mult * (x(B + 1) - 1);

mult = B;

% Multiply by a factor N(i) to compress the state vector in a scalar
for i = 1 : B
    state = state + x(i) * mult;
    mult = mult * (N(i) + 1);
end

state = state + 1;

end