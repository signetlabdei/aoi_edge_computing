function [delta] = expected_aoi_fifo(N, mu, nu, b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       function: expected_aoi_fifo                       %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the expected AoI in the FIFO system                            %
%                                                                         %
% Inputs:                                                                 %
% -N:   the number of clients in each batch [1*B]                         %
% -mu:  the service rate [scalar]                                         %
% -nu:  the period between one batch and the next [1*B]                   %
% -b:   the current batch index [scalar]                                  %
%                                                                         %
% Outputs:                                                                %
% -delta:   the expected AoI [scalar]                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the AoI using the geometric method (E[TY]/E[Y]+E[Y^2]/2E[Y])
EY = real(expected_interpacket_fifo(N, mu, nu, b));
delta = real(expected_TY_fifo(N, mu, nu, b)) / EY;
delta = delta + real(squared_interpacket_fifo(N, mu, nu, b)) / 2 / EY;


end