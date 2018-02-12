function [ k_t ] = get_transfer_rate_matrix( FRETprops, k2, r )
%GET_TRANSFER_RATE_MATRIX calculation using matrix operations
%   Detailed explanation goes here

%%%%%%% Physical constants %%%%%%%
%avogadros number                %
N_a = 6.022e23; %[mol^-1]        %
%index of refraction             %
n   = FRETprops.refIndex;%  1.4 aqueous sol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contant for Forster integral calculation
f1 = 9000 * log(10) / (128 * pi^5 * N_a); % [mol cm^3]

%%%%%%% Chromophore properties %%%
%QY of the donor                 %
Q_d = FRETprops.quantumYieldD;   % 0.68 for GFP emerald
% lifetime of donor
t_d = FRETprops.lifetimeD;  % 2.7 [ns] for GFP
% Spectral overlap
J = FRETprops.J; % [mol^-1 cm^-1 nm^4]
f2 = 1e14; % [nm^2 * cm^-2] This factor takes into account the units of J.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% transfer rate
k_t = single(f1 * f2 * ( ((Q_d * J * k2) ./ (t_d * n^4) ) ./ (r.^6) ));% [ns^-1];
            

end

