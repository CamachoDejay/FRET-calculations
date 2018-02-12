function [ y_fit, epsilon, RMSD, eps ] = SFA_fit_lsqnonneg( model, y )
%POLIM_LSQNONNEG Summary of this function goes here
%   Detailed explanation goes here

assert(isa(model,'POLIM.twoDmodel'),...
          'problems with input model, it must be a POLIM.twoDmodel object')

[ Fnoet, Fet ] = model.get_noet_et;

% Solve nonnegative least-squares constraints problem
shapes = [Fnoet(:), Fet(:)];
eps    = lsqnonneg(shapes,y);
y_fit  = shapes*eps;
y_dif  = y - y_fit;
RMSD   = (sum((y_dif).^2)/size(y_fit,1))^0.5;

% now, carefull! eps still not contain the values of epsilon:
% eps(1) = (1-epsilon); eps(2)= epsilon
% this is clear if we add eps(1)+eps(2), because the result is not always 1
% This is because model and experimental intensities are scaled/normalized
% in different ways. Thus we have an extra scaling facotr:
% eps(1) = S*(1-epsilon); eps(2) = S*epsilon
% thus to recover epsilon we need to compute:
% S = eps(1)+eps(2); epsilon = eps(2)/S;
S = sum(eps);
epsilon = eps(2)/S;        
end

