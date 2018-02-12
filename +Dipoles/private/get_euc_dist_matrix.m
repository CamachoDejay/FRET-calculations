function [ r ] = get_euc_dist_matrix( pos )
%GET_EUC_DIST_MATRIX vectorized version for euc dist calculation between a
%list of positions
%   Detailed explanation goes here
n_dipoles = size(pos,2);
r = single(NaN (n_dipoles,n_dipoles));
for i = 1:n_dipoles
   r(i,:) = single((sum((repmat(pos(:,i),1,n_dipoles) - pos).^2)).^.5);
end

            
end

