function [ p_t ] = get_transfer_prob_matrix( k_t, t_d)
%GET_TRANSFER_PROB_MATRIX translating rates into probabilities
%   Detailed explanation goes here

n_dipoles = size(k_t,2);
emission_rate = (1/t_d);
for i = 1:n_dipoles
	k_t(i,i) = emission_rate;
end
tot = sum(k_t,2);
p_t = k_t./repmat(tot,1,n_dipoles);
            
end

