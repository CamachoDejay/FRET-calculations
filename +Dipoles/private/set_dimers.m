function [ pos ] = set_dimers( pos, dimer_prob, dimer_distance )
%SET_DIMERS set dimers into monomer lattice. Dimer sites are chosen
%randomly. 
%   Detailed explanation goes here

% total number of monomer sites
n_dipoles = size(pos,2);
% total number of monomer sites that must be exchanged by dimers
dimer_n   = round(n_dipoles*dimer_prob);
% find random indices to replace
dimer_ind = randperm(n_dipoles);
dimer_ind = dimer_ind(1:dimer_n);

% to generate a dimer I have to change 1 position in the cubic
% lattice by 2 positions. The tric is that the distance between
% those points is fixed to 'dimer_distance' nm, but the
% relative orientation of the points to its center of mass
% (lattice position) is random. Thus after finding the lattice
% positions to replace I must now find by how much in x y z

r_unit_vectors   = get_random_orient_vec( dimer_n );
displacement     = r_unit_vectors.*(dimer_distance/2);

d_pos_pos = pos(:,dimer_ind)+displacement;
d_pos_neg = pos(:,dimer_ind)-displacement;

pos(:,dimer_ind) = d_pos_pos;
pos(:,end+1:end+dimer_n) = d_pos_neg;

            

end

