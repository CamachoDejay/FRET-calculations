function [ emission, et_steps, res_energy] = get_em_after_ET( ET_matrix, ind_int )
%GET_EM_AFTER_ET Summary of this function goes here
%   Detailed explanation goes here

% now what is the effect of EET. It will redistribute the
% excitation energy according to the transfer matrix.

% number of dipoles in the system
n_dipoles = size(ET_matrix,1);

% initial absorption
%I only consider the center dipole for which I created the buffer.
ex_energy = zeros(1,n_dipoles);
ex_energy(ind_int) = 1;

% amount of light emitted changes in each loop. I will keep track of it in
% the variable emission
emission  = zeros(1,n_dipoles);

% now we do the proper transfer loop
res_energy = sum(ex_energy);
et_steps = 0;

% we will iterate while there is still energy in the dipoles
while res_energy > 1e-6
    
    % we consider what happens in each dipole separately, thus we create
    % a diagonal matrix out of the excitation energy
    abs_i = diag(ex_energy);
    % Then redistribute the energy according to the ET matrix    
    redist = abs_i * ET_matrix;    
    % find how much energy was not transferrred (thus emitted) by each dipole
    em_i = diag(redist)';
    % find where the energy of the system ended up considering now all
    % dipoles
    tr_i = sum(redist,1);
    % remove from the excited state energy that which was emitted
    ex_energy = tr_i - em_i;
    % keep track of where light is being emitted
    emission = emission + em_i;
    % check how much energy is left to redistribute
    res_energy = sum(ex_energy);
    % keep track of how many ET steps Im considering
    et_steps = et_steps + 1;
    
end


end

