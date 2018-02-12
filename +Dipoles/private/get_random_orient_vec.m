function [ ori ] = get_random_orient_vec( n_dipoles )
%GET_RANDOM_ORIENTATIONS This function returns random oriented unitary
%vectors in the whole sphere. It is based on the random piking of a point
%on the surface of a unit spherethe. The direction vector of dipole i is
%given by: vi = v(:,i); in x y z
%   Detailed explanation goes here

ori = NaN(3,n_dipoles);

%  1) spherical coordinates
phi      = rand(1,n_dipoles)*2*pi;         % azimuthal angle
theta    = acos(2*rand(1,n_dipoles)-1);    % polar angle

% 2) cartesian
ori(1,:) = cos(phi).*sin(theta);  %x
ori(2,:) = sin(phi).*sin(theta);  %y
ori(3,:) = cos(theta);            %z


            

end

