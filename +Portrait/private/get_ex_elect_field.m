function [ ex_ori ] = get_ex_elect_field( ex_angles_r )
%GET_EX_ORIENTATION from the excitation angles in radians calculates
%electric field vector of the excitation light
%   Detailed explanation goes here

% Excitation light travels along the x axis and the electric
% field vector is rotated on the y-z plane.
ex_ori = zeros(size(ex_angles_r,2),3);
%  1) spherical coordinates
phi      = ones(size(ex_angles_r)).*pi/2;% azimuthal angle
theta    = ex_angles_r;%90*pi/180;    % polar angle
% 2) cartesian
ex_ori(:,1) = cos(phi).*sin(theta);  %x
ex_ori(:,2) = sin(phi).*sin(theta);  %y
ex_ori(:,3) = cos(theta);            %z
% to remove numerical noise
test = ex_ori(:,1);
assert (sum(test>1e-16) == 0)
ex_ori(:,1)=ex_ori(:,1).*0;      

end

