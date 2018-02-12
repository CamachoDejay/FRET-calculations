function [ em_ori ] = get_em_pol( em_angles_r )
%GET_EM_POL Summary of this function goes here
%   Detailed explanation goes here

% Emssion is observed from the x axis, thus the emission
% polarizer is rotated along the z-y plane;
em_ori = zeros(size(em_angles_r,2),3);
%  1) spherical coordinates
phi      = ones(size(em_angles_r)).*pi/2;% azimuthal angle
theta    = em_angles_r;%90*pi/180;    % polar angle
% 2) cartesian
em_ori(:,1) = cos(phi).*sin(theta);  %x
em_ori(:,2) = sin(phi).*sin(theta);  %y
em_ori(:,3) = cos(theta);            %z
% to remove numerical noise
test = em_ori(:,1);
assert (sum(test>1e-16) == 0)
em_ori(:,1)=em_ori(:,1).*0;

            
end

