function [ ori ] = get_orient_vec( n_dipoles, sigma )
%GET_RANDOM_ORIENTATIONS This function returns random oriented unitary
%vectors in the whole sphere. It is based on the random piking of a point
%on the surface of a unit spherethe. The direction vector of dipole i is
%given by: vi = v(:,i); in x y z
%   Detailed explanation goes here

ori = NaN(3,n_dipoles);

% building the data vector
res = 1001;
dataVec = linspace(0,1,res);

% building the weights.
mu = 0;
w1 = (1/(sigma*sqrt(2*pi)))*exp((-1.*(dataVec-mu).^2)./(2*sigma^2));
Sw1 = sum(w1);

mu = 1;
w2 = (1/(sigma*sqrt(2*pi)))*exp((-1.*(dataVec-mu).^2)./(2*sigma^2));
Sw2 = sum(w2);

assert(abs(Sw1/Sw2 - 1)/ Sw1 < 0.02, 'unexpected problems with uneven dist')

w2 = w1 + w2;



weights = w2;

weights = weights./sum(weights);

% sample the data with a weight given by the normal dist
N = n_dipoles;
y = datasample(dataVec,N,'Weights',weights);


%  1) spherical coordinates
phi      = rand(1,n_dipoles)*2*pi;         % azimuthal angle
% phi      = y*2*pi;         % azimuthal angle
% theta    = acos(2*rand(1,n_dipoles)-1);    % polar angle
theta    = acos(2*y-1);    % polar angle

% 2) cartesian
ori(1,:) = cos(phi).*sin(theta);  %x
ori(2,:) = sin(phi).*sin(theta);  %y
ori(3,:) = cos(theta);            %z

ang = pi/2;
Rx = [1 0 0; 0 cos(ang) -sin(ang); 0 sin(ang) cos(ang)];
Ry = [cos(ang) 0 sin(ang); 0 1 0; -sin(ang) 0 cos(ang)];
ori = Rx * ori;

end

