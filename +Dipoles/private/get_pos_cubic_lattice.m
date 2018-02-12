function [ pos ] = get_pos_cubic_lattice( cube, inter_dist )
%GET_POS_CUBIC_LATTICE get positions in a cubic lattice: equidistant and
%centered around 0

x_r   = 0:inter_dist:cube.x_radius;
x_l   = fliplr(0:-inter_dist:-cube.x_radius);
x_grid   = [x_l(1:end-1) x_r];

y_r   = 0:inter_dist:cube.y_radius;
y_l   = fliplr(0:-inter_dist:-cube.y_radius);
y_grid   = [y_l(1:end-1) y_r];

z_r   = 0:inter_dist:cube.z_radius;
z_l   = fliplr(0:-inter_dist:-cube.z_radius);
z_grid   = [z_l(1:end-1) z_r];

[X,Y,Z] = meshgrid(x_grid ,y_grid, z_grid);

n_dipoles = length(x_grid)*length(y_grid)*length(z_grid);
pos = zeros (3,n_dipoles);
pos(1,:) = X(:);
pos(2,:) = Y(:);
pos(3,:) = Z(:);

end

