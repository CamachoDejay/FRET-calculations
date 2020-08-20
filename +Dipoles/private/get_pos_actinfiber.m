function [ pos ] = get_pos_actinfiber( cube, inter_dist )
%GET_POS_ACTINFIBER get positions along an actin fiber: equidistant on the
%coat of a cylinder which has the origin on its coat and proceeding in both 
% directions from this specific position (marking the excited dipole)???

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

