function [ k2 ] = get_ori_factor( pos, ori )
%GET_ORIENTATION_FACTOR vectorized verion of k2 calculation

% init orientation factors

    n_dipoles = size(pos,2);
    k2   = single(NaN(n_dipoles,n_dipoles));

    for i = 1:n_dipoles
        di     = pos(:,i);
        oi     = ori(:,i);

        Ri2   = pos-repmat(di,1,n_dipoles);
        Rnorm = sum(Ri2.^2,1).^.5;
        Ri2   = Ri2./repmat(Rnorm,3,1);

        oi_rep = repmat(oi,1,n_dipoles);
        c1 = dot(oi_rep,ori);
        c2 = dot(Ri2,oi_rep);
        c3 = dot(Ri2,ori);

        ki = c1 - 3 * c2 .* c3;
        ki = ki.^2;
        k2(i,:) = single(ki);
    end


end

