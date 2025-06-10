function [xs, ys, zs] = sphere_cloud(Nt, check)
% creates Nt points on the unit sphere with better-than-random spacing
% - check is the number of test points per point in the points cloud; a
% higher value will result in a more evenly spaced point cloud
    xs = zeros(Nt, 1); % to store point cloud
    ys = xs; zs = ys;
    
    fou = 1; % counter for number of points found
    
    r = rand(1, 3)*8 - 4; % random point
    xs(1) = r(1); ys(1) = r(2); zs(1) = r(3);
    nm = sqrt(xs(1)^2 + ys(1)^2 + zs(1)^2);
    xs(1) = xs(1)/nm; % project back to sphere
    ys(1) = ys(1)/nm; zs(1) = zs(1)/nm;
    
    while fou < Nt
        % create first random point to check against existing point cloud
        r = rand(1, 3)*8 - 4;
        x = r(1); y = r(2); z = r(3);
        nm = sqrt(x^2 + y^2 + z^2);
        x = x/nm;
        y = y/nm;
        z = z/nm;
    
        % x, y, z stores potential point to add
    
        ch = 1;
        dist = min(vecnorm([xs-x, ys-y, zs-z], 2, 2));
        % make more random points to check
        while ch < check
            r = rand(1, 3)*8 - 4;
            xt = r(1);
            yt = r(2);
            zt = r(3);
            nm = sqrt(xt^2 + yt^2 + zt^2);
            xt= xt/nm;
            yt = yt/nm;
            zt = zt/nm;
            distt = min(vecnorm([xs-xt, ys-yt, zs-zt], 2, 2));
    
            % update point to add if new point is farther from existing
            if distt > dist 
                x = xt; y = yt; z = zt;
                dist = distt;
            end
            ch = ch + 1;
        end
        % add a point to the final point cloud
        fou = fou + 1;
        xs(fou) = x;
        ys(fou) = y;
        zs(fou) = z;
    end
end