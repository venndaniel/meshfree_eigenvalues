function [xf, yf, zf] = place_b_points(lsf, N, xlim, ylim, zlim, eps, per)
% creates a point cloud in 3D with N points on a level set (lsf = 0)
% xlim, ylim, zlim define the edges of a box containing the level set
% eps is the tolerance for the root-finding procedure to move the points
% to the level set
% per is the number of test points per point added to point cloud
% a higher value of per will result in a more evenly-spaced point cloud
%
% bisection is used to find a point near enough to the level set

    Ntest = N*per; % total number of points to test
    x = zeros(Ntest, 1);
    y = zeros(Ntest, 1);
    z = zeros(Ntest, 1);
    % size of box
    ellx = xlim(2) - xlim(1);
    elly = ylim(2) - ylim(1);
    ellz = zlim(2) - zlim(1);
    % find Ntest total points
    for found = 1:Ntest
        % find a point inside and outside the level set boundary
        xtp = rand()*ellx + xlim(1);
        ytp = rand()*elly + ylim(1);
        ztp = rand()*ellz + zlim(1);
        xtn = rand()*ellx + xlim(1);
        ytn = rand()*elly + ylim(1);
        ztn = rand()*ellz + zlim(1);
        while lsf(xtn, ytn, ztn) > 0 % ensure this point is inside
            xtn = rand()*ellx + xlim(1);
            ytn = rand()*elly + ylim(1);
            ztn = rand()*ellz + zlim(1);
        end
        while lsf(xtp, ytp, ztp) < 0 % ensure this point is outside
            xtp = rand()*ellx + xlim(1);
            ytp = rand()*elly + ylim(1);
            ztp = rand()*ellz + zlim(1);
        end
        % bisect until points are close
        while vecnorm([xtp - xtn, ytp - ytn, ztp - ztn], 2, 2) > eps
            xtt = (xtp + xtn)/2;
            ytt = (ytp + ytn)/2;
            ztt = (ztp + ztn)/2;
            % if test point is outside, make it the new outside point
            if lsf(xtt, ytt, ztt) > 0 
                xtp = xtt; ytp = ytt; ztp = ztt;
            else % otherwise, make it the inside point
                xtn = xtt; ytn = ytt;  ztn = ztt;               
            end
        end
        x(found) = (xtp + xtn)/2; % take average as final point
        y(found) = (ytp + ytn)/2;
        z(found) = (ztp + ztn)/2;
    end
    found = 1;
    xf = zeros(N, 1);
    yf = zeros(N, 1);
    zf = zeros(N, 1);
    xf(1) = x(1); % add first point to point cloud
    yf(1) = y(1);
    zf(1) = z(1);
    while found < N
        ran = (found*per+1):(found*per + per); % batch of points to test
        % closest point in point cloud to each point in batch
        tests = min((xf(1:found) - x(ran)').^2 +...
            (yf(1:found) - y(ran)').^2 +...
            (zf(1:found) - z(ran)').^2, [], 1); 
        % find the point farthest from existing points and add to point
        % cloud
        [~, ind] = max(tests); 
        found = found + 1;
        xf(found) = x((found-1)*per+ind);
        yf(found) = y((found-1)*per+ind);
        zf(found) = z((found-1)*per+ind);
    end
end

