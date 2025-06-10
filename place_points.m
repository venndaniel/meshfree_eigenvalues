function [x, y] = place_points(lsf, N, xlim, ylim, per, xb, yb, ew)
% creates a point cloud with N points inside a domain defined by lsf < 1
% - xlim, ylim define the edges of a box containing the level set
% - eps is the tolerance for the root-finding procedure to move the points
% to the level set
% - per is the number of test points per point added to point cloud
% a higher value of per will result in a more evenly-spaced point cloud
% - ew is the edge weighting: higher value results in more points near the
% edge of the domain. A value of 0 gives an even distribution
%
% - bisection is used to find a point near enough to the level set
    x = zeros(N, 1); % to store point clou
    y = zeros(N, 1);
    found = 0;
    % size of box
    ellx = xlim(2) - xlim(1);
    elly = ylim(2) - ylim(1);
    while found < N
        foundt = 0;
        xtt = zeros(per, 1); % to store test points
        ytt = zeros(per, 1);
        while foundt < per
            xt = rand()*ellx + xlim(1); % create test point
            yt = rand()*elly + ylim(1);
            % only add to point cloud if the test point is inside level set
            if lsf(xt, yt) < 1 
                foundt = foundt + 1;
                xtt(foundt) = xt;
                ytt(foundt) = yt;
            end
        end
        
        % compute penalty value
        tests = (ew*lsf(xtt', ytt')+1).*min([(x(1:found) - xtt').^2 + ...
            (y(1:found) - ytt').^2; ...
            ((xb - xtt').^2 + (yb - ytt').^2)], [], 1);
        [~, ind] = max(tests); % find point with largest penalty
        found = found + 1;
        x(found) = xtt(ind); % add point to point cloud
        y(found) = ytt(ind);

    end
end
