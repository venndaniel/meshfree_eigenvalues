function [x, y] = place_b_points_B(lsf, N, xlim, ylim, eps, per)
% creates a point cloud in 2D with N points on a level set (lsf = 1)
% - xlim, ylim define the edges of a box containing the level set
% - eps is the tolerance for the root-finding procedure to move the points
% to the level set
% - per is the number of test points per point added to point cloud
% a higher value of per will result in a more evenly-spaced point cloud
%
% - bisection is used to find a point near enough to the level set
    x = zeros(N, 1); % to store point cloud
    y = zeros(N, 1);
    found = 0;
    % size of box
    ellx = xlim(2) - xlim(1);
    elly = ylim(2) - ylim(1);
    while found < N % continue until N points found
        xt = zeros(per, 1); % store test points
        yt = zeros(per, 1);
        foundt = 0;
        while foundt < per % find per test point on level set
            xtp = rand()*ellx + xlim(1);
            ytp = rand()*elly + ylim(1);
            xtn = rand()*ellx + xlim(1);
            ytn = rand()*elly + ylim(1);
            while lsf(xtn, ytn) > 1 % ensure this point is inside level set
                xtn = rand()*ellx + xlim(1);
                ytn = rand()*elly + ylim(1);
            end
            while lsf(xtp, ytp) < 1 % ensure this one is outside
                xtp = rand()*ellx + xlim(1);
                ytp = rand()*elly + ylim(1);
            end
            % bisection until points are within tolerance
            while (xtp - xtn)^2 + (ytp - ytn)^2 > eps^2 
                xtt = (xtp + xtn)/2;
                ytt = (ytp + ytn)/2;
                % keep one point inside, one point outside
                if lsf(xtt, ytt) > 1
                    xtp = xtt;
                    ytp = ytt;
                else
                    xtn = xtt;
                    ytn = ytt;
                end
            end
            foundt = foundt + 1;
            xt(foundt) = (xtp + xtn)/2; % add midpoint as a test point
            yt(foundt) = (ytp + ytn)/2;
        end
        if found == 0
            found = found + 1;
            % if no points in point cloud, add first test point
            x(found) = xt(1); 
            y(found) = yt(1);
        else
            tests = min((x(1:found) - xt').^2 + (y(1:found) - yt').^2, [], 1);
            [~, ind] = max(tests);
            found = found + 1;
            % otherwise, add the test point farthest from existing points
            x(found) = xt(ind);
            y(found) = yt(ind);
        end
    end
end