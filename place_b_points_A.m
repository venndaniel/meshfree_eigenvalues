function [xf, yf] = place_b_points_A(lsf, N, xlim, ylim, eps, Ntest)
% creates a point cloud in 2D with N points on a level set (lsf = 1)
% - xlim, ylim define the edges of a box containing the level set
% - eps is the tolerance for the root-finding procedure to move the points
% to the level set
% - Ntest is the total number of test points (must be >= N)
% a higher value of per will result in a more evenly-spaced point cloud
%
% - bisection is used to find a point near enough to the level set
    x = zeros(Ntest, 1); % store Ntest points
    y = zeros(Ntest, 1);
    found = 0;
    % size of box
    ellx = xlim(2) - xlim(1); 
    elly = ylim(2) - ylim(1);
    while found < Ntest % add points until Ntest points are found
        xtp = rand()*ellx + xlim(1);
        ytp = rand()*elly + ylim(1);
        xtn = rand()*ellx + xlim(1);
        ytn = rand()*elly + ylim(1);
        while lsf(xtn, ytn) > 1 % ensure this point is inside level set
            xtn = rand()*ellx + xlim(1);
            ytn = rand()*elly + ylim(1);
        end
        while lsf(xtp, ytp) < 1 % ensure this point is outside
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
        found = found + 1;
        x(found) = (xtp + xtn)/2; % add the midpoint as a test point
        y(found) = (ytp + ytn)/2;
    end
    found = 1;
    xf = zeros(N, 1);
    yf = zeros(N, 1);
    xf(1) = x(1); % add first point to point cloud
    yf(1) = y(1);
    while found < N
        tests = min((xf(1:found) - x').^2 + (yf(1:found) - y').^2, [], 1);
        [~, ind] = max(tests);
        found = found + 1;
        % add the test point farthest from existing points until N points
        % are found
        xf(found) = x(ind); 
        yf(found) = y(ind);
    end
end