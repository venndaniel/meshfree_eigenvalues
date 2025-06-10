% Code for Subsection 4.8 (Surface Steklov)
% runs the Poisson problem test on a catenoid with a wavy edge
% intially run in MATLAB R2023a

%% Problem Setup

% set plot text to LaTex
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

% N_s in paper, half of N_\partial (ex. Nu=90/2 means N_\partial=90)
Nu = 90/2; 
Nv = Nu/3;
N = Nu*Nv; % total number of points
rng('default')
mult = 1;

% parameters
T = 5;
q = 4;
conds = 1; % number of conditioning/rescaling iterations
omax = 15; % test uses (2*omax + 1)^3) Fourier modes
ell = 5; % edge length of box (centred at origin)
wave = 0.1;

Nb = (omax*2 + 1)^3; % number of Fourier modes

% create point cloud
% u, v parameters
uo = linspace(0, 2*pi, Nu+1)';
uo = uo(1:Nu);
% v = linspace(-1, 1, Nv)';
v = sin(pi/2*linspace(-1, 1, Nv))';

% u and v for parametrization of catenoid
[u, v] = meshgrid(uo, v);
u = reshape(u, [N, 1]);
v = reshape(v, [N, 1]);
v = v + wave*sin(3*u); % make "wavy" catenoid

[x, y, z] = get_pos(u, v); % point cloud
[xb, yb, zb] = get_pos(uo, ones(Nu, 1) + wave*sin(3*uo)); % boundary points

% boundary curve tangent vectors
[tbx, tby, tbz] = get_su(uo, ones(Nu, 1) + wave*sin(3*uo));
[tbxa, tbya, tbza] = get_sv(uo, ones(Nu, 1) + wave*sin(3*uo));
tbx = tbx + tbxa.*3*wave.*cos(3*uo);
tby = tby + tbya.*3*wave.*cos(3*uo);
tbz = tbz + tbza.*3*wave.*cos(3*uo);
% normal to surface on boundary (\hat{n})
[snbx, snby, snbz] = get_n(uo, ones(Nu, 1) + wave*sin(3*uo));

% normal to boundary (nu \propto n \times t)
Nn = cross([tbx, tby, tbz], [snbx, snby, snbz], 2);
nn = vecnorm(Nn, 2, 2);
Nn = Nn ./ nn;
nbx = [-Nn(:, 1)];
nby = [-Nn(:, 2)];
nbz = [-Nn(:, 3)];

% repeat the process on the other boundary
[xb1, yb1, zb1] = get_pos(uo, -ones(Nu, 1) + wave*sin(3*uo));

[tbx, tby, tbz] = get_su(uo, -ones(Nu, 1) + wave*sin(3*uo));
[tbxa, tbya, tbza] = get_sv(uo, -ones(Nu, 1) + wave*sin(3*uo));
tbx = tbx + tbxa.*3*wave.*cos(3*uo);
tby = tby + tbya.*3*wave.*cos(3*uo);
tbz = tbz + tbza.*3*wave.*cos(3*uo);
[snbx, snby, snbz] = get_n(uo, -ones(Nu, 1) + wave*sin(3*uo));

Nn = cross([tbx, tby, tbz], [snbx, snby, snbz], 2);
nn = vecnorm(Nn, 2, 2);
Nn = Nn ./ nn;
nbx = [nbx; Nn(:, 1)];
nby = [nby; Nn(:, 2)];
nbz = [nbz; Nn(:, 3)];

[nx, ny, nz] = get_n(u, v); % normal vectors on surface

% boundary points
xb = [xb; xb1];
yb = [yb; yb1];
zb = [zb; zb1];

% scatter3(x, y, z);
% hold on;
% quiver3(xb, yb, zb, nbx, nby, nbz)
% axis('equal')

freq = (-omax:omax)*2*pi/ell;

% omega arrays
[ox, oy, oz] = meshgrid(freq);
ox = reshape(ox, [Nb, 1]);
oy = reshape(oy, [Nb, 1]);
oz = reshape(oz, [Nb, 1]);
oa = vecnorm([ox, oy, oz], 2, 2);

di = 1 ./ (exp(q*sqrt(2*pi/T)) + exp(q*sqrt(oa))); % d^{-1/2}

% interpolation matrices
V = di'.*exp(1i*(x.*ox' + y.*oy' + z.*oz')); % function values
Vb = di'.*exp(1i*(xb.*ox' + yb.*oy' + zb.*oz')); % boundary
uq = [0];
[xq, yq, zq] = get_pos(uq, 1 + wave*sin(3*uq));
Vq = di'.*exp(1i*(xq.*ox' + yq.*oy' + zq.*oz')); % evaluate at a_j
%Vq = Vb(1, :);
% boundary normal derivative
Vbn = 1i*(nbx.*ox' + nby.*oy' + nbz.*oz').*Vb;

Vlap = -oa.^2'.*V; % Laplace-Beltrami
Vn = 1i*(nx.*ox' + ny.*oy' + nz.*oz').*V; % first normal derivative
Vnn = -(nx.^2.*ox.^2' + ny.^2.*oy.^2' + nz.^2.*oz.^2' + ...
    2*nx.*ny.*(ox.*oy)' + ...
    2*nx.*nz.*(ox.*oz)' + ...
    2*nz.*ny.*(oz.*oy)').*V; % second normal derivative
 
Vs = [Vlap - Vnn; Vn; Vbn; Vq]; % full matrix (V_{N_b} in paper)
Vs2 = [Vlap*0; Vn*0; -Vb; 0*Vq];



"computing matrices"
% faster to store pre-computed blocks of matrices for larger eigenvalue
% tests
F1 = Vs*Vs'; % Phi_0
F2 = Vs*Vs2'; % Phi_1 / 2
F3 = Vs2*Vs2'; % Phi_2
'done matrices'


%% Newton's Method
f = [zeros(2*N, 1); zeros(2*Nu,1); 1];
ind = length(f);

P0 = F1; % Phi_0
P1 = F2 + F2'; % Phi_1
P2 = F3; % Phi_2

lams = []; % to store eigenvalues
starts = 0.46; % starting point for Newton's method

for s = starts
    change = 1;
    lam = s
    iter = 0;
    while abs(change) > 1e-8 && iter < 15
   
        F = lam^2*P2 + lam*P1 + P0; % Phi
        
        % compute first and second derivatives
        beta = cond_solve(F, f, conds);
        f2 = (2*lam*P2 + P1)*beta;
        beta2 = cond_solve(F, f2, conds);
        dN = real(-beta'*f2);
        ddN = real(-2*beta'*P2*beta + 2*beta2'*f2);

        if ddN < 0
            iter = 500000; % break if concave down
        end
        
        change = -dN/ddN;
        lam = lam + change; % iterate
        iter = iter + 1;
        
    end
    if (ddN > 0 && (change < 1e-5)) && (min([1000,abs(lam - lams)]) > 1e-4)
        lams = [lams lam]; % add eigenvalue
    end
end
lams = sort(lams)


%% Plot an eigenfunction (run above section first at one starting point)
scatter3(x, y, z, 30, real(V*((Vs+lam*Vs2)'*beta)), 'filled')
axis('equal')
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
colormap('jet')
colorbar()


function [xg, yg, zg] = get_pos(u, v)
    % maps parametrization of surface (u, v) to (x, y, z) coordinates
    xg = cosh(v).*cos(u);
    yg = cosh(v).*sin(u);
    zg = v;
end
    
function [ux, uy, uz] = get_su(u, v)
    % u partial derivative of parametrization
    ux = -cosh(v).*sin(u);
    uy = cosh(v).*cos(u);
    uz = zeros(length(u), 1);
end

function [vx, vy, vz] = get_sv(u, v)
    % v partial derivative of parametrization
    vx = sinh(v).*cos(u);
    vy = sinh(v).*sin(u);
    vz = ones(length(u), 1);
end

function [nx, ny, nz] = get_n(u, v)
    % gets the normal vector to the surface at (u, v)
    [ux, uy, uz] = get_su(u, v);
    [vx, vy, vz] = get_sv(u, v);
    N = cross([ux, uy, uz], [vx, vy, vz], 2);
    nn = vecnorm(N, 2, 2);
    N = N ./ nn;
    nx = N(:, 1);
    ny = N(:, 2);
    nz = N(:, 3);
end
