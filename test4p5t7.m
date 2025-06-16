% Code for Subsections 4.5, 4.6, and 4.7, except 4.5.1 (Steklov on Disk)
% intially run in MATLAB R2023a, tested again in R2024b

set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

Np = 40; % number of boundary points
Nt = round(Np^2/4); % number of interior points

% 1 for Figure 5 tests, 2 for Table 7, 3 for Table 9, 4 for Table 10
test = 4; 

% parameters
wmax = 75; 
Nb = (2*wmax + 1)^2; % number of Fourier modes
q = 4;
T = 1;
ell = 4; % size of box containing domain
mult = 1;

rng('default');
p = 2;

% Fourier frequencies
w = 2*pi/ell*(-wmax:wmax)';
[wx, wy] = meshgrid(w); 
wx = reshape(wx, [Nb, 1]);
wy = reshape(wy, [Nb, 1]);
w = sqrt(wx.^2 + wy.^2);

di = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(abs(w)))); % d^{1/2}

lsf = @(x,y) abs(x).^p + abs(y).^p;

% create point cloud
[xb, yb] = place_b_points_A(lsf, Np, [-2, 2], [-2, 2], 1e-15, 1500);
[x, y] = place_points(lsf, Nt-Np, [-2, 2], [-2, 2], 100, [xb], [yb], 4);
x = [x; xb];
y = [y; yb];

r = sqrt(x.^2 + y.^2);

% boundary normal vectors
nx = xb.^(p-1);
ny = yb.^(p-1);
nn = vecnorm([nx, ny], 2, 2);
nx = nx ./ nn;
ny = ny ./ nn;

V = di'.*exp(1i*(wx'.*x + wy'.*y)); % function values
Vlap = -w.^2'.*V; % Laplacian
Vb = di'.*exp(1i*(wx'.*xb + wy'.*yb)); % boundary values
Vbn = 1i*(wx'.*nx + wy'.*ny).*Vb; % normal derivative on boundary

% points to set function value = 1 (a_j in paper)
%inds = randsample(Np, mult);
inds = 1;
xq = xb(inds');
yq = yb(inds');

Vq = di'.*exp(1i*(wx'.*xq + wy'.*yq)); % function values at a_j

% q from Schrodinger-Steklov
Q = 0;
if test == 3
    Q = (2.404825557695773^2);
elseif test == 4
     Q = -(r/2 + 1/5*cos(5*r))./(2*r.^3 + 1);
end

Vs1 = [Vlap+Q.*V; Vbn; Vq];
Vs2 = [Vlap*0; -Vb; 0*Vq];
"computing matrices"
% faster to store pre-computed blocks of matrices for larger eigenvalue
% tests
F1 = Vs1*Vs1'; % Phi_0
F2 = Vs1*Vs2'; % Phi_1 / 2
F3 = Vs2*Vs2'; % Phi_2
'done matrices'

f = [zeros(Nt, 1); zeros(Np, 1); 1];
P0 = F1; % Phi_0
P1 = F2 + F2'; % Phi_1
P2 = F3; % Phi_2

lams = []; % list of eigenvalues
exval = 10;
starts = 10; % starting point for Newton's method
if test == 2
    starts = 0:0.3/pi:19;
elseif test == 3
    starts = 0.9;
    exval = 0.891592981473392;
elseif test == 4
    starts = 10;
    exval = 10.00807486;
end

cw = 1; % change weight
for s = starts
    change = 1;
    lam = s
    iter = 0;
    while abs(change) > 1e-8 && iter < 20
   
        F = lam^2*P2 + lam*P1 + P0; % Phi
      
        % first and second derivatives
        beta = cond_solve(F, f, 2);
        f2 = (2*lam*P2 + P1)*beta;
        beta2 = cond_solve(F, f2, 2);
        dN = real(-beta'*f2);
        ddN = real(-2*beta'*P2*beta + 2*beta2'*f2);

        if ddN < 0
            iter = 500000; % break if concave down
        end

        change = -cw*dN/(ddN);
        lam = lam + change; % iterate
        iter = iter + 1;
    end

    if (ddN > 0 && (change < 1e-4)) && (min([1000,abs(lam - lams)]) > 1e-2)
        lams = [lams lam]; % add new eigenvalue
    end
end
lams = sort(lams)

min(abs(exval - lams))/exval % relative error for tests 1, 3, 4

prev = lams(1);

