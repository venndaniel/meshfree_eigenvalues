% Code for Subsubsection 4.5.1 (Error Estimation)

%% Problem Setup

% set plot text to LaTeX
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

Np = 60; % number of boundary points
Nt = round(Np^2/4); % number of interior points

wmax = 75;
Nb = (2*wmax + 1)^2; % total number of Fourier modes
% parameters
q = 4;
T = 1;
ell = 4; % length of box containing domain
mult = 1;

rng('default')
p = 2;

% Fourier frequencies (\omega in paper)
w = 2*pi/ell*(-wmax:wmax)';
[wx, wy] = meshgrid(w);
wx = reshape(wx, [Nb, 1]);
wy = reshape(wy, [Nb, 1]);
w = sqrt(wx.^2 + wy.^2); % norm of \omega

di = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(abs(w)))); % d^{1/2}

% level set function and derivatives
lsf = @(x,y) abs(x).^p + abs(y).^p;

% create point cloud
[xb, yb] = place_b_points_A(lsf, Np, [-2, 2], [-2, 2], 1e-15, 1500);
%[xb, yb] = place_b_points_2(lsf, Np, [-2, 2], [-2, 2], 1e-15, 40);
[x, y] = place_points(lsf, Nt-Np, [-2, 2], [-2, 2], 100, xb, yb, 4);
x = [x; xb];
y = [y; yb];

r = sqrt(x.^2 + y.^2); % distance from origin

% normal vector on boundary
nx = xb.^(p-1); 
ny = yb.^(p-1);
nn = vecnorm([nx, ny], 2, 2);
nx = nx ./ nn;
ny = ny ./ nn;


V = di'.*exp(1i*(wx'.*x + wy'.*y)); % function values
Vlap = -w.^2'.*V; % Laplacian
Vb = di'.*exp(1i*(wx'.*xb + wy'.*yb)); % boundary values
Vbn = 1i*(wx'.*nx + wy'.*ny).*Vb; % normal derivative on boundary

% point(s) to set function = 1 on boundary (a_j in paper)
inds = 1;
xq = xb(inds');
yq = yb(inds');

Vq = di'.*exp(1i*(wx'.*xq + wy'.*yq)); % function value at a_j

% Q = (r/2 + 1/5*cos(5*r))./(2*r.^3 + 1); % potential
Q = 0;
Vs1 = [Vlap+Q.*V; Vbn; Vq];
Vs2 = [Vlap*0; -Vb; 0*Vq];
"computing matrices"
% faster to store pre-computed blocks of matrices for larger eigenvalue
% tests
F1 = Vs1*Vs1'; % Phi_0
F2 = Vs1*Vs2'; % Phi_1 / 2
F3 = Vs2*Vs2'; % Phi_2
'done matrices'

%% Error Estimation
f = [zeros(Nt, 1); zeros(Np, 1); 1];
P0 = F1; % Phi_0
P1 = F2 + F2'; % Phi_1
P2 = F3; % Phi_2

lams = []; % list of eigenvalues
exval = 10;
starts = 10; % starting point for Newton's method

cw = 1; % change weight
for s = starts
    change = 1;
    lam = s
    iter = 0;
    while abs(change) > 1e-8 && iter < 20
   
        F = lam^2*P2 + lam*P1 + P0; % Phi for current lambda

        % code for first and second derivative
        beta = cond_solve(F, f, 2);
        f2 = (2*lam*P2 + P1)*beta;
        beta2 = cond_solve(F, f2, 2);
        dN = real(-beta'*f2);
        ddN = real(-2*beta'*P2*beta + 2*beta2'*f2);

        if ddN < 0
            iter = 500000; % break look if concave down
        end

        change = -cw*dN/(ddN);
        lam = lam + change; % update lambda
        iter = iter + 1;
    end

    if (ddN > 0 && (change < 1e-5)) && (min([1000,abs(lam - lams)]) > 1e-4)
        lams = [lams lam]; % add to eigenvalue list
    end
end
lams = sort(lams)

test = lams(1); % minimum (assuming single eigenvalue found)

% list of lambda values near computed minimum, distanced logarithmically
pows = 1:0.1:7;
lams = test + [-10.^(-pows) 10.^(-pows)];

norms = zeros(length(lams), 1);

% indices for smaller Phi to compute norm ratio
sele = round(length(P1)*0.9);
selind = [randperm(length(P1)-1, sele) length(P1)];

% check norm ratio near minimum
for j = 1:length(lams)
    lam = lams(j);
    F = lam^2*P2 + lam*P1 + P0;
    Fs = F(selind, selind); % smaller Phi matrix

    beta = cond_solve(F, f, 2);
    fs = f(selind);
    betas = cond_solve(Fs, fs, 2);
    norms(j) = sqrt(abs(beta'*f))/sqrt(abs(betas'*fs));
    "Ratios: " + j/length(lams)*100 + "% Complete"
end


% loglog(abs(lams - test), norms, 'o');
[mn, ind] = min(norms); % minimum ratio
fails = (norms/mn)>=mn; % indices where ratio > minimum ratio squared
err = min(abs(exval - test)) % actual error
min(abs(lams(fails)-lams(ind))) % error estimate
% loglog(abs(lams-test), norms/mn - 1, 'o');
% grid on;
% hold on;


