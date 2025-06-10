% Code for Subsection 4.4 (Laplacian on an Irregular 2D Domain)
% Can run each code section separately as long as "Problem Setup" is run
% first
% intially run in MATLAB R2024b

%% Problem Setup
% set text to LaTeX for plots
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

Np = 120; % number of boundary points
Nt = round(Np^2/30); % number of interior points

wmax = 75;
Nb = (2*wmax + 1)^2; % total number of Fourier modes
% parameters
q = 2.5;
T = 0.5;
ell = 3; % length of box to contain domain

rng('default')

% create Fourier frequencies (\omega)
w = 2*pi/ell*(-wmax:wmax)'; 
[wx, wy] = meshgrid(w);
wx = reshape(wx, [Nb, 1]);
wy = reshape(wy, [Nb, 1]);
w = sqrt(wx.^2 + wy.^2); % norm of \omega

di = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(abs(w)))); % d^{-1/2}

% create level set function and compute derivatives
syms x y
phi = x.^2 + y.^2 - 0.2*sin(5*x)-0.2*sin(6*y) + 0.1./(x.^2 + y.^2);
phixs = diff(phi, x); phiys = diff(phi, y);
phix = matlabFunction(phixs, "Vars",[x y]);
phiy = matlabFunction(phiys, "Vars",[x y]);
lsf = matlabFunction(phi); % level set function

% create point clouds
[xb, yb] = place_b_points_B(lsf, Np, [-2, 2], [-2, 2], 1e-15, 2000);
[x, y] = place_points(lsf, Nt, [-2, 2], [-2, 2], 200, [20], [20], 0);
nx = phix(xb,yb);
ny = phiy(xb,yb);

% normal vectors to boundary (not used in this test)
nn = vecnorm([nx, ny], 2, 2);
nx = nx ./ nn;
ny = ny ./ nn;

% function values, Laplacian, boundary values, boundary normal derivative
V = di'.*exp(1i*(wx'.*x + wy'.*y));
Vlap = -w.^2'.*V;
Vb = di'.*exp(1i*(wx'.*xb + wy'.*yb));
Vbn = 1i*(wx'.*nx + wy'.*ny).*Vb;

% points to set function value = 1 (a_j in paper)
% use a boundary point for Steklov problem, interior otherwise
inds = 1;
% xq = xb(inds');
% yq = yb(inds');
xq = x(inds'); 
yq = y(inds');

Vq = di'.*exp(1i*(wx'.*xq + wy'.*yq)); % to evaluate function at a_j

%Q = (r/2 + 1/5*cos(5*r))./(2*r.^3 + 1);
Q = 0;

% Steklov matrices
% Vs1 = [Vlap+Q.*V; Vbn; Vq];
% Vs2 = [Vlap*0; -Vb; 0*Vq];

% Dirichlet eigenvalues matrices
Vs1 = [Vlap+Q.*V; Vb; Vq];
Vs2 = [V; -Vb*0; 0*Vq];

%% Create Phi_j matrices
f = [zeros(Nt, 1); zeros(Np, 1); 1];
F1 = Vs1*Vs1'; 'one' % Phi_0 in text
F2 = Vs1*Vs2'; 'two' % Phi_1 / 2 in text
F3 = Vs2*Vs2'; % Phi_2 in text
'done matrices'

%% Find an eigenvalue via Newton's method
P0 = F1; % Phi_0
P1 = F2 + F2'; % Phi_1
P2 = F3; % Phi_2

lams = []; % list of eigenvalues

starts = 28.5; % iteration starting point

cw = 1; % change weight0
tic;
lam = starts;
% precondition with Cholesky decomp of regularized Phi at starting lambda
% experimentally stabilizes result
% Rc = chol(lam^2*P2 + lam*P1 + P0 + 1e-15*eye(length(P0)));
for s = starts
    change = 1;
    lam = s
    iter = 0;
    while abs(change) > 1e-5 && iter < 20
   
        % compute first and second derivative
        F = lam^2*P2 + lam*P1 + P0;
        % Fp = (Rc'\(Rc'\F)')';
        % Fp = 1/2*real(Fp + Fp'); 
        beta = cond_solve(F, f, 2);
        %beta = Rc\(Fp\(Rc'\f));
        f2 = (2*lam*P2 + P1)*beta;
        %beta2 = Rc\(Fp\(Rc'\f2));
        beta2 = cond_solve(F, f2, 2);
        dN = real(-beta'*f2); % first derivative
        ddN = real(-2*beta'*P2*beta + 2*beta2'*f2); % second

        if ddN < 0
            iter = 500000; % break loop if concave down
        end

        change = -cw*dN/(ddN);
        lam = lam + change % update lambda
        iter = iter + 1;
    end

    if (ddN > 0 && (change < 1e-5)) && (min([1000,abs(lam - lams)]) > 1e-4)
        lams = [lams lam]; % add an eigenvalue
    end
end
toc
lams = sort(lams)

%% plot an eigenfunction
lam = 28.484092068881939;
f = [zeros(Nt, 1); zeros(Np, 1); 1];
a = lsqminnorm(Vs1 + lam*Vs2, f); % solves via complete orthogonal decomp

% plot scaled eigenfunction (max value = 1)
scatter(x, y, 10, real(V*a)./max(abs(V*a)), 'filled');
axis('equal')
hold on;
scatter(xb, yb, 20, 'k', 'filled')
hold off;
xlabel('$x$'); ylabel('$y$')
colorbar;