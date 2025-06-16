% code for Subsection 4.3 (Laplace--Beltrami on a Genus 2 Surface)
% intially run in MATLAB R2023a, tested again in R2024b
% Can run each code section separately as long as "Problem Setup" is run
% first

%% Problem Setup
% run this section to create the necessary matrices to run the experiments
% in the other sections

Nt = 1200; % number of random points on the surface
mul = 1; % multiplicity to test for
rng('default')

% whether to use the curvature to compute \Delta_S, as in Eq. (46)
% false is Eq. (43), true is (46)
use_curv = false;

% set plot text to LaTex
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

% level set function and derivatives
syms x y z
phi = @(x,y,z) 0.25./((x-1).^2 + y.^2) + 0.25./((x+1).^2 + y.^2)+ z.^2-...
    1 + 0.1*x.^2 + 0.25*y.^2;
phix = @(x,y,z)  -0.25*2*(x-1)./((x-1).^2 + y.^2).^2 -...
    0.25*2*(x+1)./((x+1).^2 + y.^2).^2 + 0.1*2*x;
phiy = @(x,y,z)  -0.25*2*y./((x-1).^2 + y.^2).^2 -...
    0.25*2*y./((x+1).^2 + y.^2).^2 + 0.25*2*y;
phiz = @(x,y,z)  2*z;

% norm of gradient
gnorm = @(x,y,z) sqrt(phix(x,y,z).^2 + phiy(x,y,z).^2 + phiz(x,y,z).^2);

% curvature
curv = diff(phix(x,y,z)/gnorm(x,y,z), x) + ...
diff(phiy(x,y,z)/gnorm(x,y,z), y) + diff(phiz(x,y,z)/gnorm(x,y,z), z);
curv = matlabFunction(curv);

% make point cloud
[xs, ys, zs] = place_b_points(phi, Nt, [-5,5], [-4,4], [-4,4], 1e-13, 40);

kappa = curv(xs, ys, zs); % curvature on point cloud

% normal vector
[nx, ny, nz] = nml(xs, ys, zs, phix, phiy, phiz);

% size of expanded domain (Omega)
ellx = 10; elly = 6; ellz = 3;

q = 5; % shape parameter
wm = 15; % positive freq per dim
Nb = (2*wm + 1).^3; % total Fourier basis functions
freq = (-wm):wm;
tic;
% create vectors of frequencies
[xfreq, yfreq, zfreq] = meshgrid(freq, freq, freq);
xfreq = reshape(xfreq, [Nb, 1]);
yfreq = reshape(yfreq, [Nb, 1]);
zfreq = reshape(zfreq, [Nb, 1]);

omega_x = 2*pi/ellx*xfreq;
omega_y = 2*pi/elly*yfreq;
omega_z = 2*pi/ellz*zfreq;

omega_abs = vecnorm([omega_x, omega_y, omega_z], 2, 2);

T = 12; % oscillation width
sDi = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(omega_abs))); % d^(-1/2)


% construct interpolation matrices in Fourier basis
% function value
V = exp(1i*(xs.*omega_x' + ys.*omega_y' + zs.*omega_z'));
% first normal derivative
V1 = 1i*(nx.*omega_x' + ny.*omega_y' + nz.*omega_z').*V;
% second normal derivative
V2 = - (nx.^2.*(omega_x.^2)' + ny.^2.*(omega_y.^2)' + ...
    nz.^2.*(omega_z.^2)' + ...
    2*nx.*ny.*(omega_x.*omega_y)' + ...
    2*nx.*nz.*(omega_x.*omega_z)' + ...
    2*nz.*ny.*(omega_z.*omega_y)').*V;

lap = -omega_abs.^2; % Laplacian in Fourier basis

% points to set u non-zero (query points)

xq = xs(1:mul);
yq = ys(1:mul);
zq = zs(1:mul);

if mul > 1
    valq = rand(mul, 1)*2 -1; % random non-zero value
else
    valq = 1;
end

% interpolation matrix for query points
Vq = exp(1i*(xq.*omega_x' + yq.*omega_y' + zq.*omega_z'));

if use_curv
    % using curvature
    Vs = [lap'.*V - V2 - kappa.*V1; Vq];
    Vs2 =  [V; 0*Vq];
    f = [zeros(Nt, 1); valq];
else
    % set up matrices to solve PDE
    % without using curvature
    Vs = [lap'.*V - V2; V1; Vq];
    Vs2 =  [V; 0*V1; 0*Vq];
    f = [zeros(2*Nt, 1); valq];
end


% multiply by d^(-1/2)
tVs = sDi'.*Vs;
tVs2 = sDi'.*Vs2;

"computing matrices"
% faster to store pre-computed blocks of matrices for larger eigenvalue
% tests
F1 = tVs*tVs';
F2 = tVs*tVs2';
"half"
F2 = F2;
F3 = tVs2*tVs2';
'done matrices'
ind = length(f);
toc;
% scatter3(xs, ys, zs);
% axis('equal');
% hold on;

%% linear search for eigenvalues, plot of H-norm, (Fig 2)
nc = 203; % number of lambda values to check

check = linspace(0.1, 2, nc); % create range of values to check
norms = zeros(nc, 1);
gue = ones(2*Nt + mul, 1);
I = eye(2*Nt + mul);

% kernel matrix method used here: much faster for repeated solves
% ill-conditioning of kernel matrix does not affect plot at this point
% density
for j = 1:nc
    i = check(j);
    F = F1 + i*F2 + i'*F2' + abs(i)^2*F3; % form Phi
    rd = 1./sqrt(vecnorm(F, 2, 2));

    R = sparse(diag(rd)); % simple pre-conditioner
    
    Fp = R*F*R';
    Fp = 1/2*real(Fp + Fp'); 

    beta = R*(Fp\(R*f));
    norms(j) =  (f'*beta); % compute norm squared (note Phi*beta = f)
    if mod(j, 10) == 0
        j
    end
end

% H-norm plot
plot(check,sqrt(abs(norms)), '-'); % plot norms

hold on;

%title("Eigenvalue Test for Genus 2 Surface");

xlabel("$\lambda$");
ylabel("$\mathcal{H}$-Norm of Solution");
set(gca, 'yscale', 'log');
grid on;

%% Newton's method search for eigenvalues (Table 2-5)

% Phi_0, Phi_1, Phi_2 in paper
P0 = F1;
P1 = F2 + F2';
P2 = F3;

lams = []; % store eigenvalue(s)
starts = 0.3; % points to start iteration, 0.3 for Table 2, 4
% starts = 0.63; % Table 3, 5

tic;
for s = starts
    change = 1;
    lam = s
    iter = 0;
    while abs(change) > 1e-12 && iter < 20
   
        F = lam^2*P2 + lam*P1 + P0; % Phi matrix for current lambda
        
        rd = 1./sqrt(vecnorm(F, 2, 2));
        
        R = sparse(diag(rd)); % simple pre-conditioner
        
        Fp = R*F*R';
        Fp = 1/2*real(Fp + Fp'); 
        
        beta = R*(Fp\(R*f));
        f2 = (2*lam*P2 + P1)*beta;
        beta2 = R*(Fp\(R*f2));
    
        % first and second derivative
        dN = real(-beta'*f2);
        ddN = real(-2*beta'*P2*beta + 2*beta2'*f2);
        if ddN < 0  % break if at a point with negative second derivative
            iter = 500000;
        end
        
        change = -dN/ddN;
        lam = lam + change;
        iter = iter + 1
        
    end
    % only save ending value as an eigenvalue if
    %    - second derivative positive
    %    - iteration has converged to within ~1e-4
    %    - at least 1e-4 from other eigenvalues
    if (ddN > 0 && (change < 1e-4)) && (min([1000,abs(lam - lams)]) > 1e-4)
        lams = [lams lam]; % add new eigenvalue
    end
end
toc;
lams = sort(lams)

%% plot an eigenfunction
lams = 0.6263515;
b = lsqminnorm(tVs + lams*tVs2, f); % compute eigenfunction coefficients
u = real(sDi'.*V*b); % compute eigenfunction values on point cloud

% plot
scatter3(xs, ys, zs, 10,u, 'filled')
hold on;

xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
%title("$\lambda="+lams+",\tilde{N}="+Nt+"$ Laplace-Beltrami Eigenfunction")
axis('equal')
colorbar()


function [nx, ny, nz] = nml(x, y, z, phix, phiy, phiz)
    nx = phix(x,y,z); ny = phiy(x,y,z); nz = phiz(x,y,z);
    nm = vecnorm([nx ny nz],2,2);
    nx = nx./nm; ny = ny./nm; nz = nz./nm;
end
