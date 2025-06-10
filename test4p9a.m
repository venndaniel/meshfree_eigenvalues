% Code for the first test of Subsection 4.9 (Figure 6)

%% Problem setup
% run this section to create the necessary matrices to run the experiments
% in the other section
% intially run in MATLAB R2023a

Nt = 400; % number of random points on sphere

mul = 5; % multiplicity to test for (1, 3, 5 needed for Fig. 6)

% set plot text to LaTex
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

rng('default') % default rng for replicability

% create (non-uniform) random points on sphere 
% by placing random points in an 8 by 8 box and projecting to 
% unit sphere
xs = zeros(Nt, 1);
ys = xs;
zs = ys;
fou = 1;

r = rand(1, 3)*8 - 4;
xs(1) = r(1);
ys(1) = r(2);
zs(1) = r(3);
nm = sqrt(xs(1)^2 + ys(1)^2 + zs(1)^2);
xs(1) = xs(1)/nm;
ys(1) = ys(1)/nm;
zs(1) = zs(1)/nm;


check = 40;

while fou < Nt
    r = rand(1, 3)*8 - 4;
    x = r(1);
    y = r(2);
    z = r(3);
    nm = sqrt(x^2 + y^2 + z^2);
    x = x/nm;
    y = y/nm;
    z = z/nm;

    ch = 1;
    dist = min(vecnorm([xs-x, ys-y, zs-z], 2, 2));
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
        if distt > dist
            x = xt;
            y = yt;
            z = zt;
            dist = distt;
        end
        ch = ch + 1;
    end
    fou = fou + 1;
    xs(fou) = x;
    ys(fou) = y;
    zs(fou) = z;
end

% set normal vectors
nx = xs;
ny = ys;
nz = zs;

% size of expanded domain (Omega)
ellx = 4;
ellz = 4;
elly = 4;

q = 4; % shape parameter
wm = 15; % positive freq per dim
Nb = (2*wm + 1).^3; % total Fourier basis functions
freq = (-wm):wm;

% create vectors of frequencies
[xfreq, yfreq, zfreq] = meshgrid(freq, freq, freq);
xfreq = reshape(xfreq, [Nb, 1]);
yfreq = reshape(yfreq, [Nb, 1]);
zfreq = reshape(zfreq, [Nb, 1]);

omega_x = 2*pi/ellx*xfreq;
omega_y = 2*pi/elly*yfreq;
omega_z = 2*pi/ellz*zfreq;

omega_abs = vecnorm([omega_x, omega_y, omega_z], 2, 2);

T = 4; % oscillation width
sDi = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(omega_abs))); % d^(-1/2)


% construct interpolation matrices in Fourier basis
V = exp(1i*(xs.*omega_x' + ys.*omega_y' + zs.*omega_z'));
V1 = 1i*(nx.*omega_x' + ny.*omega_y' + nz.*omega_z').*V;
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
% 
% xq = 0;
% yq = 0;
% zq = 1;

% project query points to sphere
nq = vecnorm([xq, yq, zq], 2, 2);
xq = xq ./ nq;
yq = yq ./ nq;
zq = zq ./ nq;
valq = rand(mul, 1)*2 -1; % random non-zero value
% interpolation matrix for query points
Vq = exp(1i*(xq.*omega_x' + yq.*omega_y' + zq.*omega_z'));

% set up matrices to solve PDE
Vs = [lap'.*V - V2; V1; Vq];
Vs2 =  [V; 0*V1; 0*Vq];
tVs = sDi'.*Vs;
tVs2 = sDi'.*Vs2;
f = [zeros(Nt, 1); zeros(Nt,1); valq];

"computing matrices"
% faster to store pre-computed blocks of matrices for larger eigenvalue
% tests
F1 = tVs*tVs'; % Phi_0
F2 = tVs*tVs2'; % Phi_1 / 2
F3 = tVs2*tVs2'; % Phi_2
'done matrices'
ind = length(f);

%% linear search for eigenvalues, plot
nc = 403; % number of lambda values to check

check = linspace(-1, 8, nc); % create range of values to check
norms = zeros(nc, 1);

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
    hold on;
    if mod(j, 10) == 0
        j
    end
end

cor = ((0:6).*(1:7))'; % correct eigenvalues 
semilogy(check,sqrt(abs(norms)), '-'); % plot norms
hold on;
lp = linspace(min(sqrt(norms)), max(sqrt(norms)), 10);
%onorms = sqrt(abs(norms));
%semilogy(cor.*ones(1, 10), lp, 'r-'); % plot true eigenvalues
%title("Eigenvalue Multiplicity Test");

xlabel("$\lambda$");
ylabel("$\mathcal{H}$-Norm of Solution");
set(gca, 'yscale', 'log');
grid on;

%legend(["$\tilde{N}_a=" + 1 +"$", "$\tilde{N}_a=" + 3 +"$","$\tilde{N}_a="
%+ 5 +"$"]);
