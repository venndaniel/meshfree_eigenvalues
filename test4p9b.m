% Code for the Tables 12 and 13 in Subsection 4.9 (Figure 6)

%% Problem setup
% intially run in MATLAB R2023a

Nt = 1200; % max number of points (1200 used in paper)
mul = 27; % max multiplicity to test for

% set plot text to LaTex
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

rng('default') % default rng for replicability

% create (non-uniform) random points on sphere 
% by placing random points in an 8 by 8 box and projecting to 
% unit sphere

[xs, ys, zs] = sphere_cloud(Nt, 40);


% set normal vectors
nx = xs; ny = ys; nz = zs;

% size of expanded domain (Omega)
ellx = 4; ellz = 4; elly = 4;

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



% points to set u non-zero (query points a_j)

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
valq = 1./(1:mul)';
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

% scatter3(xs, ys, zs);
% axis('equal');
% hold on;

%% ratio test for Tables 12 & 13
N1 = 700; % smaller N for ratio
N2 = round(10/9*N1); % larger N
% eigenvalue to test
%i = 56; % Table 12
i = 156; % Table 13
Ff = F1 + i*F2 + i'*F2' + abs(i)^2*F3; % form Phi

%muls = 14:17; % Table 12
muls = 24:27; % Table 13

for mul = muls
    
    % indices for matrices with N2 and N1 points
    sind = [1:N2,(Nt+1:Nt+N2),(Nt*2+1:Nt*2+mul)];
    sinds = [1:N1,(Nt+1:Nt+N1),(Nt*2+1:Nt*2+mul)];

    F = Ff(sind, sind); % larger Phi matrix
   
    Fs = Ff(sinds,sinds); % smaller Phi matrix
    rd = 1./sqrt(vecnorm(F, 2, 2));
    rds = 1./sqrt(vecnorm(Fs, 2, 2));

    R = sparse(diag(rd)); % simple pre-conditioner/rescaling
    Rs = sparse(diag(rds)); % simple pre-conditioner/rescaling
    ff = f(sind);
    Fp = R*F*R';
    Fp = 1/2*real(Fp + Fp'); 
    beta = R*(Fp\(R*ff));

    Fps = Rs*Fs*Rs';
    Fps = 1/2*real(Fps + Fps'); 
    fs = f(sinds);
    betas = Rs*(Fps\(Rs*fs));

    mul + " " + sqrt(abs((ff'*beta)./(fs'*betas))) % compute norm ratio
    %abs(((ff'*beta) - (fs'*betas))./((fs'*betas)))
end
