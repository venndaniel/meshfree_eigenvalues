% code for Subsection 4.2 (Laplace--Beltrami on a Sphere)

%% Problem setup
% run this section to create the necessary matrices to run the experiments
% intially run in MATLAB R2023a, tested again in R2024b

Nt = 650; % number of random points on sphere

mul = 1; % multiplicity to test for

% set plot text to LaTeX
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

rng('default') % default rng for replicability

[xs, ys, zs] = sphere_cloud(Nt, 40);

% set normal vectors
nx = xs; ny = ys; nz = zs;

% size of expanded domain (Omega)
ellx = 4; elly = 4; ellz = 4;

q = 4; % shape parameter
wm = 15; % positive freq per dim
Nb = (2*wm + 1).^3; % total Fourier basis functions
freq = (-wm):wm;

tic; % start timing matrix formation
% create vectors of frequencies
[xfreq, yfreq, zfreq] = meshgrid(freq, freq, freq);
xfreq = reshape(xfreq, [Nb, 1]);
yfreq = reshape(yfreq, [Nb, 1]);
zfreq = reshape(zfreq, [Nb, 1]);

% create omega vectors
omega_x = 2*pi/ellx*xfreq; 
omega_y = 2*pi/elly*yfreq; 
omega_z = 2*pi/ellz*zfreq;

% norm of omega
omega_abs = vecnorm([omega_x, omega_y, omega_z], 2, 2);

T = 4; % oscillation width for phi function
sDi = 1./(exp(q*sqrt(2*pi/T)) + exp(q*sqrt(omega_abs))); % d^(-1/2)


% construct interpolation matrices in Fourier basis
% function values
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

% points to set u non-zero (a_j points)
xq = 0;
yq = 0;
zq = 1; % use (0, 0, 1) if multiplicity is 1
if mul > 1 % otherwise, use first mul points in point cloud
    xq = xs(1:mul);
    yq = ys(1:mul);
    zq = zs(1:mul);
end

% values to set a_j points (b_j)
if mul == 1
    valq = 1;
else
    valq = rand(mul, 1)*2 -1; % random non-zero value
end

% interpolation matrix for a_j points
Vq = exp(1i*(xq.*omega_x' + yq.*omega_y' + zq.*omega_z'));

% set up matrices to solve PDE
Vs = [lap'.*V - V2; V1; Vq];
Vs2 =  [V; 0*V1; 0*Vq]; % matrix scaled by lambda
tVs = sDi'.*Vs; % multiply basis functions by d^(-1/2)
tVs2 = sDi'.*Vs2;
f = [zeros(Nt, 1); zeros(Nt,1); valq];


"computing matrices"
% faster to store pre-computed blocks of matrices for larger eigenvalue
% tests
% denoted Phi_0, Phi_1, Phi_2 in paper

P0 = tVs*tVs';
P1 = tVs*tVs2';
P1 = P1 + P1';
P2 = tVs2*tVs2';
"done matrices"
"matrix formation: " + toc

%% Newton's Method
ind = length(f);
% change = 1;
lams = []; % store eigenvalues
starts = ((1:30)/2).^2; % starting points for Newton's method
%starts = 56.25
%starts = 800;
tic;
for s = starts
    change = 1;
    lam = s % print current iteration
    iter = 0;
    while abs(change) > 1e-13 && iter < 20 % stop condition for Newton's
   
        F = lam^2*P2 + lam*P1 + P0; % Phi matrix for current lambda
        
        rd = 1./sqrt(vecnorm(F, 2, 2));
        
        R = sparse(diag(rd)); % simple pre-conditioner / rescaling matrix
        % more stable if rows/columns are of similar magnitude,
        
        Fp = R*F*R'; % ensure matrix remains symmetric after rescaling
        % Fp = 1/2*real(Fp + Fp') already algebraically, but should be
        % forced numerically
        Fp = 1/2*real(Fp + Fp'); 
        
        % compute first and second derivatives
        beta = R*(Fp\(R*f));
        f2 = (2*lam*P2 + P1)*beta;
        beta2 = R*(Fp\(R*f2));
    
        dN = real(-beta'*f2); % first derivative
        ddN = real(-2*beta'*P2*beta + 2*beta2'*f2); % second derivative

        if ddN < 0 % break if at a point with negative second derivative
            iter = 500000;
        end

        change = -dN/ddN;
        lam = lam + change; % update lambda
        iter = iter + 1;
    end
    % only save ending value as an eigenvalue if
    %    - second derivative positive
    %    - iteration has converged to within ~1e-4
    %    - at least 1e-4 from other eigenvalues
    if (ddN > 0 && (change < 1e-4)) && (min([1000,abs(lam - lams)]) > 1e-4)
        lams = [lams lam]; % add new eigenvalue
    end
end
"Newton's method time: " + toc
lams = sort(lams);

min(abs(lams - 56))/56 % relative error for lambda = 56

lams'
