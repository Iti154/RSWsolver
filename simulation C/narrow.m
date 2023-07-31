clear all; clc; close all;

%% Physical constants

Omega = 2 * pi / 86400;
R_earth = 5.447890490694152e+07;
g = 9.81;

%% Solver parameters

lat = 45;                          % Mid-latitude (deg)
Lx_norm = 120;                           % Channel length (Rd)
Ly_norm = 4;                          % Channel width (Rd)

H0 = 80;                          % Channel depth (m)
A = 1;                             % Wave amplitude
nl = 0;                            % Linear or non-linear
ah = 0.e-3;                         % Diffusion constant
bf = 0.e-10;                        % Bottom friction

save_data = 1;                     % Save data or not

gifname = 'C_narrow.gif';   % Name of gif file
dataname = 'C_narrow.mat';  % Name of data file


% Coriolis parameters
f0 = 2 * Omega * sind(lat);
beta = 2 * Omega * cosd(lat) / R_earth;

% Scales

c = sqrt(g*H0);
T = 1/f0;
u = A*c/H0;
L = c/f0;

b = beta*L/f0;

dx = L/10;                          % Spatial resolution x-axis (m)
dy = L/10;                          % Spatial resolution y-axis (m)
dt = T/100;                            % Time step

Ly = Ly_norm*L;
Lx = Lx_norm*L;

time_f = 100000;

N = fix(time_f*T/dt);                        % Number of iterations
plot_rate = T/dt/0.05;        % sampling every plot_rate iterations

% Grid properties
% Zonal and meridional vectors
x = -Lx / 2 : dx : (Lx / 2 + dx);
y = -Ly / 2 : dy : (Ly / 2 + dy);

% Grid size
nx = length(x);
ny = length(y);

fU = repmat(f0 + beta * (y - 0.5 * dy), nx, 1);
fV = repmat(f0 + beta * y, nx, 1);

% Relevant matrix indices.
i = 2 : (nx - 1);
im = i - 1;
ip = i + 1;
im(1) = nx - 1; ip(end) = 2; % Periodical boundaries in zonal direction
j = 2 : (ny - 1);
jm = j - 1;
jp = j + 1;

% Other variables
d2x = 2 * dx;
d4x = 4 * dx;
dx2 = dx^2;
d2y = 2 * dy;
d4y = 4 * dy;
dy2 = dy^2;
tm = 1; tn = 2; tp = 3; % Time indices
delt = dt;
U = zeros(nx, ny, 3);
V = zeros(nx, ny, 3);
H = zeros(nx, ny, 3);
grav = zeros(nx, ny);
vel = zeros(nx, ny);
adv = zeros(nx, ny);
xadv1 = zeros(nx, ny);
xadv2 = zeros(nx, ny);
xadv = zeros(nx, ny);
diff_x = zeros(nx, ny);
diff_y = zeros(nx, ny);
DH1 = zeros(nx, ny);
DH2 = zeros(nx, ny);



% Stability checks
staba = 1 / dx2 + 1 / dy2;
stabc = bf / H0^2;
stabd = (stabc^2 + f0^2) / stabc;
stab1 = (sqrt(stabc^2 + 16 * g * H0 * staba) - stabc) / (4 * g * H0 * staba);
stab2 = (sqrt(stabd^2 + 8 * g * H0 * staba) - stabd) / (2 * g * H0 * staba);
stab3 = 1 / stabd;
stable = dt < min([stab1 stab2 stab3]);

% For data logging
data_idx = 1;
data_logs = 1 + floor(N / plot_rate) + 1 * (mod(N, plot_rate) ~= 0);
if (save_data)
    dataf = matfile(dataname, 'writable', true);
    dataf.H = zeros(ny - 2, nx - 2, data_logs);
    dataf.U = zeros(ny - 2, nx - 2, data_logs);
    dataf.V = zeros(ny - 2, nx - 2, data_logs);
    dataf.y = y(j)/L;
    dataf.x = x(i)/L;
end

%% Initial conditions 
D = 10;

h0_y = ones(size(y));
h0_x = (heaviside(x+L*D)-heaviside(x-L*D));

H(:, :, tm) = h0_x' * h0_y + H0;
U(:, :, tm) = 0;
V(:, :, tm) = 0;

% Periodic boundaries in zonal direction.
H(1, :, tm) = H(nx - 1, :, tm); H(nx, :, tm) = H(2, :, tm);
U(1, :, tm) = U(nx - 1, :, tm); U(nx, :, tm) = U(2, :, tm);
V(1, :, tm) = V(nx - 1, :, tm); V(nx, :, tm) = V(2, :, tm);

% Rigid boundaries in meridional direction.
V(:, 1, tm) = 0; V(:, ny - 1, tm) = 0;

% Setting current and past same to be the same for first iteration.
H(:, :, tn) = H(:, :, tm);
U(:, :, tn) = U(:, :, tm);
V(:, :, tn) = V(:, :, tm);

if ~stable
    disp('Vlues do not satisfy stability conditions');
end
hold off;


%% Time stepping loop
lev1 = [-2:0.1:2];

for t = 0 : N

    % Plotting
    if (mod(t, plot_rate) == 0) || t == N
        cla;
        uu = u*H(i, j, tn);
        VV = V(i, j, tn)'./uu';
        pcolor(x(i)/L,y(j)/L,H(i,j,tn)'-H0)
        shading interp
        
        colormap(jet);
        colorbar;
        
        xlabel('x');
        ylabel('y');
        title(sprintf('t = %g',t * dt /T))
        xlim([-Lx_norm/2 Lx_norm/2])
        ylim([-Ly_norm/2 Ly_norm/2])
        caxis([-1.2 1.2])
        xlabel('x');
        xlabel('y');
       
        drawnow;
        
        if (save_data)
            % GIF creation
            frame = getframe(1);
            img = frame2im(frame);
            [imind, cm] = rgb2ind(img, 256);
            if (t == 0)
                imwrite(imind, cm, gifname, 'gif', 'Loopcount', inf, 'DelayTime', 0);
            else
                imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
            end 
            
            % Data logging
            dataf.time(1, data_idx) = t * dt /T;
            dataf.H(:,:, data_idx) = H(i, j, tn)'-H0;
            dataf.U(:,:, data_idx) = U(i, j, tn)'./uu';
            dataf.V(:,:, data_idx) = V(i, j, tn)'./uu';
            data_idx = data_idx + 1;
            
        end
    end
    
    % Advancing H.
    H(i, j, tp) = H(i, j, tm) - delt * ((U(i, j, tn) - U(im, j, tn)) / dx + (V(i, j, tn) - V(i, jm, tn)) / dy);
    H(1, :, tp) = H(nx - 1, :, tp); H(nx, :, tp) = H(2, :, tp); % Periodic boundaries in zonal direction
                                                
    % Advancing U.
    grav(i, j) = (H(ip, j, tn) - H(i, j, tn)) .* (H(ip, j, tn) + H(i, j, tn)) / d2x; % Gravitation
    vel(i, j) = (V(i, j, tn) + V(i, jm, tn) + V(ip, j, tn) + V(ip, jm, tn)) / 4; % Average velocity
    adv(i, j) = ((U(ip, j, tn) + U(i, j, tn)).^2 ./ H(ip, j, tn) - (U(i, j, tn) + U(im, j, tn)).^2 ./ H(i, j, tn)) / d4x; % Advection
    % Cross advection
    DH1(i, j) = H(i, j, tn) + H(i, jp, tn) + H(ip, j, tn) + H(ip, jp, tn);
    DH2(i, j) = H(i, j, tn) + H(i, jm, tn) + H(ip, j, tn) + H(ip, jm, tn);
    xadv1(i, j) = (U(i, jp, tn) + U(i, j, tn)) .* (V(i, j, tn) + V(ip, j, tn)) ./ DH1(i, j);
    xadv1(:, ny - 1) = 0;
    xadv2(i, j) = (U(i, j, tn) + U(i, jm, tn)) .* (V(ip, jm, tn) + V(i, jm, tn)) ./ DH2(i, j);
    xadv2(:, 2) = 0;
    xadv(i, j) = (xadv1(i, j) - xadv2(i, j)) / dy;
    % Diffusion
    diff_x(i, j) = (U(im, j, tn) - 2 * U(i, j, tn) + U(ip, j, tn)) / dx2;
    diff_y(i, j) = (U(i, jm, tn) - 2 * U(i, j, tn) + U(i, jp, tn)) / dy2;
    diff_y(:, 2) = 0; diff_y(:, ny - 1) = 0;
    % Finite difference equation
    U(i, j, tp) = U(i, j, tm) - delt *(g * grav(i, j) - fU(i, j) .* vel(i, j) + bf * U(i, j, tn) - ah * (diff_x(i, j) + diff_y(i, j))) - ...
                  nl * delt * (adv(i, j) + xadv(i, j));
    U(1, :, tp) = U(nx - 1, :, tp); U(nx, :, tp) = U(2, :, tp); % Periodic boundaries in zonal direction
    
    % Advancing V.
    grav(i, j) = (H(i, jp, tn) - H(i, j, tn)) .* (H(i, jp, tn) + H(i, j, tn)) / d2y; % Gravitation
    vel(i, j) = (U(i, j, tn) + U(im, j, tn) + U(i, jp, tn) + U(im, jp, tn)) / 4; % Average velocity
    adv(i, j) = ((V(i, jp, tn) + V(i, j, tn)).^2 ./ H(i, jp, tn) - (V(i, j, tn) + V(i, jm, tn)).^2 ./ H(i, j, tn)) / d4y; % Advection
    % Cross advection
    DH1(i, j) = H(i, j, tn) + H(i, jp, tn) + H(ip, j, tn) + H(ip, jp, tn);
    DH2(i, j) = H(i, j, tn) + H(im, j, tn) + H(i, jp, tn) + H(im, jp, tn);
    xadv1(i, j) = (U(i, jp, tn) + U(i, j, tn)) .* (V(i, j, tn) + V(ip, j, tn)) ./ DH1(i, j);
    xadv2(i, j) = (U(im, j, tn) + U(im, jp, tn)) .* (V(i, j, tn) + V(im, j, tn)) ./ DH2(i, j);
    xadv(i, j) = (xadv1(i, j) - xadv2(i, j)) / dx;
    % Diffusion
    diff_x(i, j) = (V(im, j, tn) - 2 * V(i, j, tn) + V(ip, j, tn)) / dx2;
    diff_y(i, j) = (V(i, jm, tn) - 2 * V(i, j, tn) + V(i, jp, tn)) / dy2;
    diff_y(:, 2) = 0; diff_y(:, ny - 2) = 0;
    % Finite difference equation
    V(i, j, tp) = V(i, j, tm) - delt *(g * grav(i, j) + fV(i, j) .* vel(i, j) + bf * V(i, j, tn) - ah * (diff_x(i, j) + diff_y(i, j))) - ...
                  nl * delt * (adv(i, j) + xadv(i, j));
    V(1, :, tp) = V(nx - 1, :, tp); V(nx, :, tp) = V(2, :, tp); % Periodic boundaries in zonal direction
    V(:, 1, tp) = 0; V(:, ny - 1, tp) = 0; % Rigid boundaries in meridional direction
    
    % Robert-Asselin Filter.
    if (t > 5)
        U(:, :, tn) = U(:, :, tn) + 0.05 * (U(:, :, tp) - 2 * U(:, :, tn) + U(:, :, tm));
        V(:, :, tn) = V(:, :, tn) + 0.05 * (V(:, :, tp) - 2 * V(:, :, tn) + V(:, :, tm));
        H(:, :, tn) = H(:, :, tn) + 0.05 * (H(:, :, tp) - 2 * H(:, :, tn) + H(:, :, tm));
    end
    
    % Moving on in time.
    delt = 2 * dt;
    tm = mod(tm, 3) + 1;
    tn = mod(tn, 3) + 1;
    tp = mod(tp, 3) + 1;
end
