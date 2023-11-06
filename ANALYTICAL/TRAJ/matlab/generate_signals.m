
%load trajectory from bin and compute signals
dt = 1e-4;

gwf_fn = "gwf.mat";

load(gwf_fn, 'gwf') %load the waveform (s)

trajectory_fn = "traj_intra_extra_spheres_d50_k0_gamma_dist.bin";

[Npart, T, Nt, X,Y,Z] =  load_traj_from_bin(trajectory_fn); %load trajectory

signal = zeros(size(gwf, 1), 1);
for i = 1:size(gwf, 1)
    
phase = computePhase(X, Y, Z, squeeze(gwf(i, :, :)), dt); %compute accumulated phase assuming gwf is N_time x 3
signal(i) = (1/numel(phase))*sum(cos(phase)); %compute signal
end


plot(signal, 'o-')




%%

function [Npart, T, Nt, X,Y,Z] =  load_traj_from_bin(r_fn)
if isfile(r_fn)
    fileID = fopen(r_fn);
    Npart = fread(fileID, [1,1], 'int32');
    T = fread(fileID, [1,1], 'single');
    Nt = fread(fileID, [1,1], 'int32');
    tmp_x = fread(fileID, [Npart*Nt,1], 'single');
    tmp_y = fread(fileID, [Npart*Nt,1], 'single');
    tmp_z = fread(fileID, [Npart*Nt,1], 'single');
    fclose(fileID);
else
    error("File " + r_fn + " not found.")
end

X = zeros(Npart, Nt);
Y = zeros(size(X));
Z = zeros(size(X));

for c = 1:Npart
    idx = ((c-1)*Nt+1):c*Nt;
    X(c,:) = tmp_x(idx);
    Y(c,:) = tmp_y(idx);
    Z(c,:) = tmp_z(idx);
end

end

function phase = computePhase(X, Y, Z, gwf, dt)


gammaConst = 42.576e6*2*pi;


%Compute phase evolution
Nt    = size(gwf, 1);
Npart = size(X, 1); % Number of particles

sum_g_times_r = zeros(Npart, 1); % Initialize phase integral
for i = 1:Nt
    sum_g_times_r = sum_g_times_r ...
        + squeeze(gwf(i, 1)).*squeeze(X(:, i)) ...
        + squeeze(gwf(i, 2)).*squeeze(Y(:, i)) ...
        + squeeze(gwf(i, 3)).*squeeze(Z(:, i));
end
phase = gammaConst.*sum_g_times_r.*dt; % Evaluate phase integral
end


