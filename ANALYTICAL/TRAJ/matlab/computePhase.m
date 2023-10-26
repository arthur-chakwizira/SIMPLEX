function phase = computePhase(X, Y, Z, gwf, dt)


gammaConst = 42.576e6*2*pi;


%% Compute phase evolution
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
