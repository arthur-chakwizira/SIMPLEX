function [c2, c4] = computeCumulants(X, Y, Z, gwf, dt)
%
% INPUTS
%   X, Y, Z           Particle traces
%   gwf               Gradient waveform
%   dt                Time step
%
% OUTPUTS
%   c2                Second-order phase cumulant (equivalent to -beta)
%   c4                Fourth-order phase cumulant (equivalent to zeta)
%

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


%% Compute phase cumulants
c2 = -(1/2)*mean(phase.^2);
c4 = (1/24)*(mean(phase.^4) - 3*mean(phase.^2)^2);