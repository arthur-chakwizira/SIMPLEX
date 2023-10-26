function [chivXX, chivYY, chivZZ, chivXY, chivXZ, chivYZ, t] = computeChivFromParticleTrace(X, Y, Z, dt, flagSymTime)
%
% INPUTS
%   X, Y, Z         Particle traces
%   dt              Time step
%   flagSymTime     (optional) Boolean indicating if the cross-correlation functions
%                              should be computed on a time grid symmetric around zero 
%                              (They are even in time)
%
% OUTPUTS
%   chiv...         Velocity cross-correlation functions
%   t               Time array (symmetric around zero)
%

if nargin < 5 || isempty(flagSymTime)
    flagSymTime = 1;
end


%%
% Initialize
chivXX = [];
chivYY = [];
chivZZ = [];
chivXY = [];
chivXZ = [];
chivYZ = [];
t      = (1:size(X, 2))'.*dt - dt; % Time array starting in zero

if flagSymTime
    t = [-flip(t(2:end)) ; t];     % Time array symmetric around zero
    t = t(2:end-1);
end

% Compute 
if ~isempty(X)
    vx     = (X(:, 2:end) - X(:, 1:end-1))./dt;
    chivXX = getChivFromVelocities(vx, vx, flagSymTime);
end

if ~isempty(Y)
    vy     = (Y(:, 2:end) - Y(:, 1:end-1))./dt;
    chivYY = getChivFromVelocities(vy, vy, flagSymTime);
end

if ~isempty(Z)
    vz     = (Z(:, 2:end) - Z(:, 1:end-1))./dt;
    chivZZ = getChivFromVelocities(vz, vz, flagSymTime);
end

if ~isempty(X) && ~isempty(Y) 
    chivXY = getChivFromVelocities(vx, vy, flagSymTime);
end

if ~isempty(X) && ~isempty(Z) 
    chivXZ = getChivFromVelocities(vx, vz, flagSymTime);
end

if ~isempty(Y) && ~isempty(Z) 
    chivYZ = getChivFromVelocities(vy, vz, flagSymTime);
end


%% Set as zero function if empty
if isempty(chivXX)
    chivXX = zeros(size(t));
end
if isempty(chivYY)
    chivYY = zeros(size(t));
end
if isempty(chivZZ)
    chivZZ = zeros(size(t));
end
if isempty(chivXY)
    chivXY = zeros(size(t));
end
if isempty(chivXZ)
    chivXZ = zeros(size(t));
end
if isempty(chivYZ)
    chivYZ = zeros(size(t));
end

end

%%
function chiv = getChivFromVelocities(v1, v2, flagSymTime)
    if nargin < 3 || isempty(flagSymTime)
        flagSymTime = 0;
    end

    chiv = mean(v1.*v2(:, 1), 1)';
    
    % Perform cosine transform (see J. Brabec's undulation paper)
    chiv(1) = chiv(1)/2; 
    
    if flagSymTime
        % Even in time
        chiv = [flip(chiv(2:end)) ; chiv]; 
    end
end