%load trajectories from bin
close all

r_fn = "..\output\test_traj.bin";

selection = 1:100000;
% selection = 50001:100000;
Np = numel(selection);

[Npart, T, Nt, X,Y,Z] =  load_traj_from_bin(r_fn);

X = X(selection,:);
Y = Y(selection,:);
Z = Z(selection,:);



%save X,Y,Z to where you want





%plot some statistics

dx = zeros(Np,1);
dy = zeros(Np, 1);
dz = zeros(Np, 1);
dr = zeros(Np, 1);

for c = 1:Np
    
    dx(c) = range(X(c, :));
    
    dy(c) = range(Y(c, :));
    
    dz(c) = range(Z(c, :));
    
    dr(c) = max( sqrt( (X(c,:)-X(c,1)).^2  + (Y(c,:)-Y(c,1)).^2) + sqrt( (Z(c,:)-Z(c,1)).^2));
    
end

figure('Position',[230 444 1546 420]);
subplot(1,3,1)
histogram(dx*1e6)
title("x")
xlabel("Distance [um]")
subplot(1,3,2)
histogram(dy*1e6)
title("y")
xlabel("Distance [um]")
subplot(1,3,3)
histogram(dr*1e6)
title("xy")
xlabel("Distance [um]")



figure('Color', 'w');
ax = axes(gcf);
for c = randi(min(numel(selection), numel(dx)), 1, 5000)%[1 round(Npart/3) 2*round(Npart/3)  Npart]
    plt = plot3(ax, X(c,:)*1e6, Y(c,:)*1e6, Z(c,:)*1e6, '.-', 'MarkerSize', 2);
    % set(plt, 'MarkerFaceColor', plt.Color, 'MarkerSize', 4)
    hold on
end
xlabel('x [um]')
ylabel('y [um]')
zlabel('z [um]')
axis square


function [Npart, T, Nt, X,Y,Z] = load_traj_from_bin(r_fn)
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
