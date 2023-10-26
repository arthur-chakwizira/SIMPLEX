%analytical substrate generator for the SIMPLEX simulator
%make substrates using equations

%regular
close all

regular = 1;
random_pack = 0;
gamma_dist = 0;


if regular
    clear all
    save_geo = 1;

%     mean_ds = 5;%[1 2 3 4 5 6 8 10 12 14 16 18 20];
    %     vs = 0.25e-6;
    mean_ds = 10;
    for c_d = 1:numel(mean_ds)
        
        mean_d = mean_ds(c_d);
        r = (mean_d/2)*1e-6;
        
        min_r = 2.5*r/5;%r/40;%default is 2*r/5;%30*r/2
        vs = r/10;%r/50;%min_r/2;r/5;%
        
        name = "..\substrates\cylinders_d" + num2str(mean_d) + "_regular";
        fn = name +".bin";
        
        Ncells_1D = 2;%14; %cells in one direction
        
%         x = -((2*r+min_r)*Ncells_1D):(2*r+min_r):((2*r+min_r)*Ncells_1D);
%         y = x;
%         z = y;
        
        max_x = (2*r+min_r)*Ncells_1D;
        max_y = max_x;
%         max_z = max_y;
        xc = 0:(2*r+min_r-5*vs):max_x; %want to make sure 0 is centre
        xc = [-flip(xc(2:end)) xc];
        yc = 0:(2*r+ min_r):max_y; %want to make sure 0 is centre
        yc = [-flip(yc(2:end)) yc];
        zc = yc;
        
        max_x = max(xc); min_x = min(xc);
        max_y = max(yc); min_y = min(yc);
        max_z = max(zc); min_z = min(zc);
        
        [Xc,Yc] = meshgrid(xc,yc);
        Yc(:,2:2:end) = Yc(:, 2:2:end)+ r+0.5*min_r; %for hexagonal packing

        Xc = Xc(:); %substrate centres
        Yc = Yc(:);
        
        figure;
        ax = axes();
        hold(ax, 'on')
        
        Ncells_3D = numel(Xc)^2;
        
        centres = zeros(Ncells_3D, 2);
        radii = zeros(Ncells_3D, 1);
        
        cell_counter = 0;
        
        for c_c = 1:(numel(Xc))
            tmp_x = Xc(c_c);%+min_r+r; %x(c_x) is not the centre
            tmp_y = Yc(c_c);
            xx = linspace(tmp_x-r, tmp_x+r);
                if (tmp_x<=max_x & tmp_x>=min_x & tmp_y <= max_y & tmp_y >= min_y)
                    yy = linspace(tmp_y-r, tmp_y+r);
                    o_plus = sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                    o_plus = real(o_plus);
                    %             if any(~isreal(o_plus)); error('Sir'); end
                    o_minus = -sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                    o_minus = real(o_minus);
                    fill(xx, o_plus, 'k');
                    fill(xx, o_minus, 'k');
                    
                    
                    cell_counter = cell_counter+1;
                    centres(cell_counter, :) = [tmp_x, tmp_y];
                    radii(cell_counter) = r;
                end
        end
        
        %         max_x = max(centres(:,1));
        %         min_x = min(centres(:,1));
        %         max_y = max(centres(:,2));
        %         min_y = min(centres(:,2));
        
        centres = centres(1:cell_counter, :);
        radii = radii(1:cell_counter);
        axis equal
        xlim(ax, [min_x max_x])
        ylim(ax, [min_y max_y])
        
        return
        %         movegui('center')
        
        f1 = compute_f1(centres, radii, max_x)% sum(pi*radii.^2)/((2*max(x))^2)
%         return
        [tab, cell_idx] = generate_lookup_table(centres, radii, max_x, max_y, vs);
        if numel(unique(tab)) == numel(tab); disp('All good'); end
       if save_geo; save_substrate_file(centres, radii, tab, cell_idx, f1, vs, max_x, max_y, max_z, fn, name); end
        
    end
    
end

return

 
% R = move_particles_in_substrate(centres, radii, max(Xc), max(Yc), max(max_z));




%%
if random_pack
    
    clear all
    %     vs = 0.25e-6;
    mean_ds = [1 2 3 4 5 6 8 10 12 14 16 18 20];
    for c_d = 1%:numel(mean_ds)
        
        d = mean_ds(c_d);
        r = (d/2)*1e-6;
        
        min_r = 2*r/5;%20*r/5;%2*r/5;%
        vs = r/5;%min_r/2;
%         
        name = "C:\Users\Arthur\source\repos\PaSiD\substrates\analytical_v2\cylinders_d" + num2str(d) + "_rand_pack";
%         name = "C:\\Users\\Arthur\\Documents\\LUND_UNIVERSITY\\PHD\PAPER_6\\data\\substrate\\rand_pack_d10";
        fn = name +".bin";
                figure;
                hold on
        
        
        %random pack
        N = 40; %cells in one direction
        
        min_x = -(r+min_r/2)*N;
        max_x = (r+min_r/2)*N;
        max_y = max_x;
        max_z = max_y;
        
        Ncells = 1000;
        centres = zeros(Ncells,2); %centre coordinates of existing cells
        radii = zeros(Ncells, 1);
        
        counter = 0;
        
        while counter < Ncells
            tmp_x = min_x + 2*rand*max_x;
            xx = linspace(tmp_x-r, tmp_x+r);
            tmp_y = min_x + 2*rand*max_x;
            
            out_of_world = any([abs(tmp_x-min_x), abs(max_x-tmp_x), abs(max_x-tmp_y),  abs(tmp_y-min_x)   ] <= (r+0.5*min_r));
            
            if out_of_world; continue; end
            
            if counter > 0
                min_dist = ((tmp_x-centres(1:counter, 1)).^2  + (tmp_y-centres(1:counter, 2)).^2 );
                if any(min_dist <= (2*r + min_r)^2)
                    continue
                else
                    counter = counter +1;
                                        o_plus = sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                                        o_minus = -sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                                        %                 if any(~isreal([o_plus o_minus])); error('Sir'); end
                                        fill(xx, real(o_plus), 'k');
                                        fill(xx, real(o_minus), 'k');
                    centres(counter,:) = [tmp_x, tmp_y];
                    radii(counter) = r;
                end
            else
                                o_plus = sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                                o_minus = -sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                                fill(xx, real(o_plus), 'k');
                                fill(xx, real(o_minus), 'k');
                centres(1,:) = [tmp_x, tmp_y];
                radii(1) = r;
                counter = counter +1;
            end
            
        end
        
        f1 = sum(pi*radii.^2)/((2*max_x)^2)
        
        [tab, cell_idx] = generate_lookup_table(centres, radii, max_x, max_y, vs);
        
        save_substrate_file(centres, radii, tab, cell_idx, f1, vs, max_x, max_y, max_z, fn, name);
        
                    xlim([min_x max_x])
                    ylim([min_x max_x])
%                 axis square
        %         movegui('center')
    end
end


%%

if gamma_dist
    clear all
    %     vs = 0.25e-6;
%     mean_ds = [1 2 3 4 5 6 7 8 9 10 12 14];
    mean_ds = [1 3  5 7 9 11];
    for c_d = 1:numel(mean_ds)
        
        mean_d = mean_ds(c_d);
        mean_r = (mean_d/2)*1e-6;
        
        %         min_r = 0.5e-6;
        
        %         r = (d/2)*1e-6;
        
        min_r = 2*mean_r/10;%0.5e-6;
        vs = min_r/2;
        
        name = "..\substrates\cylinders_d" + num2str(mean_d) + "_gamma_dist";
        fn = name +".bin";
        
        figure;
        hold on
        
        N = 30; %cells in one direction
        
        min_x = -(mean_r + min_r/2)*N;
        max_x = (mean_r + min_r/2)*N;
        max_y = max_x;
        max_z = max_y;
        
        Ncells = 2000;
        centres = zeros(Ncells,2); %centre coordinates of existing cells
        radii = zeros(Ncells, 1);
        
        counter = 0;
        
        while counter < Ncells
            r = gamrnd(5,mean_r/5);
            %     r = round(r, 2);
            if r < min_r; continue; end %minimum step size is 0.35e-6 m.
            
            
            success = false;
            for trial = 1:10
                
                tmp_x = min_x + 2*rand*max_x;
                xx = linspace(tmp_x-r, tmp_x+r);
                tmp_y = min_x + 2*rand*max_x;
                
                out_of_world = any([abs(tmp_x-min_x), abs(max_x-tmp_x), abs(max_x-tmp_y),  abs(tmp_y-min_x)   ] <= (r+ 0.5*min_r));
                
                if out_of_world; continue; end
                
                if counter > 0
                    %         if ismember(r, radii(1:counter)); continue; end
                    min_dist = ((tmp_x-centres(1:counter, 1)).^2  + (tmp_y-centres(1:counter, 2)).^2 ); %distance to all centres
                    if any(min_dist <= ((radii(1:counter)+ r + min_r)).^2)
                        continue
                    else
                        counter = counter +1;
                        o_plus = sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                        o_minus = -sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                        %                 if any(~isreal([o_plus o_minus])); error('Sir'); end
                        fill(xx, real(o_plus), 'k');
                        fill(xx, real(o_minus), 'k');
                        centres(counter,:) = [tmp_x, tmp_y];
                        radii(counter) = r;
                        success = true;
                    end
                else
                    o_plus = sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                    o_minus = -sqrt( r^2 - (xx - tmp_x).^2) + tmp_y;
                    fill(xx, real(o_plus), 'k');
                    fill(xx, real(o_minus), 'k');
                    centres(1,:) = [tmp_x, tmp_y];
                    counter = counter +1;
                    radii(1) = r;
                    success = true;
                end
                if success; break;else; continue; end
            end
                        if success;  disp(counter); end
        end
        
        xlim([min_x max_x])
        ylim([min_x  max_x])
        axis square
        movegui('center')
        
        2*mean(radii)
        
        ( mean((2*radii).^6)/ mean((2*radii).^2))^(1/4)
        
        f1 = sum(pi*radii.^2)/((2*max_x)^2)
        
        
        [tab, cell_idx] = generate_lookup_table(centres, radii, max_x,max_y, vs);
        
        save_substrate_file(centres, radii, tab, cell_idx, f1, vs, max_x, max_y, max_z, fn, name);
    end
end

% R = move_particles_in_substrate(centres, radii, max_x, max_y, max_z);



%% lookup table generator
function [tab, cell_idx] = generate_lookup_table(centres, radii, max_x,max_y, vs)
xs = -(max_x/vs):(max_x/vs-1); %x coordinates of voxels, corresponding to lower left corner
ys = -(max_y/vs):(max_y/vs-1);


xs = round(xs);
ys = round(ys);

[X,Y] = meshgrid(xs,ys); %lookup table
X = X(:); %voxel positions, the lookup table
Y = Y(:);


%
%     for c = 1:numel(X)
%         rectangle('Position', [X(c)*vs, Y(c)*vs, vs vs], 'EdgeColor', 'r');
%     end


% return
%find intersections
disp("Finding intersections...\n")
tic
n_int = zeros(numel(X), 1);
cell_idx = zeros(numel(X), 1);


parfor c = 1:numel(X)
    tmp_x = X(c)*vs;
    tmp_y = Y(c)*vs;
    ds1 = (centres(:,1)-tmp_x).^2 + (centres(:,2) - tmp_y).^2;
    in1 = (ds1 <= radii.^2);
    ds2 = (centres(:,1)-(tmp_x+vs)).^2 + (centres(:,2) - (tmp_y)).^2;
    in2 = (ds2 <= radii.^2);
    ds3 = (centres(:,1)-(tmp_x+vs)).^2 + (centres(:,2) - (tmp_y+vs)).^2;
    in3 = (ds3 <= radii.^2);
    ds4 = (centres(:,1)-(tmp_x)).^2 + (centres(:,2) - (tmp_y+vs)).^2;
    in4 = (ds4 <= radii.^2);
    
    id1 = find(in1);
    id2 = find(in2);
    id3 = find(in3);
    id4 = find(in4);
    
    idx = [id1;id2;id3;id4];
    
    n_int(c) = numel(unique(idx));
    
    if n_int(c) > 1; error('Multiple cells per voxel'); end
    
    if ~isempty(idx); cell_idx(c) = max(idx); else; cell_idx(c) = -1; end
    %         if ~isempty(id1); cell_idx(c) = id1; else; cell_idx(c) = -1; end
    
    %          disp(100*c/numel(X))
end
toc

disp("Generating lookup table...\n")
tic

%make a voxel lookup table
tab = zeros(size(X));
for c = 1:numel(X)
    a = X(c);
    b = Y(c);
    p = szudzik(a,b);
    tab(c) = p;
end

%We need to sort the table for binary search to work for lookup
[tab, sort_idx] = sort(tab, 'ascend');
cell_idx = cell_idx(sort_idx);

toc
% % %     %test this mechanism
%     Nsteps = 1000;
% 
%     r0 = [0 0];
% 
%     for i = 1:Nsteps
%         step = randn(1,2);
%         r0 = r0 + vs*step/norm(step);
%         a = floor(r0(1)/vs);
%         b = floor(r0(2)/vs);
%         p = szudzik(a,b);
%         c_idx = cell_idx(tab == p);
%         if c_idx > 0
%         xc = centres(c_idx,1);
%         yc = centres(c_idx,2);
%         r = radii(c_idx);
%         xx = linspace(xc-r, xc+r);
%         o_plus = sqrt( r^2 - (xx - xc).^2) + yc;
%         o_minus = -sqrt( r^2 - (xx - xc).^2) + yc;
%         fill(xx, real(o_plus), 'g');
%         fill(xx, real(o_minus), 'g');
%         plot(r0(1), r0(2), 'yo', 'MarkerFaceColor', 'y')
%         drawnow
%         end
%     end

end



%% save substrate file
function save_substrate_file(centres, radii, tab, cell_idx, f1, vox_size, max_x, max_y, max_z, fn, name)
N = size(centres, 1);
centres_x = centres(:,1);
centres_y = centres(:,2);
N_vox = numel(tab);

%convert to correct type
N = int64(N);
centres_x = single(centres_x);
centres_y = single(centres_y);
N_vox = int64(N_vox);
max_x = single(max_x);
max_y = single(max_y);
max_z = single(max_z);
vox_size = single(vox_size);
f1 = single(f1);
radii = single(radii);
tab = int64(tab);
cell_idx = int64(cell_idx);




fileID = fopen(fn, 'w');
fwrite(fileID, N, 'int64');
fwrite(fileID, N_vox, 'int64');
fwrite(fileID, max_x, 'single');
fwrite(fileID, max_y, 'single');
fwrite(fileID, max_z, 'single');
fwrite(fileID, vox_size, 'single');
fwrite(fileID, f1, 'single');
fwrite(fileID, centres_x, 'single');
fwrite(fileID, centres_y, 'single');
fwrite(fileID, radii, 'single');
fwrite(fileID, tab, 'int64');
fwrite(fileID, cell_idx, 'int64');

fclose(fileID);

%check that the writing went well
fileID = fopen(fn, 'r');
num_cells = fread(fileID, [1,1], 'int64');
tmp_N_vox = fread(fileID, [1,1], 'int64');
tmp_max_x = fread(fileID, [1,1], 'single');
tmp_max_y = fread(fileID, [1,1], 'single');
tmp_max_z = fread(fileID, [1,1], 'single');
tmp_vox_size = fread(fileID, [1,1], 'single');
tmp_f1 = fread(fileID, [1, 1], 'single');
tmp_centres_x = fread(fileID, [num_cells,1], 'single');
tmp_centres_y = fread(fileID, [num_cells,1], 'single');
tmp_radii = fread(fileID, [num_cells,1], 'single');
tmp_tab = fread(fileID, [N_vox,1], 'int64');
tmp_cell_idx = fread(fileID, [N_vox,1], 'int64');
fclose(fileID);

tests = [
    num_cells == N;
    tmp_max_x == max_x;
    tmp_max_y == max_y;
    tmp_max_z == max_z;
    tmp_vox_size == vox_size;
    tmp_f1 == f1;
    isequal(tmp_centres_x, centres_x);
    isequal(tmp_centres_y, centres_y);
    isequal(tmp_radii,  radii);
    tmp_N_vox == N_vox;
    isequal(tmp_tab, tab);
    isequal(tmp_cell_idx, cell_idx);
    ];
if ~all(tests)
    disp("FAILED.")
else
    %save information about this substrate to a mat file to be used as
    %ground truth in fitting
    info.fn = fn;
    info.num_cells = numel(radii);
    info.diameters = 2*radii;
    info.voxel_size = vox_size;
    info.surface_to_volume = 4./info.diameters;
    info.mean_d = mean(info.diameters);
    info.mr_mean_d = (mean(info.diameters.^6)/mean(info.diameters.^2))^(1/4);
    info.mr_surface_to_volume = 4/info.mr_mean_d;
    info.mean_ks = [0 1 2 3 4 5 6 8 10 12 14 16 18 20 25 30 35 40 100 200];
    info.kappas = info.mean_ks/info.mr_surface_to_volume; %this gives the correct permeability
    
    
    save(name+"_info.mat", "info", '-v7.3')
    
    
    disp(name + " =>> " + " SUCCESS.")
end

end



function R = move_particles_in_substrate(centres, radii, max_x, max_y, max_z)
%test substrate
x_length = 2*max_x;
y_length = 2*max_y;
z_length = 2*max_z;

Nsteps = 100000;
dr = sqrt(2*2e-9*1e-6);

R = zeros(2,Nsteps,3);
R(1,1,1:2) = centres(1,:);
R(2,1,1:2) = centres(1,:)+ 1.05*radii(1);

states = zeros(2,Nsteps);
states(1,1) = 1;

h = zeros(2,3);


for part = 1:2
    tmp_R = squeeze(R(part, 1,:));
    tmp_R = tmp_R';
    for step = 2:Nsteps
        
        dR = -dr + 2*rand(1,3)*dr;
        tmp_R = tmp_R + dR;
        
        %check in-out cells
        ds = (centres(:,1)-tmp_R(1)).^2 + (centres(:,2)-tmp_R(2)).^2;
        in = any(ds <= radii.^2);
        if in && (states(part, step-1) == 0); tmp_R = tmp_R-dR; states(part, step) = 0; end
        if in && (states(part, step-1) == 1); states(part, step) = 1; end
        if ~in && (states(part, step-1) == 1); tmp_R = tmp_R-dR; states(part, step) = 1;  end
        if ~in && (states(part, step-1) == 0); states(part, step) = 0;  end
        
        %check in-out substrate
        
        if tmp_R(1) < -max_x %left in negative x
            tmp_R(1) = tmp_R(1) + x_length;
            h(part,1) = h(part,1) - 1;
        end
        
        if tmp_R(1) >= max_x %left in +ve x
            tmp_R(1) = tmp_R(1) - x_length;
            h(part,1) = h(part,1) + 1;
        end
        
        if tmp_R(2) < -max_y %left in negative y
            tmp_R(2) = tmp_R(2) + y_length;
            h(part,2) = h(part,2) - 1;
        end
        
        if tmp_R(2) >= max_y %left in +ve y
            tmp_R(2) = tmp_R(2) - y_length;
            h(part,2) = h(part,2) + 1;
        end
        
        if tmp_R(3) < -max_z %left in negative z
            tmp_R(3) = tmp_R(3) + z_length;
            h(part,3) = h(part,3) - 1;
        end
        
        if tmp_R(3) >= max_z %left in +ve z
            tmp_R(3) = tmp_R(3) - z_length;
            h(part,3) = h(part,3) + 1;
        end
        
        R(part, step, :) = tmp_R + h(part, :).*[x_length, y_length, z_length];
        
    end
end

figure;
plot3(squeeze(R(1,:,1)), squeeze(R(1,:,2)), squeeze(R(1,:,3)), '.-')
hold on
plot3(squeeze(R(2,:,1)), squeeze(R(2,:,2)), squeeze(R(2,:,3)), '.-')
movegui('center')
end

function p = szudzik(a,b)
if a >= 0; a = 2*a; else; a = -2*a-1; end
if b >= 0; b = 2*b; else; b = -2*b-1; end

if a>= b; p = a^2+a+b; else; p = b^2+a; end
end


function f1 = compute_f1(centres, radii, max_x)
%compute correct f1 taking into account half and quarter cells in the
%substrate
area_total = (2*max_x)^2;
area_in = 0;
for c = 1:numel(radii)
    c_x = centres(c,1);
    c_y = centres(c,2);
    r = radii(c);
    whole = ((c_x-r > -max_x) & (c_x+r < max_x) & (c_y-r > -max_x) & (c_y+r < max_x));
    half_x_neg = ((c_x-r <= -max_x) & (c_x+r < max_x) & (c_y-r > -max_x) & (c_y+r < max_x));
    half_x_pos = ((c_x-r > -max_x) & (c_x+r >= max_x) & (c_y-r > -max_x) & (c_y+r < max_x));
    half_y_neg = ((c_x-r > -max_x) & (c_x+r < max_x) & (c_y-r <= -max_x) & (c_y+r < max_x));
    half_y_pos = ((c_x-r > -max_x) & (c_x+r < max_x) & (c_y-r > -max_x) & (c_y+r >= max_x));
    half = half_x_neg | half_x_pos | half_y_neg | half_y_pos;
    
    quarter_top_left = ((c_x-r <= -max_x) & (c_x+r < max_x) & (c_y-r > -max_x) & (c_y+r >= max_x));
    quarter_top_right = ((c_x-r > -max_x) & (c_x+r >= max_x) & (c_y-r > -max_x) & (c_y+r >= max_x));
    quarter_bottom_left = ((c_x-r <= -max_x) & (c_x+r < max_x) & (c_y-r <= -max_x) & (c_y+r < max_x));
    quarter_bottom_right = ((c_x-r > -max_x) & (c_x+r >= max_x) & (c_y-r <= -max_x) & (c_y+r < max_x));
    
    quarter = quarter_top_left | quarter_top_right | quarter_bottom_left | quarter_bottom_right;
    
    if whole; area_in = area_in + (pi*r^2);  end
    if half; area_in = area_in + 0.5*(pi*r^2);  end
    if quarter; area_in = area_in + 0.25*(pi*r^2); end
    
end

f1 = area_in/area_total;
end