
close all
run_sim = 1;
if run_sim
    %run simulation with only intracellular spins and no exchange as unit test
    main_folder = fileparts(pwd); %folder containing this matlab script
    cuda_folder = [main_folder, '\cuda\']; %folder containing cuda source and executable
    opt_path = [cuda_folder '\options.txt']; %path to options file
    exe_fn = 'simplex_ct.exe'; %executable name
    %check if options file exists
    if ~isfile(opt_path); error('Options file does not exist. Abort'); end
    %create options struct
    opt.NUMBER_OF_PARTICLES = 100000;
    opt.TOTAL_SIMULATION_TIME = 105e-3; %s
    opt.SIMULATION_TIME_STEP = 1e-5; %s
    opt.SAMPLING_TIME_INTERVAL = 1e-4; %s
    opt.INTRA_DIFFUSIVITY = 1e-9; %m^2/s
    opt.EXTRA_DIFFUSIVITY = 1e-9; %m^2/s
    opt.MEMBRANE_PERMEABILITY = 0; %m/s
    opt.INITIALISE_ALL_SPINS_INTRA = 1;  %1 = yes, 0 = no
    opt.INITIALISE_ALL_SPINS_EXTRA = 0; %1 = yes, 0 = no
    opt.ALLOW_INTRA_TO_EXTRA_TRANSITIONS = 0; %1 = yes, 0 = no
    opt.ALLOW_EXTRA_TO_INTRA_TRANSITIONS = 0; %1 = yes, 0 = no
    opt.SUBSTRATE_FILE_NAME = '..\substrates\test_cylinders.bin';
    opt.TRAJECTORY_FILE_NAME = '..\output\test_traj_ct.bin';
    opt.SAVE_STATE_HISTORY = 1; %1 = yes, 0 = no
    opt.STATE_FILE_NAME = '..\output\test_stat_ct.bin';
    
    set_options(opt_path, opt); %update options file
    STATUS = system(['cd ' cuda_folder ' & ' exe_fn], '-echo');
end

%compute cumulants and velocity autocorrelation functions
gwf_fn =  "fwf_resex.mat";
load(gwf_fn, "gwf", "dt")
G_ind = [4, 9, 5, 7, 1, 11]; %indices of increasing Gamma and Vomega
V_ind = [10, 2, 12, 3, 6,  8];
ind = [G_ind, V_ind];
gwf = gwf(ind,:); %sort gwf in increasing Gamma then increasing Vomega
dt_new = opt.SAMPLING_TIME_INTERVAL; %NOTE: using a coarse sampling resolution such as 1e-3 causes a visible bias in the predicted signals
time = (0:(size(gwf,2)*round(dt/dt_new)-1))*dt_new;
gwf_fwf = zeros(12, numel(time));
for c = 1:12
    [g, rf] = gwf_subsample_1d(gwf(c,:), dt, dt_new);
    gwf_fwf(c,1:numel(g)) = g.*rf;
end
dt = dt_new;
[gamma, vomega] = resex_mc_protocol_to_gamma_vomega(gwf_fwf, dt);

%also do sde
[gwf_sde, td] = make_sde_protocol(99.9e-3, dt);
N = size(gwf_sde,2);
gwf = zeros(24, N);
gwf(1:12, :) = gwf_sde;
for c = 13:24
    gwf(c,1:numel(gwf_fwf(c-12,:))) = gwf_fwf(c-12,:);
end
gamma = [gamma*NaN; gamma];
vomega = [vomega*NaN; vomega];
td = [td'; td'*NaN];

%make 3d
gwf_sde_fwf = gwf;
gwf = zeros(24, size(gwf_sde_fwf, 2), 3);
for c = 1:24
    gwf(c,:,1) = gwf_sde_fwf(c,:);
    gwf(c,:,2) = gwf_sde_fwf(c,:)*0;
    gwf(c,:,3) = gwf_sde_fwf(c,:)*0;
end


bvalues = [0 1 2 3 4]*1e9;
prog = 0;


r_fn = opt.TRAJECTORY_FILE_NAME;
[Npart, T, Nt, X,Y,Z] =  load_traj_from_bin(r_fn);

[chivXX, chivYY, chivZZ, chivXY, chivXZ, chivYZ, ~] = computeChivFromParticleTrace(X, Y, Z, dt, 0);
figure;
% plot(chivXX); hold on
% plot(chivYY)
% plot(chivZZ)
len = numel(chivXX)*4;
Dw = (fftshift(fft(chivXX*dt, len)));
Nw = numel(Dw);
df = 1/(Nw*dt);
if mod(Nw, 2) == 0; offset = 1; else; offset = 0; end
w = 2*pi*df*(-floor(Nw/2):floor(Nw/2)-offset); %omega

d = 5e-6;
D0 = opt.INTRA_DIFFUSIVITY;
num_k = 200;
Dw_theo = resex_mc_compute_dw(w, d, D0, num_k);

figure('Color', 'w');
plot(w, real(Dw), 'g-')
hold on
plot(w, Dw_theo, 'k-')



Ngwf = size(gwf, 1);
Nb = numel(bvalues);
c2 = zeros(Ngwf, Nb);
c4 = zeros(Ngwf, Nb);
S = zeros(Ngwf, Nb); %signal

c2_theo = zeros(Ngwf, Nb);


figure('Color', 'w');
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
ax6 = subplot(2,3,6);

sde = 0;
fwf = 0;
for c_gwf = 1:size(gwf, 1)
    for c_bv = 1:numel(bvalues)
        tmp_gwf_x = squeeze(gwf(c_gwf, :,1));
        tmp_gwf_y = squeeze(gwf(c_gwf, :,2));
        tmp_gwf_z = squeeze(gwf(c_gwf, :,3));
        tmp_q_x = msf_const_gamma*cumsum(tmp_gwf_x)*dt;
        tmp_q_y = msf_const_gamma*cumsum(tmp_gwf_y)*dt;
        tmp_q_z = msf_const_gamma*cumsum(tmp_gwf_z)*dt;
        tmp_b = sum(tmp_q_x.^2 + tmp_q_y.^2 + tmp_q_z.^2)*dt;
        if round(tmp_b)==0
            tmp_gwf_X = tmp_gwf_x*0;
            tmp_gwf_Y = tmp_gwf_y*0;
            tmp_gwf_Z = tmp_gwf_z*0;
            
            tmp_q = tmp_q_x*0;
        else
            tmp_gwf_X = tmp_gwf_x*sqrt(bvalues(c_bv)/tmp_b);
            tmp_gwf_Y = tmp_gwf_y*sqrt(bvalues(c_bv)/tmp_b);
            tmp_gwf_Z = tmp_gwf_z*sqrt(bvalues(c_bv)/tmp_b);
            
            tmp_q = tmp_q_x*sqrt(bvalues(c_bv)/tmp_b);
        end
        
        [Qw, w] = resex_mc_qt_to_qw(tmp_q, dt);
        Dw = resex_mc_compute_dw(w, d, D0, num_k);
        dw = w(2)-w(1);
        
        c2_theo(c_gwf, c_bv) = (-1/(2*pi))*trapz(w, Qw.*Dw);
        
        tmp_gwf = [tmp_gwf_X', tmp_gwf_Y', tmp_gwf_Z'];
        
        phase = computePhase(X, Y, Z, tmp_gwf, dt);
        
        c2(c_gwf, c_bv) = -(1/2)*mean(phase.^2);
        
        c4(c_gwf, c_bv) = (1/24)*(mean(phase.^4) - 3*mean(phase.^2)^2);
        
        S(c_gwf, c_bv) = (1/numel(phase))*sum(cos(phase));
        
        prog = prog+1;
        disp("Done " + num2str(prog) + " of " + num2str(size(gwf, 1)*numel(bvalues)))
        
    end
    
    if c_gwf <= 12 %sde
        sde = sde+1;
        col = [1 0 0]*(sde-1)/(12-1);
        plt = plot(ax1, bvalues, c2(c_gwf, :), 'o', 'Color', col);
        hold (ax1, 'on')
        plot(ax1, bvalues, c2_theo(c_gwf, :), '--', 'Color', plt.Color);
        plot(ax2, bvalues, c4(c_gwf, :), 'o-',  'Color', col);
        hold (ax2, 'on')
        plot(ax3, bvalues, S(c_gwf, :), 'o-',  'Color', col);
        hold (ax3, 'on')
    else
        fwf = fwf+1;
        col = [0 0 1]*(fwf-1)/(12-1);
        plt = plot(ax4, bvalues, c2(c_gwf, :), 'o', 'Color', col);
        hold (ax4, 'on')
        plot(ax4, bvalues, c2_theo(c_gwf, :), '--', 'Color', plt.Color);
        plot(ax5, bvalues, c4(c_gwf, :), 'o-',  'Color', col);
        hold (ax5, 'on')
        plot(ax6, bvalues, S(c_gwf, :), 'o-',  'Color', col);
        hold (ax6, 'on')
    end
    drawnow
end


title(ax1, 'SDE c_2')
title(ax2, 'SDE c_4')
title(ax3, 'SDE Signal')
title(ax4, 'FWF c_2')
title(ax5, 'FWF c_4')
title(ax6, 'FWF Signal')

max_dif = 100*max(abs(c2(:)-c2_theo(:))./c2_theo(:));
disp("Maximum error = " + num2str(abs(max_dif)) + " %");
if abs(max_dif) < 1; disp("Validation result: SUCCESS")
else; disp("Validation result: FAILED"); end


%compare with theory





%%
function [gwf, td] = make_sde_protocol(T, dt)
delta = 5e-3;
Gmax = 500e-3;
Smax = 500;
min_td = (2/3)*delta + 3*dt + 1e-3 + 5.2e-3;
max_td = T - delta - delta/3;
Ntd = 12;
td = linspace(min_td, max_td, Ntd);

opt.T = T; %total encoding time in seconds
opt.dt = dt; %time-step in seconds
opt.g_max = Gmax; %maximum gradient amplitude in Tesla
opt.slew_max = Smax; %maximum slew rate in Tesla/seconds
opt.t_start = 0e-3; %start-time of first gradient pulse in seconds
opt.sde.delta = delta; %gradient pulse duration in seconds
opt.time = linspace(0, opt.T,round(opt.T/opt.dt)+1); %time vector = linspace(0, T, N_steps) in seconds
opt.N_steps = length(opt.time); %number of time-steps = round(T/dt)+1
opt.gam_ma = 2*pi*42.6e6; %gamma constant for the 1H nucleus, 42.6e6 Hz/T

gwf = zeros(numel(td), (round(T/dt)+1), 1);

for c = 1:numel(td)
    tmp_td = td(c);
    [g, ~] = get_sde_gwf(tmp_td, dt, opt);
    %     gwf_to_save(c, 1:numel(g)) = g;
    %     [g, rf] = gwf_subsample_1d(g, dt, dt);
    %     g = g.*rf;
    gwf(c, 1:numel(g)) = g;
end
end


function [gwf, dt] = get_sde_gwf(td, dt, opt)
gwf = zeros(numel(td), opt.N_steps);
for c_td = 1:numel(td)
    DELTA = td(c_td) + opt.sde.delta/3;
    tmp_interval = round((DELTA-opt.sde.delta)*1e5)/1e5-dt;
    opt.sde.interval = tmp_interval; %time interval between end of pulse 1 and start of pulse 2 in seconds (interval + delta = DELTA)
    waveform_info = waveform_sde(opt);
    gwf(c_td, 1:numel(waveform_info.waveform)) = waveform_info.waveform';
end

dt = opt.dt;
end

function set_options(opt_fn, opt)
%changes the settings in the options file at opt_fn using the values stored
%in opt
f = fopen(opt_fn, 'w');
fprintf(f, '%s\n', 'NUMBER_OF_PARTICLES');
fprintf(f, '%s\n',  num2str(opt.NUMBER_OF_PARTICLES));
fprintf(f, '%s\n', 'TOTAL_SIMULATION_TIME_[s]');
fprintf(f, '%s\n',  num2str(opt.TOTAL_SIMULATION_TIME));
fprintf(f, '%s\n', 'SIMULATION_TIME-STEP_[s]');
fprintf(f, '%s\n', num2str(opt.SIMULATION_TIME_STEP));
fprintf(f, '%s\n', 'SAMPLING_TIME_INTERVAL_[s]');
fprintf(f, '%s\n', num2str(opt.SAMPLING_TIME_INTERVAL));
fprintf(f, '%s\n', 'INTRA_DIFFUSIVITY_[m^2/s]');
fprintf(f, '%s\n', num2str(opt.INTRA_DIFFUSIVITY));
fprintf(f, '%s\n', 'EXTRA_DIFFUSIVITY_[m^2/s]');
fprintf(f, '%s\n', num2str(opt.EXTRA_DIFFUSIVITY));
fprintf(f, '%s\n', 'MEMBRANE_PERMEABILITY[m/s]');
fprintf(f, '%s\n', num2str(opt.MEMBRANE_PERMEABILITY));
fprintf(f, '%s\n', 'INITIALISE_ALL_SPINS_INTRA?_[1=yes,0=no]');
fprintf(f, '%s\n', num2str(opt.INITIALISE_ALL_SPINS_INTRA));
fprintf(f, '%s\n', 'INITIALISE_ALL_SPINS_EXTRA?_[1=yes,0=no]');
fprintf(f, '%s\n', num2str(opt.INITIALISE_ALL_SPINS_EXTRA));
fprintf(f, '%s\n', 'ALLOW_INTRA_TO_EXTRA_TRANSITIONS?_[1=yes,0=no]');
fprintf(f, '%s\n', num2str(opt.ALLOW_INTRA_TO_EXTRA_TRANSITIONS));
fprintf(f, '%s\n', 'ALLOW_EXTRA_TO_INTRA_TRANSITIONS?_[1=yes,0=no]');
fprintf(f, '%s\n', num2str(opt.ALLOW_EXTRA_TO_INTRA_TRANSITIONS));
fprintf(f, '%s\n', 'SUBSTRATE_FILE_NAME');
fprintf(f, '%s\n', num2str(opt.SUBSTRATE_FILE_NAME));
fprintf(f, '%s\n', 'TRAJECTORY_FILE_NAME');
fprintf(f, '%s\n', num2str(opt.TRAJECTORY_FILE_NAME));
fprintf(f, '%s\n', 'SAVE_STATE_HISTORY?_[1=yes,0=no]');
fprintf(f, '%s\n', num2str(opt.SAVE_STATE_HISTORY));
fprintf(f, '%s\n', 'STATE_FILE_NAME');
fprintf(f, '%s\n', num2str(opt.STATE_FILE_NAME));
fclose(f);
end
