
%MACRO for running simulations with the FLEXD framework (voxelised version)

recompile = 1; %will recompile code if true

main_folder = fileparts(pwd); %folder containing this matlab script
cuda_folder = [main_folder, '\cuda\']; %folder containing cuda source and executable
opt_path = [cuda_folder '\options.txt']; %path to options file
exe_fn = 'simplex_ct.exe'; %executable name
exe_path = [cuda_folder, exe_fn]; %path to executable
cu_path = [cuda_folder, '\simplex_ct.cu']; %path to cu source file

if recompile || ~isfile(exe_path) %try to recompile if asked or if executable does not exist
    compile_cmd = ['nvcc ', cu_path, ' -o ', exe_path]; %compilation command
    STATUS = system(compile_cmd);
    if STATUS ~= 0; error ('Compilation failed. Abort.'); end
end

%check if options file exists
if ~isfile(opt_path); error('Options file does not exist. Abort'); end

%create options struct
opt.NUMBER_OF_PARTICLES = 100;
opt.TOTAL_SIMULATION_TIME = 100e-3; %s
opt.SIMULATION_TIME_STEP = 1e-6; %s
opt.SAMPLING_TIME_INTERVAL = 1e-3; %s
opt.INTRA_DIFFUSIVITY = 1e-9; %m^2/s
opt.EXTRA_DIFFUSIVITY = 1e-9; %m^2/s
opt.MEMBRANE_PERMEABILITY = 0; %m/s
opt.INITIALISE_ALL_SPINS_INTRA = 1;  %1 = yes, 0 = no
opt.INITIALISE_ALL_SPINS_EXTRA = 1; %1 = yes, 0 = no
opt.ALLOW_INTRA_TO_EXTRA_TRANSITIONS = 1; %1 = yes, 0 = no
opt.ALLOW_EXTRA_TO_INTRA_TRANSITIONS = 0; %1 = yes, 0 = no
opt.SUBSTRATE_FILE_NAME = '..\substrates\test_mesh.bin';
opt.TRAJECTORY_FILE_NAME = '..\output\traj.bin';
opt.SAVE_STATE_HISTORY = 1; %1 = yes, 0 = no
opt.STATE_FILE_NAME = '..\output\stat.bin';


%now iteratively change these options and run the simulator
%the example below will vary the membrane permeability (denoted kappa)

kappa_vec = [0 0.03 0.1]*1e-3;

for c_kappa = 1%
    opt.MEMBRANE_PERMEABILITY = kappa_vec(c_kappa); %change permeability
    opt.TRAJECTORY_FILE_NAME = ['..\output\traj', num2str(c_kappa) ,'.bin']; %change trajectory file name
    opt.STATE_FILE_NAME = ['..\output\stat', num2str(c_kappa) ,'.bin']; %change compartment/state/transition history file name
    
    set_options(opt_path, opt); %update options file

    STATUS = system(['cd ' cuda_folder ' & ' exe_fn], '-echo');
    %NOTE: set CUDA_VISIBLE_DEVICES=1 selects GPU 2 on multi-GPU computers
    %so that simulations can be run on multiple GPUs in parallel
    if STATUS~=0; warning('Failed to launch FLEXD simulator.'); end
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
fprintf(f, '%s\n', 'SAMPLE_TRAJECTORIES?_[1=yes,0=no]');
fprintf(f, '%s\n', num2str(opt.SAMPLE_TRAJECTORIES));
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

