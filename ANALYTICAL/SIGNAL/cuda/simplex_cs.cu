#include <iostream>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <random>
//This program simulates restricted diffusion with exchange in a substrate of 
//analytically defined parallel cylinders
//It computes signals at a resolution saved by the user.

//Author: 
//Arthur Chakwizira, Medical Radiation Physics, Lund University, Sweden
//arthur.chakwizira@med.lu.se


//simulation options; will be read from an options txt file
struct options
{
	long Npart;
	float T;
	float sim_dt;
	bool do_samp;
	float samp_dt;
	int n_dim;
	float D0_intra, D0_extra;
	long sim_Nt;
	long save_Nt;
	float ds_intra, ds_extra;
	long N_save; //N time points x N particles
	long N_sim;
	bool save_states; //save particle state history to file or not
	float kappa, p_12, p_21; //membrane permeability
	//option to initialise all spins in one compartment and allow one-directional transitions
	bool all_intra, all_extra, intra_to_extra_allowed, extra_to_intra_allowed;
	//additional gwf options
	long n_acq = 0; //# waveforms
	long n_gwf_points = 0; // # time points in each waveform
	float gamma; //gyromagnetic ratio
	long delay; //number of time steps to take before acquiring signals
	bool save_positions; //save final positions or not; useful for debugging
	float gwf_dt; //gradient waveform resolution
	bool save_phase; //save phase or not, might be useful for computing cumulants
};


//substrate info
struct world
{
	long long num_cells, num_voxels;
	float max_x, max_y, max_z, x_length, y_length, z_length, f1, vox_size;
};

//initialise rng with seeds based on system clock
__global__ void random_init(curandState* states)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(clock64()+index, index, 0, &states[index]);
}

//function for pairing particle coordinates using the Szudzik pairing algorithm 
__device__ void pair(long long x, long long y, long long& xy)
{
	//device function for pairing particle coordinates
	x >= 0 ? x = 2 * x : x = -2 * x - 1;
	y >= 0 ? y = 2 * y : y = -2 * y - 1;

	x >= y ? xy = x * x + x + y : xy = y * y + x;
}

//binary search to find index of element in array
__device__ long long binary_search_iter(long long* A, long long lower, long long upper, long long x)
{
	while (upper >= lower) {
		long long mid = lower + (upper - lower) / 2;
		if (A[mid] == x) return mid;
		(A[mid] > x) ? upper = mid - 1 : lower = mid + 1;
	}
	return -1;
}

//determine if particle is inside any cylinder in the substrate
__device__ bool is_particle_in_any_cell(float tmp_x, float tmp_y, world* w, float* centre_x, float* centre_y, float* radii, long long* table, long long* cell_idx)
{
	long long x_pos = floor(tmp_x / w->vox_size);
	long long y_pos = floor(tmp_y / w->vox_size);
	long long xy; //voxel identifier
	long long which_cell, which_voxel;
	bool inside = false; //zero means outside all cells
	float distance, r, cx, cy;
	
	pair(x_pos, y_pos, xy); //get the identifier

	which_voxel = binary_search_iter(table, 0, w->num_voxels - 1, xy); //iterative binary search	
	
	which_cell = cell_idx[which_voxel]-1; //index of cell containing the voxel containing the particle
	//note the minus 1 to take into account that MATLAB numbering starts at 1. Need to fix this later
    
    r = radii[which_cell];
    cx = centre_x[which_cell];
    cy = centre_y[which_cell];
	
	if (which_cell >= 0) //-1 means voxel is not in any cell
	{
		distance = (tmp_x-cx)*(tmp_x-cx) + (tmp_y-cy)*(tmp_y-cy);
		if ( distance <= (r*r)  ) //means particle is in this cell
		{
			inside = true;
		}
	}
	return inside;
}


//make diffusion step; allows for different intra/extra diffusivities
__device__ void diffuse(float& tmp_x, float& tmp_y, float& tmp_z, float& tmp_dx,
	float& tmp_dy, float& tmp_dz, curandState* states, int index, options* opt, int tmp_loc)
{
	tmp_dx = curand_normal(&states[index]);
	tmp_dy = curand_normal(&states[index]);
	tmp_dz = curand_normal(&states[index]);

	float norm;
	if (tmp_loc == 1) //particle is in mesh
	{
		norm = (*opt).ds_intra * rnorm3df(tmp_dx, tmp_dy, tmp_dz);
	}
	else  //particle is not in mesh
	{
		norm = (*opt).ds_extra * rnorm3df(tmp_dx, tmp_dy, tmp_dz);
	}
	
	tmp_dx *= norm;
	tmp_dy *= norm;
	tmp_dz *= norm;

	tmp_x += tmp_dx;
	tmp_y += tmp_dy;
	tmp_z += tmp_dz;
}


//implementation of periodic boundary conditions
__device__ void restrict_to_world(float& e_x, float& e_y, float& e_z, options* opt, world* w, float& tmp_x, float& tmp_y, float& tmp_z)
{
	if (tmp_x < -w->max_x) { tmp_x += w->x_length; e_x -= w->x_length; }
	if (tmp_x >= w->max_x) { tmp_x -= w->x_length; e_x += w->x_length; }

	if (tmp_y < -w->max_y) { tmp_y += w->y_length; e_y -= w->y_length; }
	if (tmp_y >= w->max_y) { tmp_y -= w->y_length; e_y += w->y_length; }

	if (tmp_z < -w->max_z) { tmp_z += w->z_length; e_z -= w->z_length; }
	if (tmp_z >= w->max_z) { tmp_z -= w->z_length; e_z += w->z_length; }
}

//membrane intersection checks and exchange mechanism
__device__ void check_state(world* w, options* opt, float& tmp_x, float& tmp_y, float& tmp_z, float& tmp_dx,
	float& tmp_dy, float& tmp_dz, int& tmp_loc, float* centre_x, float* centre_y, float *radii, long long* table, long long * cell_idx,  curandState* states, int index)
{
	bool reject = false; //reject move or not
	bool inside;

	inside = is_particle_in_any_cell(tmp_x, tmp_y, w, centre_x, centre_y, radii, table, cell_idx); //determine whether particle is in any cell

	if (inside) //"now intra"
	{
		if (tmp_loc == 1) { reject = false; } //was intra before
		else //was not intra before
		{
			if (curand_uniform(&states[index]) < opt->p_21) { reject = false; tmp_loc = 1; }
			else { reject = true; tmp_loc = 0; }
		}
	}

	if (!inside) //"now extra"
	{
		if (tmp_loc == 0) { reject = false; } //was extra before
		else //was not extra before
		{
			if (curand_uniform(&states[index]) < opt->p_12) { reject = false; tmp_loc = 0; }
			else { reject = true; tmp_loc = 1; }
		}
	}


	if (reject) { tmp_x -= tmp_dx; tmp_y -= tmp_dy; tmp_z -= tmp_dz; }

}


__global__ //this kernel runs the actual simulation
void engine(float* x, float* y, float* z, int* loc, float* centre_x, float* centre_y, float* radii, long long* table, long long* cell_idx,
	curandState* states, options* opt, world* w, float* phase, float * gwf_x, float *gwf_y, float* gwf_z)
{
	long index = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	//int stride = blockDim.x * gridDim.x;

	long phase_entry;
	long save;
	long save_c_t;
	float tmp_x, tmp_y, tmp_z;
	float tmp_dx, tmp_dy, tmp_dz;
	float e_x, e_y, e_z; //keep track of hyperposition
	int tmp_loc;
	long sampling_interval = (long) std::round(opt->samp_dt / opt->sim_dt);

	long c_p = index;
	if (c_p < opt->Npart) //ensure we keep within bounds
	{
		tmp_x = x[c_p];
		tmp_y = y[c_p];
		tmp_z = z[c_p];
		tmp_dx = 0;
		tmp_dy = 0;
		tmp_dz = 0;
		tmp_loc = loc[c_p];
		
		if (index == opt->Npart-1) {printf("Preparing simulation loop ... \n");}
		
		//adding a delay loop here
		e_x = 0; e_y = 0; e_z = 0;
		for (int c_t = 0; c_t < 50; c_t++)
		{
			diffuse(tmp_x, tmp_y, tmp_z, tmp_dx, //make diffusion steps
				tmp_dy, tmp_dz, states, index, opt, tmp_loc);
			restrict_to_world(e_x, e_y, e_z, opt, w, tmp_x, tmp_y, tmp_z);
			check_state(w, opt, tmp_x, tmp_y, tmp_z, tmp_dx,
				tmp_dy, tmp_dz, tmp_loc, centre_x, centre_y, radii, table, cell_idx, states, index);
			restrict_to_world(e_x, e_y, e_z, opt, w, tmp_x, tmp_y, tmp_z);
		}

		if (index == opt->Npart-1) {printf("Starting simulation loop ... \n");}
		
		e_x = 0; e_y = 0; e_z = 0;
		save_c_t = -1;
		save = 0;
		for (long c_t = 0; c_t < (*opt).sim_Nt; c_t++)
		{
			
			diffuse(tmp_x, tmp_y, tmp_z, tmp_dx, //make diffusion steps
				tmp_dy, tmp_dz, states, index, opt, tmp_loc);
			restrict_to_world(e_x, e_y, e_z, opt, w, tmp_x, tmp_y, tmp_z);
			check_state(w, opt, tmp_x, tmp_y, tmp_z, tmp_dx,
				tmp_dy, tmp_dz, tmp_loc, centre_x, centre_y, radii, table, cell_idx, states, index);
			restrict_to_world(e_x, e_y, e_z, opt, w, tmp_x, tmp_y, tmp_z);
		
			save++;

			if (save == sampling_interval)
			{
				save_c_t++;
				for (long c_a = 0; c_a < opt->n_acq; c_a++)
					{
						phase_entry = c_p * (*opt).n_acq + c_a;
						phase[phase_entry] += opt->gamma*((tmp_x+e_x)*gwf_x[c_a*(*opt).save_Nt+save_c_t] + (tmp_y+e_y)*gwf_y[c_a*(*opt).save_Nt+save_c_t] + (tmp_z+e_z)*gwf_z[c_a*(*opt).save_Nt+save_c_t])*opt->samp_dt;
					}
					save = 0;
				if (index==opt->Npart-1 && save_c_t%100==0){printf("Step %d of %ld \n", save_c_t, (*opt).save_Nt);}
			}
		}
		if (index == opt->Npart-1) {printf("Finishing simulation loop ... \n");}
	//save last position and state
	x[c_p] = tmp_x+e_x; y[c_p] = tmp_y+e_y; z[c_p] = tmp_z+e_z; loc[c_p] = tmp_loc;
	}
}


bool set_options(options* opt, char* pos_fn, char* sta_fn, char* sub_fn, char* gwf_fn,
char* sig_fn, char* phase_fn)
{
	//sets options loaded from file
	std::string opt_fn = ".\\options.txt";
	std::string dummy;

	std::ifstream tf;
	tf.open(opt_fn, std::ios::in);
	if (tf.is_open())
	{
	tf >> dummy >> opt->Npart;
	tf >> dummy >> opt->sim_dt;
	tf >> dummy >> opt->D0_intra;
	tf >> dummy >> opt->D0_extra;
	tf >> dummy >> opt->kappa;
	tf >> dummy >> opt->all_intra;
	tf >> dummy >> opt->all_extra;
	tf >> dummy >> opt->intra_to_extra_allowed;
	tf >> dummy >> opt->extra_to_intra_allowed;
	tf >> dummy >> sub_fn;
	tf >> dummy >> gwf_fn;
	tf >> dummy >> sig_fn;
	tf >> dummy >> opt->save_positions;
	tf >> dummy >> pos_fn;
	tf >> dummy >> opt->save_states;
	tf >> dummy >> sta_fn;
	tf >> dummy >> opt->save_phase;
	tf >> dummy >> phase_fn;
	tf.close();


	opt->n_dim = 3;
	opt->gamma = 2.675129e8;
	
	opt->ds_extra = (float)sqrt(2 * (*opt).n_dim * (*opt).D0_extra * (*opt).sim_dt); //diffusion step size extra
	opt->ds_intra = (float)sqrt(2 * (*opt).n_dim * (*opt).D0_intra * (*opt).sim_dt); //diffusion step size intra

	//Add transition probablities here //Lee-2021, Eq.3
	//double Cd = 2.0/3, k = opt->kappa, ds1 = (opt->ds_blood), D1 = opt->D0_blood;
	//double ds2 = (opt->ds_extra), D2 = opt->D0_extra;
	//opt->p_12 = (k*ds1*Cd/D1) / (1 + (k/2)*( (ds1/D1)  + (ds2/D2) )*Cd);
	//opt->p_21 = (k*ds2*Cd/D2) / (1 + (k/2)*( (ds2/D2)  + (ds1/D1) )*Cd);

	opt->p_21 = (float)opt->kappa * sqrt(8 * opt->sim_dt / (3 * opt->D0_extra)); 
	if (!opt->extra_to_intra_allowed) {opt->p_21 = 0;};
	opt->p_12 = (float)opt->kappa * sqrt(8 * opt->sim_dt / (3 * opt->D0_intra));
	if (!opt->intra_to_extra_allowed) {opt->p_12 = 0;};
	
	std::cout << "Loaded options from: " << opt_fn << std::endl;
	std::cout << "Substrate filename is : " << sub_fn << std::endl;
	std::cout << "Simulation time step =  " << opt->sim_dt << std::endl;
	std::cout << "Intra diffusivity = " << opt->D0_intra << std::endl;
	std::cout << "Membrane permeability = " << opt->kappa << std::endl;
	std::cout << "Transition probability intra->extra : " << opt->p_12 << std::endl;
	std::cout << "Transition probability extra->intra : " << opt->p_21 << std::endl;
	return true;
	}
	else
	{std::cout << "\n ERROR loading options from: " << opt_fn << " \n" << std::endl;
		return false;
	}
}

//saves trajectories
void save_trajectory(float* x, float* y, float* z, char* r_fn, options* opt)
{
	long save_Nt = 1;
	std::cout << "Saving final position to: " << r_fn << std::endl;
	FILE* tf;
	tf = fopen(r_fn, "wb");
	if (tf)
	{
	fwrite(&(opt->Npart), sizeof(long), 1, tf);
	fwrite(&(opt->T), sizeof(float), 1, tf);
	fwrite(&(save_Nt), sizeof(long), 1, tf);
	fwrite(x, sizeof(float), opt->Npart, tf); 
	fwrite(y, sizeof(float), opt->Npart, tf);
	fwrite(z, sizeof(float), opt->Npart, tf);
	fclose(tf);
	}
	else
	{
		std::cout << "\n ERROR: Failed to save final positions to: " << r_fn << std::endl;
	}
}

//saves history of particle identities/compartment identities
void save_state_history(int* s, char* s_fn, options* opt)
{
	std::cout << "Saving transition history to: " << s_fn << std::endl;
	FILE* tf;
	tf = fopen(s_fn, "wb");
	if (tf)
	{
	fwrite(&(opt->Npart), sizeof(long), 1, tf);
	fwrite(&(opt->T), sizeof(float), 1, tf);
	fwrite(&(opt->save_Nt), sizeof(long), 1, tf);
	fwrite(s, sizeof(int), opt->N_save, tf);
	fclose(tf);
	}
	else
	{
		std::cout << "\n ERROR: Failed to save transition history to: " << s_fn << std::endl;
	}
}

//save signal to file
void save_signal(float* h_signal, char* sig_fn, options* opt)
{
	
	std::cout << "Saving signals to: " << sig_fn << std::endl;
	FILE* tf;
	tf = fopen(sig_fn, "wb");
	if (tf)
	{
	fwrite(h_signal, sizeof(float), opt->n_acq, tf);
	fclose(tf);
	}
	else
	{
		std::cout << "\n ERROR: Failed to save signals to: " << sig_fn << std::endl;
	}
}

//open substrate file and get number of voxels
bool get_substrate_info(options* opt, world* w, char* g_fn)
{
	FILE* sf;
	sf = fopen(g_fn, "rb");
	if (sf)
	{
	fread(&(w->num_cells), sizeof(long long), 1, sf);
	fread(&(w->num_voxels), sizeof(long long), 1, sf);
	fread(&(w->max_x), sizeof(float), 1, sf);
	fread(&(w->max_y), sizeof(float), 1, sf);
	fread(&(w->max_z), sizeof(float), 1, sf);
	fread(&(w->vox_size), sizeof(float), 1, sf);
	fclose(sf);
	if (w->vox_size < std::max(opt->ds_intra, opt->ds_extra)){
		std::cout << "\n ERROR: Substrate voxel size of " << w->vox_size << 
		" is smaller than the simulation step-size of " << std::max(opt->ds_intra, opt->ds_extra) <<std::endl;
		return false;
	}
	else{
	return true;}
	}
	else 
	{std::cout << "\n ERROR: Specified substrate file: " << g_fn << "  does not exist! \n" <<std::endl;
		return false;}
}

//open gwf file and get info about the gwf
bool get_gwf_info(char* gwf_fn, options* opt)
{
	std::cout << "Getting gradient waveform info from: " << gwf_fn << std::endl;
	FILE* sf;
	sf = fopen(gwf_fn, "rb");
	if (sf)
	{
	fread(&(opt->n_acq), sizeof(long), 1, sf);
	fread(&(opt->n_gwf_points), sizeof(long), 1, sf);
	fread(&(opt->gwf_dt), sizeof(float), 1, sf);
	fclose(sf);
	
	//sampling time step is equal to gwf time-step
	opt->samp_dt = opt->gwf_dt;
	if (opt->samp_dt < opt->sim_dt){
		std::cout << "\n ERROR: Simulation time-step of " << opt->sim_dt << "is larger than the gradient waveform resolution of " << opt->samp_dt <<std::endl;
		return false;
	}
	else{
	//total simulation time is number of points in gwf*gwf_dt
	opt->sim_Nt = opt->n_gwf_points * (opt->samp_dt/opt->sim_dt) +1; //total simulation time steps
	//set this also in variable T
	opt->T = (opt->sim_Nt)*opt->sim_dt;
	
	opt->save_Nt = opt->n_gwf_points;//(long long)round(opt->T / opt->samp_dt);
	opt->N_save = (long)opt->Npart * opt->save_Nt; //N time points x N particles
	opt->N_sim = (long)opt->Npart * opt->sim_Nt;
	return true;
	}
	}
	else
	{
		std::cout << "\n ERROR: Specified gradient waveform file: " << gwf_fn << "  does not exist! \n" <<std::endl;
		return false;
	}
}

//load simulation substrate from file
bool load_substrate(float* h_centre_x, float* h_centre_y, float* h_radii, long long* h_table, long long *h_cell_idx, world* w, char* g_fn)
{
	FILE* sf;
	sf = fopen(g_fn, "rb");
	if (sf)
	{
	fread(&(w->num_cells), sizeof(long long), 1, sf);
	fread(&(w->num_voxels), sizeof(long long), 1, sf);
	fread(&(w->max_x), sizeof(float), 1, sf);
	fread(&(w->max_y), sizeof(float), 1, sf);
	fread(&(w->max_z), sizeof(float), 1, sf);
	fread(&(w->vox_size), sizeof(float), 1, sf);
	fread(&w->f1, sizeof(float), 1, sf);
	fread(h_centre_x, sizeof(float), w->num_cells, sf); 
	fread(h_centre_y, sizeof(float), w->num_cells, sf); 
	fread(h_radii, sizeof(float), w->num_cells, sf); 
	fread(h_table, sizeof(long long), w->num_voxels, sf); 
	fread(h_cell_idx, sizeof(long long), w->num_voxels, sf); 
	fclose(sf);

	w->x_length = 2 * w->max_x;
	w->y_length = 2 * w->max_y;
	w->z_length = 2 * w->max_z;

	std::cout << "Loaded substrate from: " << g_fn << std::endl;
	std::cout << "Intra volume fraction =  " <<  w->f1 << std::endl;
	std::cout << "Number of voxels =  " <<  w->num_voxels << std::endl;
	std::cout << "Voxel size =  " <<  w->vox_size << std::endl;
	return true;
	}
	{
	std::cout << "\n ERROR reading substrate from:  " << g_fn << " \n" << std::endl;		
	return false;
	}
}

//load gradient waveforms from file
bool load_gwf(char* gwf_fn, options* opt, float* h_gwf_x, float* h_gwf_y, float* h_gwf_z)
{
	std::cout << "Loading gradient waveform..."<< std::endl;
	//load substrate from file
	FILE* sf;
	sf = fopen(gwf_fn, "rb");
	if (sf)
	{
	fread(&(opt->n_acq), sizeof(long), 1, sf);
	fread(&(opt->n_gwf_points), sizeof(long), 1, sf);
	fread(&(opt->gwf_dt), sizeof(float), 1, sf);
	fread(h_gwf_x, sizeof(float), opt->n_acq*opt->n_gwf_points, sf);
	fread(h_gwf_y, sizeof(float), opt->n_acq*opt->n_gwf_points, sf);
	fread(h_gwf_z, sizeof(float), opt->n_acq*opt->n_gwf_points, sf); 
	fclose(sf);
	return true;
	}
	else
	{
		std::cout << "\n ERROR loading gradient waveforms from:  " << gwf_fn << " \n" << std::endl;
		return false;
	}

}

//convert phase to signal
void convert_phase_to_signal(float *h_phase, float *h_signal, options *opt)
{
	std::cout << "Converting phase to signal..." << std::endl;
	float sum_cos_phase;
	long phase_entry;
	for (long ca = 0; ca<opt->n_acq; ca++)
	{
		sum_cos_phase = 0;
	for (long c_p = 0; c_p < opt->Npart; c_p++)
		{
			phase_entry = c_p * (*opt).n_acq + ca;
			sum_cos_phase += cos(h_phase[phase_entry]);
		}
		h_signal[ca] = sum_cos_phase/opt->Npart;
	}
	 std::cout << "Done." << std::endl;
}

//save phase to file
void save_phase(float* h_phase, options *opt, char* sta_fn)
{
		
	long N = opt->Npart*opt->n_acq;
	std::cout << "Saving phase to: " << sta_fn << std::endl;
	FILE* tf;
	tf = fopen(sta_fn, "wb");
	if (tf)
	{
	fwrite(&N, sizeof(long), 1, tf);
	fwrite(h_phase, sizeof(float), opt->n_acq*opt->Npart, tf);
	fclose(tf);
	std::cout << "Done." << std::endl;
	}
	else
	{
		std::cout << "\n ERROR saving phase to: " << sta_fn << std::endl;
	}
}



//generate starting positions for all particles
__global__ void generate_initial_distribution(float* x, float* y, float* z, int* loc, float* centre_x,
	float* centre_y, float* radii, long long* table, long long* cell_idx, curandState* states, options* opt, world* w)
{
	long index = blockIdx.x * blockDim.x + threadIdx.x; //get thread idx
	float tmp_x, tmp_y, tmp_z; //frac defines intra-extra split of initial populations
	long idx;
	bool inside;
	bool success = false;
	float intra_frac; //fraction of spins to put intra
	if (opt->all_intra && !opt->all_extra)  {intra_frac = 1;};
	if (opt->all_extra && !opt->all_intra)  {intra_frac = 0;};
	if (opt->all_intra && opt->all_extra)  {intra_frac = 0.5;}; //if both are true, divide by half
	if (!opt->all_intra && !opt->all_extra)  {intra_frac = w->f1;};
	
	if (index < opt->Npart)
	{
		//places particles in initial positions all over substrate

		if (curand_uniform(&states[index]) < intra_frac) //intra
		{
			while (!success)
			{
			tmp_x = -w->max_x + 2 * curand_uniform(&states[index]) * w->max_x; //suggest initial position
			tmp_y = -w->max_y + 2 * curand_uniform(&states[index]) * w->max_y; //suggest initial position
			tmp_z = -w->max_z + 2 * curand_uniform(&states[index]) * w->max_z; //suggest initial position
			
			inside = is_particle_in_any_cell(tmp_x, tmp_y, w, centre_x, centre_y, radii, table, cell_idx); //determine whether particle is in any cell
			
			if (inside)
			{
				idx = index;
				x[idx] = tmp_x;
				y[idx] = tmp_y;
				z[idx] = tmp_z;
				loc[idx] = 1;
				success = true;
			}
			}
		}
		else
		{
			while (!success)
			{
			tmp_x = -w->max_x + 2 * curand_uniform(&states[index]) * w->max_x; //suggest initial position
			tmp_y = -w->max_y + 2 * curand_uniform(&states[index]) * w->max_y; //suggest initial position
			tmp_z = -w->max_z + 2 * curand_uniform(&states[index]) * w->max_z; //suggest initial position
			
			inside = is_particle_in_any_cell(tmp_x, tmp_y, w, centre_x, centre_y, radii, table, cell_idx); //determine whether particle is in any cell
			
			if (!inside)
			{
				idx = index;
				x[idx] = tmp_x;
				y[idx] = tmp_y;
				z[idx] = tmp_z;
				loc[idx] = 0;
				success = true;
			}
			}
		}
	}
}



int main(void)
{
	//OPENING STATEMENTS
	//__________________________________________________________________
	std::clock_t start;
	float duration;
	start = std::clock();
	cudaError error = cudaSuccess;
	int nDevices;
	cudaGetDeviceCount(&nDevices);
	printf("\n PARALLEL SIMULATIONS OF DIFFUSION WITH EXCHANGE (SIMPLEX) \n");
	printf(" --- Hardware information ---\n");
	printf("Number of GPUs: %d\n", nDevices);
	int activeDevice;
	cudaGetDevice(&activeDevice);
	printf("Active GPU index: %d\n", activeDevice);
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, activeDevice);
	printf("GPU name: %s\n", prop.name);
	printf(" --------------------  \n");
	//__________________________________________________________________
	
	//DECLARATIONS
	//__________________________________________________________________
	bool PROCEED = true; //If false, simulation will abort
	//non-numerical options, host-only
	char pos_fn[500], sta_fn[500], sub_fn[500], gwf_fn[500], sig_fn[500], phase_fn[500]; //filenames for trajectories, states, substrate, gradient, signal and phase
	//options structure for host
	options* opt, * dev_opt;
	//declare gradient waveform arrays
	float *h_gwf_x, *h_gwf_y, *h_gwf_z; //waveform in x y z
	float *dev_gwf_x, *dev_gwf_y, *dev_gwf_z;
	//declare substrate structure for device and host
	world* h_w, * dev_w;
	//declare substrate arrays on device and host
	float* h_centre_x, *h_centre_y, *h_radii;
	float* dev_centre_x, *dev_centre_y, *dev_radii;
	long long *h_table, *h_cell_idx;
	long long *dev_table, *dev_cell_idx;
	//declare traj arrays and particle location (compartment id)
	int* h_loc, * dev_loc;
	float* h_x, * h_y, * h_z; //for the host
	float* dev_x, * dev_y, * dev_z; //for the device
	//cuRand states
	curandState* dev_states;
	//declare signal and phase arrays
	float* h_signal, *h_phase; //for the host
	float * dev_phase; //for the device
	//__________________________________________________________________
	
	//OPTIONS
	//__________________________________________________________________
	opt = (options*)malloc(sizeof(options));
	PROCEED = PROCEED && set_options(opt, pos_fn, sta_fn, sub_fn, gwf_fn, sig_fn, phase_fn);
	if (PROCEED){
	error = cudaMalloc(&dev_opt, sizeof(options));
	if (error != cudaSuccess){
		std::cout << "\n ERROR: FAILED TO ALLOCATE OPTIONS MEMORY ON GPU! \n" << std::endl;
		PROCEED = false;}
	}
	//__________________________________________________________________
	
	//GRADIENT WAVEFORMS
	//__________________________________________________________________
	//get info about gradient waveforms (how many and how many time points in each)
	if (PROCEED){
		PROCEED = get_gwf_info(gwf_fn, opt);}
		
	if (PROCEED){
	//allocate gradient waveform memory on host
	h_gwf_x = (float*)malloc(opt->n_acq*opt->n_gwf_points* sizeof(float));
	h_gwf_y = (float*)malloc(opt->n_acq*opt->n_gwf_points* sizeof(float));
	h_gwf_z = (float*)malloc(opt->n_acq*opt->n_gwf_points* sizeof(float));
	//load gradient waveform
	PROCEED = load_gwf(gwf_fn, opt, h_gwf_x, h_gwf_y, h_gwf_z);}
	if (PROCEED){
	//allocate gradient waveform memory on device
	cudaMalloc(&dev_gwf_x, opt->n_acq*opt->n_gwf_points * sizeof(float));
	cudaMalloc(&dev_gwf_y, opt->n_acq*opt->n_gwf_points * sizeof(float));
	error = cudaMalloc(&dev_gwf_z, opt->n_acq*opt->n_gwf_points * sizeof(float));
	if (error != cudaSuccess){
		std::cout << "\n ERROR: FAILED TO ALLOCATE GRADIENT WAVEFORM MEMORY ON GPU! \n" << std::endl;
		PROCEED = false;}
	else{
		cudaMemcpy(dev_gwf_x, h_gwf_x, opt->n_acq*opt->n_gwf_points * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_gwf_y, h_gwf_y, opt->n_acq*opt->n_gwf_points * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_gwf_z, h_gwf_z, opt->n_acq*opt->n_gwf_points * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_opt, opt, sizeof(options), cudaMemcpyHostToDevice); //copy updated opt to device
		}
	}
	//__________________________________________________________________
	
	//SIMULATION SUBSTRATE
	//__________________________________________________________________
	//substrate
	if (PROCEED){
	h_w = (world*)malloc(sizeof(world));
	//Load info about the substrate so we know how much memory to allocate for it
	PROCEED = get_substrate_info(opt, h_w, sub_fn);}
	
	//allocate substrate memory on host
	if (PROCEED)
	{
	h_centre_x = (float*)malloc(h_w->num_cells * sizeof(float));
	h_centre_y = (float*)malloc(h_w->num_cells * sizeof(float));
	h_radii = (float*)malloc(h_w->num_cells * sizeof(float));
	h_table = (long long*)malloc(h_w->num_voxels * sizeof(long long));
	h_cell_idx = (long long*)malloc(h_w->num_voxels * sizeof(long long));
	
	//load the substrate
	PROCEED = load_substrate(h_centre_x, h_centre_y, h_radii, h_table, h_cell_idx, h_w, sub_fn);
	}

	if (PROCEED)
	{
	//allocate substrate memory on device and check if we stil have space
	cudaMalloc(&dev_centre_x, h_w->num_cells * sizeof(float));
	cudaMalloc(&dev_centre_y, h_w->num_cells * sizeof(float));
	cudaMalloc(&dev_radii, h_w->num_cells * sizeof(float));
	cudaMalloc(&dev_cell_idx, h_w->num_voxels * sizeof(long long));
	cudaMalloc(&dev_w, sizeof(world));
	error = cudaMalloc(&dev_table, h_w->num_voxels * sizeof(long long));
	if (error != cudaSuccess)
	{
		std::cout << "\n ERROR: FAILED TO ALLOCATE SUBSTRATE MEMORY ON GPU! \n" << std::endl;
		PROCEED = false;//throw error;
	}
	else {	
		//copy the substrate to the GPU 
		cudaMemcpy(dev_table, h_table, h_w->num_voxels * sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_w, h_w, sizeof(world), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_centre_x, h_centre_x, h_w->num_cells * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_centre_y, h_centre_y, h_w->num_cells * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_radii, h_radii, h_w->num_cells * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_cell_idx, h_cell_idx, h_w->num_voxels * sizeof(long long), cudaMemcpyHostToDevice);
		}
	}
	//__________________________________________________________________
	
	//SIGNALS AND PHASES
	//__________________________________________________________________
	//allocate them on host
	if (PROCEED)
	{
	h_signal = (float*)malloc(opt->n_acq * sizeof(float));
	h_phase = (float*)malloc(opt->n_acq*opt->Npart * sizeof(float));
	//initialise phase array, it's important that it contains only zeros
	for (long c = 0; c < opt->n_acq*opt->Npart; c++) {h_phase[c] = 0;}
	//allocate memory for arrays on device
	error = cudaMalloc(&dev_phase, opt->n_acq*opt->Npart * sizeof(float));
	
	if (error != cudaSuccess)
	{
		std::cout << "\n ERROR: FAILED TO ALLOCATE PHASE MEMORY ON GPU! \n" << std::endl;
		PROCEED = false;//throw error;
	}
	else {	
		//copy phase array to device
		cudaMemcpy(dev_phase, h_phase, opt->n_acq*opt->Npart * sizeof(float), cudaMemcpyHostToDevice);
		}
	}
	//__________________________________________________________________


	//FINAL POSITIONS AND TRANSITION HISTORY
	//__________________________________________________________________
	if (PROCEED)
	{
	//allocate them on host
	h_loc = (int*)malloc(opt->Npart* sizeof(int));
	h_x = (float*)malloc(opt->Npart * sizeof(float));
	h_y = (float*)malloc(opt->Npart * sizeof(float));
	h_z = (float*)malloc(opt->Npart * sizeof(float));

	//allocate memory for arrays on device
	cudaMalloc(&dev_loc, opt->Npart * sizeof(int));
	cudaMalloc(&dev_x, opt->Npart * sizeof(float));
	cudaMalloc(&dev_y, opt->Npart * sizeof(float));
	error = cudaMalloc(&dev_z, opt->Npart * sizeof(float));

	if (error != cudaSuccess)
	{
		std::cout << "\n ERROR: FAILED TO ALLOCATE TRAJECTORY MEMORY ON GPU. \n" << std::endl;
		PROCEED = false;//throw error;
	}
	else {
			//copy x,y,z  and id arrays to device
		cudaMemcpy(dev_loc, h_loc, opt->Npart * sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_x, h_x, opt->Npart * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_y, h_y, opt->Npart * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_z, h_z, opt->Npart * sizeof(float), cudaMemcpyHostToDevice);
		}	
	}
	//__________________________________________________________________


	//DEFINE GRID TEXTURE
	//__________________________________________________________________
	int blockSize;
	64 > opt->Npart ? blockSize = (int)opt->Npart : blockSize = 64;
	int numBlocks = (int)(opt->Npart + blockSize - 1) / blockSize; //make sure to round up in case N is not an integer multiple of blockSize
	
	//allocate curandState for every CUDA thread on the host
	error = cudaMalloc(&dev_states, blockSize * numBlocks * sizeof(curandState));
	if (error != cudaSuccess)
	{
		std::cout << "\n ERROR: FAILED TO ALLOCATE RNG STATES MEMORY ON GPU. \n" << std::endl;
		PROCEED = false;//throw error;
	}
	//__________________________________________________________________
	
	
	//BEGIN SIMULATION
	//__________________________________________________________________
	if (PROCEED)
	{
	std::cout << "Generating GPU grid texture.." << std::endl;
	std::cout << "Number of blocks = " << numBlocks << " AND block size = " << blockSize << std::endl;
	
	std::cout << "Generating initial spin distribution..." << std::endl;
	
	//initialise RNG for all threads
	random_init << < numBlocks, blockSize >> > (dev_states);
	//generate initial particle distribution
	generate_initial_distribution << < numBlocks, blockSize >> > (dev_x, dev_y, dev_z, dev_loc, dev_centre_x, dev_centre_y, dev_radii, dev_table, dev_cell_idx, dev_states, dev_opt, dev_w);
	cudaDeviceSynchronize(); //Tell CPU to wait until kernel is done before accessing results. This is necessary because
							//cuda kernel launches do not block the calling CPU thread.
	//launch simulation engine
	std::cout << "Running simulation..." << std::endl;
	engine << < numBlocks, blockSize >> > (dev_x, dev_y, dev_z, dev_loc, dev_centre_x, dev_centre_y, dev_radii, dev_table, dev_cell_idx, dev_states, dev_opt, dev_w, dev_phase,
	dev_gwf_x, dev_gwf_y, dev_gwf_z);
	cudaDeviceSynchronize(); 
	//__________________________________________________________________
	
	//DOWNLOAD AND SAVE RESULTS
	//__________________________________________________________________
	//copy phase array back to host machine
	cudaMemcpy(h_phase, dev_phase, opt->n_acq*opt->Npart * sizeof(float), cudaMemcpyDeviceToHost);
	//copy final particle positions for diagnostics
	cudaMemcpy(h_x, dev_x, opt->Npart * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_y, dev_y, opt->Npart * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_z, dev_z, opt->Npart * sizeof(float), cudaMemcpyDeviceToHost);

	convert_phase_to_signal(h_phase, h_signal, opt);
	//write results to binary files
	if (opt->save_phase) {save_phase(h_phase, opt, phase_fn);}
	save_signal(h_signal, sig_fn, opt);
	if (opt->save_positions) {save_trajectory(h_x, h_y, h_z, pos_fn, opt);};
	if (opt->save_states) {save_state_history(h_loc, sta_fn, opt);};
	//__________________________________________________________________
	
	duration = (std::clock() - start) / (float)CLOCKS_PER_SEC;
	std::cout << "Elapsed time is: " << duration << " seconds." << std::endl;
	}
	else
	{
		std::cout << "\n FAILED TO LAUNCH SIMULATION ENGINE. ABORTING. \n" << std::endl;
	}

	//FREE RESOURCES
	//__________________________________________________________________
	// Free memory on host
	free(h_signal);
	free(h_phase);
	free(h_gwf_x);
	free(h_gwf_y);
	free(h_gwf_z);
	free(h_x);
	free(h_y);
	free(h_z);
	free(h_w);
	free(opt);
	free(h_loc);
	free(h_centre_x);
	free(h_centre_y);
	free(h_radii);
	free(h_table);
	free(h_cell_idx);
	//free memory on device
	cudaFree(dev_phase);
	cudaFree(dev_gwf_x);
	cudaFree(dev_gwf_y);
	cudaFree(dev_gwf_z);
	cudaFree(dev_x);
	cudaFree(dev_y);
	cudaFree(dev_z);
	cudaFree(dev_states);
	cudaFree(dev_loc);
	cudaFree(dev_w);
	cudaFree(dev_opt);
	cudaFree(dev_centre_x);
	cudaFree(dev_centre_y);
	cudaFree(dev_radii);
	cudaFree(dev_table);
	cudaFree(dev_cell_idx);
	//__________________________________________________________________
	return 0;
}
