/* 
MIT License Copyright 2024 Francesco Savelli - see LICENSE file.

Code for reading and handling the simulation parameters provided in the file
contained in the simulation folder.

*/

#ifndef SIMIO_HPP 
#define SIMIO_HPP

#include <string>
#include "network.hpp"

using namespace std;

typedef unsigned long long int  timestamp_t;
typedef vector<vector<string> > string_matrix_t;

plasticity_kind_t str2plasticity(string from);
string plasticity2str(plasticity_kind_t from);

void read_pos_file(const string& namefile, 
                 vector<timestamp_t>& timestamps,
                 vector<double>& posxs, vector<double>& posys);
 
string_matrix_t& read_params_file(const string& namefile);

struct parameters_record_t{

	string pos_file_relpath; // relative path of the trajectory data

	timestamp_t session_start_time; // first useful timestamp
	timestamp_t session_end_time; // last useful timestamp
	unsigned int start; // skip these many minutes after session_start_time 
	unsigned int duration; // length of simulaton in minutes
  	bool phase2; // whether to concatenate a repeat of the session
	// specifies cell_tags for which cells need weighs saved periodically 
	vector<string> weights_to_save; 
	double camera_resolution; // size (cm) of camera-tracked pos file pixels

	double grid_smallest_scale; // smallest scale for the grid input (cm)
	double grid_largest_scale; // smallest scale for the grid input (cm)
  	double grid_stretch; // magnitude of compression (<1) or expansion (>1)
	double grid_shear; // magnitude of shearing (0.0 is no shearing)
	// how many discrete scales in the simulations
	// ranging from the smallest to the largest scale	
	double grid_min_fraction_peak; // used for intervertex peak variability
	unsigned int grid_scale_steps; 
	unsigned int modules_per_place_cell;
	double module_scale_ratio;
	vector<double> module_orientations; // orientation for each module
	vector<double> module_stretch_directions; 
	vector<double> module_stretch_magnitudes; 
	vector<double> module_shear_directions; 
	vector<double> module_shear_magnitudes; 
	// max rate of a grid cell in at the center of the vertex	
	double grid_peak_rate;
	double grid_rise; // regulates how much of the field "emerges"
  	// how many grid cells are created; overall or per module in modular simulations
	neural_idx_t grid_cells_number;
 	// how many place cells are created; overall or per module in modular simulations
	neural_idx_t place_cells_number; 
  	// how many grid cells are in input to any single place cell
	neural_idx_t grids_per_cell_number; 
	
	double gDelta;
	double gtau;
	plasticity_kind_t plasticity_rule;
  	//inter plasticity interval (how often plasticity rule is applied)
	delay_t IPI; 
	double plasticity_rate; // multiplicative constant
	double theta_p;
	double theta_d;
	double cut_off;
	double fr_tau;
	double prefr_tau;

	unsigned long int spiking_rng_seed;
	
	//static const unsigned int how_many_parameters=26;

	parameters_record_t(){

		session_start_time=0;
		session_end_time=0; // last useful timestamp
		start=0; // minutes start
		duration=0; //minutes duration
    	phase2=false;
		
		camera_resolution=0; // size (cm) of camera-tracked pos file pixels

		grid_smallest_scale=0; // smallest scale for the grid input (cm)
		grid_largest_scale=0; // largest scale for the grid input (cm)
		grid_scale_steps=0;
		modules_per_place_cell=0;
		module_scale_ratio=1.5;
		grid_stretch=1;
		grid_shear=0;
		grid_peak_rate=0;
		grid_min_fraction_peak=1.0;
		grid_rise=0;
		grid_cells_number=0;
		place_cells_number=0;
		grids_per_cell_number=0;

		gDelta=0;
		gtau=0;
		plasticity_rule=plasticity_off;
		plasticity_rate=0;
		IPI=100;
		
		theta_p=0;
		theta_d=0;
		cut_off=0;

		fr_tau=100;
		prefr_tau=100;
		

		spiking_rng_seed=1971;
	};

	void from_file(const string& namefile);
  	
};

void print_params(const parameters_record_t& pars);

#endif

