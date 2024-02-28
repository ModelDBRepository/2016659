/* 
MIT License Copyright 2024 Francesco Savelli - see LICENSE file.

Procedures to assemble networks of grid inputs and place cells in the desired
configuration for the simulations. Mostly relies on the templates and engine
provided in grid_input.hpp and network.hpp.
*/

#ifndef MODELS_HPP
#define MODELS_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

/* Older implementation depended on the Gnu Scientific Library for random
 * number generation. Current implementation uses equivalent rng facilities now
 * available in the C++ standard library. Older GSL calls are commented out for
 * record. New implementation produces the same sequences of random numbers as
 * the older one would.

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
*/
#include <random>

#include "grid_input.hpp"
#include "network.hpp"
#include "util.hpp"
#include "simIO.hpp"

#define PI 3.14159265

using namespace std;

// to choose randomly from a set (C++17 has std::sample which should be equivalent 
// coded after gsl_ran_choose implementation, assumes source.size() >= dest.size()
inline void cpp_ran_choose(mt19937* cpp_rng, 
    vector<neural_idx_t>& dest, const vector<neural_idx_t>& source){

  int source_idx=0;
  int dest_idx=0;
  int source_n=source.size();
  int dest_n=dest.size();

  uniform_real_distribution<float> URD(0.0, 1.0);
  while (source_idx<source_n && dest_idx<dest_n){
    if ((source_n-source_idx)*URD(*cpp_rng) < dest_n-dest_idx){
      dest[dest_idx]=source[source_idx];
      dest_idx++;
    }
    source_idx++;
  }
}


struct grid_module_t{

  double scale;
  double orientation;
  double stretch_direction;
  double stretch_magnitude;
  double shear_direction;
  double shear_magnitude;
  
  // indexes to the vector of grid objects produced externally to keep 
  // track of partnership to this module
  vector<neural_idx_t> grids_indices; 
};

void build_fixed_modular_grid_cell_pool(const parameters_record_t& pars, 
                mt19937* connect_rng, mt19937* grid_rng, mt19937* spike_rng,
                vector<grid_cell_t>& grid_cells, vector<string>& grid_tags, 
                vector<bool>& grids_to_save, vector<grid_module_t>& modules,
                FILE* grid_parameters_file, FILE* log_file){

  
  if (!(pars.grid_scale_steps==pars.module_orientations.size() &&
      pars.module_orientations.size()==pars.module_stretch_directions.size() &&
      pars.module_stretch_directions.size()==pars.module_stretch_magnitudes.size() &&
      pars.module_stretch_magnitudes.size()==pars.module_shear_directions.size() &&
      pars.module_shear_directions.size()==pars.module_shear_magnitudes.size())) {
    
    cerr<<endl<<"grid_scale_steps different from number of geometrical properties: "
      <<"scale steps "<<pars.grid_scale_steps<<endl
      <<"orientations "<<pars.module_orientations.size()<<endl
      <<"stretch directions "<<pars.module_stretch_directions.size()<<endl
      <<"stretch magnitudes "<<pars.module_stretch_magnitudes.size()<<endl
      <<"shear directions "<<pars.module_shear_directions.size()<<endl
      <<"shear magnitudes "<<pars.module_shear_magnitudes.size()<<endl;
    exit(1);
  }

  for (unsigned int scale_step=0; scale_step<pars.grid_scale_steps; 
                                                          scale_step++){
    grid_module_t module;
    // eventually make 1.5 scale ratio a parameter read from file?
    module.scale=pars.grid_smallest_scale*
                (pow(pars.module_scale_ratio, scale_step)); // cm
    double scale_in_pixels=module.scale/pars.camera_resolution; // pixels
    /* old code to sample orientation randomly
    uniform_int_distribution<unsigned long int> UID(0, 59);
    module.orientation=UID(*connect_rng); // degrees */
    module.orientation=pars.module_orientations[scale_step];
    module.stretch_direction=pars.module_stretch_directions[scale_step];
    module.stretch_magnitude=pars.module_stretch_magnitudes[scale_step];
    module.shear_direction=pars.module_shear_directions[scale_step];
    module.shear_magnitude=pars.module_shear_magnitudes[scale_step];
    printf("MODULE %u \n", scale_step);
    printf("Scale (cm): %f Orientation (deg): %f \n", 
        module.scale, module.orientation);
    printf("Stretch direction (deg): %f, magnitude: %f \n",
        module.stretch_direction, module.stretch_magnitude);
    printf("Shear direction (deg): %f, magnitude: %f \n",
        module.shear_direction, module.shear_magnitude);
    fprintf(log_file, "MODULE %u \n", scale_step);
    fprintf(log_file, "Scale (cm): %f Orientation (deg): %f \n", 
        module.scale, module.orientation);
    fprintf(log_file, "Stretch direction (deg): %f, magnitude: %f \n",
        module.stretch_direction, module.stretch_magnitude);
    fprintf(log_file, "Shear direction (deg): %f, magnitude: %f \n",
        module.shear_direction, module.shear_magnitude);
    
    // make a reference grid for the module at phase 0, 0 
    // will be needed to keep random phases within the hexagon
    peaks_discounts_t peaks_discounts(100, 100, pars.grid_min_fraction_peak, grid_rng);
    grid_t ref_grid(0.0, 0.0, module.orientation, scale_in_pixels, pars.grid_rise, 
      module.stretch_direction, module.stretch_magnitude, 
      module.shear_direction, module.shear_magnitude, peaks_discounts);

    unsigned int distinct_poses=0; 
    while(distinct_poses<pars.grid_cells_number){ 
      /* legacy GSL implementation 
      (this is also a different distribution than evolved 
      below to sample from within the hexagon)
      double x=gsl_rng_uniform(connect_rng)*module.scale;
      double y=gsl_rng_uniform(connect_rng)*module.scale;
      */
      uniform_real_distribution<float> URD(0.0, 1.0);
      double x=2*URD(*connect_rng)*scale_in_pixels - scale_in_pixels;
      double y=2*URD(*connect_rng)*scale_in_pixels - scale_in_pixels;
      // keep only if it falls within the hexagonal perimeter
      sq_dist_cycles_t sdc = ref_grid.sq_dist_cycles_at(x, y);
      if (-1<sdc.nonround_cycles1 && sdc.nonround_cycles1<1 &&
          -1<sdc.nonround_cycles2 && sdc.nonround_cycles2<1 &&
          -1<sdc.nonround_cycles3 && sdc.nonround_cycles3<1) {
        // make a grid cell with this phase
        peaks_discounts_t peaks_discounts(100, 100, pars.grid_min_fraction_peak, grid_rng);
        // add to the population of grid cells
        tag_t gridtag = string("GC")+itos(scale_step)+
                          string("_")+itos(distinct_poses);
        grid_cells.push_back(
            grid_cell_t(x, y, module.orientation, scale_in_pixels, pars.grid_rise, 
              module.stretch_direction, module.stretch_magnitude, 
              module.shear_direction, module.shear_magnitude,
              pars.grid_peak_rate, peaks_discounts, spike_rng, gridtag));
        // keep track of the index of this grid cell in the module
        module.grids_indices.push_back(grid_cells.size()-1);
        // keep track on whether the spikes from this grid are to be recorded
        grids_to_save.push_back(true);
        grid_tags.push_back(gridtag);
        fprintf(grid_parameters_file, 
            "%s X:%f Y:%f O:%f S:%f R:%f \n", gridtag.c_str(), 
            x, y, module.orientation, scale_in_pixels, pars.grid_rise);
        distinct_poses++;
      }
    }

    printf("%u grid cells in this module \n\n", module.grids_indices.size());
    fprintf(log_file, "%u grid cells in this module \n\n", module.grids_indices.size());
    modules.push_back(module); // add to the module collection
  }
}

// grid_cells below is used only for a sanity check, erase from function parameters
// or fuse this function with build grid cell, eventually
void build_modular_network(const parameters_record_t& pars, mt19937* connect_rng,  
                    const vector<grid_cell_t>& grid_cells, 
                    const vector<string>& grid_tags,
                    const vector<grid_module_t>& modules, 
                    network_t& net,
                    vector<string>& cell_tags, 
                    vector<bool>& cells_to_save,
                    ofstream& connectivity_file,
                    FILE* log_file){

  unsigned int number_grid_modules=pars.grid_scale_steps;
  neural_idx_t number_grid_cells_per_module=pars.grid_cells_number;
  neural_idx_t total_number_grid_inputs=grid_cells.size();
  if (total_number_grid_inputs!=number_grid_modules*number_grid_cells_per_module){
    cerr<<endl<<"total_number_grid_inputs is inconsistent";
    exit(1);
  }
  neural_idx_t number_place_cells_per_group=pars.place_cells_number;
  unsigned int number_modules_per_group=pars.modules_per_place_cell;
  neural_idx_t number_grids_per_place=pars.grids_per_cell_number;
  
  // add grid inputs to the network 
  // these inputs are assumed to have a one-to-one relationship to the grids in 
  // the vector grid_cells built elsewhere to which the indices in each entry of 
  // the vector modules refer to 
  for (neural_idx_t input_index=0; 
          input_index<total_number_grid_inputs; input_index++)
    net.add_input(input_t(pars.fr_tau, grid_tags[input_index]));

  // specify generic place cell
  double Cm=2.0;
  double Vthr=-50.0;
  double Vres=-70.0;
  double Vap=30.0;
  double initial_V=-60.0;
  double timestep=1.0;
  delay_t refractory_steps=3;
  double fr_tau=pars.fr_tau;
  
  // add place cells to the network, one group at a time, each group 
  // taking inputs from a subset of consecutive grid cell modules
  neural_idx_t global_pc_idx=0; // unique place cell index across all groups
  for (int group_idx=0; 
        group_idx <= number_grid_modules-number_modules_per_group; 
        group_idx++){

    printf("Building place cell group %u from grid modules (scales): ", group_idx); 
    fprintf(log_file, 
        "Building place cell group %u from grid modules (scales): ", group_idx);
    // prepare admissible grid inputs for this place cell group
    vector<neural_idx_t> grids_for_this_group;
    for (int input_module=group_idx; 
            input_module<group_idx+number_modules_per_group; input_module++){
      grids_for_this_group.insert(grids_for_this_group.end(), 
                                  modules[input_module].grids_indices.begin(),
                                  modules[input_module].grids_indices.end());
      printf(" %f ", modules[input_module].scale);
      fprintf(log_file, " %f ", modules[input_module].scale);
    }
    printf(" (%u grids in input) \n", grids_for_this_group.size());
    fprintf(log_file, " (%u grids in input) \n", grids_for_this_group.size());
     
    // add each place cell from this group to the network, 
    // plus its synapses from random grid inputs
    for (neural_idx_t pc_idx=0; pc_idx<number_place_cells_per_group; pc_idx++){
      tag_t place_cell_tag=string("PC")+itos(group_idx)+string("_")+itos(pc_idx);
      net.add_cell(
          cell_t(Cm, Vthr, Vres, Vap, initial_V, 
            timestep, refractory_steps, fr_tau, place_cell_tag));
      cells_to_save.push_back(true);
      //weights_to_save.push_back(false);
      cell_tags.push_back(place_cell_tag);
      
      // start sampling from the pool of admissible grid inputs
      
      /* legacy GSL implementation, variables proabably now renamed
      neural_idx_t grid_cells_range[number_grid_cells];
      for (neural_idx_t range_ind=0; range_ind<number_grid_cells; range_ind++)
        grid_cells_range[range_ind]=range_ind;
      neural_idx_t random_afferents[number_grids_per_place];
      gsl_ran_choose(connect_rng, random_afferents, number_grids_per_place, 
                    grid_cells_range, number_grid_cells, sizeof(neural_idx_t));
      */
      
      // vector where the indices for the grid inputs sampled for 
      // this place cell will go 
      vector<neural_idx_t> random_afferents(number_grids_per_place); 
      // sample
      cpp_ran_choose(connect_rng, random_afferents, grids_for_this_group);

      // end of sampling process, the grid inputs for this place cell
      // have been determined and their indices are in "random afferents"

      // specify a generic synapse
      double gDelta=pars.gDelta;
      double gtau=pars.gtau;
      double E=0.0;
      double timestep=1.0;
      delay_t transmission_delay=2;
      
      // add one synapse for each random input
      //connectivity_file<<"Cell #"<<global_pc_idx<<" "<<place_cell_tag<<" <--- ";
      connectivity_file<<place_cell_tag<<" <--- ";
      for (neural_idx_t ag_id=0; ag_id<number_grids_per_place; ag_id++){
        synapse_t new_synapse=synapse_t(gDelta, gtau, E, timestep);
        new_synapse.set_connections(random_afferents[ag_id], 
                              global_pc_idx, transmission_delay);
        new_synapse.set_rate_plasticity(pars.plasticity_rule,
                      pars.theta_d*0.001, pars.theta_p*0.001, pars.cut_off*0.001,
                      pars.plasticity_rate, pars.IPI);
        net.add_input_to_cell(new_synapse);
        //connectivity_file<<random_afferents[ag_id]<<" ";
        connectivity_file<<grid_tags[random_afferents[ag_id]]<<" ";
      }
      connectivity_file<<endl;
      global_pc_idx++;
    }
  }
  
}

#endif
