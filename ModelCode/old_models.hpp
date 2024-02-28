/* 
MIT License Copyright 2024 Francesco Savelli - see LICENSE file.

This is some of the code currently used in models.hpp plus legacy code used to
setup different simulations that were not used in Savelli JNeurosci 2024. Some
of this code was adapted from simulations that were used in Savelli&Knierim
JNeurophy 2010, but some issues of compatibility may arise with the newer
version of the codebase. 
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


// Assemble various network architectures

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
                FILE* log_file){

  
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
    for (unsigned int distinct_poses=0; 
                      distinct_poses<pars.grid_cells_number; distinct_poses++){ 
      /* legacy GSL implementation
      double x=gsl_rng_uniform(connect_rng)*module.scale;
      double y=gsl_rng_uniform(connect_rng)*module.scale;
      */
      uniform_real_distribution<float> URD(0.0, 1.0);
      double x=URD(*connect_rng)*scale_in_pixels;
      double y=URD(*connect_rng)*scale_in_pixels;
      peaks_discounts_t peaks_discounts(100, 100, pars.grid_min_fraction_peak, grid_rng);
      // add to the population of grid cells
      grid_cells.push_back(
          grid_cell_t(x, y, module.orientation, scale_in_pixels, pars.grid_rise, 
            module.stretch_direction, module.stretch_magnitude, 
            module.shear_direction, module.shear_magnitude,
            pars.grid_peak_rate, peaks_discounts, spike_rng));
      // keep track of the index of this grid cell in the module
      module.grids_indices.push_back(grid_cells.size()-1);
      // keep track on whether the spikes from this grid are to be recorded
      grids_to_save.push_back(true);
      grid_tags.push_back(string("GC")+itos(scale_step)+
                          string("_")+itos(distinct_poses));
    }
    printf("%u grid cells in this module \n\n", module.grids_indices.size());
    fprintf(log_file, "%u grid cells in this module \n\n", module.grids_indices.size());
    modules.push_back(module); // add to the module collection
  }
}

void build_nonmodular_grid_cell_pool(const parameters_record_t& pars, 
                          mt19937* connect_rng, mt19937* spike_rng,
                          const vector<double>& posxs, const vector<double>& posys,
                          vector<grid_cell_t>& grid_cells, vector<string>& grid_tags, 
                          vector<bool>& grids_to_save){

  double peak_rate=pars.grid_peak_rate;
  double rise=pars.grid_rise;
  double smallest_scale=pars.grid_smallest_scale;
  double scale_increment=(pars.grid_largest_scale-pars.grid_smallest_scale)/
                                                    (pars.grid_scale_steps - 1);

  double max_x=*(max_element(posxs.begin(), posxs.end()));
  double min_x=*(min_element(posxs.begin(), posxs.end()));
  double max_y=*(max_element(posys.begin(), posys.end()));
  double min_y=*(min_element(posys.begin(), posys.end()));
  double enclosure_center_x=(max_x-min_x)/2;
  double enclosure_center_y=(max_y-min_y)/2;
  cerr<<"Enclosure Center "<<enclosure_center_x<<" "<<enclosure_center_y<<endl;

  // specify grid inputs
  neural_idx_t grid_name_counter=0;
  for (unsigned int scale_step=0; scale_step<pars.grid_scale_steps; scale_step++){
    double scale=(smallest_scale+(scale_step*scale_increment))/
                                        pars.camera_resolution; // pixels
    cerr<<"Scale: "<<scale<<endl;
    // legacy GSL implementation
    // unsigned int first_orientation=gsl_rng_uniform_int(connect_rng, 6); // degrees
    //
    uniform_int_distribution<unsigned long int> UID(0,5);
    unsigned int first_orientation=UID(*connect_rng);
    for (unsigned int orientation_step=0; orientation_step<10; orientation_step++){
      unsigned int orientation=first_orientation+orientation_step*6; // degrees
      for (unsigned int distinct_poses=0; distinct_poses<10; distinct_poses++){ 
        // legacy GSL implementation
        //double x=min_x+gsl_rng_uniform(connect_rng)*(max_x-min_x);
        //double y=min_y+gsl_rng_uniform(connect_rng)*(max_y-min_y);
        //
        uniform_real_distribution<float> URD(0.0, 1.0);
        double x=min_x+URD(*connect_rng)*(max_x-min_x);
        double y=min_y+URD(*connect_rng)*(max_y-min_y);
        double arbitrary_stretch_direction=45.0, stretch_magnitude=pars.grid_stretch;
        double arbitrary_shear_direction=0.0, shear_magnitude=pars.grid_shear;
        grid_cells.push_back(
            grid_cell_t(x, y, orientation, scale, rise, 
              arbitrary_stretch_direction, stretch_magnitude, 
              arbitrary_shear_direction, shear_magnitude,
              peak_rate, spike_rng));
        grids_to_save.push_back(true);
        grid_tags.push_back(string("GC_")+itos(grid_name_counter));
        grid_name_counter++;
      }
    }
  }
}

void insert_fb_interneurons(mt19937* connect_rng,
                              network_t& net,
                              vector<string>& cell_tags,
                              vector<bool>& cells_to_save,
                              vector<bool>& weights_to_save){

  // specify interneurons

  double Cm=2.0;
  double Vthr=-50.0;
  double Vres=-70.0;
  double Vap=30.0;
  double initial_V=-60.0;
  double timestep=1.0;
  delay_t refractory_steps=3;

  
  for (neural_idx_t in_id=0; in_id<50; in_id++){
    net.add_cell(
        cell_t(Cm, Vthr, Vres, Vap, initial_V, timestep, refractory_steps, 100));
    cells_to_save.push_back(true);
    weights_to_save.push_back(false);
    cell_tags.push_back(string("IN_")+itos(in_id));
  }

  // specify synapses to interneurons
  for (neural_idx_t pc_id=0; pc_id<500; pc_id++){
    
    // start sampling from the pool
    /* legacy GSL implementation
    neural_idx_t interneurons_range[50];
    for (neural_idx_t range_ind=0; range_ind<50; range_ind++)
      interneurons_range[range_ind]=range_ind;
    neural_idx_t random_afferents[40];
    gsl_ran_choose(connect_rng, random_afferents, 40, interneurons_range, 50, 
                                              sizeof(neural_idx_t));
    */
    vector<neural_idx_t> random_afferents(40), interneurons_range(50);
    for (neural_idx_t range_ind=0; range_ind<interneurons_range.size(); range_ind++) 
      interneurons_range[range_ind]=range_ind;
    cpp_ran_choose(connect_rng, random_afferents, interneurons_range);
    // end sampling process, the random input have been determined and are in
    // "random afferents"

    for (neural_idx_t ra_id=0; ra_id<40; ra_id++){

      //double gDelta=0.4; ORIGINAL
      double gDelta=0.8;
      double gtau=2.0;
      double E=0.0;
      double timestep=1.0;
      delay_t transmission_delay=2;

      synapse_t new_synapse=synapse_t(gDelta, gtau, E, timestep);//, 100);
      new_synapse.set_connections(pc_id, 500+random_afferents[ra_id], transmission_delay);
      net.add_cell_to_cell(new_synapse);
    }
  }

  // specify synapses from interneurons
  for (neural_idx_t in_id=0; in_id<50; in_id++){

    // start sampling from the pool
    /* legacy GSL implementation
    neural_idx_t pcells_range[500];
    for (neural_idx_t range_ind=0; range_ind<500; range_ind++)
      pcells_range[range_ind]=range_ind;
    neural_idx_t random_afferents[300];
    gsl_ran_choose(connect_rng, random_afferents, 300, pcells_range, 500, 
                                          sizeof(neural_idx_t));
    */
    vector<neural_idx_t> random_afferents(300), pcells_range(500);
    for (neural_idx_t range_ind=0; range_ind<pcells_range.size(); range_ind++) 
      pcells_range[range_ind]=range_ind;
    cpp_ran_choose(connect_rng, random_afferents, pcells_range);
    // end sampling process the random input have been determined and are in 
    // "random afferents"

    for (neural_idx_t ra_id=0; ra_id<300; ra_id++){

      double gDelta=0.2;
      double gtau=6.0;
      double E=-70.0;
      double timestep=1.0;
      delay_t transmission_delay=2;

      synapse_t new_synapse=synapse_t(gDelta, gtau, E, timestep);//, 100);
      new_synapse.set_connections(500+in_id, random_afferents[ra_id], transmission_delay);
      net.add_cell_to_cell(new_synapse);
    }
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
                    vector<bool>& weights_to_save, 
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
    net.add_input(input_t(pars.fr_tau));

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
      net.add_cell(
          cell_t(Cm, Vthr, Vres, Vap, initial_V, timestep, refractory_steps, fr_tau));
      cells_to_save.push_back(true);
      weights_to_save.push_back(false);
      string place_cell_tag=string("PC")+itos(group_idx)+string("_")+itos(pc_idx);
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

  
  //insert_fb_interneurons(connect_rng, net, cell_tags, cells_to_save, weights_to_save);
}

void build_nonmodular_network(const parameters_record_t& pars, mt19937* connect_rng,  
                    const vector<double>& posxs, const vector<double>& posys,
                    const vector<grid_cell_t>& grid_cells, network_t& net,
                    vector<string>& cell_tags, vector<bool>& cells_to_save,
                    vector<bool>& weights_to_save, ofstream& connectivity_file){

  neural_idx_t number_grid_cells=1000;
  neural_idx_t number_place_cells=pars.place_cells_number;
  neural_idx_t number_grids_per_place=pars.grids_per_cell_number;

  for (neural_idx_t input_index=0; input_index<grid_cells.size(); input_index++)
    net.add_input(input_t(pars.fr_tau));

  // specify place cells
  double Cm=2.0;
  double Vthr=-50.0;
  double Vres=-70.0;
  double Vap=30.0;
  double initial_V=-60.0;
  double timestep=1.0;
  delay_t refractory_steps=3;
  double fr_tau=pars.fr_tau;

  for (neural_idx_t pc_id=0; pc_id<number_place_cells; pc_id++){
    net.add_cell(
        cell_t(Cm, Vthr, Vres, Vap, initial_V, timestep, refractory_steps, fr_tau));
    cells_to_save.push_back(true);
    weights_to_save.push_back(false);
    cell_tags.push_back(string("PC_")+itos(pc_id));
  }

  for (neural_idx_t pc_id=0; pc_id<number_place_cells; pc_id++){
    
    // start sampling from the pool
    
    /* legacy GSL implementation
    neural_idx_t grid_cells_range[number_grid_cells];
    for (neural_idx_t range_ind=0; range_ind<number_grid_cells; range_ind++)
      grid_cells_range[range_ind]=range_ind;
    neural_idx_t random_afferents[number_grids_per_place];
    gsl_ran_choose(connect_rng, random_afferents, number_grids_per_place, 
                  grid_cells_range, number_grid_cells, sizeof(neural_idx_t));
    */
    vector<neural_idx_t> random_afferents(number_grids_per_place); 
    vector<neural_idx_t> grid_cells_range(number_grid_cells);
    for (neural_idx_t range_ind=0; range_ind<grid_cells_range.size(); range_ind++) 
      grid_cells_range[range_ind]=range_ind;
    cpp_ran_choose(connect_rng, random_afferents, grid_cells_range);
    // end sampling process the random input have been determined 
    // and are in "random afferents"
    connectivity_file<<"Cell #"<<pc_id<<" <--- ";
    for (neural_idx_t ag_id=0; ag_id<number_grids_per_place; ag_id++){
      double gDelta=pars.gDelta;
      double gtau=pars.gtau;
      double E=0.0;
      double timestep=1.0;
      delay_t transmission_delay=2;
      
      synapse_t new_synapse=synapse_t(gDelta, gtau, E, timestep);
      new_synapse.set_connections(random_afferents[ag_id], pc_id, transmission_delay);
      new_synapse.set_rate_plasticity(pars.plasticity_rule,
                    pars.theta_d*0.001, pars.theta_p*0.001, pars.cut_off*0.001,
                    pars.plasticity_rate, pars.IPI);
      net.add_input_to_cell(new_synapse);
      connectivity_file<<random_afferents[ag_id]<<" ";
    }
    connectivity_file<<endl;
  }
  //insert_fb_interneurons(connect_rng, net, cell_tags, cells_to_save, weights_to_save);
}

// older for comparison/debug ---------------------------------------------
void old_build_network(const parameters_record_t& pars,
                    mt19937* connect_rng, 
                    const vector<grid_cell_t>& grid_cells, network_t& net,
                    vector<string>& cell_tags, vector<bool>& cells_to_save,
                    vector<bool>& weights_to_save, ofstream& connectivity_file){

  cerr<<endl<<endl<<" grid_cells.size "<<grid_cells.size()<<endl<<endl;
  neural_idx_t number_grid_cells=1000;
  neural_idx_t number_place_cells=pars.place_cells_number;
  neural_idx_t number_grids_per_place=pars.grids_per_cell_number;

  for (neural_idx_t input_index=0; input_index<grid_cells.size(); input_index++)
    net.add_input(input_t(pars.fr_tau));

  // specify place cells
  double Cm=2.0;
  double Vthr=-50.0;
  double Vres=-70.0;
  double Vap=30.0;
  double initial_V=-60.0;
  double timestep=1.0;
  delay_t refractory_steps=3;
  double fr_tau=pars.fr_tau;

  for (neural_idx_t pc_id=0; pc_id<number_place_cells; pc_id++){
    net.add_cell(
        cell_t(Cm, Vthr, Vres, Vap, initial_V, timestep, refractory_steps, fr_tau));
    cells_to_save.push_back(true);
    weights_to_save.push_back(false);
    cell_tags.push_back(string("PC_")+itos(pc_id));
  }

  for (neural_idx_t pc_id=0; pc_id<number_place_cells; pc_id++){
    
    // start sampling from the pool
    
    /* legacy GSL implementation
    neural_idx_t grid_cells_range[number_grid_cells];
    for (neural_idx_t range_ind=0; range_ind<number_grid_cells; range_ind++)
      grid_cells_range[range_ind]=range_ind;
    neural_idx_t random_afferents[number_grids_per_place];
    gsl_ran_choose(connect_rng, random_afferents, number_grids_per_place, 
                  grid_cells_range, number_grid_cells, sizeof(neural_idx_t));
    */
    vector<neural_idx_t> random_afferents(number_grids_per_place); 
    vector<neural_idx_t> grid_cells_range(number_grid_cells);
    for (neural_idx_t range_ind=0; range_ind<grid_cells_range.size(); range_ind++) 
      grid_cells_range[range_ind]=range_ind;
    cpp_ran_choose(connect_rng, random_afferents, grid_cells_range);
    // end sampling process the random input have been determined 
    // and are in "random afferents"
    connectivity_file<<"Cell #"<<pc_id<<" <--- ";
    for (neural_idx_t ag_id=0; ag_id<number_grids_per_place; ag_id++){
      double gDelta=pars.gDelta;
      double gtau=pars.gtau;
      double E=0.0;
      double timestep=1.0;
      delay_t transmission_delay=2;
      
      synapse_t new_synapse=synapse_t(gDelta, gtau, E, timestep);
      new_synapse.set_connections(random_afferents[ag_id], pc_id, transmission_delay);
      new_synapse.set_rate_plasticity(pars.plasticity_rule,
                    pars.theta_d*0.001, pars.theta_p*0.001, pars.cut_off*0.001,
                    pars.plasticity_rate, pars.IPI);
      net.add_input_to_cell(new_synapse);
      connectivity_file<<random_afferents[ag_id]<<" ";
    }
    connectivity_file<<endl;
  }
  //insert_fb_interneurons(connect_rng, net, cell_tags, cells_to_save, weights_to_save);
}

#endif
