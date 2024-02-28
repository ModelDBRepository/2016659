/* 
MIT License Copyright 2024 Francesco Savelli - see LICENSE file.

This is the source file producing the executable for running the simulation.
See usage, command line params, and input folder structure in README file.
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>
#include <algorithm>
#include <unistd.h>

/* Older implementation depended on the Gnu Scientific Library for random
 * number generation. Current implementation uses equivalent rng facilities now
 * available in the C++ standard library. Older GSL calls are commented out for
 * record. New implementation produces the same sequences of random numbers as
 * the older one would.

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
*/

#include <random>

#include "util.hpp"
#include "simIO.hpp"
#include "grid_input.hpp"
#include "network.hpp"
#include "models.hpp"


//# include "IO_network.hpp"

#define PI 3.14159265

using namespace std;

void simulate(const vector<timestamp_t>& timestamps, 
      const vector<double>& posxs, const vector<double>& posys,
      timestamp_t start_time, timestamp_t end_time,
      vector<grid_cell_t>& grid_cells, network_t& network,
      const vector<bool>& grids_to_save, const vector<bool>& cells_to_save,
      const vector<tag_t>& grid_tags, const vector<tag_t>& cell_tags,
      const vector<tag_t>& weights_to_save, 
      FILE* spike_file, FILE* log_file, FILE** weight_file_array){  

  timestamp_t max_pos_index=timestamps.size()-1;
  timestamp_t time=start_time;
  timestamp_t pos_index=0;

  double posx; 
  double posy;
  // make auxiliary boolean vector for momentary input spike tracking 
  vector<bool> inputs_that_spiked(grid_cells.size(), false);
  
  clock_t t0= clock();
  clock_t t1=t0;
  if (t0==clock_t(-1))
    fprintf(log_file, "NO TIME 1\n");
  while (time<end_time){
    if (((time-start_time)%10000000)==0){
      clock_t t2= clock();
      if (t2==clock_t(-1))
        fprintf(log_file, "NO TIME 1\n");
      fprintf(log_file, "Simulated %llu seconds in %f real seconds\n",
                  (time-start_time)/1000000, double(t2-t0)/CLOCKS_PER_SEC );

      double realtime_since_start=double(t2-t0)/CLOCKS_PER_SEC;
      timestamp_t simtime_since_start=(time-start_time)/1000000;
      double current_simtorealtime_ratio=(double(t2-t1)/CLOCKS_PER_SEC)/10;
      timestamp_t estimated_remaining_minutes=
            int(current_simtorealtime_ratio*((end_time-time)/1000000) / 60)+1;
      printf("Current simulated-to-real time ratio (over last 10 simulated seconds): %f ", 
                  current_simtorealtime_ratio);
      printf("Estimated remaining time: %llu minutes \n", estimated_remaining_minutes);
      printf("(Simulated %llu seconds in %f real seconds.)\n", 
                  simtime_since_start, realtime_since_start);
      
      fprintf(log_file, 
          "Current simulated-to-real time ratio (over last 10 simulated seconds): %f ", 
           current_simtorealtime_ratio);
      fprintf(log_file, "Estimated remaining time: %llu minutes \n",
          estimated_remaining_minutes);
      fprintf(log_file, "(Simulated %llu seconds in %f real seconds.)\n", 
                  simtime_since_start, realtime_since_start);

      fprintf(log_file, "timestamp %llu\n", time);
      t1=t2;
    }
    // increment the current index in the positional data as needed
    while ((time > timestamps[pos_index+1]) and (pos_index+2<=max_pos_index))
      pos_index=pos_index+1;
    
    // decide which position is the current one
    if (time<= timestamps[pos_index+1])
      if (time-timestamps[pos_index] < timestamps[pos_index+1]-time){
        posx=posxs[pos_index];
        posy=posys[pos_index];
      }
      else{
        posx=posxs[pos_index+1];
        posy=posys[pos_index+1];
      }
    
    
    // simulation step for inputs (grid cells) 
    // it must be: grid_cells.size() <= (network->_inputs).size()
    for(neural_idx_t gc_index=0; gc_index<grid_cells.size(); gc_index++){
      grid_cells[gc_index].step(posx, posy);
      if (grid_cells[gc_index].just_spiked())
        inputs_that_spiked[gc_index]=true;
      else
        inputs_that_spiked[gc_index]=false;
    }
    // save time and position
    fprintf(spike_file, "%llu", time);
    fprintf(spike_file, " (%f", posx);
    fprintf(spike_file, " %f) ", posy);
    // save grid cells spikes 
    for (neural_idx_t gindex=0; gindex<grid_cells.size(); gindex++)
      if (inputs_that_spiked[gindex] && grids_to_save[gindex]) 
        fprintf(spike_file, ", %s", grid_tags[gindex].c_str());
    
    // simulation step throughout the network excepts inputs
    vector<bool> cells_that_spiked=network.step(inputs_that_spiked);
    // save spikes from all other cells
    for (neural_idx_t cindex=0; cindex<cells_that_spiked.size(); cindex++){
      if (cells_that_spiked[cindex] && cells_to_save[cindex])
        fprintf(spike_file, ", %s", cell_tags[cindex].c_str());
    }
    fprintf(spike_file, "\n");
   
    // save synaptic weights every 500 ms of simulated time
    if (((time-start_time)%(1000*500))==0){ 
      for (neural_idx_t cell_idx=0; cell_idx<weights_to_save.size(); cell_idx++){
        tag_t cell_tag=weights_to_save[cell_idx];
        vector<tag_t> afferent_tags; 
        vector<double> afferent_weights;
        network.get_labeled_synweights_for(cell_tag, 
                                    afferent_tags, afferent_weights);
        //fprintf(weight_file_array[cell_idx], "%s %llu ", cell_tag.c_str(), time);
        fprintf(weight_file_array[cell_idx], "%llu ", time);
        for (neural_idx_t ag_id=0; ag_id<afferent_tags.size(); ag_id++){
          fprintf(weight_file_array[cell_idx], "[%s %f] ", 
                    afferent_tags[ag_id].c_str(), afferent_weights[ag_id]);
        }
        fprintf(weight_file_array[cell_idx],"\n");
      }
    } 
    time+=1000; // timestamps in microseconds
  }

}    



int main (int argc, char* argv[]){

  string path(argv[1]);
  string par_path=path+"parameters.txt";
  //string pos_path=path+"Pos.p.ascii";
  cout<<"Path: "<<path<<endl<<"Parameters path: "<<par_path<<endl;

  // initialize output files
  string spike_file_name=path+"spike_file";
  FILE* spike_file=fopen(spike_file_name.c_str(), "w");
  string log_file_name=path+"log_file";
  FILE* log_file=fopen(log_file_name.c_str(), "w");
  string grid_parameters_file_name=path+"grid_parameters_file";
  FILE* grid_parameters_file=fopen(grid_parameters_file_name.c_str(), "w");
  cout<<spike_file_name<<" "<<log_file_name<<" "<<grid_parameters_file_name<<endl;
  string connectivity_file_name=path+"connectivity";
  ofstream connectivity_file(connectivity_file_name.c_str());
  string initial_weights_file_name=path+"initial_weights";
  ofstream initial_weights_file(initial_weights_file_name.c_str());
  string final_weights_file_name=path+"final_weights";
  ofstream final_weights_file(final_weights_file_name.c_str());

  // read simulation parameters
  parameters_record_t pars;
  pars.from_file(par_path);
  print_params(pars);

  string pos_path=pars.pos_file_relpath;
  cout<<"Positional file path: "<<pos_path<<endl;

  // read animal trajectory
  vector<timestamp_t> timestamps;
  vector<double> posxs, posys;
  read_pos_file(pos_path, timestamps, posxs, posys);
  /* legacy GSL implementation
  // c-style allocation and initialization of the two Mersenne-Twister
  // random number generators from the Gnu Scientific Library
  gsl_rng* connect_rng = gsl_rng_alloc(gsl_rng_mt19937); // for connectivity 
  gsl_rng_set(connect_rng, 2005);
  gsl_rng* spike_rng = gsl_rng_alloc(gsl_rng_mt19937); // for Poisson spike trains
  gsl_rng_set(spike_rng, pars.spiking_rng_seed); 
  */

  // initialization of the Mersenne-Twister random number generators 
  // all randomness in the simulation is provided by these two rngs alone
  mt19937* connect_rng=new mt19937{2005}; 
  mt19937* grid_rng=new mt19937{2005}; 
  mt19937* spike_rng=new mt19937{pars.spiking_rng_seed};

  // create neural network
  vector<grid_cell_t> grid_cells;
  network_t network;
  vector<string> grid_tags;
  vector<string> cell_tags;
  vector<bool> grids_to_save;
  vector<bool> cells_to_save;
  vector<string> weights_to_save;
  vector<grid_module_t> grid_modules;
  
  printf("Building network\n");
  fprintf(log_file, "Building network\n");
  
  build_fixed_modular_grid_cell_pool(pars, connect_rng, grid_rng, spike_rng, 
                              grid_cells, grid_tags, grids_to_save, grid_modules,
                              grid_parameters_file, log_file);
  fclose(grid_parameters_file);

  printf("%lu grid cells in input \n", grid_cells.size());
  fprintf(log_file, "%lu grid cells in input \n", grid_cells.size());
  
  build_modular_network(pars, connect_rng, grid_cells, grid_tags, grid_modules,
    network, cell_tags, cells_to_save, connectivity_file, log_file);
  

  printf("Network ready\n");
  fprintf(log_file, "Network ready\n");
  
  connectivity_file.close();

  if (pars.weights_to_save.size()>0 && pars.weights_to_save[0]=="all"){
    weights_to_save=cell_tags;
  }
  else{
    weights_to_save=pars.weights_to_save;
  }
  FILE** weight_file_array=(FILE**)calloc(weights_to_save.size(), sizeof(FILE*));

  printf("Weights to save: ");
  fprintf(log_file, "Weights to save: ");
  for (neural_idx_t file_ind=0; file_ind<weights_to_save.size(); file_ind++){
    printf("%lu - %s; ", file_ind, weights_to_save[file_ind].c_str());
    fprintf(log_file, "%lu - %s; ", file_ind, weights_to_save[file_ind].c_str());
    string file_name=path+"weights_"+weights_to_save[file_ind];
    weight_file_array[file_ind]=fopen(file_name.c_str(), "w");
  }
  printf("\n");
  fprintf(log_file, "\n");

  // establish exact starting/ending point of the simulation
  timestamp_t start_time, end_time;
  start_time=pars.session_start_time+
        (static_cast<unsigned long long int>(pars.start)*60*1000000);
  end_time=start_time+
        (static_cast<unsigned long long int>(pars.duration)*60*1000000);
  printf("start_time %llu\n", start_time);
  fprintf(log_file, "start_time %llu\n", start_time);
  printf("end_time %llu\n", end_time);
  fprintf(log_file, "end_time %llu\n", end_time);
  printf("end_time within session time: %s\n\n", 
      (end_time<=pars.session_end_time)? "True" : "False");
  fprintf(log_file, "end_time within session time: %s\n\n", 
      (end_time<=pars.session_end_time)? "True" : "False");
  printf("Positional and parameters files loaded \n");
  fprintf(log_file, "Positional and parameters files loaded \n");

  
  // save initial weights
  for (neural_idx_t pc_idx=0; pc_idx<cell_tags.size(); pc_idx++){
    initial_weights_file<<cell_tags[pc_idx]<<" <--- ";
    vector<tag_t> afferent_tags; 
    vector<double> afferent_weights;
    network.get_labeled_synweights_for(cell_tags[pc_idx], 
                                afferent_tags, afferent_weights);
    for (neural_idx_t ag_id=0; ag_id<afferent_tags.size(); ag_id++){
      initial_weights_file<<"["<<afferent_tags[ag_id]<<" "<<afferent_weights[ag_id]<<"] ";
    }
    initial_weights_file<<endl;
  }
  initial_weights_file.close();

  // actual simulation
  simulate(timestamps, posxs, posys, start_time, end_time, grid_cells, network,
      grids_to_save, cells_to_save, grid_tags, cell_tags,
      weights_to_save, spike_file, log_file, weight_file_array);
    
  // save final weights
  for (neural_idx_t pc_idx=0; pc_idx<cell_tags.size(); pc_idx++){
    final_weights_file<<cell_tags[pc_idx]<<" <--- ";
    vector<tag_t> afferent_tags; 
    vector<double> afferent_weights;
    network.get_labeled_synweights_for(cell_tags[pc_idx], 
                                afferent_tags, afferent_weights);
    for (neural_idx_t ag_id=0; ag_id<afferent_tags.size(); ag_id++){
      final_weights_file<<"["<<afferent_tags[ag_id]<<" "<<afferent_weights[ag_id]<<"] ";
    }
    final_weights_file<<endl;
  }
  final_weights_file.close();
  
  printf("Done.\n");  
  fprintf(log_file, "Done.\n");  
  
  // releasing memory, closing files
  fclose(spike_file);
  fclose(log_file);
  for (neural_idx_t file_ind=0; file_ind<weights_to_save.size(); file_ind++){
    fclose(weight_file_array[file_ind]);
  }
  free(weight_file_array);

  /* legacy GSL implementation
  gsl_rng_free(connect_rng);
  gsl_rng_free(spike_rng);
  */

  
  // dummy file created to signal that the simulation is complete 
  string done_signal_file_name=path+"DONE";
  ofstream done_signal_file(done_signal_file_name.c_str());
  done_signal_file.close();

  return 0;
  
}


