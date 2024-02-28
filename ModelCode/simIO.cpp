/* 
MIT License Copyright 2024 Francesco Savelli - see LICENSE file.

Code for reading and handling the simulation parameters provided in the file
contained in the simulation folder.

*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstring>
#include <string>
#include <algorithm>

#include "util.hpp"
#include "simIO.hpp"

using namespace std;

plasticity_kind_t str2plasticity(string from){

	if (from=="off")
		return plasticity_off;
	else if (from=="presynaptically_gated")
		return presynaptically_gated;
	else if (from=="postsynaptically_gated")
		return postsynaptically_gated;
	else{
		cerr<<endl<<"string of plasticity kind not recognized: "<<from<<endl;
		return plasticity_off;
	}
}

string plasticity2str(plasticity_kind_t from){

	if (from==plasticity_off)
		return "off";
	else if (from==presynaptically_gated)
		return "presynaptically_gated";
	else if (from==postsynaptically_gated)
		return "postsynaptically_gated";
	else{
		cerr<<endl<<" plasticity kind not recognized"<<endl;
		return "not recognized";
	}
}


void read_pos_file(const string& namefile, 
                 vector<timestamp_t>& timestamps,
                 vector<double>& posxs, vector<double>& posys){

  vector<string> lines = readlines(namefile);

  char* remainder;

  timestamps.clear();
  posxs.clear();
  posys.clear();

  for (unsigned long int frame_ind=0; frame_ind<lines.size(); frame_ind++){

    if (lines[frame_ind][0]!='%'){

      vector<string> params = splitstring(lines[frame_ind], ',');
      timestamp_t timestamp=strtoull(params[0].c_str(), &remainder, 0);
      double posx=strtod(params[1].c_str(), &remainder);
      double posy=strtod(params[2].c_str(), &remainder);
      timestamps.push_back(timestamp);
      posxs.push_back(posx);
      posys.push_back(posy);
    }
  }

  // transform the image coordinates into euclidean coordinates
  // in the generic frame of reference with origin (min_x, max_y)
  // of the original image frame

  vector<double> posxs_wo_nans(posxs.size());
  auto itx = copy_if(posxs.begin(), posxs.end(), posxs_wo_nans.begin(), 
      [](double x){return !(std::isnan(x));});
  posxs_wo_nans.resize(std::distance(posxs_wo_nans.begin(),itx));

  vector<double> posys_wo_nans(posys.size());
  auto ity = copy_if(posys.begin(), posys.end(), posys_wo_nans.begin(), 
      [](double x){return !(std::isnan(x));});
  posys_wo_nans.resize(std::distance(posys_wo_nans.begin(),ity));

  double min_x= *(min_element(posxs_wo_nans.begin(), posxs_wo_nans.end()));
  double max_y= *(max_element(posys_wo_nans.begin(), posys_wo_nans.end()));

  for (unsigned long int ind=0; ind<posxs.size(); ind++)
    posxs[ind]=posxs[ind]-min_x;
  for (unsigned long int ind=0; ind<posys.size(); ind++)
    posys[ind]=max_y-posys[ind];
}


void parameters_record_t::from_file(const string& namefile){

  
  vector<string> lines = readlines(namefile);
  
 	parameters_record_t pars;
  unsigned int pars_count=0;
  
 	for (unsigned int lidx=0; lidx<lines.size(); lidx++){

    vector<string> words = splitstring(lines[lidx], '=');

		if (words[0]== "pos_file_relpath"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			pos_file_relpath=words[1];
			pars_count++;
   	}
		if (words[0]== "session_start_time"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			session_start_time=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
		if (words[0]== "session_end_time"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			session_end_time=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
		if (words[0]== "start"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			start=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
		if (words[0]== "duration"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			duration=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
    if (words[0]== "phase2"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
      if (words[1]== "true")
        phase2=true;
			pars_count++;
   	}
		if (words[0]== "camera_resolution"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			camera_resolution=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
		if (words[0]== "grid_smallest_scale"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_smallest_scale=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
		if (words[0]== "grid_largest_scale"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_largest_scale=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "weights_to_save"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
      vector<string> subwords = splitstring(words[1], ' ');
      for (int m=0; m<subwords.size(); m++){
			  weights_to_save.push_back(subwords[m]);
      }
			pars_count++;
   	}
	if (words[0]== "module_orientations"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
      vector<string> subwords = splitstring(words[1], ' ');
      for (int m=0; m<subwords.size(); m++){
			  module_orientations.push_back(
          strtod(subwords[m].c_str(), &remainder));
      }
			pars_count++;
   	}    
    if (words[0]== "module_stretch_directions"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
      vector<string> subwords = splitstring(words[1], ' ');
      for (int m=0; m<subwords.size(); m++){
			  module_stretch_directions.push_back(
          strtod(subwords[m].c_str(), &remainder));
      }
			pars_count++;
   	}
    if (words[0]== "module_stretch_magnitudes"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
      vector<string> subwords = splitstring(words[1], ' ');
      for (int m=0; m<subwords.size(); m++){
			  module_stretch_magnitudes.push_back(
          strtod(subwords[m].c_str(), &remainder));
      }
			pars_count++;
   	}
    if (words[0]== "module_shear_directions"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
      vector<string> subwords = splitstring(words[1], ' ');
      for (int m=0; m<subwords.size(); m++){
			  module_shear_directions.push_back(
          strtod(subwords[m].c_str(), &remainder));
      }
			pars_count++;
   	}
    if (words[0]== "module_shear_magnitudes"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
      vector<string> subwords = splitstring(words[1], ' ');
      for (int m=0; m<subwords.size(); m++){
			  module_shear_magnitudes.push_back(
          strtod(subwords[m].c_str(), &remainder));
      }
			pars_count++;
   	}
    if (words[0]== "grid_scale_steps"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_scale_steps=strtoul(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
    if (words[0]== "modules_per_place_cell"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			modules_per_place_cell=strtoul(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
    if (words[0]== "module_scale_ratio"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			module_scale_ratio=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "grid_stretch"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_stretch=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "grid_shear"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_shear=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
		if (words[0]== "grid_peak_rate"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_peak_rate=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "grid_min_fraction_peak"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_min_fraction_peak=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "grid_rise"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_rise=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "grid_cells_number"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grid_cells_number=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
    if (words[0]== "place_cells_number"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			place_cells_number=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
    if (words[0]== "grids_per_cell_number"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			grids_per_cell_number=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
		if (words[0]== "gDelta"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			gDelta=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    	if (words[0]== "gtau"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			gtau=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
		if (words[0]== "plasticity_rule"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			plasticity_rule=str2plasticity(words[1]);
			pars_count++;
   	}
    if (words[0]== "IPI"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			IPI=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
   		 if (words[0]== "plasticity_rate"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			plasticity_rate=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
		if (words[0]== "theta_p"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			theta_p=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
		if (words[0]== "theta_d"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			theta_d=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
 		if (words[0]== "cut_off"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			cut_off=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "fr_tau"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			fr_tau=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "prefr_tau"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			prefr_tau=strtod(words[1].c_str(), &remainder);
			pars_count++;
   	}
    if (words[0]== "spiking_rng_seed"){
  		//cerr<<words[0]<<" "<<words[1]<<endl;
			char* remainder;
			spiking_rng_seed=strtoull(words[1].c_str(), &remainder, 0);
			pars_count++;
   	}
	}

  /*
	if (pars_count<how_many_parameters){
		cerr<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
		cerr<<"Fewer parameters read than expected!"<<endl;
		cerr<<"Only "<<pars_count<<" read."<<endl;
		cerr<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	}
  */
}

void print_params(const parameters_record_t& pars){

  printf("pos_file_relpath %s\n", pars.pos_file_relpath.c_str());
  printf("session_start_time %llu\n", pars.session_start_time);
  printf("session_end_time %llu\n", pars.session_end_time);
  printf("start %u\n", pars.start);
  printf("duration %u\n", pars.duration);
  printf("camera_resolution %f\n", pars.camera_resolution);
  for (int cidx=0; cidx<pars.weights_to_save.size(); cidx++)
    printf("weights_to_save \"%s\"\n", pars.weights_to_save[cidx].c_str());
  printf("grid_smallest_scale %f\n", pars.grid_smallest_scale);
  printf("grid_largest_scale %f\n", pars.grid_largest_scale);
  printf("grid_scale_steps %u\n", pars.grid_scale_steps);
  printf("module_scale_ratio %f\n", pars.module_scale_ratio);
  for (int midx=0; midx<pars.module_orientations.size(); midx++)
    printf("module_orientations %u - %f\n", midx, pars.module_orientations[midx]);
  for (int midx=0; midx<pars.module_stretch_directions.size(); midx++)
    printf("module_stretch_directions %u - %f\n", midx, 
                                                  pars.module_stretch_directions[midx]);
  for (int midx=0; midx<pars.module_stretch_magnitudes.size(); midx++)
    printf("module_stretch_magnitudes %u - %f\n", midx, 
                                                  pars.module_stretch_magnitudes[midx]);
  for (int midx=0; midx<pars.module_shear_directions.size(); midx++)
    printf("module_shear_directions %u - %f\n", midx, 
                                                  pars.module_shear_directions[midx]);
  for (int midx=0; midx<pars.module_shear_magnitudes.size(); midx++)
    printf("module_shear_magnitudes %u - %f\n", midx, 
                                                  pars.module_shear_magnitudes[midx]);
  printf("grid_stretch %f\n", pars.grid_stretch);
  printf("grid_shear %f\n", pars.grid_shear);
  printf("grid_peak_rate %f\n", pars.grid_peak_rate);
  printf("grid_min_fraction_peak %f\n", pars.grid_min_fraction_peak);
  printf("grid_rise %f\n", pars.grid_rise);
  printf("grid_cells_number %u\n", pars.grid_cells_number);
  printf("place_cells_number %u\n", pars.place_cells_number);
  printf("grids_per_cell_number %u\n", pars.grids_per_cell_number);
  printf("modules_per_place_cell %u\n", pars.modules_per_place_cell);
  printf("gDelta %f\n", pars.gDelta);
  printf("gtau %f\n", pars.gtau);
  printf("plasticity_rule %s\n", (plasticity2str(pars.plasticity_rule)).c_str());
  printf("IPI %u\n", pars.IPI);
  printf("plasticity_rate %f\n", pars.plasticity_rate);
  printf("theta_p %f\n", pars.theta_p);
  printf("theta_d %f\n", pars.theta_d);
  printf("cut_off %f\n", pars.cut_off);
  printf("fr_tau %f\n", pars.fr_tau);
  printf("prefr_tau %f\n", pars.prefr_tau);
  printf("spiking_rng_seed %lu\n", pars.spiking_rng_seed);
  
}
