/* 
MIT License Copyright 2024 Francesco Savelli - see LICENSE file.

Basic types of spiking units (input and integrate-and-fire) that can be combined
into a neural network. External code will use these classes to build the network
to be used in the simulation. Each component has a ".step()" method performing
the relevant operations at each simulation timestep.

*/

#ifndef NETWORK_HPP
#define NETWORK_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

enum plasticity_kind_t {plasticity_off, presynaptically_gated, postsynaptically_gated};

typedef unsigned long int neural_idx_t;
typedef unsigned int delay_t;
typedef std::string tag_t;

class input_t{
  
  public:
    // general constructor for this type
    input_t(double fr_tau, tag_t tag){
      this->tag=tag;
      firing_rate=0.0;
      this->fr_tau=fr_tau;
      is_spiking=false;
    }
  
    // what to do at any step of the overall simulation for this type
    inline void step(bool make_spike){
      // decay factor with tau=100 ms
      firing_rate= firing_rate*exp(-1.0/fr_tau); 
      is_spiking=make_spike;
      if (make_spike==true)
        firing_rate+= (1.0/fr_tau); // 1/tau with tau=100 ms, max (1.0) corresponds to 1000Hz
    }
    
    inline bool just_spiked(){return is_spiking;}
    inline double get_firing_rate(){return firing_rate;}
    inline tag_t get_tag(){return tag;}

  private:
    tag_t tag;
    double firing_rate;
    double fr_tau;
    bool is_spiking;   
    
};


class cell_t{

  public:
    // general constructor for this type
    cell_t(double Cm, double Vthr, double Vres, double Vap, 
           double initial_V, double timestep, 
           delay_t refractory_steps, double fr_tau, tag_t tag) : 

           Cm(Cm), Vthr(Vthr), Vres(Vres), Vap(Vap),
           V(initial_V), delta_t(timestep), 
           refractory_steps(refractory_steps), fr_tau(fr_tau), tag(tag) {
          
      refractory_to_go=0;
      is_spiking=false;
      // leaky components
      sum_of_conductances=0.2;
      sum_of_currents=-65.0*0.2;

      firing_rate=0.0;
    }

    inline tag_t get_tag(){return tag;}
    inline bool just_spiked(){return is_spiking;}
    inline double get_firing_rate(){return firing_rate;}
    inline void add_conductance(double g){sum_of_conductances+=g;}
    inline void add_current(double I){sum_of_currents+=I;}
  
    // what to do at any step of the overall simulation for this type
    inline void step(double I);
    
  private:

    tag_t tag; 
    // Capacitance in nF, voltage in mV, current in nA, time in ms
    
    double Cm; 
    double Vthr; 
    double Vres;
    double Vap;
    double V;
    
    double delta_t;

    delay_t refractory_steps;
    delay_t refractory_to_go;

    double sum_of_conductances;
    double sum_of_currents;

    bool is_spiking;
    double firing_rate;
    double fr_tau;
  
};

inline void cell_t::step(double I=0.0){

  firing_rate= firing_rate*exp(-1.0/fr_tau);
  // decay factor with tau=100 ms

  if (refractory_to_go>0)
    refractory_to_go -= 1;

  if (is_spiking){
    V=Vres;
    is_spiking=false;
    refractory_to_go=refractory_steps;
  }
  else if (V>=Vthr && refractory_to_go==0){
    V=Vap;
    is_spiking=true;
    // 1/tau with tau=100 ms, max (1.0) corresponds to 1000Hz
    firing_rate+= (1.0/fr_tau);
  }
  else{
    // sum_of_currents and sum_of_conductances are supposed to have 
    // been already updated (see synapse.step())
    double VBase=(sum_of_currents+I)/sum_of_conductances;
    double tau=Cm/sum_of_conductances; // results in msec

    V=VBase+(V-VBase)*exp(-(delta_t/tau)); // exponential Euler
    //V=V+(VBase-V)*(delta_t/tau);

    if (V<-100.0) V=-100.0;
    if (V>100.0) V=100.0;

  }

  // leaky components
  sum_of_conductances=0.2;
  sum_of_currents=-65.0*0.2;

}


class synapse_t{

  public:

    // General constructor for this type. 
    // Plasticity is independently set by other member functions.
    synapse_t(double gDelta, double gtau, double E, double timestep){

      this->gDelta=gDelta;
      this->E=E;
      g=0.0;
      gDecay=exp(-(timestep/gtau));
      plasticity_rule=plasticity_off;
      connections_are_set=false;
    }
    
    void set_connections(neural_idx_t presyn_idx, 
                         neural_idx_t postsyn_idx,
                         delay_t axonal_delay) {

      this->presyn_idx=presyn_idx;
      this->postsyn_idx=postsyn_idx;

      this->axonal_delay=axonal_delay;
      axonal_spike_buffer=vector<int>(axonal_delay,0);
      axonal_fr_buffer=vector<double>(axonal_delay,0.0);
      postsyn_just_spiked=false;

      connections_are_set=true;
    }
    
    void set_rate_plasticity(plasticity_kind_t plasticity_rule, 
                             double theta_d, double theta_p, double fr_cutoff,
                             double plasticity_rate, delay_t IPI){

      if (this->plasticity_rule!=plasticity_off){
        cout<<"WARNING: non-OFF plasticity rule changed to a new kind";
      }
      
      this->plasticity_rule=plasticity_rule;
      this->theta_d=theta_d;
      this->theta_p=theta_p;
      this->fr_cutoff=fr_cutoff;
      this->plasticity_rate=plasticity_rate;
      this->IPI=IPI;
      steps_to_next_syn_change=IPI;
    }

    void set_plasticity_off(){plasticity_rule=plasticity_off;} 
   
    neural_idx_t get_presyn_idx() {return presyn_idx;}
    neural_idx_t get_postsyn_idx() {return postsyn_idx;}
    double get_syn_weight() {return gDelta;}

    inline void presyn_cell_just_spiked(){axonal_spike_buffer[axonal_delay-1]=1;}
    inline void postsyn_cell_just_spiked(){postsyn_just_spiked=true;}
        
    // simulation step for this synapse
    inline void step(bool pre_just_spiked, double pre_fr, double post_fr,
                      double &post_conductance, double &post_current);

  private:

    // specific to synaptic transmission - set by constructor
    double g;
    double gDelta;
    double E;
    double gDecay;

    
    // specific to connectivity to pre and post cells in the network
    // set by set_connections(...)
    neural_idx_t presyn_idx, postsyn_idx;
    delay_t axonal_delay; //size of axonal spike and fr buffer
    vector<int> axonal_spike_buffer; // delays spikes 
    vector<double> axonal_fr_buffer; // delays presyn firing rate 
    bool postsyn_just_spiked;
    bool connections_are_set;

    // plasticity parameters - set by set_rate_plasticity(...)
    plasticity_kind_t plasticity_rule; 
    double theta_d, theta_p, fr_cutoff; // params of synaptic rule
    delay_t IPI; // how often the rule is applied
    double plasticity_rate; // learning rate, takes into account IPI
    delay_t steps_to_next_syn_change; // step counter to next rule application

};

inline void synapse_t::step(bool pre_just_spiked, double pre_fr, double post_fr,
                              double &post_conductance, double &post_current){

  if (connections_are_set) {
    // decay factor with tau=100 ms
    //g=gBase+(g-gBase)*gDecay;
    g=g*gDecay;
    // if a spike has arrived (delayed by axonal buffer)
    if (axonal_spike_buffer[0]==1) 
      g+=gDelta; // 1/tau with tau=100 ms, max (1.0) corresponds to 1000Hz
    // update
    post_conductance = g;
    post_current = g*E;
    // the presynaptic firing rate to use in the synaptic rule 
    // is also delayed by the buffer
    double delayed_pre_fr=axonal_fr_buffer[0];
    // shift the spike and firing rate axonal buffers
    for (delay_t axonal_step=1; axonal_step<axonal_delay; axonal_step++){
      axonal_spike_buffer[axonal_step-1]=axonal_spike_buffer[axonal_step];
      axonal_fr_buffer[axonal_step-1]=axonal_fr_buffer[axonal_step];
    }
    axonal_spike_buffer[axonal_delay-1]= pre_just_spiked ? 1 : 0;
    axonal_fr_buffer[axonal_delay-1]=pre_fr;

    //PLASTICITY
    
    if (plasticity_rule==presynaptically_gated)
      if (steps_to_next_syn_change==0){
        // synaptic plasticity rule
        double clipped_pre_fr=delayed_pre_fr-fr_cutoff;
        if (clipped_pre_fr<0.0)
          clipped_pre_fr=0.0;
        if (post_fr>theta_d){
          gDelta+=(post_fr-theta_p) * clipped_pre_fr * plasticity_rate;
        }

        if (gDelta>0.1)
          gDelta=0.1;
        if (gDelta<0.0)
          gDelta=0.0;

        steps_to_next_syn_change=IPI;
      }
      else
        steps_to_next_syn_change-=1;

    if (plasticity_rule==postsynaptically_gated)
      if (steps_to_next_syn_change==0){
        // synaptic plasticity rule
        double clipped_post_fr=post_fr-fr_cutoff;
        if (clipped_post_fr<0.0)
          clipped_post_fr=0.0;
        if (delayed_pre_fr>theta_d){
          gDelta+=(delayed_pre_fr-theta_p) * clipped_post_fr * plasticity_rate;  
          /*if (clipped_post_fr>0.0){
            cout<<"("<<get_presyn_idx()<<" "<<get_postsyn_idx()<<")";
            cout<<"; pre_fr "<<delayed_pre_fr<<"; post_fr "<<clipped_post_fr<<"; gDelta "<<gDelta;
            cout<<"; cutoff "<<fr_cutoff<<"; IPI "<<IPI<<"; theta_p "<<theta_p<<"; theta_d "<<theta_d<<endl; 
          }*/
        }
          
        if (gDelta>0.1)
          gDelta=0.1;
        if (gDelta<0.0)
          gDelta=0.0;

        steps_to_next_syn_change=IPI;
      }
      else{
        steps_to_next_syn_change-=1;
      }

    // assume no postsyn spike is about to occur, if not
    // postsyn_cell_just_spiked() will be externally invoked
    postsyn_just_spiked=false;

  }
}


class network_t{

  private:
      
    vector<input_t> inputs;
    vector<cell_t> cells;
    
    vector<synapse_t> inputs_to_cells;
    vector<synapse_t> cells_to_cells;

  public:

    inline void add_input(const input_t &new_input){
      inputs.push_back(new_input);
    }

    inline void add_cell(const cell_t &new_cell){
      cells.push_back(new_cell);
    }

    inline void add_input_to_cell(const synapse_t &new_synapse){
      inputs_to_cells.push_back(new_synapse);
    }

    inline void add_cell_to_cell(const synapse_t &new_synapse){
      cells_to_cells.push_back(new_synapse);
    }

    inline tag_t get_tag_of_cell_idx(const neural_idx_t &cellidx){
      return cells[cellidx].get_tag();
    }

    inline tag_t get_tag_of_input_idx(const neural_idx_t &inputdx){
      return inputs[inputdx].get_tag();
    }

    inline neural_idx_t get_idx_of_cell_tag(const tag_t &celltag){
      neural_idx_t cellidx=0;
      bool found=false; 
      while(cellidx<cells.size() && !found)
        if (cells[cellidx].get_tag()==celltag)
          found=true;
        else
          cellidx++;
      if (!found){
        cerr<<"Did not find cell in get_idx_of_celltag"<<celltag;
        exit(-1);
      }
      return cellidx;
    }
  
    inline neural_idx_t get_idx_of_input_tag(const tag_t &inputtag){
      neural_idx_t inputidx=0;
      bool found=false; 
      while(inputidx<inputs.size() && !found)
        if (inputs[inputidx].get_tag()==inputtag)
          found=true;
        else
          inputidx++;
      if (!found){
        cerr<<"Did not find cell in get_idx_of_inputtag"<<inputtag;
        exit(-1);
      }
      return inputidx;
    }

    inline void get_labeled_synweights_for(const tag_t &cell_tag,
                            vector<tag_t> &afferent_tags, 
                            vector<double> &afferent_weights){

      neural_idx_t cellidx=get_idx_of_cell_tag(cell_tag);
      for (vector<synapse_t>::iterator syn_it=inputs_to_cells.begin(); 
                                syn_it<inputs_to_cells.end(); syn_it++){
          if (syn_it->get_postsyn_idx()==cellidx){
            afferent_tags.push_back(get_tag_of_input_idx(syn_it->get_presyn_idx()));
            afferent_weights.push_back(syn_it->get_syn_weight());
          }
      }
      for (vector<synapse_t>::iterator syn_it=cells_to_cells.begin(); 
                                syn_it<cells_to_cells.end(); syn_it++){
          if (syn_it->get_postsyn_idx()==cellidx){
            afferent_tags.push_back(get_tag_of_input_idx(syn_it->get_presyn_idx()));
            afferent_weights.push_back(syn_it->get_syn_weight());
          }
      }
                            
    }        

    // simulation step for this synapse
    inline vector<bool> step(const vector<bool> &inputs_that_spiked, double I=0.0);

};


inline vector<bool> network_t::step(const vector<bool> &inputs_that_spiked, double I){

  // First, update inputs
  for (neural_idx_t i=0; i<inputs.size(); i++){
    inputs[i].step(inputs_that_spiked[i]);
  }

  // then update synapses
  vector<synapse_t>::iterator synapse_it;
  
  for (synapse_it=inputs_to_cells.begin(); 
       synapse_it<inputs_to_cells.end(); synapse_it++){
   neural_idx_t pre_idx=synapse_it->get_presyn_idx(); 
   neural_idx_t post_idx=synapse_it->get_postsyn_idx(); 
   double post_conductance=0.0;
   double post_current=0.0;
   synapse_it->step(inputs[pre_idx].just_spiked(),
                    inputs[pre_idx].get_firing_rate(),
                    cells[post_idx].get_firing_rate(),
                    post_conductance, post_current);
   cells[post_idx].add_conductance(post_conductance);
   cells[post_idx].add_current(post_current);
  }

  for (synapse_it=cells_to_cells.begin(); 
       synapse_it<cells_to_cells.end(); synapse_it++){
   neural_idx_t pre_idx=synapse_it->get_presyn_idx(); 
   neural_idx_t post_idx=synapse_it->get_postsyn_idx(); 
   double post_conductance=0.0;
   double post_current=0.0;
   synapse_it->step(cells[pre_idx].just_spiked(),
                    cells[pre_idx].get_firing_rate(),
                    cells[post_idx].get_firing_rate(),
                    post_conductance, post_current);
   cells[post_idx].add_conductance(post_conductance);
   cells[post_idx].add_current(post_current);
  } 
  
  // finally update cells and return if they spiked
  vector<bool> cells_that_spiked(cells.size(), false);
  for (neural_idx_t cidx=0; cidx<cells.size(); cidx++){ 
    cells[cidx].step(I);
    cells_that_spiked[cidx]=cells[cidx].just_spiked();
  }
  return cells_that_spiked;
}


#endif

