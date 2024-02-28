#ifndef GRID_INPUT
#define GRID_INPUT

#include <iostream>
#include <cmath>

/* Older implementation depended on the Gnu Scientific Library for random
 * number generation. Current implementation uses equivalent rng facilities now
 * available in the C++ standard library. Older GSL calls are commented out for
 * record. New implementation produces the same sequences of random numbers as
 * the older one would.

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
*/
#include <random>

#include "network.hpp"

#define PI 3.14159265

#define RESTING_ACTIVITY 0.0

// Geometric definitions and transformation are w.r.t. a conventionally 
// oriented xy coordinate system. Trajectory data was transformed from 
// the camera frame coordinate system with y pointing downward at the 
// time it was read in.

struct sq_dist_cycles_t{
  double sq_dist; 
  int cycles1; 
  int cycles2; 
  int cycles3;
  double nonround_cycles1; 
  double nonround_cycles2; 
  double nonround_cycles3;
};

struct peaks_discounts_t{
  
  unsigned int rows, columns;
  vector<vector<double>> mat;

  peaks_discounts_t(unsigned int rows, unsigned int columns, double min_fraction_peak,
                    mt19937* grid_rng) : rows(rows), columns(columns) {

    if (min_fraction_peak<1.0){
      uniform_real_distribution<float> peaks_distribution (min_fraction_peak, 1.0);
      for (unsigned int r=0; r<rows; r++){
        vector<double> row; 
        for (unsigned int c=0; c<columns; c++)
          row.push_back(peaks_distribution(*grid_rng));
        mat.push_back(row);
      }
    }
    else{
      for (unsigned int r=0; r<rows; r++){
        vector<double> row; 
        for (unsigned int c=0; c<columns; c++)
          row.push_back(1.0);
        mat.push_back(row);
      }
    }
  }

};

class grid_t{

  public:

    grid_t(double x, double y, double orientation, double scale, double rise, 
                      double stretch_direction, double stretch_magnitude, 
                      double shear_direction, double shear_magnitude,
                      peaks_discounts_t peaks_discounts);


    inline double activation_at(double x, double y){

      if (isnan(x)||isnan(y)){ // invalid position
        return 0.0; // zero activation, makes it impossible to spike
      }
      else{

        if (stretch_magnitude!=1.0){
          double angle = stretch_direction * PI/180;
          rotate(angle, x, y);
          stretchx(stretch_magnitude, x, y);
          rotate(-angle, x, y);
        }

        if (shear_magnitude!=0.0){
          double angle = shear_direction * PI/180;
          rotate(angle, x, y);
          shearx(shear_magnitude, x, y);
          rotate(-angle, x, y);
        }

        sq_dist_cycles_t res = sq_dist_from_closest_vertex(x, y);
        // determine the lookup table row/col, based on vertex cycles, 
        // attempt to keep between bounds
        int row=(peaks_discounts.rows/2+res.cycles1) % peaks_discounts.rows; 
        int col=(peaks_discounts.columns/2+res.cycles2) % peaks_discounts.columns; 
        double peak_discount = peaks_discounts.mat[row][col];
        return exp(- res.sq_dist / (rise*scale*scale)) * peak_discount;
      }
    }

    inline sq_dist_cycles_t sq_dist_cycles_at(double x, double y){

      // primarily thougth to be used for trimming
      // random phases into an hexagon
      // doesn't handle nans as activation_at 

      if (stretch_magnitude!=1.0){
        double angle = stretch_direction * PI/180;
        rotate(angle, x, y);
        stretchx(stretch_magnitude, x, y);
        rotate(-angle, x, y);
      }

      if (shear_magnitude!=0.0){
        double angle = shear_direction * PI/180;
        rotate(angle, x, y);
        shearx(shear_magnitude, x, y);
        rotate(-angle, x, y);
      }

      sq_dist_cycles_t res = sq_dist_from_closest_vertex(x, y);

      return res;

    }

  private:

    double m1, q1, delta_q1;
    double m2, q2, delta_q2;
    double m3, q3, delta_q3;

    double x, y, orientation, scale, rise;
    double stretch_direction, stretch_magnitude;
    double shear_direction, shear_magnitude;

    peaks_discounts_t peaks_discounts;

    inline void rotate(double angle, double& x, double& y){
      double tx, ty;
      tx = x*cos(angle) - y*sin(angle);
      ty = x*sin(angle) + y*cos(angle);
      x=tx; y=ty;
    }

    inline void stretchx(double scale, double& x, double& y){
      double tx, ty;
      tx = scale*x;
      ty = y;
      x=tx; y=ty;
    }

    inline void shearx(double shear, double& x, double& y){
      double tx, ty;
      tx = x + shear*y;
      ty = y;
      x=tx; y=ty;
    }

    inline sq_dist_cycles_t sq_dist_from_closest_vertex(double xp, double yp);
  
};

class grid_cell_t{

  
  public:

    grid_cell_t(double x, double y, double orientation, double scale, double rise,
        double stretch_direction, double stretch_magnitude,
        double shear_direction, double shear_magnitude,
        double peak_rate, peaks_discounts_t peaks_discounts, 
        mt19937* spike_rng, tag_t tag) :
                                
                                  grid(x, y, orientation, scale, rise, 
                                      stretch_direction, stretch_magnitude, 
                                      shear_direction, shear_magnitude, 
                                      peaks_discounts){
       
      
      this->tag=tag;
      this->peak_rate=peak_rate;

      is_spiking=false;
      is_mute=false;

      this->spike_rng=spike_rng;

      URD=uniform_real_distribution<float>(0.0, 1.0);
      ED=exponential_distribution<float>(peak_rate);
      generate_next_tentative_ISI();
    }

    inline void reset_peak_rate(double new_peak_rate){
      peak_rate=new_peak_rate;
      ED=exponential_distribution<float>(peak_rate);
    }

    inline void set_mute(bool flag){is_mute=flag;}

    inline bool just_spiked(){return !is_mute && is_spiking;}

    inline void step(double posx, double posy){

      steps_to_next_tentative_spike -= 1;
        
      if(steps_to_next_tentative_spike==0){
        // thin the homogeneous Poisson process (at peak rate) based on the 
        // grid firing rate an the current position and draw the next ISI of 
        // the homogeneous process
        
        /* legacy GSL implementation
        if (gsl_rng_uniform(spike_rng)<=grid.activation_at(posx, posy))
        */
        if (URD(*spike_rng)<=grid.activation_at(posx, posy))
          is_spiking=true;// spike!
        else
          is_spiking=false; // this potential spike is thinned off
        // regardless, it is time to restart a new ISI
        generate_next_tentative_ISI();
      }
      else // still waiting for this ISI to be over
        is_spiking=false;
    }

    tag_t get_tag() {return tag;}

  private:

    tag_t tag;
    grid_t grid;
    
    double peak_rate;

    delay_t steps_to_next_tentative_spike;
    bool is_spiking;
    bool is_mute;
    
    mt19937* spike_rng;
    uniform_real_distribution<float> URD;
    exponential_distribution<float> ED;
   
    inline void generate_next_tentative_ISI(){
      // homogeneous Poisson process based on max firing rate
      /* legacy GSL implementation
      steps_to_next_tentative_spike= 
            int(round(gsl_ran_exponential(spike_rng, 1.0/peak_rate)*1000)); // result in ms
      */
      steps_to_next_tentative_spike= 
            int(round(ED(*spike_rng)*1000)); // result in ms
      if (steps_to_next_tentative_spike<3) 
        // refractory period: never spike twice inside of 3ms
        steps_to_next_tentative_spike=3;
    }

};


// definitions

grid_t::grid_t(double x, double y, double orientation, double scale, double rise, 
                            double stretch_direction, double stretch_magnitude,
                            double shear_direction, double shear_magnitude,
                            peaks_discounts_t peaks_discounts) :

          x(x), y(y), orientation(orientation), scale(scale), rise(rise),
          stretch_direction(stretch_direction), stretch_magnitude(stretch_magnitude),
          shear_direction(shear_direction), shear_magnitude(shear_magnitude),
          peaks_discounts(peaks_discounts) {

  // from degrees to radians
  double headings[3]={orientation*PI/180, 
			(orientation+60.0)*PI/180,
			(orientation+120.0)*PI/180};
  
  double abs_cos_headings[3]={fabs(cos(headings[0])), 
				fabs(cos(headings[1])), 
				fabs(cos(headings[2]))};
  
  unsigned int best_h1, best_h2, worst_h;
  if (abs_cos_headings[0]<=abs_cos_headings[1] &&
      abs_cos_headings[0]<=abs_cos_headings[2]){
    best_h1=1; 
    best_h2=2; 
    worst_h=0;
  }
  else if (abs_cos_headings[1]<=abs_cos_headings[0] &&
	   abs_cos_headings[1]<=abs_cos_headings[2]){
    best_h1=0; 
    best_h2=2;
    worst_h=1;
  }
  else{
    best_h1=0;
    best_h2=1; 
    worst_h=2;
  }

  m1=tan(headings[best_h1]);
  m2=tan(headings[best_h2]);
  // the worst heading is accounted as x=my+q instead of the usual y=mx+q
  headings[worst_h]=PI/2 - headings[worst_h];
  m3=tan(headings[worst_h]);

  q1= y - m1*x;
  q2= y - m2*x;
  // the worst heading... 
  q3= x - m3*y;
  
  double delta_x1, delta_y1, x1, y1;
  double delta_x2, delta_y2, x2, y2;
  double delta_x3, delta_y3, x3, y3;
  
  delta_x1 = scale*cos(headings[best_h1]);
  delta_y1 = scale*sin(headings[best_h1]);
  
  delta_x2 = scale*cos(headings[best_h2]);
  delta_y2 = scale*sin(headings[best_h2]);
  
  // ...again...
  delta_y3 = scale*cos(headings[worst_h]);
  delta_x3 = scale*sin(headings[worst_h]);

  x1 = x+delta_x1;
  y1 = y+delta_y1;
  
  x2 = x+delta_x2; 
  y2 = y+delta_y2;
  
  x3 = x+delta_x3; 
  y3 = y+delta_y3;
  
  delta_q1 = fabs(q1-(y2-m1*x2));
  delta_q2 = fabs(q2-(y1-m2*x1));
  // ... and again ...
  delta_q3= fabs(q3-(x2-m3*y2));
  
}

inline sq_dist_cycles_t grid_t::sq_dist_from_closest_vertex(double xp, double yp){
		
  double qp1= yp-m1*xp;
  double qp2= yp-m2*xp;
  // this line is in the form x=m3*y+q3
  double qp3= xp-m3*yp;
 
  double nonround_cycles1= (qp1-q1)/delta_q1;
  double nonround_cycles2= (qp2-q2)/delta_q2;
  double nonround_cycles3= (qp3-q3)/delta_q3;

  double cycles1= round(nonround_cycles1);
  double cycles2= round(nonround_cycles2);
  double cycles3= round(nonround_cycles3);

  double qg1= q1 + delta_q1 * cycles1;
  double qg2= q2 + delta_q2 * cycles2;
  double qg3= q3 + delta_q3 * cycles3;
  
  //intersect line 1 and 2
  double xg_12= (qg2-qg1)/(m1-m2);
  double yg_12= m1*xg_12+qg1;
  
  //intersect line 1 and 3
  double yg_13= (m1*qg3+qg1)/(1-m1*m3);
  double xg_13= m3*yg_13+qg3;
  
  //intersect line 2 and 3
  double yg_23= (m2*qg3+qg2)/(1-m2*m3);
  double xg_23= m3*yg_23+qg3;
  
  double sq_dist_12= (xp-xg_12)*(xp-xg_12) + (yp-yg_12)*(yp-yg_12);
  double sq_dist_13= (xp-xg_13)*(xp-xg_13) + (yp-yg_13)*(yp-yg_13);
  double sq_dist_23= (xp-xg_23)*(xp-xg_23) + (yp-yg_23)*(yp-yg_23);
  
  double sq_dist;
  if (sq_dist_12<=sq_dist_13 and sq_dist_12<=sq_dist_23)
    sq_dist=sq_dist_12;
  else if (sq_dist_13<=sq_dist_12 and sq_dist_13<=sq_dist_23)
    sq_dist=sq_dist_13;
  else
    sq_dist=sq_dist_23;

  sq_dist_cycles_t res;
  res.sq_dist=sq_dist; 
  res.cycles1=cycles1; 
  res.cycles2=cycles2; 
  res.cycles3=cycles3;
  res.nonround_cycles1=nonround_cycles1; 
  res.nonround_cycles2=nonround_cycles2; 
  res.nonround_cycles3=nonround_cycles3;

  return res;
  
}

#endif
