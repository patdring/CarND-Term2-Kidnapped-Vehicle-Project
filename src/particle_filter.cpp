/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "map.h" 

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	cout << "[0] Initialization" << endl;
  
	num_particles = 100;   
	default_random_engine gen;
    
	// Create normal (Gaussian) distributions for x,y, theta with  std[[standard deviation of x [m], 
	// 																	standard deviation of y [m], 
	// 																	standard deviation of yaw [rad]]
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
    
    // Sample and add them from normal distrubtions 
	for (int i = 0; i < num_particles; i++) {
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;
		// Print particle to the terminal.
		cout << "[0] Particle" << particle.id << " " << particle.x << " " << particle.y << " " << particle.theta << " " << particle.weight << endl;
        particles.push_back(particle);
	}
	is_initialized = true;    
	return;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine
  
	cout << "[1] Prediction" << endl;
    default_random_engine gen;
    
    for(Particle particle : particles) {
		particle.x = particle.x + velocity/yaw_rate * abs( sin(particle.theta+(yaw_rate*delta_t)) - sin(particle.theta) );
		particle.y = particle.y + velocity/yaw_rate * abs( cos(particle.theta) - cos(particle.theta+(yaw_rate*delta_t)) );
		particle.theta = particle.theta + yaw_rate*delta_t;
      
		normal_distribution<double> dist_x(particle.x, std_pos[0]);
		normal_distribution<double> dist_y(particle.y, std_pos[1]);
		normal_distribution<double> dist_theta(particle.theta, std_pos[2]);
      
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen); 
        cout << "[1] Particle" << particle.id << " " << particle.x << " " << particle.y << " " << particle.theta << " " << particle.weight << endl;
    }                                                                       
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    cout << "[2] Update Weights" << endl;
    // Transform Observations to map coord.
    std::vector<LandmarkObs> transf_observations;
  
    cout << "[2.1] Transform Observations" << endl;
  
    for( Particle particle : particles) {
        for ( LandmarkObs obs : observations) {
            LandmarkObs tobs;
    		tobs.x = particle.x + cos(particle.theta * obs.x) - sin(particle.theta * obs.y);
        	tobs.y = particle.y + sin(particle.theta * obs.x) + cos(particle.theta * obs.y); 
            transf_observations.push_back(tobs);
        }
          
        // calc. dist between tobs and map_landmarks, founded min. -> set id of landmark as tobs id
        cout << "[2.2] Find closest landmark" << endl;              
        for ( LandmarkObs tobs : transf_observations) {   
            // landmark sensor range as inital value for min. distance
            double min_dist = sensor_range;
            int min_id = -1;
          
        	for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
            	double d = dist(tobs.x, tobs.y, (double)map_landmarks.landmark_list[i].x_f, \
                                   (double)map_landmarks.landmark_list[i].y_f);
            	
                // reresh id if a new min. distance is detected
              	if (d < min_dist) {
                	min_dist = d;
                  	min_id = map_landmarks.landmark_list[i].id_i;
                }
            } 
          	tobs.id = min_id;
        }
     
        // calc. weight value
        cout << "[2.3] Calculate particle weight value" << endl;
      	// calculate normalization term
		double gauss_norm= (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
      
        for ( LandmarkObs tobs : transf_observations) {           
          	if (tobs.id == -1) {
            	continue;
            }
          
            double mu_x = 0.0;
          	double mu_y = 0.0;
            
            for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
            	if (tobs.id == map_landmarks.landmark_list[i].id_i) {
                    mu_x = (double)map_landmarks.landmark_list[i].x_f;
                    mu_y = (double)map_landmarks.landmark_list[i].y_f;
                	break;
                }
            }
          
			// calculate exponent
			double exponent= ((tobs.x - mu_x)*(tobs.x - mu_x))/(2 * std_landmark[0]*std_landmark[0] ) \
          						+ ((tobs.y - mu_y)*(tobs.y - mu_y))/(2 * std_landmark[1]*std_landmark[1]);
			// calculate weight using normalization terms and exponent
			particle.weight *= gauss_norm * exp(-exponent);
        }
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
  	discrete_distribution<int> dist_i(0, num_particles);
    int index = dist_i(gen);
    double beta = 0.0;
    double mv = 0.0;
  
    for( Particle particle : particles) {
    	if (particle.weight > mv) {
        	mv = particle.weight;
        }
    }
  
    std::vector<Particle> new_particles;
  
    for (int i = 0; i < particles.size(); i++) {
        normal_distribution<double> dist_b(0, 2*mv);
        beta +=  dist_b(gen);
      	while (beta > particles[index].weight) {
          	beta -= particles[index].weight;
          	index = (index + 1) % particles.size();
        }
      	new_particles.push_back(particles[index]);
    } 
  	// Copying vector by copy function 
    copy(new_particles.begin(), new_particles.end(), back_inserter(particles));
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
