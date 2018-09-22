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
  
	num_particles = 10;   
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
    
  	// Calculate new state.
	for (int i = 0; i < num_particles; i++) {
     
        // use different motion model when there is now yaw rate
        if (fabs(yaw_rate) < 0.00001) {  
      		particles[i].x += velocity * delta_t * cos(particles[i].theta);
      		particles[i].y += velocity * delta_t * sin(particles[i].theta);
    	} else {
     		particles[i].x += velocity/yaw_rate * (sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
      		particles[i].y += velocity/yaw_rate * (cos(particles[i].theta)-cos(particles[i].theta+yaw_rate * delta_t));
      		particles[i].theta += yaw_rate * delta_t;
        }
      
    	normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
  		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
  		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
      
    	particles[i].x = dist_x(gen);
    	particles[i].y = dist_y(gen);
    	particles[i].theta = dist_theta(gen);
      
    	cout << "[1] Particle" << particles[i].id << " " << particles[i].x << " " \
        	 << particles[i].y << " " << particles[i].theta << " " << particles[i].weight << endl;
  	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  	cout << "[2] Data Association" << endl;
	for (int i = 0; i < observations.size(); i++) {
    	double min_dist = numeric_limits<double>::max();
    	int min_id = -1;
    
    	for (unsigned int j = 0; j < predicted.size(); j++) {
      		double d = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

      		if (d < min_dist) {
        		min_dist = d;
        		min_id = predicted[j].id;
      		}
    	}	
    	observations[i].id = min_id;
      	cout << "[2.1] Associat obs. with id "<< min_id << endl;
  	}   
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
  	cout << "[3] Update weights" << endl;
  
    // calculate normalization term (multivariant gaussian distribution)
	double gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
    double support_term1 = 2 * std_landmark[0]*std_landmark[0];
    double support_term2 = 2 * std_landmark[1]*std_landmark[1];
  
	for (int i = 0; i < num_particles; i++) {       
    	vector<LandmarkObs> landmarks;
        // extract landmark positions and ids and store them in a new vector obs
    	for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {	  
      		LandmarkObs landmark; // temporary placeholder
      
      		landmark.x = (double) map_landmarks.landmark_list[j].x_f;
      		landmark.y = (double) map_landmarks.landmark_list[j].y_f;
      		landmark.id = map_landmarks.landmark_list[j].id_i;
      
            if (dist(landmark.x, landmark.y, particles[i].x, particles[i].y) <= sensor_range) {
      			landmarks.push_back(landmark);
            }            
    	}

    	vector<LandmarkObs> tobs;
        // transform positions from car to map coord. and store them in a new vector tobs
    	for (int j = 0; j < observations.size(); j++) {
      		LandmarkObs lobs; // temporary placeholder
      		lobs.x = particles[i].x + cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y;
      		lobs.y = particles[i].y + sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;
      		tobs.push_back(lobs);
    	}

      	// find associated landmark for given observations
    	dataAssociation(landmarks, tobs);
    	// reset weight
    	particles[i].weight = 1.0;
  
    	for (unsigned int j = 0; j < tobs.size(); j++) {
            double mu_x = 0.0;
      		double mu_y = 0.0;

            // do following only if there is a associated landmarks
          	if (tobs[j].id != -1) {
              	// extract landmark positions for give (transf.) observation
      			for (int k = 0; k < landmarks.size(); k++) {
      				if (landmarks[k].id == tobs[j].id) {
          				mu_x = landmarks[k].x;
          				mu_y = landmarks[k].y;
        			}
      			}
            }
          
        	// calculate exponent (multivariant gaussian distribution)
			double exponent = ((tobs[j].x - mu_x)*(tobs[j].x - mu_x))/support_term1 \
          					+ ((tobs[j].y - mu_y)*(tobs[j].y - mu_y))/support_term2;
      
          	double w = gauss_norm * exp(-exponent);
      		// calculate/update weight using normalization terms and exponent
            if (w != 0) {
      			particles[i].weight *= w;
            } else {
              	particles[i].weight *= 0.00001; // update with very low weight/probability if calc. weight is 0
            }
    	}
  	}  
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution  
    default_random_engine gen;
    cout << "[3] Resample values" << endl;
	vector<Particle> new_particles;
    double mw = 0.0;
    
    for (int i = 0; i < num_particles; i++) {
    	if(particles[i].weight > mw) {
          	mw = particles[i].weight;
        }
    }

    // generate random starting index
    uniform_int_distribution<int> uniintdist(0, num_particles-1);
    int index = uniintdist(gen);
  
    normal_distribution<double> dist_b(0, mw);
    double beta = 0.0;

    // resample wheel
    for (int i = 0; i < num_particles; i++) {
        beta += dist_b(gen) * 2.0;
      
        while (beta > particles[index].weight) {
            beta -= particles[index].weight;
            index = (index + 1) % num_particles;
        }
        new_particles.push_back(particles[index]);
    }
    particles = new_particles;  
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
