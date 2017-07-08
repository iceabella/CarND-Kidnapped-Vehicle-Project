/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *
 *  Updated on: June 26, 2017
 *  	Author: Isabella Johansson
 */

#include <random> // Need this for sampling from distributions
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	//  Sets the number of particles. Initializes all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 10; //100 // Number of particles, tunable parameter

	// initalize a random engine generator
	std::default_random_engine gen;

	// Normal (Gaussian) distribution for x, y and theta
	std::normal_distribution<double> distribution_x(x, std[0]);
	std::normal_distribution<double> distribution_y(y, std[1]);
	std::normal_distribution<double> distribution_theta(theta, std[2]);

	// Create num_particlese Gaussian distributed particles
	Particle p;
	for(int i=0; i<num_particles; i++){
		p.id = i;
		p.x = distribution_x(gen); 
		p.y = distribution_y(gen);
		p.theta = distribution_theta(gen);
		p.weight = 1.0; // uniform distribution
		particles.push_back(p);
	}
	// done initializing
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Adds measurements to each particle and random Gaussian noise.
	// NOTE: When adding noise std::normal_distribution and std::default_random_engine are useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// initalize a random engine generator
	std::default_random_engine gen;

	double xf,yf;
	double thetaf;

	// Normal (Gaussian) distribution for 0 mean noise
	std::normal_distribution<double> distribution_x(0, std_pos[0]);
	std::normal_distribution<double> distribution_y(0, std_pos[1]);
	std::normal_distribution<double> distribution_theta(0, std_pos[2]);

	for(auto& p:particles){
		// avoid division with zero (if zero yaw_rate)
		if(fabs(yaw_rate) == 0){
			xf = p.x + velocity*cos(p.theta)*delta_t;
			yf = p.y + velocity*sin(p.theta)*delta_t;      
        	} else{
			xf = p.x + velocity/yaw_rate*(sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
			yf = p.y + velocity/yaw_rate*(cos(p.theta) - cos(p.theta + yaw_rate*delta_t)); 
       		}
		
		thetaf = p.theta + yaw_rate*delta_t;

		// Save predicted state with Gaussian noise
		p.x = xf + distribution_x(gen); 
		p.y = yf + distribution_y(gen);
		p.theta = thetaf + distribution_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Finds the predicted measurement that is closest to each observed measurement and assigns the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// Updates the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. The particles are located
	//   according to the MAP'S coordinate system. Transformation is needed between the two systems.
	//   This transformation requires both rotation AND translation (but no scaling).
	//   Theory: https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the actual implemented equation (equation 3.33) http://planning.cs.uiuc.edu/node99.html

	std::vector<double> senseGlobal_x;
	std::vector<double> senseGlobal_y;
	std::vector<int> associated_ids;
	// reset weights vector
	weights.clear();
	for(auto& p:particles){
		// reset vectors
		senseGlobal_x.clear();
		senseGlobal_y.clear();
		associated_ids.clear();
		p.weight = 1.0;

		for(auto const& o:observations){
			// 1. Transform observation from vehicle coordinates to global (map) coordinates
			double x_observation = o.x*cos(p.theta) - o.y*sin(p.theta) + p.x;
			double y_observation = o.x*sin(p.theta) + o.y*cos(p.theta) + p.y; 

			senseGlobal_x.push_back(x_observation);
			senseGlobal_y.push_back(y_observation);

			double x_landmark;
			double y_landmark;
			float min_dist = -1.0;
			int id_min;
			float dist_p;
			// 2. Associate global observation to landmark
			for(auto const& landmark:map_landmarks.landmark_list){
				dist_p = dist(landmark.x_f, landmark.y_f, x_observation, y_observation);
			
				if(min_dist > dist_p || min_dist == -1.0){
					// save data for closest landmark to observation
					min_dist = dist_p;
					id_min = landmark.id_i;
					x_landmark = landmark.x_f;
					y_landmark = landmark.y_f; 
				}
			}
			// associate observation with found closest landmark
			associated_ids.push_back(id_min);	

			// 3. Update particle's weight update based on global observation and associated landmark (non normalized)
			double w_obsi = 1/(2*M_PI*std_landmark[0]*std_landmark[1])*exp(-0.5* 
				    ((x_observation - x_landmark)*(x_observation - x_landmark)/(std_landmark[0]*std_landmark[0]) + 
				    (y_observation - y_landmark)*(y_observation - y_landmark)/(std_landmark[1]*std_landmark[1])));
			p.weight *= w_obsi;
		}
		// Add p's weight to weight vector (def in particle_filter.h)
		weights.push_back(p.weight);
		// set associations of observations and corresponding closest landmark IDs for particle p
		p = SetAssociations(p, associated_ids, senseGlobal_x, senseGlobal_y);	
	}
}

void ParticleFilter::resample() {
	// Resamples particles with replacement with probability proportional to their weight. 
	// NOTE: std::discrete_distribution is helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// Create discrete distribution and generator based on weights 
	std::default_random_engine gen;
  	std::discrete_distribution<int> dist_w(weights.begin(), weights.end());

	// create temporary particle vector, copy from particles 
	std::vector<Particle> tempParticles;
	// Sample num_particles particles from particles vector
	for(int i=0; i<num_particles; i++){
		int pout = dist_w(gen);
		tempParticles.push_back(particles[pout]);	
	}
	// set particles vector
	particles = tempParticles;	
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
