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
	// TODO: Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	int num_particles = 1000; // Number of particles, tunable parameter

	// initalize a random engine generator
	std::default_random_engine gen;

	// Normal (Gaussian) distribution for x, y and theta
	std::normal_distribution<double> distribution_x(x, std[0]);
	std::normal_distribution<double> distribution_y(y, std[1]);
	std::normal_distribution<double> distribution_theta(theta, std[2]);

	//std::vector<Particle> particles;
	for(int i=0; i<num_particles; i++){
		particles[i].id = i;
		particles[i].x = distribution_x(gen); 
		particles[i].y = distribution_y(gen);
		particles[i].theta = distribution_theta(gen);
		particles[i].weight = 1; // uniform distribution
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// initalize a random engine generator
	std::default_random_engine gen;

	for(auto p:particles){
		double xf = p.x + velocity/yaw_rate*(sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
		double yf = p.y + velocity/yaw_rate*(cos(p.theta) - cos(p.theta + yaw_rate*delta_t)); 	
		double thetaf = p.theta + yaw_rate*delta_t;
		
		// Normal (Gaussian) distribution for noise
		std::normal_distribution<double> distribution_x(xf, std_pos[0]);
		std::normal_distribution<double> distribution_y(yf, std_pos[1]);
		std::normal_distribution<double> distribution_theta(thetaf, std_pos[2]);

		// Save predicted state with Gaussian noise
		p.x = distribution_x(gen); 
		p.y = distribution_y(gen);
		p.theta = distribution_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// estimate distance between each landmark and particle position

	// calculate difference in distance between estimated distance (predicted observation) and the actual observation
	/*for(pred:predicted){ // DO WE KNOW WHAT LANDMARK THE OBSERVED LANEMARKS ARE CONNECTED TO??? THEN WE CAN JUST COMPARE??? SHOULD NOT BE POSSIBLE???

		float min_dist = -1.0;
		int id_min;
		for(obs:observations){
			dist_p = helpers.dist(pred.x, pred.y, obs.x, obs.y); //make sure in correct coordinate frame
			if(min_dist > dist_p || min_dist == -1.0){
				min_dist = dist_p;
				id_min = obs.id;
			}
		}
		// assign the observation to this landmark
		SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
	}*/


	// START WITH NEAREST NEIGHBOUR THEREAFTER TEST WITH JOINT CUMULATIVE DISTRIBUTION???


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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

	std::vector<double> senseGlobal_x;
	std::vector<double> senseGlobal_y;
	std::vector<int> associated_ids;
	for(auto& p:particles){
		// reset vectors
		senseGlobal_x.clear();
		senseGlobal_y.clear();
		associated_ids.clear();
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

			// Add weight to weight vector (def in particle_filter.h)
			weights.push_back(w_obsi);
		}
		// set associations of observations and corresponding closest landmarking IDs for particle p
		p = SetAssociations(p, associated_ids, senseGlobal_x, senseGlobal_y);		
	}

	// transform the car's measurements from its local coordinate system to the map's coordinate system

	// associate each measurement with a landmark identifier (use closest land mark to each transformed observation
	// Calculate the particles weight value

	//for(auto p:particles){ // will this need to be reference so we are sure we are updating p...
		/*w.append(p[i].measurement_prob(Z))
		prob = 1.0;
        	for i in range(len(landmarks)):
            		dist = sqrt((self.x - landmarks[i][0]) ** 2 + (self.y - landmarks[i][1]) ** 2)
            		prob *= self.Gaussian(dist, self.sense_noise, measurement[i])
        	return prob*/


		//p.weight *= 1/(2*M_PI*std_landmark[0]*std_landmark[1])*exp( -(p.sense_x -   
	//}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// Create discrete distribution and generator based on weights 
	std::default_random_engine gen;
  	std::discrete_distribution<int> dist_w(weights.begin(), weights.end());

	// create temporary particle vector, copy from particles 
	std::vector<Particle> tempParticles;
	// Sample num_particles particles from particles vector
	for(int i=0; i<num_particles; i++){
		tempParticles.push_back(particles[dist_w(gen)]);	
	}
	// set particles vector
	particles.clear();
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
