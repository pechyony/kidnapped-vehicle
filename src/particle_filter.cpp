/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *
 * Modified on: Jan 2, 2018
 *      Author: Dmitry Pechyony
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

using namespace std;

/**
* init Initialize particle filter by sampling particles from Gaussian distribution around first position.
* @param x Initial x position [m] (simulated estimate from GPS)
* @param y Initial y position [m] (simulated estimate from GPS)
* @param theta Initial orientation [rad] (simulated estimate from GPS)
* @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
*              standard deviation of yaw [rad]]
*/
void ParticleFilter::init(double x, double y, double theta, double std[]) 
{	
	Particle p;
	default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; i++) {
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		particles.push_back(p);
	}
	is_initialized = true;
}

/**
* prediction Predict the state (position and yaw) for the next time step using CTRV process model. 
*            Sample particle from the Gaussian around the predicted state.
* @param delta_t Time between time step t and t+1 in measurements [s]
* @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
*                  standard deviation of yaw [rad]]
* @param velocity Velocity of car from t to t+1 [m/s]
* @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
*/
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
	default_random_engine gen;

	for (int i = 0; i < num_particles; i++) {

		// Predict particle's position and yaw rate using CTRV model
		if (abs(yaw_rate) >= 1e-10) {
			// non-zero yaw rate
			particles[i].x += velocity * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) / yaw_rate;
			particles[i].y += velocity * (-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta)) / yaw_rate;
			particles[i].theta += yaw_rate * delta_t;
		}
		else {
			// yaw rate is close to zero
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}

		// Sample particle from the Gaussian around the predicted state
		normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
		
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
}

/**
* dataAssociation Find a mapping between observations and landmarks by using a nearest-neighbors data association.
* @param predicted Vector of predicted landmark observations
* @param observations Vector of landmark observations
*/
void ParticleFilter::dataAssociation(std::vector<Map::single_landmark_s> predicted, std::vector<LandmarkObs>& observations) 
{	
	vector<bool> assigned(predicted.size(), false); // indicator if a landmark is already assigned to some observation
	double min_distance, distance;
	int closest_landmark;

	for (int i = 0; i < observations.size(); i++) {

		// find closest unassigned landmark
		min_distance = 100000000;
		closest_landmark = -1;
		for (int j = 0; j < predicted.size(); j++) {
			if (assigned[j])
				continue;
			distance = dist(observations[i].x, observations[i].y, predicted[j].x_f, predicted[j].y_f);
			if (distance < min_distance) {
				min_distance = distance;
				closest_landmark = j;
			}
		}
		if (closest_landmark != -1) {
			// found unassigned landmark, create a mapping between observation and this landmark 
			observations[i].id = closest_landmark;
			assigned[closest_landmark] = true;
		}
	}
}

/**
* updateWeights Updates the weights for each particle based on the likelihood of the
*               observed measurements.
* @param sensor_range Range [m] of sensor
* @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
* @param observations Vector of landmark observations
* @param map Map class containing map landmarks
*/
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const vector<LandmarkObs> &observations, const Map &map_landmarks) 
{	
	for (int i = 0; i < num_particles; i++) {

		// transform observations from particle coordinate system into world coordinate system
		vector<LandmarkObs> observations_map;
		LandmarkObs observation_map;
		observation_map.id = -1;

		for (int j = 0; j < observations.size(); j++) {
			observation_map.x = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
			observation_map.y = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
			observations_map.push_back(observation_map);
		}

		// create list of landmarks that are approximately within sensor_range of the particle
		vector<Map::single_landmark_s> close_landmarks;
		for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
			if (abs(particles[i].x - map_landmarks.landmark_list[j].x_f) <= sensor_range &&
				abs(particles[i].y - map_landmarks.landmark_list[j].y_f) <= sensor_range)
				close_landmarks.push_back(map_landmarks.landmark_list[j]);

		// associate observations with the landmarks
		dataAssociation(close_landmarks, observations_map);

		// update weight of particle and prepare associations data
		vector<int> associations;
		vector<double> sense_x, sense_y;
		particles[i].weight = 1;
		for (int j = 0; j < observations_map.size(); j++) {
			if (observations_map[j].id != -1) {
				// observation is associated with some landmark, update the weight of the particle
				double new_weight;
				new_weight = normpdf(observations_map[j].x, observations_map[j].y,
					                 close_landmarks[observations_map[j].id].x_f, close_landmarks[observations_map[j].id].y_f,
					                 std_landmark[0], std_landmark[1]);
				particles[i].weight *= new_weight;

				// update associations data
				associations.push_back(close_landmarks[observations_map[j].id].id_i);
				sense_x.push_back(observations_map[j].x);
				sense_y.push_back(observations_map[j].y);
			}
		}
		
		// prevent the case when weights of all particles are zero (this might happen when the number of particles is small)
		if (particles[i].weight < 1e-20)
			particles[i].weight = 1e-20;

		// update global vector of particle weights
		weights[i] = particles[i].weight;

		// set associations of the particle
		setAssociations(particles[i], associations, sense_x, sense_y);
	}
}

/**
* resample Resample particles with replacement with probability proportional to their weight. 
*/
void ParticleFilter::resample() 
{	
	default_random_engine gen;
	discrete_distribution<> distribution(weights.begin(), weights.end());
	vector<Particle> new_particles;

	int ind;
	for (int i = 0; i < num_particles; i++) {
		ind = distribution(gen);
		new_particles.push_back(particles[ind]);
	}
	particles = new_particles;
}

/**
* setAssociations Set particle's list of associations of measurements to landmarks,
*                 along with the measurements calculated in world x,y coordinates
* @param particle Particle
* @param associations List of associations of measurements to landmarks
* @param sense_x, sense_y  Associated measurements, calculated in world x,y coordinates
*/
Particle ParticleFilter::setAssociations(Particle& particle, const vector<int>& associations, 
                                         const vector<double>& sense_x, const vector<double>& sense_y)
{
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

	return particle;
}

/**
* getAssociations Get a string with particle's associations of measurements to landmarks
* @param particle Particle
*/
string ParticleFilter::getAssociations(Particle best) const
{
	vector<int> v = best.associations;
	stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

/**
* getSenseX Get a string with X coordinates of particle's measurements that are associated with landmarks.
*           The returned coordinates are in the world coordinate system.
* @param Particle
*/
string ParticleFilter::getSenseX(Particle best) const
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

/**
* getSenseY Get a string with Y coordinates of particle's measurements that are associated with landmarks.
*           The returned coordinates are in the world coordinate system.
* @param Particle
*/
string ParticleFilter::getSenseY(Particle best) const
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
