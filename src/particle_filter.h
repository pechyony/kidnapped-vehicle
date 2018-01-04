/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

struct Particle {
	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};

class ParticleFilter {
public:
	
	// Set of current particles
	std::vector<Particle> particles;

	// Constructor
	// @param num_particles Number of particles
	ParticleFilter(int num_particles_) : num_particles(num_particles_), is_initialized(false), weights(num_particles_) {}

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initialize particle filter by sampling particles from Gaussian distribution around first position.
     * @param x Initial x position [m] (simulated estimate from GPS)
     * @param y Initial y position [m] (simulated estimate from GPS)
	 * @param theta Initial orientation [rad] (simulated estimate from GPS)
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
     *              standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta, double std[]);

	/**
	 * prediction Predict the state (position and yaw) for the next time step using CTRV process model.
	 *            Sample particle from the Gaussian around the predicted state.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *                  standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);
	
	/**
	 * dataAssociation Find a mapping between observations and landmarks by using a nearest-neighbors data association.
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	void dataAssociation(std::vector<Map::single_landmark_s> predicted, std::vector<LandmarkObs>& observations);
	
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *               observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations,
			           const Map &map_landmarks);
	
	/**
	 * resample Resample particles with replacement with probability proportional to their weight. 
	 */
	void resample();

	/**
	 * setAssociations Set particle's list of associations of measurements to landmarks, 
	 *                 along with the measurements calculated in world x,y coordinates
	 * @param particle Particle
	 * @param associations List of associations of measurements to landmarks  
	 * @param sense_x, sense_y  Associated measurements, calculated in world x,y coordinates
	 */
	Particle setAssociations(Particle& particle, const std::vector<int>& associations,
		                     const std::vector<double>& sense_x, const std::vector<double>& sense_y);

	/**
	 * getAssociations Get a string with particle's associations of measurements to landmarks
	 * @param particle Particle
	 */
	std::string getAssociations(Particle particle) const;

	/**
	* getSenseX Get a string with X coordinates of particle's measurements that are associated with landmarks.
	*           The returned coordinates are in the world coordinate system.
	* @param Particle
	*/
	std::string getSenseX(Particle particle) const;

	/**
	 * getSenseY Get a string with Y coordinates of particle's measurements that are associated with landmarks.
	 *           The returned coordinates are in the world coordinate system.
	 * @param Particle
	 */
	std::string getSenseY(Particle particle) const;

	/**
	 * initialized Returns whether particle filter is initialized yet or not.
	 */
	const bool initialized() const {
		return is_initialized;
	}

private:
	// Number of particles to draw
	int num_particles;

	// Flag, if filter is initialized
	bool is_initialized;

	// Vector of weights of all particles
	std::vector<double> weights;
};

#endif /* PARTICLE_FILTER_H_ */
