//
// Created by Andrei Vasilev on 1/13/17.
//

#include "../include/ParticleSystem.h"
#include <cmath>
#include <iostream>
#include <iterator>
#include <math.h>
#include <random>
#include <stdlib.h>

/*Class Particle definitions ------------*/
ParSim::Particle::Particle() { // default constructor
  x = 0;
  y = 0;
}

ParSim::Particle::Particle(int N, double phi, double L) {
  Lattice_initialize(N, phi, L);
}

ParSim::Particle::Particle(double x_cor, double y_cor, double v_x, double v_y,
                           double orientation) {
  x = x_cor;
  y = y_cor;
  vx = v_x;
  vy = v_y;
  alpha = orientation;
}

void ParSim::Particle::random_initialize(int N, double phi, double L) {

  // std:: cout << "Random init called with: (N,phi,L)" << N<< " "<<phi <<" " <<
  // L << std::endl;
  std::random_device rd;
  std::uniform_real_distribution<double> x_coordinate(-1, 1);
  std::uniform_real_distribution<double> y_coordinate(-1, 1);
  std::uniform_real_distribution<double> vx_dist(-1, 1);
  std::uniform_real_distribution<double> vy_dist(-1, 1);
  std::uniform_real_distribution<double> alpha_dist(-1, 1);
  std::uniform_real_distribution<double> theta_dist(-1, 1);
  std::uniform_real_distribution<double> omega_dist(-1, 1);

  x = 1 * (L / 2) * x_coordinate(rd);
  y = 1 * (L / 2) * y_coordinate(rd); // random distribution

  // Generate random particle speed. Speed is squared causing
  // particle distribution to be exponential instead of linear.
  vx = 0 * vx_dist(rd);
  vy = 0 * vy_dist(rd);

  // Generate random particle orientation (0 to 2pi) and omegas
  alpha = alpha_dist(rd);
  omega = 0 * M_PI * omega_dist(rd);

  // Generate random V0
  vx_activity = 0 * vx_dist(rd);
  vy_activity = 0 * vy_dist(rd);

  // Generate random theta
  theta = 0 * M_PI * theta_dist(rd);

  // Generatoe random omega
  omega_activity = 0 * M_PI * omega_dist(rd);
}

void ParSim::Particle::Lattice_initialize(int N, double phi, double L) {

  // std:: cout << "Random init called with: (N,phi,L)" << N<< " "<<phi <<" " <<
  // L << std::endl;

  double spacing =
      (L - 2.2) / std::sqrt(N); // 2.2 should be replaced with 1.1xsigma
  static double x_cor = 0;      // executed only once
  static double y_cor = 0;
  std::random_device rd;
  std::uniform_real_distribution<double> vx_dist(-1, 1);
  std::uniform_real_distribution<double> vy_dist(-1, 1);
  std::uniform_real_distribution<double> alpha_dist(-1, 1);
  std::uniform_real_distribution<double> theta_dist(-1, 1);
  std::uniform_real_distribution<double> omega_dist(-1, 1);

  // std:: cout << lattice_grid_ << ' ' << lattice_grid[1] << std::endl;
  if (x_cor > sqrt(N) - 1) {
    x_cor = 0;
    y_cor += 1;
  }

  // std::cout << x_cor << ' ' << y_cor << std::endl;

  // lattice grid distribution
  x = (spacing)*x_cor - 0.5 * (L - 2.2);
  y = (spacing)*y_cor - 0.5 * (L - 2.2);

  // Generate random particle speed.
  vx = 0 * vx_dist(rd);
  vy = 0 * vy_dist(rd);

  // Generate random particle orientation (0 to 2pi) and omegas
  alpha = 2 * M_PI * alpha_dist(rd);
  omega = 0 * M_PI * omega_dist(rd);

  // Generate random V0
  vx_activity = 0 * vx_dist(rd);
  vy_activity = 0 * vy_dist(rd);

  // Generate random theta -- active velocity director
  theta = 2 * M_PI * theta_dist(rd);

  // Generate random omega0
  omega_activity = (2 * M_PI / 628.0); // 0.01--anticlockwise rotation

  x_cor = (x_cor + 1);
}

/*Class Particle System definitions----------------*/
ParSim::ParticleSystem::ParticleSystem(int N, double phi, double dim) {
  this->no_of_particles = N;
  this->particle_array = new Particle[no_of_particles];
  this->L = dim;

  for (int i = 0; i < N; i++) {
    Particle temp(N, phi, L);
    particle_array[i] = temp;
  }
}

ParSim::ParticleSystem::~ParticleSystem() { delete[] particle_array; }

namespace ParSim {
Particle *const ParticleSystem::get_particles() { return particle_array; }
} // namespace ParSim

double ParSim::ParticleSystem::distance(Particle par1, Particle par2) {
  double distance =
      pow(pow((par1.x - par2.x), 2) + pow((par1.y - par2.y), 2), 0.5);
  return distance;
}

double ParSim::ParticleSystem::min_sep(double x1, double x2) {
  double dx = x1 - x2;
  if (dx < (-L / 2)) {
    dx += L;
  } else if (dx > (L / 2)) {
    dx -= L;
  }

  return dx;
}

double ParSim::ParticleSystem::nearest_img_dist(Particle par1, Particle par2) {

  double dist;

  dist = pow(
      (pow(min_sep(par1.x, par2.x), 2) + pow(min_sep(par1.y, par2.y), 2)), 0.5);

  return dist;
}

double ParSim::ParticleSystem::nearest_img_dist_wall_y(Particle par1,
                                                       Particle par2) {

  double dist;

  dist = pow((pow(min_sep(par1.x, par2.x), 2) + pow(par1.y - par2.y, 2)), 0.5);

  return dist;
}

double ParSim::ParticleSystem::dist_from_origin(Particle par) {
  double distance;
  distance = pow((pow(par.x, 2) + pow(par.y, 2)), 0.5);
  return distance;
}