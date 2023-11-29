#include "../include/Physics.h"
#include <bits/types/time_t.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <time.h>

void state_before_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                             ParSim ::Physics &physics, int steps,
                             double dimension, double phi, double Pecr);
void state_after_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                            ParSim ::Physics &physics);

int main() {

  // 1)Creating and initializaing particle system

  /*Parameters*/
  /*Try to stick to S.I units to make sense out of numbers*/
  int Number_of_particles = 784;
  int Number_of_time_steps = 2000;

  // Mips parameters
  double phi = 0.5; // packing fraction
  double L =
      sqrt(M_PI * Number_of_particles / (phi)); // periodic boundary length
  double Pecr = 50;                             // rotational Peclet number

  ParSim::ParticleSystem parsym(Number_of_particles, phi,
                                L); // create a simple system
  ParSim::Physics physics;

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to particles

  /*Setting physics parameters -- all game to be played here */
  physics.parameters[8] = 0.1;   // time step
  physics.parameters[0] = 1;     // k
  physics.parameters[1] = 2;     // interaction_diameter sigma
  physics.parameters[2] = 0.1;   // mass
  physics.parameters[5] = 1;     // gamma
  physics.parameters[10] = 0.01; // V0 --- active velocity
  physics.parameters[11] =
      physics.parameters[10] / Pecr; // Dr  -- rotational diffusion

  physics.parameters[4] = 0.8; // mu


  /*---obsolete params---*/

  physics.parameters[3] = 1; // radius

  physics.parameters[6] = 0.00000001;      // epsilon1  -- softening length
  physics.parameters[7] = M_PI / 10000000; // epsilon2 -- softening omega
  physics.parameters[9] = 0.5 * physics.parameters[5] /
                          pow((physics.parameters[2] * physics.parameters[0]),
                              0.5); // zeta

  /* 2)Reading initial state*/

  // std::ifstream input_state("init_state.txt");
  // int i = 0;

  // if (!input_state) { // file couldn't be opened
  //   std::cerr << "Error: Init_state file could not be opened" << std::endl;
  //   exit(1);
  // }

  // while (input_state >> particle[i].x >> particle[i].y >> particle[i].alpha
  // >>
  //        particle[i].vx >> particle[i].vy) {
  //   i++;
  // }

  // std::cout << "Input state read for " << i << " particles ..." << std
  // ::endl;

  // input_state.close();

  // 3)Creating a data file for storage and log-----------

  std::ofstream data_output;
  std::ofstream log;

  data_output.open("data.xyz");
  log.open("log.txt");

  // Print the state before the simulation in log
  state_before_simulation(log, parsym, physics, Number_of_time_steps, L, phi,
                          Pecr);

  log << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;
  std::cout << "-x-x-x-x-x-Simulation initiated-x-x-x-x-x- " << std::endl;
  std::cout << "No. of time steps: " << Number_of_time_steps << std::endl;

  time_t start = time(&start); // for measuring total runtime

  // 4) Main simulation loop--------------
  for (int step = 0; step < Number_of_time_steps; step++) {

    // writing data of this state to file (will be used for rendering the
    // system in ovito), write every nth state

    if (step % 10 == 0) {
      data_output << Number_of_particles << std::endl;
      data_output << "Lattice="
                  << "\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 0.0\"" << std::endl;
      // first store current configuration
      for (int i = 0; i < parsym.no_of_particles; ++i) {
        data_output << particle[i].x << ' ' << particle[i].y << ' ' << 0 << ' '
                    << cos(particle[i].alpha) << ' ' << sin(particle[i].alpha)
                    << ' ' << 0 << ' ' << particle[i].alpha << ' '
                    << particle[i].vx << ' ' << particle[i].vy << ' '
                    << particle[i].omega << ' ' << std::endl;
      }
    }

    if (step % 100 == 0) {
      std ::cout << "----------Step count: " << step << std::endl;
      log << "----------Step count: " << step << std::endl;
    }

    // Manipulate particle positions for next iteration.
    physics.evolve_system_ERM(parsym, step);
  }

  time_t end = time(&end);

  data_output.close();

  std::cout << "-x-x-x-x-x-Simulation ended-x-x-x-x-x-" << std::endl;
  log << "-x-x-x-x-x-Simulation ended-x-x-x-x-x-" << std::endl;
  log << "Runtime: " << end - start << " seconds" << std::endl;

  /*----------------------------------*/

  // Print the state before the simulation in log
  state_after_simulation(log, parsym, physics);

  return 0;
}

void state_before_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                             ParSim ::Physics &physics, int steps,
                             double dimension, double phi, double Pecr) {

  ParSim::Particle *const particle =
      parsym.get_particles(); // get access to paticles

  log << "-------Parameters and state------" << std::endl;
  log << "Number of particles: " << parsym.no_of_particles << std::endl
      << "Time step: " << physics.parameters[8] << std::endl
      << "Number of time steps: " << steps << std::endl
      << "Dimension: " << dimension << std::endl
      << "phi: " << phi << std::endl
      << "Pecr: " << Pecr << std::endl
      << "k: " << physics.parameters[0] << std::endl
      << "Interaction diameter (sigma): " << physics.parameters[1] << std::endl
      << "Mass (m): " << physics.parameters[2] << std::endl
      << "gamma: " << physics.parameters[5] << std::endl
      << "VO: " << physics.parameters[10] << std::endl
      << "Dr: " << physics.parameters[11] << std::endl
      << "mu: " << physics.parameters[4] << std::endl;


  log << "Energy-momentum before the collision: " << std::endl;
  log << "Total Energy: "
      << physics.EnergyMomentum(parsym)[0] + physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Translational K.Energy: " << physics.EnergyMomentum(parsym)[0]
      << std::endl;
  log << "Rotational K.Energy: " << physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Momentum: "
      << "(" << physics.EnergyMomentum(parsym)[2] << ", "
      << physics.EnergyMomentum(parsym)[3] << ")" << std::endl;

  log << "-------Initial conditions------" << std::endl;

  for (int i = 0; i < parsym.no_of_particles; ++i) {
    log << "Particle: " << i << " ---------" << std::endl;
    log << "x, y = " << particle[i].x << ", " << particle[i].y << std::endl;
    log << "V = " << particle[i].vx << ", " << particle[i].vy << std::endl;
    log << "Omega = " << particle[i].omega << std::endl;
    log << "V0 = " << particle[i].vx_activity << ", " << particle[1].vy_activity
        << std::endl;
    log << "Omega0 = " << particle[i].omega_activity << std::endl;
  }
};

void state_after_simulation(std::ofstream &log, ParSim::ParticleSystem &parsym,
                            ParSim ::Physics &physics) {
  log << "Energy-momentum After the collision: " << std::endl;
  log << "Total Energy: "
      << physics.EnergyMomentum(parsym)[0] + physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Translational K.Energy: " << physics.EnergyMomentum(parsym)[0]
      << std::endl;
  log << "Rotational K.Energy: " << physics.EnergyMomentum(parsym)[1]
      << std::endl;
  log << "Momentum: "
      << "(" << physics.EnergyMomentum(parsym)[2] << ", "
      << physics.EnergyMomentum(parsym)[3] << ")" << std::endl;
}
