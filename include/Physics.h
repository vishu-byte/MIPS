//
// Created by Vishu Saini on 31/08/23
//

#ifndef PARTICLE_SIMULATION_PHYSICS_H
#define PARTICLE_SIMULATION_PHYSICS_H

#include "ParticleSystem.h"
#include <math.h>
#include <random>

#include <fstream>
#include <iostream> //for writing data
#include <vector>

namespace ParSim { // for particle simulation

class Physics {
public:
  // class responsible for handling all physics behind the simulation
  std::vector<double> parameters{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  std::random_device generator;

  std::vector<double> ERM_pararms{0, 0, 0}; // e, delta, del t
  Physics();                                // constructor

  virtual ~Physics(){}; // destructor
public:
  /*Force linker + integrators-- */
  void Force_PP(ParticleSystem &, std::ofstream &);        // Main force linker
  void Force_PP_PBC(ParticleSystem &, std::ofstream &);    // Main force linker
  void Integrator(ParticleSystem &, int, std::ofstream &); // Main integrator
  void Euler_Integrator(Particle &, int, std::ofstream &);
  void Vel_Verlet_Integrator(Particle &, int, std::ofstream &);
  void ERM_Integrator1(Particle &, int, std::ofstream &);
  void ERM_Integrator2(Particle &, double, int, std::ofstream &);

  /*Conserved quantities*/
  std::vector<double> EnergyMomentum(ParticleSystem &);

  void evolve_system(ParticleSystem &, int,
                     std::ofstream &file); // takes a particle system and moves
                                           // it forward in time

  void evolve_system_ERM(ParticleSystem &, int, std::ofstream &file);
};


} // namespace ParSim

#endif // PARTICLE_SIMULATION_PHYSICS_H
