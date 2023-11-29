//
// Created by Andrei Vasilev on 1/13/17.
//

#include "../include/Physics.h"
#include "../include/ParticleSystem.h"
#include <cmath>
#include <fstream>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <vector>

// void ParSim::Physics::default_physics(){
//  move(elapsed);
// }

ParSim::Physics::Physics() { int sample = 0; }

/*Force linker + integrators-- */

void ParSim::Physics::Force_PP(ParSim::ParticleSystem &ps) {

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // log << "1st loop -- " << std::endl;

    // Store previous step forces
    ps.particle_array[i].force_radial_prev[0] =
        ps.particle_array[i].force_radial[0];
    ps.particle_array[i].force_radial_prev[1] =
        ps.particle_array[i].force_radial[1];
    ps.particle_array[i].force_tangential_prev[0] =
        ps.particle_array[i].force_tangential[0];
    ps.particle_array[i].force_tangential_prev[1] =
        ps.particle_array[i].force_tangential[1];
    ps.particle_array[i].torque_prev = ps.particle_array[i].torque;

    // Reset the current forces
    ps.particle_array[i].force_radial[0] = 0;
    ps.particle_array[i].force_radial[1] = 0;
    ps.particle_array[i].force_tangential[0] = 0;
    ps.particle_array[i].force_tangential[1] = 0;
    ps.particle_array[i].torque = 0;

    // Unary force of damping. Always there. Translational and Rotational
    // activities added.
    ps.particle_array[i].force_radial[0] +=
        -1 * (this->parameters[5]) * ps.particle_array[i].vx +
        ps.particle_array[i].vx_activity;
    ps.particle_array[i].force_radial[1] +=
        -1 * (this->parameters[5]) * ps.particle_array[i].vy +
        ps.particle_array[i].vy_activity;
    ps.particle_array[i].torque +=
        -1 * (this->parameters[5]) * ps.particle_array[i].omega +
        ps.particle_array[i].omega_activity;

    // Knary force calculation --- Loop2: through all particles
    for (int j = 0; j < ps.no_of_particles; ++j) {
      if (j == i) { // no self coupling
        continue;
      }

      double d = ps.distance(ps.particle_array[i], ps.particle_array[j]);

      // U

      if (d <= (this->parameters[1])) {
        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction
        // radial interaction force
        ps.particle_array[i].force_radial[0] +=
            N * (ps.particle_array[i].x - ps.particle_array[j].x) /
            (d + (this->parameters[6]));

        ps.particle_array[i].force_radial[1] +=
            N * (ps.particle_array[i].y - ps.particle_array[j].y) /
            (d + (this->parameters[6]));

        // tangential friction force
        double omega_sum =
            (ps.particle_array[i].omega + ps.particle_array[j].omega);

        ps.particle_array[i].force_tangential[0] +=
            -(this->parameters[4]) *
            ((this->parameters[0]) * ((this->parameters[1]) - d)) *
            (omega_sum / (abs(omega_sum + (this->parameters[7])))) *
            ((ps.particle_array[i].y - ps.particle_array[j].y) /
             (d + (this->parameters[6])));

        ps.particle_array[i].force_tangential[1] +=
            -(this->parameters[4]) *
            ((this->parameters[0]) * ((this->parameters[1]) - d)) *
            (omega_sum / (abs(omega_sum + (this->parameters[7])))) *
            (-(ps.particle_array[i].x - ps.particle_array[j].x) /
             (d + (this->parameters[6])));

        // torque on particle

        if (omega_sum != 0) {
          ps.particle_array[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
        }

        if (omega_sum == 0) {
          ps.particle_array[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
        }
      }
    }
  }
}

void ParSim::Physics::Force_PP_PBC(ParSim::ParticleSystem &ps) {

  // Loop1: through all particles
  for (int i = 0; i < ps.no_of_particles; ++i) {

    // Store previous step forces
    ps.particle_array[i].force_radial_prev[0] =
        ps.particle_array[i].force_radial[0];
    ps.particle_array[i].force_radial_prev[1] =
        ps.particle_array[i].force_radial[1];
    ps.particle_array[i].force_tangential_prev[0] =
        ps.particle_array[i].force_tangential[0];
    ps.particle_array[i].force_tangential_prev[1] =
        ps.particle_array[i].force_tangential[1];
    ps.particle_array[i].torque_prev = ps.particle_array[i].torque;

    // Reset the current forces
    ps.particle_array[i].force_radial[0] = 0;
    ps.particle_array[i].force_radial[1] = 0;
    ps.particle_array[i].force_tangential[0] = 0;
    ps.particle_array[i].force_tangential[1] = 0;
    ps.particle_array[i].torque = 0;

    // wall forces and torque
    double force_wall_y = 0.0;
    double force_wall_fric_x = 0.0;
    double torque_wall = 0.0;
    double yp = 0.0;

    if (ps.particle_array[i].y > 0) {
      yp = ps.particle_array[i].y + (this->parameters[1] / 2);

    } else {
      yp = ps.particle_array[i].y - (this->parameters[1] / 2);
    }

    if (abs(yp) > ps.L / 2) {
      double norm = 2 * (abs(yp) - (ps.L / 2));     //kwall = 2
      force_wall_y =
          -norm * ps.particle_array[i].y / (abs(ps.particle_array[i].y));

      force_wall_fric_x =
          -this->parameters[4] * abs(norm) *
          (ps.particle_array[i].omega /
           (abs(ps.particle_array[i].omega + (this->parameters[7])))) *
          (-(ps.particle_array[i].y) / abs(ps.particle_array[i].y));

      if (ps.particle_array[i].y > 0) {
        torque_wall =
            -this->parameters[4] * abs(norm) *
            (ps.particle_array[i].omega /
             (abs(ps.particle_array[i].omega + (this->parameters[7])))) *
            (abs((ps.L / 2) - (ps.particle_array[i].y)));

      } else {
        torque_wall =
            -this->parameters[4] * abs(norm) *
            (ps.particle_array[i].omega /
             (abs(ps.particle_array[i].omega + (this->parameters[7])))) *
            (abs((ps.L / 2) + (ps.particle_array[i].y)));
      }

    } else {
      force_wall_y = 0.0;
      force_wall_fric_x = 0.0;
      torque_wall = 0.0;
    }

    // Unary forces.

    ps.particle_array[i].force_radial[0] +=
        -1 * (this->parameters[5]) * ps.particle_array[i].vx +
        ps.particle_array[i].vx_activity +
        this->parameters[10] * (this->parameters[5]) *
            cos(ps.particle_array[i].theta)+ force_wall_fric_x;

    ps.particle_array[i].force_radial[1] +=
        -1 * (this->parameters[5]) * ps.particle_array[i].vy +
        ps.particle_array[i].vy_activity +
        (this->parameters[10] * (this->parameters[5]) *
         sin(ps.particle_array[i].theta)) +
        force_wall_y;

    ps.particle_array[i].torque +=
        -1 * (this->parameters[5]) * ps.particle_array[i].omega +
        (ps.particle_array[i].omega_activity * (this->parameters[5])) +
        torque_wall;

    // Nnary force calculation --- Loop2: through all particles
    for (int j = 0; j < ps.no_of_particles; ++j) {
      if (j == i) { // no self coupling
        continue;
      }

      // distance from nearest image of jth particle
      double d = ps.nearest_img_dist_wall_y(ps.particle_array[i],
                                            ps.particle_array[j]);

      // U

      if (d <= (this->parameters[1])) {
        // mag. of radial force
        double N = (this->parameters[0]) *
                   ((this->parameters[1]) -
                    d); // magnitude of radial force used as normal reaction
        // radial interaction force

        double dx = ps.min_sep(ps.particle_array[i].x, ps.particle_array[j].x);
        double dy = ps.particle_array[i].y - ps.particle_array[j].y;

        ps.particle_array[i].force_radial[0] +=
            N * (dx) / (d + (this->parameters[6]));

        ps.particle_array[i].force_radial[1] +=
            N * (dy) / (d + (this->parameters[6]));

        // tangential friction force
        double omega_sum =
            (ps.particle_array[i].omega + ps.particle_array[j].omega);

        ps.particle_array[i].force_tangential[0] +=
            -(this->parameters[4]) *
            ((this->parameters[0]) * ((this->parameters[1]) - d)) *
            (omega_sum / (abs(omega_sum + (this->parameters[7])))) *
            ((dy) / (d + (this->parameters[6])));

        ps.particle_array[i].force_tangential[1] +=
            -(this->parameters[4]) *
            ((this->parameters[0]) * ((this->parameters[1]) - d)) *
            (omega_sum / (abs(omega_sum + (this->parameters[7])))) *
            (-(dx) / (d + (this->parameters[6])));

        // torque on particle

        if (omega_sum != 0) {
          ps.particle_array[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
        }

        if (omega_sum == 0) {
          ps.particle_array[i].torque +=
              -(this->parameters[4]) *
              ((this->parameters[0]) * ((this->parameters[1]) - d)) *
              (omega_sum / (abs(omega_sum + (this->parameters[7])))) * d;
        }
      }
    }
  }
}

void ParSim::Physics::Euler_Integrator(ParSim::Particle &par, int step) {

  double m = this->parameters[2];
  double time_step = this->parameters[8];
  double Fx;
  double Fy;
  double Tau;
  double dvx;
  double dvy;
  double dw;
  // Total force and torque on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1]; // -g newton, m = 1
  // log << "Ftangential: " << par.force_tangential[1] << std::endl;
  Tau = par.torque;

  dvx = (Fx / m) * time_step;
  dvy = (Fy / m) * time_step;
  // log << "dvy: " << dvy << std::endl;
  dw = (Tau)*time_step;

  // update the attributes

  par.x += par.vx * time_step;
  par.y += par.vy * time_step;
  par.vx += dvx;
  par.vy += dvy;

  // log << "vy: " << par.vy << std::endl;

  par.alpha += par.omega * time_step;
  par.omega += dw;
}

void ParSim::Physics::Vel_Verlet_Integrator(ParSim::Particle &par, int step) {

  double m = this->parameters[2];
  double time_step = this->parameters[8];

  // update the attributes

  par.x += (time_step * par.vx) +
           (pow(time_step, 2)) *
               (par.force_radial[0] + par.force_tangential[0]) / (2 * m);
  par.y += (time_step * par.vy) +
           (pow(time_step, 2)) *
               (par.force_radial[1] + par.force_tangential[1]) / (2 * m);

  par.alpha +=
      (time_step * par.omega) + (pow(time_step, 2)) * (par.torque) / (2 * m);

  if (step != 0) {
    par.vx += time_step *
              ((par.force_radial[0] + par.force_tangential[0]) +
               (par.force_radial_prev[0] + par.force_tangential_prev[0])) /
              (2 * m);
    par.vy += time_step *
              ((par.force_radial[1] + par.force_tangential[1]) +
               (par.force_radial_prev[1] + par.force_tangential_prev[1])) /
              (2 * m);
    par.omega += time_step * (par.torque + par.torque_prev) / (2 * m);
  }

  // log << "vy: " << par.vy << std::endl;
}

void ParSim::Physics::ERM_Integrator1(ParSim::Particle &par, int step) {

  double Fx;
  double Fy;
  double Tau;
  // Total force and torque on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1]; //

  Tau = par.torque;

  // Save present attributes
  par.position_prev[0] = par.x;
  par.position_prev[1] = par.y;

  par.velocity_prev[0] = par.vx;
  par.velocity_prev[1] = par.vy;

  par.alpha_prev = par.alpha;
  par.omega_prev = par.omega;

  // update the attributes upto midpoint (x', v')

  par.x += par.vx * (this->parameters[8]) / 2; // x'
  par.y += par.vy * (this->parameters[8]) / 2;
  par.vx += (Fx / (this->parameters[2])) * (this->parameters[8]) / 2; // v'
  par.vy += (Fy / (this->parameters[2])) * (this->parameters[8]) / 2;

  par.alpha += par.omega * (this->parameters[8]) / 2;
  par.omega += (Tau / (this->parameters[2])) * (this->parameters[8]) / 2;

  // For random active force -- use Fundamental Euler integration to update by
  // half step
  std::mt19937 mt(this->generator());
  std::normal_distribution<double> distribution(0.0, 1.0);

  par.theta += pow(this->parameters[11], 0.5) * distribution(mt) *
               (this->parameters[8] / 2); // rotational diffusion
}

void ParSim::Physics::ERM_Integrator2(ParSim::Particle &par, double L,
                                      int step) {

  double Fx;
  double Fy;
  double Tau;
  // Total force (F') and torque (Tau') on this particle
  Fx = par.force_radial[0] + par.force_tangential[0];
  Fy = par.force_radial[1] + par.force_tangential[1];

  Tau = par.torque;

  // update the attributes to next time step (x1,v1)

  par.x = par.position_prev[0] + par.vx * (this->parameters[8]); // x1
  par.y = par.position_prev[1] + par.vy * (this->parameters[8]);
  par.vx = par.velocity_prev[0] +
           (Fx / (this->parameters[2])) * (this->parameters[8]); // v1
  par.vy = par.velocity_prev[1] +
           (Fy / (this->parameters[2])) * (this->parameters[8]);

  par.alpha = par.alpha_prev + par.omega * (this->parameters[8]);
  par.omega =
      par.omega_prev + (Tau / (this->parameters[2])) * (this->parameters[8]);

  // For random active force -- use Fundamental Euler integration to update by
  // half step
  std::mt19937 mt(this->generator());
  std::normal_distribution<double> distribution(0.0, 1.0);

  par.theta += pow(this->parameters[11], 0.5) * distribution(mt) *
               (this->parameters[8] / 2); // rotational diffusion

  // Periodic boudary condition

  if (par.x > L / 2) {
    par.x -= L;
  } else if (par.x < -L / 2) {
    par.x += L;
  }

  // if (par.y > L / 2) {
  //   par.y -= L;
  // } else if (par.y < -L / 2) {
  //   par.y += L;
  // }
}

void ParSim ::Physics::Integrator(ParSim::ParticleSystem &parsym, int step) {

  for (int i = 0; i < parsym.no_of_particles; ++i) {
    // Vel_Verlet_Integrator(parsym.particle_array[i], step, log);
    ERM_Integrator1(parsym.particle_array[i], step);
    //  boundary conditions

    // if (parsym.particle_array[i].x < -1000 ||
    //     parsym.particle_array[i].x > 1000 ||
    //     parsym.particle_array[i].y < -800 || parsym.particle_array[i].y >
    //     800)
    //   parsym.particle_array[i].random_initialize();
  }
}

void ParSim::Physics::evolve_system(ParticleSystem &parsym, int step) {

  // 1)Force-linking--------
  ParSim::Physics::Force_PP(parsym); // links forces on each object
                                     // at this stage
  // 2)Integrating----------
  ParSim::Physics::Integrator(
      parsym, step); // applies forces on each object as determined above
}

void ParSim::Physics::evolve_system_ERM(ParticleSystem &parsym, int step) {

  // i) Calculate force (F) from positions and velocities (x,v) --
  Force_PP_PBC(parsym); // links forces on each object

  // ii) Update x,v to x', v'  ----
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    ERM_Integrator1(parsym.particle_array[i], step);
  }

  // iii) Again calculate force (F') from x',v' ------
  Force_PP_PBC(parsym);

  //  iv) Update x',v' to xnew, vnew --------
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    ERM_Integrator2(parsym.particle_array[i], parsym.L, step);
  }
}

std::vector<double> ParSim::Physics::EnergyMomentum(ParticleSystem &parsym) {
  std::vector<double> p{0, 0, 0, 0}; // kinetic, rotational, momenta
  double m = this->parameters[2];
  double I = m;
  for (int i = 0; i < parsym.no_of_particles; ++i) {
    p[0] += 0.5 * m *
            (pow(parsym.particle_array[i].vx, 2) +
             pow(parsym.particle_array[i].vy, 2)); // KE
    p[1] +=
        0.5 * I * (pow(parsym.particle_array[i].omega, 2)); // Rotational K.E
    p[2] += m * parsym.particle_array[i].vx;                // px
    p[3] += m * parsym.particle_array[i].vy;                // py
  }
  return p;
}
