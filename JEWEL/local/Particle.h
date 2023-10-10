#ifndef PARTICLE_H
#define PARTICLE_H

struct Particle {
  Int_t id;
  Float_t px, py, pz, energy;
  Int_t status;
};

#endif