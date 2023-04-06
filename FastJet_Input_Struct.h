#ifndef FastJet_Input_Struct_h
#define FastJet_Input_Struct_h

#include "Math/Vector4D.h"

using namespace ROOT::Math;

const int MaxParticles = 10000;

struct FastJet_Input{

int EventId;
float PtHat;
float SigmaGen;
float Weight;
int NFinal;

vector<int> ParticleId;
vector<int> ParticlePID;
vector<short> ParticleCharge;
vector<float> ParticleMass;
vector<PxPyPzEVector> Particle4Mom; 

};


#endif
