#ifndef PATHINTEGRALS_HPP
#define PATHINTEGRALS_HPP

#include "UtilityFunctions.hpp"


enum struct BC_Type { fixed, free, periodic };

template<class System> struct Action;
template<class System> struct Metropolis;
template<class System> struct PathIntegrand;

struct Particle1d_Base;
struct EuclidParticle1d;
struct EuclidFreeParticle1d;
struct EuclidHarmonicOscillator1d;


#endif   /* PATHINTEGRALS_HPP */




