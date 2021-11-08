#include <chrono>
#include <cmath>  // std::exp, std::sqrt, std::acos, std::pow
#include "DynamicalSystem.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <iostream>  // std::cout
#include <numeric>   // std::accumulate
#include <random> // std::uniform_real_distribution, std::normal_distribution
#include <utility>  // std::pair
#include <vector>

inline double ElementaryAction
 (EuclidFreeParticle1D* part, const double& t0, const double& x0,
  const double& t1, const double& x1) {
   const double mass = part->Mass(),
                deltaT = t1 - t0,
                deltaX = x1 - x0;
   return mass/2.*deltaX*deltaX/deltaT;
}

inline double ElementaryAction
 (EuclidHarmonicOscillator1D* osci, const double& t0, const double& x0,
  const double& t1, const double& x1) {
   const double mass = osci->Mass(),
                freq = osci->Frequency(),
                deltaT = t1 - t0,
                deltaX = x1 - x0;
   return mass/2.*deltaX*deltaX/deltaT
             + deltaT*mass*freq*freq/4.*( x1*x1 + x0*x0 );
}

template<class PotentialFunc> inline double ElementaryAction
 (EuclidParticle1D<PotentialFunc>* part, const double& t0, const double& x0,
  const double& t1, const double& x1) {
   const double mass = part->Mass(),
                deltaT = t1 - t0,
                deltaX = x1 - x0;
   PotentialFunc potential = part->Potential();
   return mass/2.*deltaX*deltaX/deltaT
             + deltaT/2.*( part.Potential(x1) + part.Potential(x0) );
}

template<class Obj> inline double Action(Particle1D* part) {
   const auto pathSize = part->PathSize();
   if(pathSize < 2)
      part->Throw("in function Action: need at least two points in the path");
   auto t = part->PathCoords();
   auto x = part->PathValues();
   double action = 0.;
   for(std::vector<double>::size_type k = 1; k != pathSize; ++k)
      action += ElementaryAction(dynamic_cast<Obj*>(part),
                                 t[k-1],x[k-1],t[k],x[k]);
   return action;
}

template<class Obj> inline double PathIntegrand(Particle1D* part) {
   const auto pathSize = part->PathSize();
   if(pathSize < 2)
      part->Throw("in function PathIntegrand: "
                  "need at least two points in the path");
   auto t = part->PathCoords();
   auto x = part->PathValues();
   std::vector<double>::size_type const N = pathSize - 1;
   const double mass = part->Mass();
   double integrand = 1.;
   for(std::vector<double>::size_type k = 1; k != N+1; ++k)
      integrand *= std::sqrt(mass/(2.*M_PI*(t[k] - t[k-1])))
                   *std::exp( -ElementaryAction(dynamic_cast<Obj*>(part),
                              t[k-1],x[k-1],t[k],x[k]) );
   return integrand; // exp(-S)*{normalization}
}

inline double PathIntegrand(double* x, std::size_t dim, void* params) { // to be used with gsl VEGAS routine
   const double mass = ((double*)params)[0],
                freq = ((double*)params)[1],
                step = ((double*)params)[2];
   double ret = 1.;
   for(std::size_t k = 1; k != dim + 2; ++k) {
      double deltaX = x[k] - x[k-1],
             deltaS = mass/2.*deltaX*deltaX/step
                      + step*mass*freq*freq/4.*( x[k]*x[k] + x[k-1]*x[k-1] );
      ret *= std::sqrt(mass/(2.*M_PI*step)) * std::exp(-deltaS);
   }
   return ret;
}

template<class Obj, class Generator> inline std::pair<double, double>
 CrudeMCAmplitude (Particle1D* part, std::size_t nSteps,
  double xmin, double xmax, Generator& gen, std::size_t nEvals) {
   std::uniform_real_distribution<double> distr(xmin, xmax);
   std::vector<double> evals;
   gen.seed( std::chrono::high_resolution_clock::
               now().time_since_epoch().count() );
   for(size_t k = 0; k != nEvals; ++k) {
      if( !(k % 10'000) )
         std::cout << (float)k << '\r' << std::flush;
      part->SetRandomPath(nSteps+1, distr, gen);
      evals.push_back( PathIntegrand<Obj>(part)
                        *std::pow(xmax - xmin, nSteps-1) );
   }
   double mean = std::accumulate(evals.begin(), evals.end(), 0.) / nEvals;
   double error = 0.;
   for(std::size_t k = 0; k != nEvals; ++k)
      error += std::pow(evals[k] - mean, 2.);
   error = std::sqrt( error/(nEvals*(nEvals-1)) );   // standard deviation of the mean
   return std::make_pair(mean,error);
}

// wrapper to simplify the call to CrudeMCAmplitude:
template<class Obj, class Generator> inline std::pair<double, double>
 CrudeMCAmplitude (Obj* part, std::size_t nSteps,
  double xmin, double xmax, Generator& gen, std::size_t nEvals) {
   return CrudeMCAmplitude<Obj>((Particle1D*)part, nSteps, xmin, xmax, gen, nEvals);
}

inline std::pair<double, double> VegasMCAmplitude
 (EuclidHarmonicOscillator1D& osci, std::size_t nSteps,
  double xmin, double xmax, std::size_t nEvals) {
   auto bounds = GetKeys( osci.BoundaryConditions() );
   if( bounds.size() != 2)
      osci.Throw("in function VEGASMCAmplitude: need two boundary conditions");
   const double step = (bounds[1] - bounds[0])/nSteps;
   const std::size_t dim = nSteps - 1;
   double* params = new double[3];
   params[0] = osci.Mass(),
   params[1] = osci.Frequency(),
   params[2] = step;
   gsl_monte_function func {PathIntegrand, dim, params};
   double *xl = new double[dim],
          *xu = new double[dim];
   for(std::size_t k = 0; k != dim; ++k) xl[k] = xmin,
                                         xu[k] = xmax;
   gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937_1999);
   gsl_rng_set(gen, std::chrono::high_resolution_clock::
                        now().time_since_epoch().count() );
   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
   double result, stddev;
   int err = gsl_monte_vegas_integrate(&func, xl, xu, dim, nEvals,
                                       gen, s, &result, &stddev);
   gsl_monte_vegas_free(s);
   gsl_rng_free(gen);
   delete[] params,
   delete[] xl,
   delete[] xu;
   if( err != 0 ) throw err;
   return std::make_pair(result, stddev);
}




