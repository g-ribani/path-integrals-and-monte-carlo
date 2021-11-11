#include <chrono>
#include <cmath>  // std::exp, std::sqrt, std::acos, std::pow
#include "DynamicalSystem.hpp"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
// #include <iostream>  // std::cout
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
   return mass/2.*deltaX*deltaX/deltaT
             + deltaT/2.*( part->Potential(x1) + part->Potential(x0) );
}

template<class Obj> inline double Action(Basic_Particle1D* part) {
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

template<class Obj> inline double PathIntegrand(Basic_Particle1D* part) {
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

// to be used with gsl Vegas routine:
inline double PathIntegrand(double* x, std::size_t dim, void* params) {
   const double mass = ((double*)params)[0],
                freq = ((double*)params)[1],
                step = ((double*)params)[2],
                xin = ((double*)params)[3],
                xfin = ((double*)params)[4];
   double ret = 1., deltaX, deltaS;
   deltaX = x[0] - xin,
   deltaS = mass/2.*deltaX*deltaX/step
                   + step*mass*freq*freq/4.*( x[0]*x[0] + xin*xin );
   ret *= std::sqrt(mass/(2.*M_PI*step)) * std::exp(-deltaS);
   for(std::size_t k = 1; k != dim; ++k) {
      deltaX = x[k] - x[k-1],
      deltaS = mass/2.*deltaX*deltaX/step
                      + step*mass*freq*freq/4.*( x[k]*x[k] + x[k-1]*x[k-1] );
      ret *= std::sqrt(mass/(2.*M_PI*step)) * std::exp(-deltaS);
   }
   deltaX = xfin - x[dim-1],
   deltaS = mass/2.*deltaX*deltaX/step
                   + step*mass*freq*freq/4.*( xfin*xfin + x[dim-1]*x[dim-1] );
   ret *= std::sqrt(mass/(2.*M_PI*step)) * std::exp(-deltaS);
   return ret;
}

struct MCResult {
   double res;
   double err;
};

template<class Obj, class Generator> inline MCResult
 CrudeMCAmplitude(Basic_Particle1D* part, std::size_t nSteps,
  double xmin, double xmax, Generator& gen, std::size_t nEvals) {
   std::uniform_real_distribution<double> distr(xmin, xmax);
   std::vector<double> evals;
   gen.seed( std::chrono::high_resolution_clock::
               now().time_since_epoch().count() );
   for(size_t k = 0; k != nEvals; ++k) {
      // if( !(k % 10'000) )
      //    std::cout << (float)k << '\r' << std::flush;
      part->SetRandomPath(nSteps+1, distr, gen);
      evals.push_back( PathIntegrand<Obj>(part)
                        *std::pow(xmax - xmin, nSteps-1) );
   }
   double mean = std::accumulate(evals.begin(), evals.end(), 0.) / nEvals;
   double error = 0.;
   for(std::size_t k = 0; k != nEvals; ++k)
      error += std::pow(evals[k] - mean, 2.);
   error = std::sqrt( error/(nEvals*(nEvals-1)) );   // standard deviation of the mean
   return {mean,error};
}

// wrapper to simplify the call to CrudeMCAmplitude:
template<class Obj, class Generator> inline MCResult
 CrudeMCAmplitude (Obj* part, std::size_t nSteps,
  double xmin, double xmax, Generator& gen, std::size_t nEvals) {
   return CrudeMCAmplitude<Obj>
            ((Basic_Particle1D*)part, nSteps, xmin, xmax, gen, nEvals);
}

inline MCResult VegasMCAmplitude
 (EuclidHarmonicOscillator1D* osci, std::size_t nSteps,
  double xmin, double xmax, std::size_t nEvals) {
   auto BCs = osci->BoundaryConditions();
   auto t_bounds = GetKeys(BCs);
   if( t_bounds.size() != 2)
      osci->Throw("in function VegasMCAmplitude: need two boundary conditions");
   const double step = (t_bounds[1] - t_bounds[0])/nSteps;
   const std::size_t dim = nSteps - 1;
   double* params = new double[5];
   auto x_bounds = GetValues(BCs);
   params[0] = osci->Mass(),
   params[1] = osci->Frequency(),
   params[2] = step,
   params[3] = x_bounds[0],
   params[4] = x_bounds[1];
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
   return {result, stddev};
}

// this template integrates a generic function accepting a double* as an argument
template<class Function> inline MCResult VegasIntegrate
 (Function& func, std::size_t dim, void *params, double xmin,
  double xmax, std::size_t nEvals) {
   gsl_monte_function integrand {func, dim, params};
   double *xl = new double[dim],
          *xu = new double[dim];
   for(std::size_t k = 0; k != dim; ++k) xl[k] = xmin,
                                         xu[k] = xmax;
   gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937_1999);
   gsl_rng_set(gen, std::chrono::high_resolution_clock::
                        now().time_since_epoch().count() );
   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
   double result, stddev;
   int err = gsl_monte_vegas_integrate(&integrand, xl, xu, dim, nEvals,
                                       gen, s, &result, &stddev);
   gsl_monte_vegas_free(s);
   gsl_rng_free(gen);
   delete[] xl,
   delete[] xu;
   if( err != 0 ) throw err;
   return {result, stddev};
}




