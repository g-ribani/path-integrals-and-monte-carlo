#ifndef QMPATHINTEGRALS_HPP
#define QMPATHINTEGRALS_HPP

#include <cmath>
#include <vector>
#include "UtilityFunctions.hpp"

template<class T> struct PathIntegrand;

//// one-dimensional free particle in euclidean time

struct EuclidFreeParticle1D {
   struct Point {
      double t;
      double x;
   };
   EuclidFreeParticle1D(double const &m) : mass(m) {}
   Point initial;
   Point final;
   double mass;
   EuclidFreeParticle1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      time_step = (final.t - initial.t)/N;
      return *this;
   }
   double ExactAmplitude() const {
      double const deltaT = final.t - initial.t,
                   deltaX = final.x - initial.x;
      return std::sqrt( mass/(2.*M_PI*deltaT) )
         *std::exp( -mass*SQUARE(deltaX)/(2.*deltaT) );
   }
   protected:
   std::size_t n_steps;
   double time_step;
};

template<> struct PathIntegrand<EuclidFreeParticle1D>
 : EuclidFreeParticle1D {
   double operator() (const double *x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = n_steps - 1;
      std::vector<double> deltaX(n_steps);
      deltaX[0] = x[0] - initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = final.x - x[dim-1];
      double ret = 1.;
      for(size_type k = 0; k != dim + 1; ++k) {
         double deltaS = mass/2.*SQUARE(deltaX[k])/time_step;
         ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      }
      return ret;
   }
   double operator() (std::vector<double> const &x) const {
      return operator()( x.data() );
   }
   PathIntegrand(const EuclidFreeParticle1D &part)
      : EuclidFreeParticle1D(part) {}
};

//// one-dimensional harmonic oscillator in euclidean time

struct EuclidHarmonicOscillator1D {
   struct Point {
      double t;
      double x;
   };
   EuclidHarmonicOscillator1D(double const &m, double const &f)
    : frequency(f), mass(m) {}
   double ExactAmplitude() const {
      double const deltaT = final.t - initial.t,
                   x0 = initial.x,
                   x1 = final.x;
      return std::sqrt( mass*frequency/(2.*M_PI*std::sinh(frequency*deltaT)) )
         *std::exp( -mass*frequency/2.
                     *( (x0*x0 + x1*x1)/std::tanh(frequency*deltaT)
                        - 2.*x0*x1/std::sinh(frequency*deltaT) ) );
   }
   Point initial;
   Point final;
   double frequency;
   double mass;
   EuclidHarmonicOscillator1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      time_step = (final.t - initial.t)/N;
      return *this;
   }
   protected:
   std::size_t n_steps;
   double time_step;
};

template<> struct PathIntegrand<EuclidHarmonicOscillator1D>
 : EuclidHarmonicOscillator1D {
   double operator() (const double *x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = n_steps - 1;
      std::vector<double> deltaX(n_steps);
      deltaX[0] = x[0] - initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = final.x - x[dim-1];
      double ret = 1., deltaS;
      deltaS = mass/2.*SQUARE(deltaX[0])/time_step
               + time_step*mass*SQUARE(frequency)/4.
                  *( SQUARE(x[0]) + SQUARE(initial.x) );
      ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      for(size_type k = 1; k != dim; ++k) {
         deltaS = mass/2.*SQUARE(deltaX[k])/time_step
                  + time_step*mass*SQUARE(frequency)/4.
                     *( SQUARE(x[k]) + SQUARE(x[k-1]) );
         ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      }
      deltaS = mass/2.*SQUARE(deltaX[dim])/time_step
               + time_step*mass*SQUARE(frequency)/4.
                  *( SQUARE(final.x) + SQUARE(x[dim-1]) );
      ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      return ret;
   }
   double operator() (std::vector<double> const &x) const {
      return operator()( x.data() );
   }
   PathIntegrand(const EuclidHarmonicOscillator1D &osci)
      : EuclidHarmonicOscillator1D(osci) {}
};

//// one-dimensional particle in a potential in euclidean time

template<class PotentialFunc> struct EuclidParticle1D {
   struct Point {
      double t;
      double x;
   };
   EuclidParticle1D(double const &m,
                    PotentialFunc const &V = [](double const&) { return 0.; })
    : mass(m), potential(V) {}
   Point initial;
   Point final;
   double mass;
   PotentialFunc potential;
   EuclidParticle1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      time_step = (final.t - initial.t)/N;
      return *this;
   }
   protected:
   std::size_t n_steps;
   double time_step;
};

template<class PotentialFunc>
 struct PathIntegrand<EuclidParticle1D<PotentialFunc>>
 : EuclidParticle1D<PotentialFunc> {
   double operator() (const double *x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = this->n_steps - 1;
      std::vector<double> deltaX(this->n_steps);
      deltaX[0] = x[0] - this->initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = this->final.x - x[dim-1];
      double const &mass = this->mass;
      double &time_step = this->time_step;
      PotentialFunc &potential = this->potential;
      double ret = 1., deltaS;
      deltaS = mass/2.*SQUARE(deltaX[0])/time_step
               + time_step*( potential(x[0]) + potential(this->initial.x) )/2.;
      ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      for(size_type k = 1; k != dim; ++k) {
         double deltaS = mass/2.*SQUARE(deltaX[k])/time_step
                         + time_step*( potential(x[k]) + potential(x[k-1]) )/2.;
         ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      }
      deltaS = mass/2.*SQUARE(deltaX[dim])/time_step
              + time_step*( potential(this->final.x) + potential(x[dim-1]) )/2.;
      ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      return ret;
   }
   double operator() (std::vector<double> const &x) const {
      return operator()( x.data() );
   }
   PathIntegrand(const EuclidParticle1D<PotentialFunc> &part)
      : EuclidParticle1D<PotentialFunc>(part) {}
};

#endif   /* QMPATHINTEGRALS_HPP */




