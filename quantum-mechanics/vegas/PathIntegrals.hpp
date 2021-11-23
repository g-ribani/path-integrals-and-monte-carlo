#ifndef PATHINTEGRALS_HPP
#define PATHINTEGRALS_HPP

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include "UtilityFunctions.hpp"

template<class T> struct PathIntegrand;


//// base class for all 1d particles

struct Particle1D_Base {
   Particle1D_Base(double const mass) : _mass(mass) {}
   struct Point {
      double t;
      double x;
   };
   Point Initial() const { return _initial; }
   void Initial(double const t, double const x) { _initial = {t, x}; }
   void Initial(Point const &i) { _initial = i; }
   Point Final() const { return _final; }
   void Final(double const t, double const x) { _final = {t, x}; }
   void Final(Point const &f) { _final = f; }
   double Mass() const { return _mass; }
   void Mass(double const m) { _mass = m; }
   protected:
   Point _initial;
   Point _final;
   double _mass;
};


////  1d free particle in euclidean time

struct EuclidFreeParticle1D : Particle1D_Base {
   double ClassicalAction() const {
      double const deltaT = _final.t - _initial.t,
                   deltaX = _final.x - _initial.x;
      return _mass/2.*pow_2(deltaX)/deltaT;
   }
   double ClassicalPath(double const t) const {
      double const deltaT = _final.t - _initial.t,
                   deltaX = _final.x - _initial.x;
      return _initial.x + (t - _initial.t)*deltaX/deltaT;
   }
   double ExactAmplitude() const {
      using namespace boost::math::double_constants;
      using std::exp, std::sqrt;
      double const deltaT = _final.t - _initial.t;
      return sqrt(_mass/(two_pi*deltaT))*exp(-ClassicalAction());
   }
};

template<> struct PathIntegrand<EuclidFreeParticle1D>
: EuclidFreeParticle1D {
   PathIntegrand(EuclidFreeParticle1D const &part, size_t const n_steps)
   : EuclidFreeParticle1D(part), _n_steps(n_steps) {
      _time_step = (_final.t - _initial.t)/n_steps;
   }
   double operator() (double const *x) const {
      using namespace boost::math::double_constants;
      using std::exp, std::sqrt;
      using size_type = std::vector<double>::size_type;
      auto const dim = _n_steps - 1;
      std::vector<double> deltaX(_n_steps);
      deltaX[0] = x[0] - _initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = _final.x - x[dim-1];
      double ret = 1.;
      for(size_type k = 0; k != dim + 1; ++k) {
         double deltaS = _mass/2.*pow_2(deltaX[k])/_time_step;
         ret *= sqrt(_mass/(two_pi*_time_step)) * exp(-deltaS);
      }
      return ret;
   }
   double operator() (std::vector<double> const &x) const {
      return operator()( x.data() );
   }
   size_t NSteps() const { return _n_steps; }
   void NSteps(size_t const n_steps) { _n_steps = n_steps; }
   double TimeStep() const { return _time_step; }
   void TimeStep(double const time_step) { _time_step = time_step; }
   protected:
   size_t _n_steps;
   double _time_step;
};


////  1d harmonic oscillator in euclidean time

struct EuclidHarmonicOscillator1D : Particle1D_Base {
   EuclidHarmonicOscillator1D(double const mass, double const freq)
   : Particle1D_Base(mass), _freq(freq) {}
   double ClassicalAction() const {
      using std::sinh, std::tanh;
      double const deltaT = _final.t - _initial.t,
                   x0 = _initial.x,
                   x1 = _final.x;
      return _mass*_freq/2.*( (x0*x0 + x1*x1)/tanh(_freq*deltaT)
                                    - 2.*x0*x1/sinh(_freq*deltaT) );
   }
   double ClassicalPath(double const t) const {
      using std::sinh;
      return ( _initial.x*sinh(_freq*(_final.t - t))
               + _final.x*sinh(_freq*(t - _initial.t)) )
             /sinh(_freq*(_final.t - _initial.t));
   }
   double ExactAmplitude() const {
      using namespace boost::math::double_constants;
      using std::exp, std::sqrt;
      double const deltaT = _final.t - _initial.t;
      return sqrt( _mass*_freq/(two_pi*sinh(_freq*deltaT)) )
               *exp(-ClassicalAction());
   }
   double Frequency() const { return _freq; }
   void Frequency(double const freq) { _freq = freq; }
   protected:
   double _freq;
};

template<> struct PathIntegrand<EuclidHarmonicOscillator1D>
: EuclidHarmonicOscillator1D {
   PathIntegrand(EuclidHarmonicOscillator1D const &osci,
                 size_t const n_steps)
   : EuclidHarmonicOscillator1D(osci), _n_steps(n_steps) {
      _time_step = (_final.t - _initial.t)/n_steps;
   }
   double operator() (double const *x) const {
      using namespace boost::math::double_constants;
      using std::exp, std::sqrt;
      using size_type = std::vector<double>::size_type;
      auto const dim = _n_steps - 1;
      std::vector<double> deltaX(_n_steps);
      deltaX[0] = x[0] - _initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = _final.x - x[dim-1];
      double ret = 1., deltaS;
      deltaS = _mass/2.*pow_2(deltaX[0])/_time_step
               + _time_step*_mass*pow_2(_freq)/4.
                  *( pow_2(x[0]) + pow_2(_initial.x) );
      ret *= sqrt(_mass/(two_pi*_time_step)) * exp(-deltaS);
      for(size_type k = 1; k != dim; ++k) {
         deltaS = _mass/2.*pow_2(deltaX[k])/_time_step
                  + _time_step*_mass*pow_2(_freq)/4.
                     *( pow_2(x[k]) + pow_2(x[k-1]) );
         ret *= sqrt(_mass/(two_pi*_time_step)) * exp(-deltaS);
      }
      deltaS = _mass/2.*pow_2(deltaX[dim])/_time_step
               + _time_step*_mass*pow_2(_freq)/4.
                  *( pow_2(_final.x) + pow_2(x[dim-1]) );
      ret *= sqrt(_mass/(two_pi*_time_step)) * exp(-deltaS);
      return ret;
   }
   double operator() (std::vector<double> const &x) const {
      return operator()( x.data() );
   }
   size_t NSteps() const { return _n_steps; }
   void NSteps(size_t const n_steps) { _n_steps = n_steps; }
   double TimeStep() const { return _time_step; }
   void TimeStep(double const time_step) { _time_step = time_step; }
   protected:
   size_t _n_steps;
   double _time_step;
};


////  1d particle in a potential in euclidean time

template<class PotentialFunc> struct EuclidParticle1D : Particle1D_Base {
   EuclidParticle1D(double const mass,
                    PotentialFunc const &V = [](double const) { return 0.; })
   : Particle1D_Base(mass), _potential(V) {}
   PotentialFunc Potential() const { return _potential; }
   void Potential(PotentialFunc const &pot) { _potential = pot; }
   protected:
   PotentialFunc _potential;
};

template<class PotentialFunc>
 struct PathIntegrand<EuclidParticle1D<PotentialFunc>>
 : EuclidParticle1D<PotentialFunc> {
   PathIntegrand(EuclidParticle1D<PotentialFunc> const &part,
                 size_t const n_steps)
   : EuclidParticle1D<PotentialFunc>(part), _n_steps(n_steps) {
      _time_step = (this->_final.t - this->_initial.t)/n_steps;
   }
   double operator() (double const *x) const {
      using namespace boost::math::double_constants;
      using std::exp, std::sqrt;
      using size_type = std::vector<double>::size_type;
      auto const dim = _n_steps - 1;
      std::vector<double> deltaX(_n_steps);
      deltaX[0] = x[0] - this->_initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = this->_final.x - x[dim-1];
      double const mass = this->_mass;
      PotentialFunc const &potential = this->_potential;
      double ret = 1., deltaS;
      deltaS = mass/2.*pow_2(deltaX[0])/_time_step + _time_step
                  *( potential(x[0]) + potential(this->_initial.x) )/2.;
      ret *= sqrt(mass/(two_pi*_time_step)) * exp(-deltaS);
      for(size_type k = 1; k != dim; ++k) {
         double deltaS = mass/2.*pow_2(deltaX[k])/_time_step + _time_step
                           *( potential(x[k]) + potential(x[k-1]) )/2.;
         ret *= sqrt(mass/(two_pi*_time_step)) * exp(-deltaS);
      }
      deltaS = mass/2.*pow_2(deltaX[dim])/_time_step + _time_step
                  *( potential(this->_final.x) + potential(x[dim-1]) )/2.;
      ret *= sqrt(mass/(two_pi*_time_step)) * exp(-deltaS);
      return ret;
   }
   double operator() (std::vector<double> const &x) const {
      return operator()( x.data() );
   }
   size_t NSteps() const { return _n_steps; }
   void NSteps(size_t const n_steps) { _n_steps = n_steps; }
   double TimeStep() const { return _time_step; }
   void TimeStep(double const time_step) { _time_step = time_step; }
   protected:
   size_t _n_steps;
   double _time_step;
};

#endif   /* PATHINTEGRALS_HPP */




