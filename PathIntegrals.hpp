#ifndef PATHINTEGRALS_HPP
#define PATHINTEGRALS_HPP

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <functional>
#include <random>
#include <boost/math/statistics/univariate_statistics.hpp>
// #include <algorithm> // std::shuffle
#include <stdexcept> // std::domain_error
#include "UtilityFunctions.hpp"


template<class System> struct Action;
template<class System> struct Metropolis;
template<class System> struct PathIntegrand;


//// Base class for all 1d particles:

struct Particle1d_Base {
   Particle1d_Base(double const mass) : _mass(mass) {}
   struct Point {
      double t;
      double x;
   };
   Point Initial() const { return _initial; }
   void Initial(Point const &i) { _initial = i; }
   Point Final() const { return _final; }
   void Final(Point const &f) { _final = f; }
   double Mass() const { return _mass; }
   void Mass(double const m) { _mass = m; }
   protected:
   Point _initial;
   Point _final;
   double _mass;
};


////  Possible boundary conditions on paths,
////  for use, e.g. in the Metropolis algorithm:

enum struct BC_Type { fixed, free, periodic };


////  1d particle in a potential in euclidean time:

struct EuclidParticle1d : Particle1d_Base {
   typedef std::function<double(double)> PotentialFunc;
   EuclidParticle1d(double const mass) : Particle1d_Base(mass) {}
   EuclidParticle1d(double const mass, PotentialFunc const &potential)
   : Particle1d_Base(mass), _potential(potential) {}
   PotentialFunc Potential() const { return _potential; }
   void Potential(PotentialFunc const &potential) { _potential = potential; }
   protected:
   PotentialFunc _potential;
};

template<> struct Action <EuclidParticle1d> {
   Action(EuclidParticle1d const *part, size_t const n_steps)
   : _n_steps(n_steps), _particle(part) {}
   template<class Vector> double operator() (Vector const &x) const {
      size_t const dim = _n_steps - 1;
      Particle1d_Base::Point const initial = _particle->Initial(),
                                   final = _particle->Final();
      double const mass = _particle->Mass(),
                   time_step = (final.t - initial.t)/_n_steps;
      EuclidParticle1d::PotentialFunc const potential = _particle->Potential();
      double action = 0.;
      action += 0.5*mass*POW_2(x[0] - initial.x)/time_step + time_step
                  *0.5*( potential(x[0]) + potential(initial.x) );
      for(size_t k = 1; k != dim; ++k)
         action += 0.5*mass*POW_2(x[k] - x[k-1])/time_step + time_step
                     *0.5*( potential(x[k]) + potential(x[k-1]) );
      action += 0.5*mass*POW_2(final.x - x[dim-1])/time_step + time_step
                  *0.5*( potential(final.x) + potential(x[dim-1]) );
      return action;
   }
   protected:
   size_t const _n_steps;
   EuclidParticle1d const *const _particle;
};

template<> struct Metropolis <EuclidParticle1d> {
   Metropolis(EuclidParticle1d const *part, size_t const n_steps)
   : _n_steps(n_steps), _particle(part) {}
   std::vector<double> AbsErrors() const { return _stddev; }
   double AbsError() const { return _stddev[0]; }
   double AcceptRate() const { return _accept_rate; }
   Metropolis& operator()
   (std::function<double(std::vector<double>)> kernel,
    size_t const n_conf, size_t const n_corr, size_t const n_therm,
    double const epsi,
    BC_Type const bc_type) {
      std::vector<std::function<double(std::vector<double>)>>
         kernels(1, kernel);
      return (*this)(kernels, n_conf, n_corr, n_therm, epsi, bc_type);
   }
   Metropolis& operator()
   (std::vector<std::function<double(std::vector<double>)>> kernels,
    size_t const n_conf, size_t const n_corr, size_t const n_therm,
    double const epsi,
    BC_Type const bc_type) {
      // Initialization:
      _mass = _particle->Mass(),
      _potential = _particle->Potential(),
      _time_step = (_particle->Final().t - _particle->Initial().t)/_n_steps;
      switch(bc_type) {
         case BC_Type::free:
         _path = LinearRange(_particle->Initial().x,
                             _particle->Final().x,
                             _n_steps+1);
         _start = 0, _stop = _n_steps;
         break;
         case BC_Type::fixed:
         _path = LinearRange(_particle->Initial().x,
                             _particle->Final().x,
                             _n_steps+1);
         _start = 1, _stop = _n_steps - 1;
         break;
         case BC_Type::periodic:
         _path.assign(_n_steps+1,
                      0.5*(_particle->Initial().x + _particle->Final().x));
         _start = 0, _stop = _n_steps - 1;
      }
      std::mt19937_64 gen(std::random_device{}());
      std::vector<std::vector<double>>
         results(kernels.size(), std::vector<double>(n_conf));
      // Thermalization:
      for(size_t k = 0; k != n_therm; ++k) _Update(epsi, gen, bc_type);
      // Run:
      _accepted = 0;
      for(size_t c = 0; c != n_conf; ++c) {
         for(size_t k = 0; k != n_corr; ++k) _Update(epsi, gen, bc_type);
         for(size_t k = 0; k != kernels.size(); ++k)
            results[k][c] = kernels[k](_path);
      }
      // Store results:
      using namespace boost::math::statistics;
      _accept_rate =
         ((double)_accepted)/((double)(_stop-_start+1)*n_conf*n_corr);
      _average.resize(kernels.size()),
      _stddev.resize(kernels.size());
      for(size_t k = 0; k != kernels.size(); ++k)
         _average[k] = mean(results[k]),
         _stddev[k] = std::sqrt(sample_variance(results[k])/n_conf);
      return *this;
   }
   std::vector<double> Results() const { return _average; }
   double Result() const { return _average[0]; }
   protected:
   template<class Vector> double
   _ActionPiece(Vector const &x, size_t const k, BC_Type const bc_type) const {
      switch(bc_type) {
      case BC_Type::fixed: {
         if(0 < k and k < _n_steps)
            return _mass*x[k]*(x[k] - x[k+1] - x[k-1])/_time_step
                     + _time_step*_potential(x[k]);
      } break;
      case BC_Type::free: {
         if(0 < k and k < _n_steps)
            return _mass*x[k]*(x[k] - x[k+1] - x[k-1])/_time_step
                     + _time_step*_potential(x[k]);
         else if (k == 0)
            return _mass*x[k]*(0.5*x[k] - x[k+1])/_time_step
                     + _time_step*0.5*_potential(x[k]);
         else if (k == _n_steps)
            return _mass*x[k]*(0.5*x[k] - x[k-1])/_time_step
                     + _time_step*0.5*_potential(x[k]);
      } break;
      case BC_Type::periodic: {
            size_t prec{}, seq{};
            if(0 < k and k < _n_steps) prec = k-1, seq = k+1;
            else if(k == 0) prec = _n_steps-1, seq = 1;
            else break;
            return _mass*x[k]*(x[k] - x[seq] - x[prec])/_time_step
                     + _time_step*_potential(x[k]);
         }
      }
      throw std::domain_error("in method _ActionPiece: "
                  "parameter k is outside the valid range "
                  "with the given boundary conditions");
   }
   template<class Generator>
   void _Update(double const epsi, Generator &&gen, BC_Type const bc_type) {
      std::uniform_real_distribution unif(0., 1.);
      for(size_t k = _start; k <= _stop; ++k) {
         double old_x = _path[k],
                old_S = _ActionPiece(_path, k, bc_type);
         _path[k] = _path[k] + (2.*unif(gen)-1.)*epsi;
         double delta_S = _ActionPiece(_path, k, bc_type) - old_S;
         if(delta_S > 0. and std::exp(-delta_S) < unif(gen))
            _path[k] = old_x;
         else {
            _accepted++;
            if(k == 0 and bc_type == BC_Type::periodic)
               _path[_n_steps] = _path[0];
         }
      }
   }
   unsigned long _accepted;
   double _accept_rate;
   std::vector<double> _average;
   double _mass;
   size_t const _n_steps;
   EuclidParticle1d const *const _particle;
   std::vector<double> _path;
   EuclidParticle1d::PotentialFunc _potential;
   size_t _start;
   std::vector<double> _stddev;
   size_t _stop;
   double _time_step;
};

template<> struct PathIntegrand <EuclidParticle1d> {
   PathIntegrand(EuclidParticle1d const *part, size_t const n_steps)
   : _n_steps(n_steps), _particle(part) {}
   template<class Vector> double operator() (Vector const &x) const {
      using namespace boost::math::double_constants;
      double const time_step = (_particle->Final().t
                                - _particle->Initial().t)/_n_steps;
      return std::pow(_particle->Mass()/(two_pi*time_step), half*_n_steps)
               *std::exp(-Action<EuclidParticle1d>
                          (_particle, _n_steps)(x));
   }
   protected:
   size_t const _n_steps;
   EuclidParticle1d const *const _particle;
};


////  1d free particle in euclidean time:

struct EuclidFreeParticle1d : EuclidParticle1d {
   EuclidFreeParticle1d(double const mass) : EuclidParticle1d(mass) {
      Potential( [](double) { return 0.; } );
   }
   double ClassicalAction() const {
      return 0.5*_mass*POW_2(_final.x - _initial.x)/(_final.t - _initial.t);
   }
   double ClassicalPath(double const t) const {
      return _initial.x + (t - _initial.t)*(_final.x - _initial.x)
                                          /(_final.t - _initial.t);
   }
   double ExactAmplitude() const {
      using namespace boost::math::double_constants;
      return std::sqrt(_mass/(two_pi*(_final.t - _initial.t)))
               *std::exp(-ClassicalAction());
   }
   auto Potential() const { return EuclidParticle1d::Potential(); }
   protected:
   using EuclidParticle1d::Potential;
};


////  1d harmonic oscillator in euclidean time:

struct EuclidHarmonicOscillator1d : EuclidParticle1d {
   EuclidHarmonicOscillator1d(double const mass, double const freq)
   : EuclidParticle1d(mass), _frequency(freq) {
      Potential( [=](double x) { return 0.5*mass*POW_2(freq*x); } );
   }
   double ClassicalAction() const {
      using std::sinh, std::tanh;
      double const delta_T = _final.t - _initial.t;
      return 0.5*_mass*_frequency
             *( (POW_2(_initial.x) + POW_2(_final.x))/tanh(_frequency*(delta_T))
                - 2.*_initial.x*_final.x/sinh(_frequency*delta_T) );
   }
   double ClassicalPath(double const t) const {
      using std::sinh;
      return ( _initial.x*sinh(_frequency*(_final.t - t))
               + _final.x*sinh(_frequency*(t - _initial.t)) )
             /sinh(_frequency*(_final.t - _initial.t));
   }
   double ExactAmplitude() const {
      using namespace boost::math::double_constants;
      return
         std::sqrt( _mass*_frequency
                           /(two_pi*sinh(_frequency*(_final.t - _initial.t))) )
         *std::exp(-ClassicalAction());
   }
   double Frequency() const { return _frequency; }
   void Frequency(double const freq) { _frequency = freq; }
   auto Potential() const { return EuclidParticle1d::Potential(); }
   protected:
   double _frequency;
   using EuclidParticle1d::Potential;
};


#endif   /* PATHINTEGRALS_HPP */




