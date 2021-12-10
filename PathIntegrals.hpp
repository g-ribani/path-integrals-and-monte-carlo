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
   explicit Particle1d_Base(double const mass) : _mass(mass) {}
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
   explicit EuclidParticle1d(double const mass) : Particle1d_Base(mass) {}
   EuclidParticle1d(double const mass, PotentialFunc const &potential)
   : Particle1d_Base(mass), _potential(potential) {}
   PotentialFunc Potential() const { return _potential; }
   void Potential(PotentialFunc const &potential) { _potential = potential; }
   protected:
   PotentialFunc _potential;
};

template<> struct Action <EuclidParticle1d> {
   Action(EuclidParticle1d const *part) : _particle(part) {}
   template<class Vector> double operator() (Vector const &x) const {
      //^ x contains the internal spatial points of the path, while
      // the boundary spacetime points are stored inside the _particle itself.
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
   double operator() (std::vector<Particle1d_Base::Point> const &path) const {
      double const mass = _particle->Mass();
      EuclidParticle1d::PotentialFunc const potential = _particle->Potential();
      Particle1d_Base::Point p1 = path[0], p2;
      double action = 0., time_step;
      for(size_t k = 1; k != path.size(); ++k) {
         p2 = path[k];
         time_step = p2.t - p1.t;
         action += 0.5*mass*POW_2(p2.x - p1.x)/time_step + time_step
                     *0.5*( potential(p2.x) + potential(p1.x) );
         p1 = p2;
      }
      return action;
   }
   size_t NSteps() const { return _n_steps; }
   void NSteps(size_t const n_steps) { _n_steps = n_steps; }
   protected:
   size_t _n_steps;
   EuclidParticle1d const *const _particle;
};

template<> struct Metropolis <EuclidParticle1d> {
   Metropolis(EuclidParticle1d const *part) : _particle(part) {}
   std::vector<double> AbsErrors() const { return _stddev; }
   double AbsError() const { return _stddev[0]; }
   double AcceptRate() const { return _accept_rate; }
   std::vector<double> Averages() const { return _average; }
   double Average() const { return _average[0]; }
   BC_Type BoundaryConditions() const { return _bc_type; }
   void BoundaryConditions(BC_Type const bc_type) { _bc_type = bc_type; }
   std::vector<double> operator()
   (std::function<double(std::vector<double>)> kernel,
    size_t const n_conf, size_t const n_corr, size_t const n_therm,
    double const epsi, size_t const bin_size = 0) {
      std::vector<std::function<double(std::vector<double>)>>
         kernels(1, kernel);
      auto configs = (*this)(kernels, n_conf, n_corr, n_therm, epsi, bin_size);
      return configs[0];
   }
   std::vector<std::vector<double>> operator()
   (std::vector<std::function<double(std::vector<double>)>> kernels,
    size_t const n_conf, size_t const n_corr, size_t const n_therm,
    double const epsi, size_t const bin_size = 0) {
//todo: implement binning
      // Initialization:
      _mass = _particle->Mass(),
      _potential = _particle->Potential(),
      _time_step = (_particle->Final().t - _particle->Initial().t)/_n_steps;
      switch(_bc_type) {
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
      _gen = new std::mt19937_64(std::random_device{}());
      std::vector<std::vector<double>> configs
         (kernels.size(), std::vector<double>(n_conf));
      // Thermalization:
      for(size_t k = 0; k != n_therm; ++k) _Update(epsi);
      // Run:
      _accepted = 0;
      for(size_t c = 0; c != n_conf; ++c) {
         for(size_t k = 0; k != n_corr; ++k) _Update(epsi);
         for(size_t k = 0; k != kernels.size(); ++k)
            configs[k][c] = kernels[k](_path);
      }
      // Store results:
      using namespace boost::math::statistics;
      _accept_rate =
         ((double)_accepted)/((double)(_stop-_start+1)*n_conf*n_corr);
      _average.resize(kernels.size()),
      _stddev.resize(kernels.size());
      for(size_t k = 0; k != kernels.size(); ++k)
         _average[k] = mean(configs[k]),
         _stddev[k] = std::sqrt(sample_variance(configs[k])/n_conf);
      // Clear and return:
      _path.clear(),
      _path.shrink_to_fit();
      delete _gen;
      return configs;
   }
   size_t NSteps() const { return _n_steps; }
   void NSteps(size_t const n_steps) { _n_steps = n_steps; }
   protected:
   double _ActionPiece(size_t const k) const {
      switch(_bc_type) {
      case BC_Type::fixed: {
         if(0 < k and k < _n_steps)
            return _mass*_path[k]*(_path[k]-_path[k+1]-_path[k-1])/_time_step
                     + _time_step*_potential(_path[k]);
      } break;
      case BC_Type::free: {
         if(0 < k and k < _n_steps)
            return _mass*_path[k]*(_path[k]-_path[k+1]-_path[k-1])/_time_step
                     + _time_step*_potential(_path[k]);
         else if (k == 0)
            return _mass*_path[k]*(0.5*_path[k]-_path[k+1])/_time_step
                     + _time_step*0.5*_potential(_path[k]);
         else if (k == _n_steps)
            return _mass*_path[k]*(0.5*_path[k] - _path[k-1])/_time_step
                     + _time_step*0.5*_potential(_path[k]);
      } break;
      case BC_Type::periodic: {
            size_t prec{}, seq{};
            if(0 < k and k < _n_steps) prec = k-1, seq = k+1;
            else if(k == 0) prec = _n_steps-1, seq = 1;
            else break;
            return _mass*_path[k]*(_path[k]-_path[seq]-_path[prec])/_time_step
                     + _time_step*_potential(_path[k]);
         }
      }
      throw std::domain_error("in method _ActionPiece: "
                  "parameter k is outside the valid range "
                  "with the given boundary conditions");
   }
   void _Update(double const epsi) {
      std::uniform_real_distribution unif(0., 1.);
      double old_x, old_S, delta_S;
      for(size_t k = _start; k <= _stop; ++k) {
         old_x = _path[k],
         old_S = _ActionPiece(k);
         _path[k] = _path[k] + (2.*unif(*_gen)-1.)*epsi;
         delta_S = _ActionPiece(k) - old_S;
         if(delta_S > 0. and std::exp(-delta_S) < unif(*_gen))
            _path[k] = old_x;
         else {
            _accepted++;
            if(k == 0 and _bc_type == BC_Type::periodic)
               _path[_n_steps] = _path[0];
         }
      }
   }
   unsigned long _accepted;
   double _accept_rate;
   std::vector<double> _average;
   BC_Type _bc_type;
   std::mt19937_64 *_gen;
   double _mass;
   size_t _n_steps;
   EuclidParticle1d const *const _particle;
   std::vector<double> _path;
   EuclidParticle1d::PotentialFunc _potential;
   size_t _start;
   std::vector<double> _stddev;
   size_t _stop;
   double _time_step;
};

template<> struct PathIntegrand <EuclidParticle1d> {
   PathIntegrand(EuclidParticle1d const *part) : _particle(part) {}
   template<class Vector> double operator() (Vector const &x) const {
      //^ x contains the internal spatial points of the path, while
      // the boundary spacetime points are stored inside the _particle itself.
      using namespace boost::math::double_constants;
      double const time_step = (_particle->Final().t
                                - _particle->Initial().t)/_n_steps;
      Action<EuclidParticle1d> action(_particle);
      action.NSteps(_n_steps);
      return std::pow(_particle->Mass()/(two_pi*time_step), half*_n_steps)
               *std::exp(-action(x));
   }
   size_t NSteps() const { return _n_steps; }
   void NSteps(size_t const n_steps) { _n_steps = n_steps; }
   protected:
   size_t _n_steps;
   EuclidParticle1d const *const _particle;
};


////  1d free particle in euclidean time:

struct EuclidFreeParticle1d : EuclidParticle1d {
   explicit EuclidFreeParticle1d(double const mass) : EuclidParticle1d(mass) {
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




