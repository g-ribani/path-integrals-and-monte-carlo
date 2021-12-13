#ifndef EUCLIDPARTICLE1D_HPP
#define EUCLIDPARTICLE1D_HPP

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <functional>
#include <random>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <stdexcept>
#include "Particle1d_Base.hpp"

// todo: prevent damages which could happen by copying mamber pointers

struct EuclidParticle1d : Particle1d_Base {
   typedef std::function<double(double)> PotentialFunc;
   explicit EuclidParticle1d(double const mass) : Particle1d_Base(mass) {}
   EuclidParticle1d(double const mass, PotentialFunc const &potential)
   : Particle1d_Base(mass), _potential(potential) {}
   template<class Vector> double Action(Vector const &x) const {
      //^ x contains the internal spatial points of the path, while
      // the boundary spacetime points are stored inside the _particle itself.
      size_t const dim = _n_steps - 1;
      double action = 0.;
      action += 0.5*_mass*POW_2(x[0] - _initial.x)/_time_step + _time_step
                  *0.5*( _potential(x[0]) + _potential(_initial.x) );
      for(size_t k = 1; k != dim; ++k)
         action += 0.5*_mass*POW_2(x[k] - x[k-1])/_time_step + _time_step
                     *0.5*( _potential(x[k]) + _potential(x[k-1]) );
      action += 0.5*_mass*POW_2(_final.x - x[dim-1])/_time_step + _time_step
                  *0.5*( _potential(_final.x) + _potential(x[dim-1]) );
      return action;
   }
   double Action(std::vector<Particle1d_Base::Point> const &path) const {
      size_t const path_size = path.size();
      if(path_size == 0)
         throw std::invalid_argument("zero size path passed to method "
                                     "EuclidParticle1d::Action");
      Particle1d_Base::Point p1 = path[0], p2;
      double action = 0., time_step;
      for(size_t k = 1; k != path_size; ++k) {
         p2 = path[k];
         time_step = p2.t - p1.t;
         action += 0.5*_mass*POW_2(p2.x - p1.x)/time_step + time_step
                     *0.5*( _potential(p2.x) + _potential(p1.x) );
         p1 = p2;
      }
      return action;
   }
   template<class Vector> double PathIntegrand(Vector const &x) const {
      //^ x contains the internal spatial points of the path, while
      // the boundary spacetime points are stored inside the _particle itself.
      using namespace boost::math::double_constants;
      return std::pow(_mass/(two_pi*_time_step), half*_n_steps)
               *std::exp(-Action(x));
   }
   PotentialFunc Potential() const { return _potential; }
   void Potential(PotentialFunc const &potential) { _potential = potential; }
   protected:
   PotentialFunc _potential;
};

template<> struct Action <EuclidParticle1d> {
   Action(EuclidParticle1d const *part) : _particle(part) {}
   template<class Vector> double operator() (Vector const &x) const {
      return _particle->Action<Vector>(x);
   }
   double operator() (std::vector<Particle1d_Base::Point> const &path) const {
      return _particle->Action(path);
   }
   protected:
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
   (std::function<double(double*)> kernel,
    size_t const n_bins, size_t const n_corr, size_t const n_therm,
    double const epsi, size_t const bin_size = 1) {
      std::vector<std::function<double(double*)>>
         kernels(1, kernel);
      auto configs = (*this)(kernels, n_bins, n_corr, n_therm, epsi, bin_size);
      return configs[0];
   }
   std::vector<std::vector<double>> operator()
   (std::vector<std::function<double(double*)>> kernels,
    size_t const n_bins, size_t const n_corr, size_t const n_therm,
    double const epsi, size_t const bin_size = 1) {
      if(bin_size == 0)
         throw std::invalid_argument("zero bin size passed to "
                                     "Metropolis<EuclidParticle1d>");
      // Initialization:
      _n_steps = _particle->NSteps(),
      _mass = _particle->Mass(),
      _potential = _particle->Potential(),
      _time_step = _particle->TimeStep();
      switch(_bc_type) {
         case BC_Type::free:
         _path = NewPointerToLinearRange(_n_steps+1,
                                         _particle->Initial().x,
                                         _particle->Final().x);
         _start = 0, _stop = _n_steps;
         break;
         case BC_Type::fixed:
         _path = NewPointerToLinearRange(_n_steps+1,
                                         _particle->Initial().x,
                                         _particle->Final().x);
         _start = 1, _stop = _n_steps - 1;
         break;
         case BC_Type::periodic:
         _path = NewPointer
            (_n_steps+1, 0.5*(_particle->Initial().x + _particle->Final().x));
         _start = 0, _stop = _n_steps - 1;
      }
      _gen = new std::mt19937_64(std::random_device{}());
      size_t const n_kernels = kernels.size();
      std::vector<std::vector<double>> configs
         (n_kernels, std::vector<double>(n_bins));
      // Thermalization:
      for(size_t k = 0; k != n_therm; ++k) _Update(epsi);
      // Run:
      double *bin_avg = new double[n_kernels];
      _accepted = 0;
      for(size_t b = 0; b != n_bins; ++b) {
         for(size_t k = 0; k != n_kernels; ++k) bin_avg[k] = 0.;
         for(size_t c = 0; c != bin_size; ++c) {
            for(size_t i = 0; i != n_corr; ++i) _Update(epsi);
            for(size_t k = 0; k != n_kernels; ++k)
               bin_avg[k] += kernels[k](_path);
         }
         for(size_t k = 0; k != n_kernels; ++k)
            configs[k][b] = bin_avg[k]/bin_size;
      }
      delete[] bin_avg;
      // Store results:
      using namespace boost::math::statistics;
      _accept_rate =
         ((double)_accepted)/((double)(_stop-_start+1)*n_bins*bin_size*n_corr);
      _average.resize(n_kernels),
      _stddev.resize(n_kernels);
      for(size_t k = 0; k != n_kernels; ++k)
         _average[k] = mean(configs[k]),
         _stddev[k] = std::sqrt(sample_variance(configs[k])/n_bins);
      // Clear and return:
      delete[] _path,
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
      throw std::domain_error("in method Metropolis<EuclidParticle1d>::"
               "_ActionPiece: parameter k is outside the valid range "
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
   double *_path;
   EuclidParticle1d::PotentialFunc _potential;
   size_t _start;
   std::vector<double> _stddev;
   size_t _stop;
   double _time_step;
};

template<> struct PathIntegrand <EuclidParticle1d> {
   PathIntegrand(EuclidParticle1d const *part) : _particle(part) {}
   template<class Vector> double operator() (Vector const &x) const {
      return _particle->PathIntegrand<Vector>(x);
   }
   protected:
   EuclidParticle1d const *const _particle;
};


#endif   /* EUCLIDPARTICLE1D_HPP */




