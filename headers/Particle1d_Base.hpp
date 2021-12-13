#ifndef PARTICLE1D_BASE_HPP
#define PARTICLE1D_BASE_HPP

#include <stdexcept>
#include "PathIntegrals.hpp"


struct Particle1d_Base {
   explicit Particle1d_Base(double const mass) : _mass(mass) {}
   struct Point {
      double t;
      double x;
   };
   Point Initial() const { return _initial; }
   void Initial(Point const &i) {
      _initial = i;
      _AdjustTimeStep();
   }
   Point Final() const { return _final; }
   void Final(Point const &f) {
      _final = f;
      _AdjustTimeStep();
   }
   double Mass() const { return _mass; }
   void Mass(double const m) { _mass = m; }
   size_t NSteps() const { return _n_steps; }
   void NSteps(size_t const n_steps) {
      if(n_steps == 0)
         throw std::invalid_argument("zero argument passed to method "
                                     "Particle1d_Base::NSteps");
      _n_steps = n_steps;
      _AdjustTimeStep();
   }
   double TimeStep() const { return _time_step; }
   void TimeStep(double const time_step) { _time_step = time_step; }
   protected:
   double _AdjustTimeStep() {
      return _time_step = (_final.t - _initial.t)/_n_steps;
   }
   Point _initial;
   Point _final;
   double _mass;
   size_t _n_steps = 1;
   double _time_step;
};


#endif   /* PARTICLE1D_BASE_HPP */




