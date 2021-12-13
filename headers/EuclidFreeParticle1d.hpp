#ifndef EUCLIDFREEPARTICLE1D_HPP
#define EUCLIDFREEPARTICLE1D_HPP

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "EuclidParticle1d.hpp"


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


#endif   /* EUCLIDFREEPARTICLE1D_HPP */




