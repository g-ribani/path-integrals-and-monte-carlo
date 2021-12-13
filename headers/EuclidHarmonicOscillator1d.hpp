#ifndef EUCLIDHARMONICOSCILLATOR1D_HPP
#define EUCLIDHARMONICOSCILLATOR1D_HPP

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "EuclidParticle1d.hpp"


struct EuclidHarmonicOscillator1d : EuclidParticle1d {
   EuclidHarmonicOscillator1d(double const mass, double const freq)
   : EuclidParticle1d(mass), _frequency(freq) {
      Potential( [&](double x) { return 0.5*_mass*POW_2(_frequency*x); } );
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


#endif   /* EUCLIDHARMONICOSCILLATOR1D_HPP */




