#ifndef DYNAMICAL_SYSTEMS_HPP
#define DYNAMICAL_SYSTEMS_HPP
#include <cmath>
#include <boost/core/demangle.hpp>
#include <exception>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <UtilityFunctions.hpp>

template<class T> class DynamicalException : public std::exception {
   public:
   DynamicalException(T* sys, const char* s = "") noexcept : obj(sys) {
      (_explain = ClassName(*sys)).append(": ").append(s);
   }
   T* const obj;
   const char* what() const noexcept { return _explain.data(); }
   private:
   std::string _explain;
};

template<class CT, class VT> class DynamicalSystem {
   public:
   typedef CT _CoordType;
   typedef VT _ValueType;
   void ClearPath() { _Path.clear(); }
   virtual ~DynamicalSystem() = 0;
   std::map<_CoordType, _ValueType> Path() const { return _Path; }
   void PrintPath(std::ostream& os = std::cout) const { os << _Path << '\n'; }
   void SetPath(const std::map<_CoordType, _ValueType>& P) { _Path = P; }
   protected:
   std::map<_CoordType, _ValueType> _Path;
};

template<class CT, class VT> DynamicalSystem<CT, VT>::~DynamicalSystem() {}

// euclidean one-dimensional free particle:
class FreeParticle1D : public DynamicalSystem<double, double> {
   public:
   double Action() const  {
      auto pathSize = _Path.size();
      if(pathSize < 2) Throw("need at least 2 points in the path "
                                  "to compute the action");
      auto t = GetKeys(_Path);
      auto x = GetValues(_Path);
      // x[0] = x_in, x[N] = x_fin
      std::vector<double>::size_type N = pathSize-1;
      double action = 0.;
      for(std::vector<double>::size_type k = 1; k != N+1; ++k) {
         double deltaT = t[k] - t[k-1],
                deltaX = x[k] - x[k-1];
         action += _mass/2.*deltaX*deltaX/deltaT;
      }
      return action;
   }
   void AddBoundaryCondition(double t, double x) {
      _BCs.insert_or_assign(t, x);
   }
   double At(double tau, bool scheckIfSolvable = true) const {
      if(scheckIfSolvable and _BCs.size() != 2)
         Throw("need exactly 2 boundary conditions "
                "to compute the classical path");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return ( (t[1] - tau)*x[0] + (tau - t[0])*x[1] )/( t[1] - t[0] );
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   double ExactAmplitude() const {
      if(_BCs.size() != 2) Throw("need exactly 2 boundary conditions "
                                  "to compute the amplitude");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0],
             deltaX = x[1] - x[0];
      static constexpr double pi = std::acos(-1.);
      return std::sqrt( _mass/(2.*pi*deltaT) )
         *std::exp( -_mass*deltaX*deltaX/(2.*deltaT) );
   }
   FreeParticle1D(double mass = 1.) : _mass(mass) {
      if(mass <= 0.) Throw("must have positive mass");
   }
   auto BoundaryConditions() const { return _BCs; }
   double Mass() const { return _mass; }
   void SetClassicalPath(double step) {
      if(_BCs.size() != 2) Throw("need exactly 2 boundary conditions "
                                  "to compute the classical path");
      // boundary conditions are included in the classical path:
      std::map<_CoordType, _ValueType> path = _BCs;
      auto t = GetKeys(_BCs);
      for(double tau = t[0] + step; tau < t[1]; tau += step)
         path.insert_or_assign(tau, this->At(tau, false));
      SetPath(path);
   }
   private:
   void Throw(const char* why) const { throw DynamicalException(this, why); }
   std::map<double, double> _BCs;
   const double _mass;
};

// // euclidean nDim-dimensional free particle
// template<size_t nDim> class EuclidFreeParticleND
//  : public DynamicalSystem<double, std::array<nDim, double>> {
// };

// euclidean one-dimensional harmonic oscillator:
class HarmonicOscillator1D : public DynamicalSystem<double, double> {
   public:
   double Action() const  {
      auto pathSize = _Path.size();
      if(pathSize < 2) Throw("need at least 2 points in the path "
                              "to compute the action");
      auto t = GetKeys(_Path);
      auto x = GetValues(_Path);
      // x[0] = x_in, x[N] = x_fin
      std::vector<double>::size_type N = pathSize - 1;
      double action = 0.;
      for(std::vector<double>::size_type k = 1; k != N+1; ++k) {
         double deltaT = t[k] - t[k-1],
                deltaX = x[k] - x[k-1];
         action += _mass/2.*deltaX*deltaX/deltaT
                   + deltaT*_mass*_freq*_freq/4.*( x[k]*x[k] + x[k-1]*x[k-1] );
      }
      return action;
   }
   void AddBoundaryCondition(double t, double x) {
      _BCs.insert_or_assign(t, x);
   }
   double At(double tau, bool scheckIfSolvable = true) const {
      if(scheckIfSolvable and _BCs.size() != 2)
         Throw("need exactly 2 boundary conditions "
                "to compute the classical path");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return ( std::sinh( _freq*(t[1] - tau) )*x[0]
                              + std::sinh( _freq*(tau - t[0]) )*x[1] )
                              / std::sinh( _freq*(t[1] - t[0]) );
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   HarmonicOscillator1D(double mass = 1., double freq = 1.)
    : _freq(freq), _mass(mass) {
      if(mass <= 0. or freq <= 0.) Throw("must have positive mass "
                                        "and positive frequency");
   }
   double ExactAmplitude() const {
      if(_BCs.size() != 2) Throw("need exactly 2 boundary conditions "
                                  "to compute the amplitude");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0];
      static constexpr double pi = std::acos(-1.);
      return std::sqrt( _mass/(2.*pi*deltaT) )
         *std::exp( -_mass*_freq/2.
                     *( (x[0]*x[0] + x[1]*x[1])/std::tanh(_freq*deltaT)
                        - 2*x[0]*x[1]/std::sinh(_freq*deltaT) ) );
   }
   double Frequency() const { return _freq; }
   auto BoundaryConditions() const { return _BCs; }
   double Mass() const { return _mass; }
   void SetClassicalPath(double step) {
      if(_BCs.size() != 2) Throw("need exactly 2 boundary conditions "
                                  "to compute the classical path");
      // boundary conditions are included in the classical path:
      std::map<_CoordType, _ValueType> path = _BCs;
      auto t = GetKeys(_BCs);
      for(double tau = t[0] + step; tau < t[1]; tau += step)
         path.insert_or_assign(tau, this->At(tau, false));
      SetPath(path);
   }
   private:
   void Throw(const char* why) const { throw DynamicalException(this, why); }
   std::map<double, double> _BCs;
   const double _freq;
   const double _mass;
};

// euclidean one-dimensional particle in a potential
template<class Potential> class Particle1D :
 public DynamicalSystem<double, double> {
   public:
   double Action() const {
      auto pathSize = this->_Path.size();
      if(pathSize < 2) Throw("need at least 2 points in the path "
                              "to compute the action");
      auto t = GetKeys(this->_Path);
      auto x = GetValues(this->_Path);
      // x[0] = x_in, x[N] = x_fin
      std::vector<double>::size_type N = pathSize - 1;
      double action = 0.;
      for(std::vector<double>::size_type k = 1; k != N+1; ++k) {
         double deltaT = t[k] - t[k-1],
                deltaX = x[k] - x[k-1];
         action += _mass/2.*deltaX*deltaX/deltaT
                   + deltaT/2.*( _potential(x[k]) + _potential(x[k-1]) );
      }
      return action;
   }
   void AddBoundaryCondition(double t, double x) {
      _BCs.insert_or_assign(t, x);
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   // void SetClassicalPath(double epsi); // to be implemented
   auto BoundaryConditions() const { return _BCs; }
   double Mass() const { return _mass; }
   Particle1D(double mass = 1.,
    const Potential& pot = [](double){ return 0.; }) :
     _mass(mass), _potential(pot) {
            if(mass <= 0.) Throw("must have positive mass");
   }
   private:
   void Throw(const char* why) const { throw DynamicalException(this, why); }
   std::map<double, double> _BCs;
   const double _mass;
   const Potential _potential;
};

#endif // DYNAMICAL_SYSTEMS_HPP




