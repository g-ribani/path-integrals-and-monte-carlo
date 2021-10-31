#ifndef DYNAMICAL_SYSTEMS_HPP
#define DYNAMICAL_SYSTEMS_HPP
#include <cmath>
#include <boost/core/demangle.hpp>
#include <exception>
#include <iostream>
#include <random>
#include <string>
#include <UtilityFunctions.hpp>

template<typename, typename> class DynamicalSystem;

template<class T> class DynamicalException : public std::exception {
   public:
   // DynamicalException() {};
   DynamicalException(T* sys, const char* s = "") noexcept : obj(sys) {
      (_explain = ClassName(*obj)).append(": ").append(s);
   }
   // DynamicalException(const DynamicalException&) = default;
   // DynamicalException& operator=(const DynamicalException&) = default;
   // DynamicalException(DynamicalException&&) = default;
   // DynamicalException& operator=(DynamicalException&&) = default;
   // ~DynamicalException() {};
   T* const obj;
   const char* what() const noexcept { return _explain.data(); }
   private:
   std::string _explain;
};

template<class CT, class VT> class DynamicalSystem {
   public:
   typedef CT _CoordType;
   typedef VT _ValueType;
   void AddToPath(_CoordType c, _ValueType v) {
      _Path.insert_or_assign(c, v);
   }
   void ClearPath() { _Path.clear(); }
   virtual ~DynamicalSystem() = 0;
   // virtual bool IsExactlySolvable() const = 0;
   std::map<_CoordType, _ValueType> Path() const { return _Path; }
   void PrintPath
    (const std::string& s = "", std::ostream& o = std::cout) const {
      o << s << _Path << std::endl;
   }
   void SetPath(const std::map<_CoordType, _ValueType>& P) {
      ClearPath();
      _Path = P;
   }
   protected:
   std::map<_CoordType, _ValueType> _Path;
};

template<class CT, class VT> DynamicalSystem<CT, VT>::~DynamicalSystem() {}

// this is for dynamical systems for which a known exact (classical) solution
// exists, given the appropriate boundary conditions
template<class CT, class VT> class SolvableDynamicalSystem :
 public DynamicalSystem<CT, VT> {
   public:
   using typename DynamicalSystem<CT, VT>::_CoordType;
   using typename DynamicalSystem<CT, VT>::_ValueType;
   void AddToClassicalPath(_CoordType c) {
      this->_Path.insert_or_assign( c, operator()(c) );
   }
   std::map<_CoordType, _ValueType> ClassicalPath
    (const std::vector<_CoordType>& C) const {
      // IsSolvable_Assert();
      std::map<_CoordType, _ValueType> P;
      for(auto c : C) P.insert_or_assign(c, operator()(c) );
      return P;
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   auto GetBoundaryConditions() const { return _BCs; }
   // checks if boundary conditions uniquely determine the classical path
   virtual bool IsSolvable() const = 0;
   // computes the classical path
   virtual _ValueType operator()(_CoordType) const = 0;
   void SetBoundaryCondition(_CoordType c, _ValueType v) {
      _BCs.insert_or_assign(c, v);
   }
   void SetClassicalPath(const std::vector<_CoordType>& C) {
      this->SetPath( ClassicalPath(C) );
   }
   protected:
   std::map<_CoordType, _ValueType> _BCs;
   void IsSolvable_Assert() const {
      if(!IsSolvable())
         throw DynamicalException(this, "cannot compute the classical path");
   }
};

// euclidean one-dimensional harmonic oscillator
class HarmonicOscillator1D :
 public SolvableDynamicalSystem<double, double> {
   public:
   double Action() const  {
      if(_Path.size() < 2)
         throw DynamicalException(this, "cannot compute the action");
      auto t = GetKeys(_Path);
      auto x = GetValues(_Path);
      // x[0] = x_in, x[N] = x_fin
      std::size_t N = _Path.size()-1;
      double action = 0.;
      for(std::size_t k = 1; k != N+1; ++k) {
         double deltaT = t[k] - t[k-1],
                deltaX = x[k] - x[k-1];
         action += _mass/2.*deltaX*deltaX/deltaT
                   + deltaT*_mass*_freq*_freq/4.*( x[k]*x[k] + x[k-1]*x[k-1] );
      }
      return action;
   }
   HarmonicOscillator1D(double m = 1., double omega = 1.)
    : _freq(omega), _mass(m) {
      if(m <= 0. or omega <= 0.)
         throw DynamicalException(this, "must have positive mass "
                                          "and positive frequency");
   }
   // ~HarmonicOscillator1D() {};
   double ExactAmplitude() const {
      IsSolvable_Assert();
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0], deltaX = x[1] - x[0];
      static constexpr double pi = std::acos(-1.);
      return std::sqrt( _mass/(2.*pi*deltaT) )
         *std::exp( -_mass*deltaX*deltaX/(2.*deltaT) );
   }
   double Frequency() const { return _freq; }
   bool IsSolvable() const { return _BCs.size() == 2; }
   // bool IsExactlySolvable() const { return true; }
   double Mass() const { return _mass; }
   double operator()(double tau) const {
      IsSolvable_Assert();
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      if(_freq != 0.) return ( std::sinh( _freq*(t[1] - tau) )*x[0]
                              + std::sinh( _freq*(tau - t[0]) )*x[1] )
                              / std::sinh( _freq*(t[1] - t[0]) );
      // free particle if _freq == 0
      // else return ( (t[1] - tau)*x[0] + (tau - t[0])*x[1] )/( t[1] - t[0] );
   }
   private:
   const double _freq;
   const double _mass;
};

// euclidean one-dimensional free particle
class FreeParticle1D :
 public SolvableDynamicalSystem<double, double> {
   public:
   double Action() const  {
      if(_Path.size() < 2)
         throw DynamicalException(this, "cannot compute the action");
      auto t = GetKeys(_Path);
      auto x = GetValues(_Path);
      // x[0] = x_in, x[N] = x_fin
      std::size_t N = _Path.size()-1;
      double action = 0.;
      for(std::size_t k = 1; k != N+1; ++k) {
         double deltaT = t[k] - t[k-1],
                deltaX = x[k] - x[k-1];
         action += _mass/2.*deltaX*deltaX/deltaT;
      }
      return action;
   }
   FreeParticle1D(double m = 1.) : _mass(m) {
      if(m <= 0.) throw DynamicalException(this, "must have positive mass");
   }
   // ~FreeParticle1D() {};
   double ExactAmplitude() const {
      IsSolvable_Assert();
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0], deltaX = x[1] - x[0];
      static constexpr double pi = std::acos(-1.);
      return std::sqrt( _mass/(2.*pi*deltaT) )
         *std::exp( -_mass*deltaX*deltaX/(2.*deltaT) );
   }
   bool IsSolvable() const { return _BCs.size() == 2; }
   // bool IsExactlySolvable() const { return true; }
   double Mass() const { return _mass; }
   double operator()(double tau) const {
      IsSolvable_Assert();
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return ( (t[1] - tau)*x[0] + (tau - t[0])*x[1] )/( t[1] - t[0] );
   }
   private:
   const double _mass;
};

// // euclidean nDim-dimensional free particle
// template<size_t nDim> class EuclidFreeParticleND
//  : public SolvableDynamicalSystem<double, std::array<nDim, double>> {
// };

// template<class V> class Particle1D : public DynamicalSystem<double, double> {
//    public:
//    Particle1D(double m = 1., const V& pot = [](double x){ return 0.; }) :
//       _mass(m), _Potential(pot) {}
//    double Action() {
//       if(this->_Path.size() < 2)
//          throw DynamicalException(this, "cannot compute the action");
//       auto t = GetKeys(this->_Path);
//       auto x = GetValues(this->_Path);
//       // x[0] = x_in, x[N] = x_fin
//       std::size_t N = this->_Path.size()-1;
//       double action = 0.;
//       for(std::size_t k = 1; k != N+1; ++k) {
//          double deltaT = t[k] - t[k-1],
//                 deltaX = x[k] - x[k-1];
//          action += _mass/2.*deltaX*deltaX/deltaT
//                    + deltaT/2.*( _Potential(x[k]) + _Potential(x[k-1]) );
//       }
//       return action;
//    }
//    void SetMass(double m) { _mass = m; }
//    void SetPotential(const V& pot) { _Potential = pot; }
//    protected:
//    double _mass;
//    V _Potential;
// };
//
// template<class V> class EuclideanHarmonicOscillator : public Particle1D<V> {
//    public:
//    EuclideanHarmonicOscillator(double mass =1., double freq = 1.)
//     : Particle1D(mass, [mass, freq](double x){ return mass*freq*freq*x*x/2.; } )
//      {}
// };

#endif // DYNAMICAL_SYSTEMS_HPP




