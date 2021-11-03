#ifndef DYNAMICAL_SYSTEM_HPP
#define DYNAMICAL_SYSTEM_HPP
#include <cmath>  // std::sinh, std::tanh
#include <boost/core/demangle.hpp>
#include <exception>
#include <iostream>
#include <map>
#include <string>
#include <UtilityFunctions.hpp>
#include <vector>

template<class T> class DynamicalException : public std::exception {
   public:
   DynamicalException(T* sys, const char* s = "") noexcept : obj(sys) {
      (_explain = ClassName(*sys)).append(": ").append(s);
   }
   T* const obj;
   const char* what() const noexcept override { return _explain.data(); }
   private:
   std::string _explain;
};

template<class CT, class VT> class DynamicalSystem {
   public:
   typedef CT _CoordType;
   typedef VT _ValueType;
   void AddBoundaryCondition(_CoordType t, _ValueType x) {
      _BCs.insert_or_assign(t, x);
   }
   void AddToPath(_CoordType t, _ValueType x) {
      _Path.insert_or_assign(t, x);
   }
   std::map<_CoordType, _ValueType> BoundaryConditions() const { return _BCs; }
   void ClearBoundaryConditions() { _BCs.clear(); }
   void ClearPath() { _Path.clear(); }
   virtual ~DynamicalSystem() {}
   std::map<_CoordType, _ValueType> Path() const { return _Path; }
   void PrintPath(std::ostream& os = std::cout) const { os << _Path << '\n'; }
   std::vector<_CoordType> PathCoords() const { return GetKeys(_Path); }
   typename std::map<_CoordType, _ValueType>::size_type PathSize() const {
      return _Path.size();
   }
   std::vector<_ValueType> PathValues() const { return GetValues(_Path); }
   template<class Generator, class...Args> void SetUserPath
    (const std::vector<_CoordType>& coords, Generator& gen, Args...args) {
      ClearPath();
      for(auto t : coords) AddToPath(t, gen(args...));
   }
   void SetPath(const std::map<_CoordType, _ValueType>& P) { _Path = P; }
   virtual void Throw(const char* why) const {
      throw DynamicalException(this, why);
   }
   protected:
   std::map<_CoordType, _ValueType> _BCs;
   std::map<_CoordType, _ValueType> _Path;
};

// specialization for _CoordType = _ValueType = double
template<> class DynamicalSystem<double, double> {
   public:
   typedef double _CoordType;
   typedef double _ValueType;
   void AddBoundaryCondition(double t, double x) {
      _BCs.insert_or_assign(t, x);
   }
   void AddToPath(double t, double x) {
      _Path.insert_or_assign(t, x);
   }
   std::map<double, double> BoundaryConditions() const { return _BCs; }
   virtual bool IsSolvable() const { return false; }
   virtual double ClassicalValue
    (double t, double epsi, bool checkIfSolvable = true) const {
      Throw("don't know how to solve this system");
      return 0.;
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   void ClearPath() { _Path.clear(); }
   virtual ~DynamicalSystem() {}
   std::map<double, double> Path() const { return _Path; }
   void PrintPath(std::ostream& os = std::cout) const { os << _Path << '\n'; }
   std::map<double, double>::size_type PathSize() const { return _Path.size(); }
   void SetPath(const std::map<double, double>& P) { _Path = P; }
   virtual void SetClassicalPath(double step, double epsi) {
      if(!IsSolvable()) Throw("cannot compute the classical path");
      if(step <= 0.) Throw("path step must be positive");
      // boundary conditions are included in the classical path:
      _Path = _BCs;
      std::vector<double> t = GetKeys(_Path);
      for(double tau = t[0] + step; tau < t[1]; tau += step)
         AddToPath(tau, ClassicalValue(tau, epsi, false));
   }
   template<class Generator, class...Args> void SetUserPath
    (double step, Generator& gen, Args...args) {
      if(step <= 0.) Throw("path step must be positive");
      // boundary conditions are included:
      _Path = _BCs;
      std::vector<double> t = GetKeys(_Path);
      for(double tau = t[0] + step; tau < t[1]; tau += step)
         AddToPath(tau, gen(args...));
   }
   template<class Generator, class...Args> void SetUserPath
    (const std::vector<double>& coords, Generator& gen, Args...args) {
      ClearPath();
      for(auto t : coords) AddToPath(t, gen(args...));
   }
   virtual void Throw(const char* why) const {
      throw DynamicalException(this, why);
   }
   protected:
   std::map<double, double> _BCs;
   std::map<double, double> _Path;
};

// euclidean one-dimensional free particle:
class EuclidFreeParticle1D : public DynamicalSystem<double, double> {
   public:
   double ClassicalValue
    (double tau, double epsi = 0., bool checkIfSolvable = true) const override {
      if(checkIfSolvable and !IsSolvable())
         Throw("need exactly 2 boundary conditions "
                "to solve the classical motion");
      std::vector<double> t = GetKeys(_BCs);
      std::vector<double> x = GetValues(_BCs);
      return ( (t[1] - tau)*x[0] + (tau - t[0])*x[1] )/( t[1] - t[0] );
   }
   EuclidFreeParticle1D(double mass = 1.) : _mass(mass) {
      if(mass <= 0.) Throw("must have positive mass");
   }
   double ExactAmplitude() const {
      if(!IsSolvable()) Throw("need exactly 2 boundary conditions "
                                 "to compute the amplitude");
      std::vector<double> t = GetKeys(_BCs);
      std::vector<double> x = GetValues(_BCs);
      double deltaT = t[1] - t[0],
             deltaX = x[1] - x[0];
      static constexpr double pi = std::acos(-1.);
      return std::sqrt( _mass/(2.*pi*deltaT) )
         *std::exp( -_mass*deltaX*deltaX/(2.*deltaT) );
   }
   bool IsSolvable() const override { return _BCs.size() == 2; }
   double Mass() const { return _mass; }
   void SetClassicalPath(double step, double epsi = 0.) override {
      DynamicalSystem::SetClassicalPath(step,epsi);
   }
   void Throw(const char* why) const override {
      throw DynamicalException(this, why);
   }
   private:
   const double _mass;
};

// euclidean one-dimensional harmonic oscillator:
class EuclidHarmonicOscillator1D : public DynamicalSystem<double, double> {
   public:
   double ClassicalValue
    (double tau, double epsi = 0., bool checkIfSolvable = true) const override {
      if(checkIfSolvable and !IsSolvable())
         Throw("need exactly 2 boundary conditions "
                "to solve the classical motion");
      std::vector<double> t = GetKeys(_BCs);
      std::vector<double> x = GetValues(_BCs);
      return ( std::sinh( _freq*(t[1] - tau) )*x[0]
                              + std::sinh( _freq*(tau - t[0]) )*x[1] )
                              / std::sinh( _freq*(t[1] - t[0]) );
   }
   EuclidHarmonicOscillator1D(double mass = 1., double freq = 1.)
    : _freq(freq), _mass(mass) {
      if(mass <= 0. or freq <= 0.) Throw("must have positive mass "
                                          "and positive frequency");
   }
   double ExactAmplitude() const {
      if(!IsSolvable()) Throw("need exactly 2 boundary conditions "
                                 "to compute the amplitude");
      std::vector<double> t = GetKeys(_BCs);
      std::vector<double> x = GetValues(_BCs);
      double deltaT = t[1] - t[0];
      static constexpr double pi = std::acos(-1.);
      return std::sqrt( _mass/(2.*pi*deltaT) )
         *std::exp( -_mass*_freq/2.
                     *( (x[0]*x[0] + x[1]*x[1])/std::tanh(_freq*deltaT)
                        - 2*x[0]*x[1]/std::sinh(_freq*deltaT) ) );
   }
   double Frequency() const { return _freq; }
   bool IsSolvable() const override { return _BCs.size() == 2; }
   double Mass() const { return _mass; }
   void SetClassicalPath(double step, double epsi = 0.) override {
      DynamicalSystem::SetClassicalPath(step,epsi);
   }
   void Throw(const char* why) const override { throw DynamicalException(this, why); }
   private:
   const double _freq;
   const double _mass;
};

// euclidean one-dimensional particle in a potential
template<class PotentialFunc> class EuclidParticle1D :
 public DynamicalSystem<double, double> {
   public:
   EuclidParticle1D(double mass = 1.,
    const PotentialFunc& pot = [](double){ return 0.; }) :
     _mass(mass), _potential(pot) {
            if(mass <= 0.) Throw("must have positive mass");
   }
   double Mass() const { return _mass; }
   PotentialFunc Potential() const { return _potential; }
   void Throw(const char* why) const override {
      throw DynamicalException(this, why);
   }
   private:
   const double _mass;
   const PotentialFunc _potential;
};

#endif // DYNAMICAL_SYSTEM_HPP




