#ifndef DYNAMICAL_SYSTEM_HPP
#define DYNAMICAL_SYSTEM_HPP
#include <cmath>  // std::sinh, std::tanh, std::acos
#include <exception>
#include <iostream>  // std::cout
#include <map>
#include <random> // std::normal_distribution
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
   void AddBoundaryCondition(const _CoordType& t, const _ValueType& x) {
      _BCs.insert_or_assign(t, x);
   }
   void AddBoundaryCondition(const std::pair<_CoordType, _ValueType>& point) {
      AddBoundaryCondition(std::get<0>(point), std::get<1>(point));
   }
   void AddToPath(const _CoordType& t, const _ValueType& x) {
      _Path.insert_or_assign(t, x);
   }
   void AddToPath(const std::pair<_CoordType, _ValueType>& point) {
      AddToPath(std::get<0>(point), std::get<1>(point));
   }
   void AddToPath(const std::map<_CoordType, _ValueType>& subpath) {
      for(const auto& point : subpath) AddToPath(point);
   }
   std::map<_CoordType, _ValueType> BoundaryConditions() const { return _BCs; }
   virtual _ValueType ClassicalValue (_CoordType t, double epsi) const {
      Throw("don't know how to solve this system");
      return _ValueType();
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   void ClearPath() { _Path.clear(); }
   virtual ~DynamicalSystem() {}
   virtual bool IsSolvable() const { return false; }
   std::map<_CoordType, _ValueType> Path() const { return _Path; }
   std::vector<_CoordType> PathCoords() const { return GetKeys(_Path); }
   typename std::map<_CoordType, _ValueType>::size_type PathSize() const {
      return _Path.size();
   }
   std::vector<_ValueType> PathValues() const { return GetValues(_Path); }
   void PrintPath(std::ostream& os = std::cout) const { os << _Path << '\n'; }
   void SetClassicalPath
    (const std::vector<_CoordType>& coords, double epsi) {
      ClearPath();
      for(const auto& t : coords) AddToPath( t, ClassicalValue(t, epsi) );
   }
   void SetPath(const std::map<_CoordType, _ValueType>& P) { _Path = P; }
   template<class Distribution, class Generator> void SetRandomPath
    (const std::vector<_CoordType>& coords,
     Distribution& distr, Generator& gen) {
      auto func = [&distr, &gen](){ return distr(gen); };
      SetUserPath(coords, func);
   }
   template<class Function> void SetUserPath
    (const std::vector<_CoordType>& coords, Function& func) {
      ClearPath();
      for(const auto& t : coords) AddToPath(t, func());
   }
   virtual void Throw(const char* why) const {
      throw DynamicalException(this, why);
   }
   protected:
   std::map<_CoordType, _ValueType> _BCs, _Path;
};

// specialization for _CoordType = _ValueType = double
template<> class DynamicalSystem<double, double> {
   public:
   typedef double _CoordType;
   typedef double _ValueType;
   void AddBoundaryCondition(const double& t, const double& x) {
      _BCs.insert_or_assign(t, x);
   }
   void AddBoundaryCondition(const std::pair<double, double>& point) {
      AddBoundaryCondition(std::get<0>(point), std::get<1>(point));
   }
   void AddToPath(const double& t, const double& x) {
      _Path.insert_or_assign(t, x);
   }
   void AddToPath(const std::pair<double, double>& point) {
      AddToPath(std::get<0>(point), std::get<1>(point));
   }
   void AddToPath(const std::map<double, double>& subpath) {
      for(const auto& point : subpath) AddToPath(point);
   }
   std::map<double, double> BoundaryConditions() const { return _BCs; }
   virtual double ClassicalValue (double t, double epsi) const {
      Throw("in function ClassicalValue: don't know how to solve this system");
      return 0.;
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   void ClearPath() { _Path.clear(); }
   virtual ~DynamicalSystem() {}
   virtual bool IsSolvable() const { return false; }
   std::map<double, double> Path() const { return _Path; }
   std::vector<double> PathCoords() const { return GetKeys(_Path); }
   std::map<double, double>::size_type PathSize() const { return _Path.size(); }
   std::vector<double> PathValues() const { return GetValues(_Path); }
   void PrintPath(std::ostream& os = std::cout) const { os << _Path << '\n'; }
   void SetPath(const std::map<double, double>& P) { _Path = P; }
   void SetClassicalPath(std::size_t nPoints, double epsi) {
      if( _BCs.size() != 2 )
         Throw("in function SetClassicalPath: "
               "need two boundary conditions");
      auto bounds = GetKeys(_BCs);
      SetClassicalPath( LinearRange(bounds[0], bounds[1], nPoints), epsi );
   }
   void SetClassicalPath
    (const std::vector<double>& coords, double epsi) {
      ClearPath();
      for(const auto& t : coords) AddToPath(t, ClassicalValue(t, epsi));
   }
   template<class Distribution, class Generator> void SetRandomPath
    (std::size_t nPoints, Distribution& distr, Generator& gen) {
      if( _BCs.size() != 2 )
         Throw("in function SetRandomPath: "
               "need two boundary conditions");
      auto bounds = GetKeys(_BCs);
      auto coords = LinearRange(bounds[0], bounds[1], nPoints);
      coords.erase(coords.begin()),
      coords.erase(coords.end()-1);
      SetRandomPath(coords, distr, gen);
      AddToPath(_BCs);
   }
   template<class Distribution, class Generator> void SetRandomPath
    (const std::vector<double>& coords, Distribution& distr, Generator& gen) {
      auto func = [&distr, &gen](){ return distr(gen); };
      SetUserPath(coords, func);
   }
   template<class Function> void SetUserPath
    (std::size_t nPoints, Function& func) {
      if( _BCs.size() != 2 )
         Throw("in function SetUserPath: "
               "need two boundary conditions");
      auto bounds = GetKeys(_BCs);
      auto coords = LinearRange(bounds[0], bounds[1], nPoints);
      coords.erase(coords.begin()),
      coords.erase(coords.end()-1);
      SetUserPath(coords, func);
      AddToPath(_BCs);
   }
   template<class Function> void SetUserPath
    (const std::vector<double>& coords, Function& func) {
      ClearPath();
      for(const auto& t : coords) AddToPath(t, func());
   }
   virtual void Throw(const char* why) const {
      throw DynamicalException(this, why);
   }
   protected:
   std::map<double, double> _BCs, _Path;
};

// euclidean one-dimensional free particle:
class EuclidFreeParticle1D : public DynamicalSystem<double, double> {
   public:
   double ClassicalAction() const {
      if( _BCs.size() != 2 )
         Throw("in function ClassicalAction: "
               "need two boundary conditions");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return _mass/2.* (x[1] - x[0]) / (t[1] - t[0]);
   }
   double ClassicalValue
    (double tau, double epsi = 0.) const override {
      if(!IsSolvable())
         Throw("in function ClassicalValue: cannot solve the classical motion");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return ( (t[1] - tau)*x[0] + (tau - t[0])*x[1] )/( t[1] - t[0] );
   }
   EuclidFreeParticle1D(double mass = 1.) : _mass(mass) {
      if(mass <= 0.) Throw("in constructor: must have positive mass");
   }
   double ExactAmplitude() const {
      static constexpr double pi = std::acos(-1.);
      if( _BCs.size() != 2 )
         Throw("in function ExactAmplitude: "
               "need two boundary conditions");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0],
             deltaX = x[1] - x[0];
      return std::sqrt( _mass/(2.*pi*deltaT) )
         *std::exp( -_mass*deltaX*deltaX/(2.*deltaT) );
   }
   bool IsSolvable() const override { return _BCs.size() == 2; }
   double Mass() const { return _mass; }
   void Throw(const char* why) const override {
      throw DynamicalException(this, why);
   }
   private:
   const double _mass;
};

// euclidean one-dimensional harmonic oscillator:
class EuclidHarmonicOscillator1D : public DynamicalSystem<double, double> {
   public:
   double ClassicalAction() const {
      if( _BCs.size() != 2 )
         Throw("in function ClassicalAction: "
               "need two boundary conditions");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0];
      return _mass*_freq/2.
                       *( (x[0]*x[0] + x[1]*x[1])/std::tanh(_freq*deltaT)
                        - 2*x[0]*x[1]/std::sinh(_freq*deltaT) );
   }
   double ClassicalValue
    (double tau, double epsi = 0.) const override {
      if(!IsSolvable())
         Throw("in function ClassicalValue: cannot solve the classical motion");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return ( std::sinh( _freq*(t[1] - tau) )*x[0]
                           + std::sinh( _freq*(tau - t[0]) )*x[1] )
                           / std::sinh( _freq*(t[1] - t[0]) );
   }
   EuclidHarmonicOscillator1D(double mass = 1., double freq = 1.)
    : _freq(freq), _mass(mass) {
      if(mass <= 0. or freq <= 0.)
         Throw("in constructor: must have positive mass "
               "and positive frequency");
   }
   double ExactAmplitude() const {
      static constexpr double pi = std::acos(-1.);
      if( _BCs.size() != 2 )
         Throw("in function ExactAmplitude: "
               "need two boundary conditions");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0];
      return std::sqrt( _mass/(2.*pi*deltaT) )
         *std::exp( -_mass*_freq/2.
                     *( (x[0]*x[0] + x[1]*x[1])/std::tanh(_freq*deltaT)
                        - 2*x[0]*x[1]/std::sinh(_freq*deltaT) ) );
   }
   double Frequency() const { return _freq; }
   bool IsSolvable() const override { return _BCs.size() == 2; }
   double Mass() const { return _mass; }
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
            if(mass <= 0.) Throw("in constructor: must have positive mass");
   }
   double Mass() const { return _mass; }
   auto Potential(double x) const { return _potential(x); }
   void Throw(const char* why) const override {
      throw DynamicalException(this, why);
   }
   private:
   const double _mass;
   const PotentialFunc _potential;
};

#endif // DYNAMICAL_SYSTEM_HPP




