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
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


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
   typedef CT CoordType;
   typedef VT ValueType;
   void AddBoundaryCondition(const CoordType& t, const ValueType& x) {
      _BCs.insert_or_assign(t, x);
   }
   void AddBoundaryCondition(const std::pair<CoordType, ValueType>& point) {
      AddBoundaryCondition(std::get<0>(point), std::get<1>(point));
   }
   void AddToPath(const CoordType& t, const ValueType& x) {
      _Path.insert_or_assign(t, x);
   }
   void AddToPath(const std::pair<CoordType, ValueType>& point) {
      AddToPath(std::get<0>(point), std::get<1>(point));
   }
   void AddToPath(const std::map<CoordType, ValueType>& subpath) {
      for(const auto& point : subpath) AddToPath(point);
   }
   std::map<CoordType, ValueType> BoundaryConditions() const { return _BCs; }
   virtual ValueType ClassicalValue (CoordType t, double epsi) const {
      Throw("don't know how to solve this system");
      return ValueType();
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   void ClearPath() { _Path.clear(); }
   CoordType& Coord
    (typename std::map<CoordType, ValueType>::iterator it) const {
      return std::get<0>(*it);
   }
   virtual ~DynamicalSystem() {}
   virtual bool IsSolvable() const { return false; }
   std::map<CoordType, ValueType> Path() const { return _Path; }
   // typename std::map<CoordType, ValueType>::iterator PathBegin() {
   //    return _Path.begin();
   // }
   std::vector<CoordType> PathCoords() const { return GetKeys(_Path); }
   // typename std::map<CoordType, ValueType>::iterator PathEnd() {
   //    return _Path.end();
   // }
   typename std::map<CoordType, ValueType>::size_type PathSize() const {
      return _Path.size();
   }
   std::vector<ValueType> PathValues() const { return GetValues(_Path); }
   void PrintPath(std::ostream& os = std::cout) const { os << _Path << '\n'; }
   void SetClassicalPath
    (const std::vector<CoordType>& coords, double epsi) {
      ClearPath();
      for(const auto& t : coords) AddToPath( t, ClassicalValue(t, epsi) );
   }
   void SetPath(const std::map<CoordType, ValueType>& P) { _Path = P; }
   template<class Distribution, class Generator> void SetRandomPath
    (const std::vector<CoordType>& coords,
     Distribution& distr, Generator& gen) {
      auto func = [&distr, &gen](){ return distr(gen); };
      SetUserPath(coords, func);
   }
   template<class Function> void SetUserPath
    (const std::vector<CoordType>& coords, Function& func) {
      ClearPath();
      for(const auto& t : coords) AddToPath(t, func());
   }
   virtual void Throw(const char* why) const {
      throw DynamicalException(this, why);
   }
   // ValueType& Value
   //  (typename std::map<CoordType, ValueType>::iterator it) const {
   //    return std::get<1>(*it);
   // }
   protected:
   std::map<CoordType, ValueType> _BCs, _Path;
};

// specialization for CoordType = ValueType = double
template<> class DynamicalSystem<double, double> {
   public:
   typedef double CoordType, ValueType;
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
      Throw("in method ClassicalValue: don't know how to solve this system");
      return 0.;
   }
   void ClearBoundaryConditions() { _BCs.clear(); }
   void ClearPath() { _Path.clear(); }
   // double& Coord (typename std::map<double, double>::iterator it) const {
   //    return std::get<0>(*it);
   // }
   virtual ~DynamicalSystem() {}
   virtual bool IsSolvable() const { return false; }
   std::map<double, double> Path() const { return _Path; }
   // typename std::map<double, double>::iterator PathBegin() {
   //    return _Path.begin();
   // }
   std::vector<double> PathCoords() const { return GetKeys(_Path); }
   // typename std::map<double, double>::iterator PathEnd() {
   //    return _Path.end();
   // }
   typename std::map<double, double>::size_type PathSize() const {
      return _Path.size();
   }
   std::vector<double> PathValues() const { return GetValues(_Path); }
   void PrintPath(std::ostream& os = std::cout) const { os << _Path << '\n'; }
   void SetPath(const std::map<double, double>& P) { _Path = P; }
   void SetClassicalPath(std::size_t nPoints, double epsi) {
      if( _BCs.size() != 2 )
         Throw("in method SetClassicalPath: "
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
         Throw("in method SetRandomPath: "
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
         Throw("in method SetUserPath: "
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
   // double& Value (typename std::map<double, double>::iterator it) const {
   //    return std::get<1>(*it);
   // }
   protected:
   std::map<double, double> _BCs, _Path;
};

class Basic_Particle1D : public DynamicalSystem<double, double> {
   public:
   double Mass() const { return _mass; }
   Basic_Particle1D(double mass) : _mass(mass) {
      if(mass <= 0.) Throw("in constructor: must have positive mass");
   }
   void Throw(const char* why) const override {
      throw DynamicalException(this, why);
   }
   protected:
   const double _mass;
};

// euclidean one-dimensional free particle:
class EuclidFreeParticle1D : public Basic_Particle1D {
   public:
   double ClassicalAction() const {
      if( _BCs.size() != 2 )
         Throw("in method ClassicalAction: "
               "need two boundary conditions");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return _mass/2.* (x[1] - x[0]) / (t[1] - t[0]);
   }
   double ClassicalValue
    (double tau, double epsi = 0.) const override {
      if(!IsSolvable())
         Throw("in method ClassicalValue: cannot solve the classical motion");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return ( (t[1] - tau)*x[0] + (tau - t[0])*x[1] )/( t[1] - t[0] );
   }
   using Basic_Particle1D::Basic_Particle1D;
   double ExactAmplitude() const {
      if( _BCs.size() != 2 )
         Throw("in method ExactAmplitude: "
               "need two boundary conditions");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0],
             deltaX = x[1] - x[0];
      return std::sqrt( _mass/(2.*M_PI*deltaT) )
         *std::exp( -_mass*deltaX*deltaX/(2.*deltaT) );
   }
   bool IsSolvable() const override { return _BCs.size() == 2; }
   void Throw(const char* why) const override {
      throw DynamicalException(this, why);
   }
};

// euclidean one-dimensional harmonic oscillator:
class EuclidHarmonicOscillator1D : public Basic_Particle1D {
   public:
   double ClassicalAction() const {
      if( _BCs.size() != 2 )
         Throw("in method ClassicalAction: "
               "need two boundary conditions");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0];
      return _mass*_freq/2.
                       *( (x[0]*x[0] + x[1]*x[1])/std::tanh(_freq*deltaT)
                        - 2.*x[0]*x[1]/std::sinh(_freq*deltaT) );
   }
   double ClassicalValue
    (double tau, double epsi = 0.) const override {
      if(!IsSolvable())
         Throw("in method ClassicalValue: cannot solve the classical motion");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      return ( std::sinh( _freq*(t[1] - tau) )*x[0]
                           + std::sinh( _freq*(tau - t[0]) )*x[1] )
                           / std::sinh( _freq*(t[1] - t[0]) );
   }
   EuclidHarmonicOscillator1D(double mass, double freq)
    :  Basic_Particle1D(mass), _freq(freq) {
      if(mass <= 0. or freq <= 0.)
         Throw("in constructor: must have positive mass "
               "and positive frequency");
   }
   double ExactAmplitude() const {
      if( _BCs.size() != 2 )
         Throw("in method ExactAmplitude: "
               "need two boundary conditions");
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      double deltaT = t[1] - t[0];
      return std::sqrt( _mass/(2.*M_PI*deltaT) )
         *std::exp( -_mass*_freq/2.
                     *( (x[0]*x[0] + x[1]*x[1])/std::tanh(_freq*deltaT)
                        - 2.*x[0]*x[1]/std::sinh(_freq*deltaT) ) );
   }
   double Frequency() const { return _freq; }
   bool IsSolvable() const override { return _BCs.size() == 2; }
   void Throw(const char* why) const override {
      throw DynamicalException(this, why);
   }
   private:
   const double _freq;
};

// euclidean one-dimensional particle in a potential
template<class PotentialFunc> class EuclidParticle1D :
 public Basic_Particle1D {
   public:
   EuclidParticle1D(double mass, const PotentialFunc& pot) :
     Basic_Particle1D(mass), _potential(pot) {
            if(mass <= 0.) Throw("in constructor: must have positive mass");
   }
   auto Potential(double x) const { return _potential(x); }
   void Throw(const char* why) const override {
      throw DynamicalException(this, why);
   }
   private:
   const PotentialFunc _potential;
};

#endif // DYNAMICAL_SYSTEM_HPP




