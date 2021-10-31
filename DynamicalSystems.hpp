#ifndef DYNAMICAL_SYSTEMS_HPP
#define DYNAMICAL_SYSTEMS_HPP
#include <cmath>
#include <boost/core/demangle.hpp>
#include <exception>
#include <iostream>
#include <random>
#include <string>
#include <UtilityFunctions.hpp>

class DynamicalException : public std::exception {
   public:
   DynamicalException() {};
   DynamicalException(const std::string& s) : _explain(s) {}
   DynamicalException(const DynamicalException&) = default;
   DynamicalException& operator=(const DynamicalException&) = default;
   DynamicalException(DynamicalException&&) = default;
   DynamicalException& operator=(DynamicalException&&) = default;
   ~DynamicalException() {};
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

template<class CT, class VT> class ExactDynamicalSystem :
 public DynamicalSystem<CT, VT> {
   public:
   using typename DynamicalSystem<CT, VT>::_CoordType;
   using typename DynamicalSystem<CT, VT>::_ValueType;
   void AddToClassicalPath(_CoordType c) {
      this->_Path.insert_or_assign( c, operator()(c) );
   }
   std::map<_CoordType, _ValueType> ClassicalPath
    (const std::vector<_CoordType>& C) const {
      IsDeterministic_Assert();
      std::map<_CoordType, _ValueType> P;
      for(auto c : C) P.insert_or_assign(c, operator()(c) );
      return P;
   }
   auto GetBoundaryConditions() { return _BCs; }
   // checks if boundary conditions uniquely determine the classical path
   virtual bool IsDeterministic() const = 0;
   // computes the classical path
   virtual _ValueType operator()(_CoordType) const = 0;
   void SetBoundaryCondition(_CoordType c, _ValueType v) {
      _BCs.insert_or_assign(c, v);
   }
   void SetClassicalPath(const std::vector<_CoordType>& C) {
      SetPath( ClassicalPath(C) );
   }
   protected:
   std::map<_CoordType, _ValueType> _BCs;
   void IsDeterministic_Assert() const {
      if(!IsDeterministic()) {
         std::string s("The motion of the object of type ");
         s.append( ClassName(*this) )
          .append(" cannot be determined");
         throw DynamicalException(s);
      }
   }
};

class EuclideanHarmonicOscillator :
 public ExactDynamicalSystem<double, double> {
   public:
   double Action() const {
      if(_Path.size() < 2) {
         std::string s(ClassName(*this));
         s.append(" does not have a discretized path "
                  "suitable to compute the action");
         throw DynamicalException(s);
      }
      auto t = GetKeys(_Path);
      auto x = GetValues(_Path);
      // x[0] = x_in, x[N] = x_fin
      std::size_t N = _Path.size()-1;
      double deltaT;
      double kineticPiece, potentialPiece, action = 0.;
      for(std::size_t k = 1; k != N+1; ++k)
         deltaT = t[k] - t[k-1],
         kineticPiece =
            _mass/(2.*deltaT)*std::pow(x[k] - x[k-1], 2.),
         potentialPiece = deltaT*_mass*_freq*_freq/4.
                           *( x[k]*x[k] + x[k-1]*x[k-1] ),
         action += kineticPiece + potentialPiece;
      return action;
   }
   EuclideanHarmonicOscillator(double omega = 1., double m = 1.)
    : _freq(omega), _mass(m) {
      if(m <= 0.) {
         std::string s;
         s.append("Cannot istantiate an object of type ")
          .append( ClassName(*this) )
          .append(" with non positive mass");
         throw DynamicalException(s);
      }
      if(omega == 0.)
         std::clog << "Warning: istantiating an object of type "
                     << ClassName(*this)
                     << " with zero frequency\n";
   }
   ~EuclideanHarmonicOscillator() {};
   double ExactAmplitude(double xin, double xfin, double beta) const {
      static constexpr double pi = std::acos(-1.);
      return std::sqrt( _mass/(2.*pi*beta) )
               *std::exp( -_mass*std::pow( xfin - xin, 2. )/(2.*beta) );
   }
   double Frequency() const { return _freq; }
   bool IsDeterministic() const { return _BCs.size() == 2; }
   // bool IsExactlySolvable() const { return true; }
   double Mass() const { return _mass; }
   double operator()(double tau) const {
      IsDeterministic_Assert();
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      if(_freq != 0.) return ( std::sinh( _freq*(t[1] - tau) )*x[0]
                              + std::sinh( _freq*(tau - t[0]) )*x[1] )
                              / std::sinh( _freq*(t[1] - t[0]) );
      // free particle if _freq == 0
      else return ( (t[1] - tau)*x[0] + (tau - t[0])*x[1] )/( t[1] - t[0] );
   }
   private:
   const double _freq;
   const double _mass;
};

#endif // DYNAMICAL_SYSTEMS_HPP




