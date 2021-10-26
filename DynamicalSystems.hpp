#ifndef DYNAMICAL_SYSTEMS_HPP
#define DYNAMICAL_SYSTEMS_HPP
#include <boost/core/demangle.hpp>
#include <cmath>
#include <exception>
#include <iostream>
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
   typedef std::pair<_CoordType, _ValueType> _DynPair;
   void AddToClassicalPath(_CoordType c) {
      IsDeterministic_Assert();
      _Path.insert_or_assign( c, operator()(c) );
   }
   void AddToPath(_CoordType c, _ValueType v) {
      _Path.insert_or_assign(c, v);
   }
   std::map<_CoordType, _ValueType> ClassicalPath
      (const std::vector<_CoordType>& C) {
      std::map<_CoordType, _ValueType> P;
      for(auto c : C) P.insert_or_assign(c, operator()(c) );
      return P;
   }
   std::string ClassName() {
      return std::string( boost::core::demangle( typeid(*this).name() ) );
   }
   void ClearPath() { if( HasPath() ) _Path.clear(); }
   virtual ~DynamicalSystem() = 0;
   bool HasPath() const { return !_Path.empty(); }
   virtual bool IsDeterministic() = 0;
   virtual bool IsExactlySolvable() = 0;
   virtual _ValueType operator() (_CoordType) = 0;
   void PrintPath
      (const std::string& s = "", std::ostream& o = std::cout) const {
         o << s << _Path << std::endl;
   }
   void SetClassicalPath (const std::vector<_CoordType>& C) {
      IsDeterministic_Assert();
      ClearPath();
      for(auto c : C) _Path.insert_or_assign( c, operator()(c) );
   }
   virtual void SetBoundaryCondition (_CoordType, _ValueType) = 0;
   void SetPath(const std::map<_CoordType, _ValueType>& P) {
      ClearPath();
      _Path = P;
   }
   protected:
   std::map<_CoordType, _ValueType> _Path;
   void IsDeterministic_Assert() {
      if(!IsDeterministic()) {
         std::string s("The motion of the object of type ");
         s.append( ClassName() )
          .append(" cannot be determined");
         throw DynamicalException(s);
      }
   }
};

template<class CT, class VT> DynamicalSystem<CT,VT>::~DynamicalSystem() {}

class EuclideanHarmonicOscillator : public DynamicalSystem<double, double> {
   public:
   double Action() {
   }
   EuclideanHarmonicOscillator(double omega = 1., double m = 1.)
      : _freq(omega), _mass(m) {
         if(m <= 0.) {
            std::string s;
            s.append("Cannot istantiate an object of type ")
             .append( ClassName() )
             .append(" with non positive mass");
            throw DynamicalException(s);
         }
         if(omega == 0.)
            std::clog << "Warning: istantiating an object of type "
                        << ClassName()
                        << " with zero frequency\n";
   }
   ~EuclideanHarmonicOscillator() {};
   double Frequency() const { return _freq; }
   bool IsDeterministic() { return _BCs.size() == 2; }
   bool IsExactlySolvable() { return true; }
   double Mass() const { return _mass; }
   _ValueType operator () (_CoordType tau) {
      IsDeterministic_Assert();
      auto t = GetKeys(_BCs);
      auto x = GetValues(_BCs);
      if(_freq != 0.) return ( std::sinh( _freq*(t[1] - tau) )*x[0]
                              + std::sinh( _freq*(tau - t[0]) )*x[1] )
                              / std::sinh( _freq*(t[1] - t[0]) );
      // free particle if _freq == 0
      else return ( (t[1] - tau)*x[0] + (tau - t[0])*x[1] )/( t[1] - t[0] );
   }
   void SetBoundaryCondition(_CoordType tau, _ValueType x) {
      _BCs.insert_or_assign(tau,x);
   }
   private:
   std::map<_CoordType, _ValueType> _BCs;
   const double _freq;
   const double _mass;
};

#endif // DYNAMICAL_SYSTEMS_HPP




