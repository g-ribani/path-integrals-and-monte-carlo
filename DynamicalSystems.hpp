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

template<class C, class V> class DynamicalSystem {
   public:
   typedef C _CoordType;
   typedef V _ValueType;
   typedef std::pair<_CoordType, _ValueType> _DynPair;

   virtual ~DynamicalSystem() = 0;
   virtual bool IsExactlySolvable() = 0;
   virtual bool IsDeterministic() = 0;
   virtual _ValueType operator() (_CoordType) = 0;

   bool HasDiscretisedPath() const { return _discPath.empty(); }
   void SetDiscretisedPath(const std::vector<_CoordType>& T) {
      IsDeterministic_Assert();
      _discPath.clear();
      for(auto t : T) _discPath.insert_or_assign( t, operator()(t) );
   }
   void AddPointToDiscretisedPath(double t) {
      _discPath.insert_or_assign( t, operator()(t) );
   }
   void PrintDiscretisedPath
      (const std::string& s = "", std::ostream& o = std::cout) const {
         o << s << _discPath << std::endl;
   }

   protected:
   std::map<_CoordType, _ValueType> _discPath;

   void IsDeterministic_Assert() {
      if(!IsDeterministic()) {
         std::string s("The motion of the ");
         s.append( boost::core::demangle( typeid(*this).name() ) )
            .append(" cannot be determined");
         throw DynamicalException(s);
      }
   }
};

template<class C, class V> DynamicalSystem<C,V>::~DynamicalSystem() {}

class EuclideanHarmonicOscillator : public DynamicalSystem<double, double> {
   public:
   using _CoordType = typename EuclideanHarmonicOscillator::_CoordType;
   using _ValueType = typename EuclideanHarmonicOscillator::_ValueType;
   using _DynPair = typename EuclideanHarmonicOscillator::_DynPair;

   EuclideanHarmonicOscillator(double omega = 1., double m = 1.)
      : _freq(omega), _mass(m) {}
   ~EuclideanHarmonicOscillator() {};

   double Frequency() const { return _freq; }
   double Mass() const { return _mass; }

   bool IsExactlySolvable() { return true; }

   void SetCondition(double tau, double x) {
      _conds.insert_or_assign(tau,x);
   }

   bool IsDeterministic() { return _conds.size() == 2; }

   double operator () (double tau) {
      IsDeterministic_Assert();
      auto t = GetKeys(_conds);
      auto x = GetValues(_conds);
      return ( std::sinh( _freq*(t[1] - tau) )*x[0]
               + std::sinh( _freq*(tau - t[0]) )*x[1] )
               / std::sinh( _freq*(t[1] - t[0]) );
   }
// implement this if useful
   double ClassicalPathAction();

   double DiscretizedAction() {

   }

   private:
   const double _freq;
   const double _mass;
   std::map<_CoordType, _ValueType> _conds;
};

#endif // DYNAMICAL_SYSTEMS_HPP




