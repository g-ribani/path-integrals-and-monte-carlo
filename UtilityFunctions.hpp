#ifndef UTILITY_FUNCTIONS_HPP
#define UTILITY_FUNCTIONS_HPP

#include <cstdlib>
#include <string>
#include <ostream>
#include <vector>
#include <array>
#include <boost/core/demangle.hpp>
#include <map>
#include <utility>   // std::pair


// Powers with small integer exponents.
template<class T> constexpr T pow_2(T const x) {
   return x*x;
}
template<class T> constexpr T pow_3(T const x) {
   return x*x*x;
}
template<class T> constexpr T pow_4(T const x) {
   return x*x*x*x;
}

// Returns a human-readable string identifying the type of the argument.
template<class T> inline std::string ClassName(T const &t) {
   return std::string( boost::core::demangle( typeid(t).name() ) );
}

// Obtains (by copy) the vector of keys from an std::map
template<class K, class V> inline std::vector<K> GetKeys
 (std::map<K,V> const &m) {
   std::vector<K> k;
   for (auto p : m) k.push_back(std::get<0>(p));
   return k;
}

// Obtains (by copy) the vector of values from an std::map.
template<class K, class V> inline std::vector<V> GetValues
 (std::map<K,V> const &m) {
   std::vector<V> v;
   for (auto p : m) v.push_back(std::get<1>(p));
   return v;
}

// Returns a linearly spaced vector of nPoints of type Real from in to fin,
// ends included.
template<class Real> inline std::vector<Real> LinearRange
 (Real const in, Real const fin, size_t const nPoints) {
      std::vector<Real> ret;
      // empty if nPoints == 0:
      if(nPoints == 0);
      // midpoint if nPoints == 1:
      else if(nPoints == 1) ret.push_back((in + fin)/2.);
      else {
         const Real step = (fin - in)/(nPoints - 1);
         for(size_t k = 0; k != nPoints; ++k)
            ret.push_back(in+k*step);
      }
      return ret;
}

// Writes an std::array on an std::ostream.
template<class T, size_t N> inline std::ostream& operator <<
   (std::ostream &os, std::array<T,N> const &a) {
      auto it = a.begin();
      os << '(';
      if(it != a.end()) {
         for(; it != a.end()-1; ++it) os << *it << ", ";
         os << *it;
      }
      return os << ')';
}

// Writes an std::vector on an std::ostream.
template<class T> inline std::ostream& operator <<
   (std::ostream &os, std::vector<T> const &v) {
      auto it = v.begin();
      os << '(';
      if(it != v.end()) {
         for(; it != v.end()-1; ++it) os << *it << ", ";
         os << *it;
      }
      return os << ')';
}

// Writes an std::pair on an std::ostream.
template<class T1, class T2> inline std::ostream& operator <<
   (std::ostream& os, std::pair<T1,T2> const &p) {
      return os << '(' << std::get<0>(p) << ", " << std::get<1>(p) << ')';
}

// Writes an std::map on an std::ostream.
template<class K, class T> inline std::ostream& operator <<
   (std::ostream& os, std::map<K,T> const &m) {
      os << "{ ";
      for(auto p : m) os << p << ' ';
      return os << '}';
}

// todo:
// template<Types...> std::ostream& operator <<
// (std::ostream&, std::tuple<Types...> const &);

// Avoids compiler warnings on unused variables etc.
template<class...T> void Unused(T const &...) {}

#endif // UTILITY_FUNCTIONS_HPP




