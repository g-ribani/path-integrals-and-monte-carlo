#ifndef UTILITY_FUNCTIONS_HPP
#define UTILITY_FUNCTIONS_HPP

#include <vector>
#include <array>
#include <utility>   // std::pair
#include <map>
#include <string>
#include <typeinfo>
#include <boost/core/demangle.hpp>
#include <ostream>


// Small powers:
template<class T> constexpr T pow_2(T const x) { return x*x; }
template<class T> constexpr T pow_3(T const x) { return x*x*x; }
template<class T> constexpr T pow_4(T const x) { return x*x*x*x; }

// Corresponding macros:
#define POW_2(x) ((x)*(x))
#define POW_3(x) ((x)*(x)*(x))
#define POW_4(x) ((x)*(x)*(x)*(x))

// Returns a human-readable string identifying the type of the argument:
template<class T> inline std::string ClassName(T const &t) {
   return std::string( boost::core::demangle( typeid(t).name() ) );
}

// Obtains (by copy) the vector of keys from an std::map:
template<class K, class V> inline std::vector<K> GetKeys
 (std::map<K,V> const &m) {
   std::vector<K> k;
   for (auto p : m) k.push_back(std::get<0>(p));
   return k;
}

// Obtains (by copy) the vector of values from an std::map:
template<class K, class V> inline std::vector<V> GetValues
 (std::map<K,V> const &m) {
   std::vector<V> v;
   for (auto p : m) v.push_back(std::get<1>(p));
   return v;
}

// Returns a linearly spaced vector of nPoints of type Real from in to fin,
// ends included:
template<class Real> inline std::vector<Real> LinearRange
 (Real const in, Real const fin, size_t const n_points) {
      std::vector<Real> ret(n_points);
      if(n_points > 1) {
         Real const step = (fin - in)/(n_points - 1);
         for(size_t k = 0; k != n_points; ++k)
            ret[k] = in + k*step;
      }
      // midpoint if n_points == 1, empty if n_points == 0:
      else if(n_points == 1) ret[0] = (in + fin)/2.;
      return ret;
}

// Writes an std::array on an std::ostream:
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

// Writes an std::vector on an std::ostream:
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

// Writes an std::pair on an std::ostream:
template<class T1, class T2> inline std::ostream& operator <<
   (std::ostream& os, std::pair<T1,T2> const &p) {
      return os << '(' << std::get<0>(p) << ", " << std::get<1>(p) << ')';
}

// Writes an std::map on an std::ostream:
template<class K, class T> inline std::ostream& operator <<
   (std::ostream& os, std::map<K,T> const &m) {
      os << "{ ";
      for(auto const &p : m) os << p << ' ';
      return os << '}';
}

// todo:
// template<Types...> std::ostream& operator <<
// (std::ostream&, std::tuple<Types...> const &);

// Avoids compiler warnings on unused variables etc:
template<class...T> void Unused(T const &...) {}

#endif // UTILITY_FUNCTIONS_HPP




