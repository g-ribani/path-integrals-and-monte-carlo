#ifndef UTILITY_FUNCTIONS_HPP
#define UTILITY_FUNCTIONS_HPP
#include <array>
#include <boost/core/demangle.hpp>
#include <map>
#include <ostream>
#include <string>
#include <utility>   // std::pair
#include <vector>

// returns a linearly spaced vector of nPoints of type Real from in to fin, included
template<class Real> inline std::vector<Real> LinearRange
 (const Real &in, const Real &fin, const std::size_t &nPoints) {
      std::vector<Real> ret;
      if(nPoints == 0); // empty if nPoints == 0
      else if(nPoints == 1) ret.push_back((in + fin)/2.); // midpoint if nPoints == 1
      else {
         const Real step = (fin - in)/(nPoints - 1);
         for(std::size_t k = 0; k != nPoints; ++k)
            ret.push_back(in+k*step);
      }
      return ret;
}

// writes an std::array on an std::ostream
template<class T, std::size_t N> inline std::ostream& operator <<
   (std::ostream& os, const std::array<T,N>& a) {
      auto it = a.begin();
      os << '(';
      if(it != a.end()) {
         for(; it != a.end()-1; ++it) os << *it << ", ";
         os << *it;
      }
      return os << ')';
}

// writes an std::vector on an std::ostream
template<class T> inline std::ostream& operator <<
   (std::ostream& os, const std::vector<T>& v) {
      auto it = v.begin();
      os << '(';
      if(it != v.end()) {
         for(; it != v.end()-1; ++it) os << *it << ", ";
         os << *it;
      }
      return os << ')';
}

// writes an std::pair on an std::ostream
template<class T1, class T2> inline std::ostream& operator <<
   (std::ostream& os, const std::pair<T1,T2>& p) {
      return os << '(' << std::get<0>(p) << ", " << std::get<1>(p) << ')';
}

// writes an std::map on an std::ostream
template<class K, class T> inline std::ostream& operator <<
   (std::ostream& os, const std::map<K,T>& m) {
      os << "{ ";
      for(auto p : m) os << p << ' ';
      return os << '}';
}

// template<Types...> std::ostream& operator << (std::ostream&, const std::tuple<Types...>&);   // implement this

// obtains (by copy) the vector of keys from an std::map
template<class K, class V> inline std::vector<K> GetKeys
 (const std::map<K,V>& m) {
   std::vector<K> k;
   for (auto p : m) k.push_back(std::get<0>(p));
   return k;
}

// obtains (by copy) the vector of values from an std::map
template<class K, class V> inline std::vector<V> GetValues
 (const std::map<K,V>& m) {
   std::vector<V> v;
   for (auto p : m) v.push_back(std::get<1>(p));
   return v;
}

// returns a human-readable string identifying the type of the argument
template<class T> inline std::string ClassName(const T& t) {
   return std::string( boost::core::demangle( typeid(t).name() ) );
}

#endif // UTILITY_FUNCTIONS_HPP




