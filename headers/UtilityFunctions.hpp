#ifndef UTILITY_FUNCTIONS_HPP
#define UTILITY_FUNCTIONS_HPP

#include <vector>
#include <array>
#include <tuple>
#include <utility>   // std::pair, std::forward
#include <map>
#include <string>
#include <initializer_list>
#include <typeinfo>
#include <type_traits>  // std::invoke_result
#include <boost/core/demangle.hpp>
#include <ostream>
#include <stdexcept>


// Small powers:
template<class T> constexpr T pow_2(T const x) { return x*x; }
template<class T> constexpr T pow_3(T const x) { return x*x*x; }
template<class T> constexpr T pow_4(T const x) { return x*x*x*x; }

// Corresponding macros:
#define POW_2(x) ((x)*(x))
#define POW_3(x) ((x)*(x)*(x))
#define POW_4(x) ((x)*(x)*(x)*(x))

// Set of classes and a function which allocate a pointer of arbitary dimension
// to arbitrary type:
// template<class T, size_t dim> struct PointerTo {
//    typedef typename PointerTo<T, dim-1>::Type* Type;
// };
// template<class T> struct PointerTo<T, 0> {
//    typedef T Type;
// };
// template<class T> inline void* Allocate(std::initializer_list<size_t> dims) {
//    PointerTo<T, size(dims)> ret;
//    size_t k = 0;
//    for(auto d : dims);
//    return ret;
// };

// Returns a human-readable string identifying the type of the argument:
template<class T> inline std::string ClassName(T const &t) {
   return std::string( boost::core::demangle( typeid(t).name() ) );
}

template<class Vector> inline size_t FindMinIndex(Vector const& v) {
   size_t const size = std::size(v);
   if(size == 0)  throw std::length_error("zero size container passed "
                                          "to function FindMinIndex");
   size_t ret = 0;
   for(size_t k = 1; k != size; ++k) if(v[k] < v[ret]) ret = k;
   return ret;
}

template<class Vector> inline size_t FindMaxIndex(Vector const& v) {
   size_t const size = std::size(v);
   if(size == 0)  throw std::length_error("zero size container passed "
                                          "to function FindMaxIndex");
   size_t ret = 0;
   for(size_t k = 1; k != size; ++k) if(v[k] > v[ret]) ret = k;
   return ret;
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

// Avoids compiler warnings on unused variables etc:
template<class...T> inline void Ignore(T&&...) {}

template<class Real> inline Real* NewPointer(size_t const dim, Real const val) {
   Real *ret = new Real[dim];
   for(size_t i = 0; i != dim; ++i) ret[i] = val;
   return ret;
}

// Set of template classes and a function which invokes every functor
// in an std::tuple with the given arguments and returns the tuple of results:
template<size_t num, class FunctorTuple, class ResultTuple, class...Args>
struct InvokeAndStore {
   InvokeAndStore(FunctorTuple &funcs, ResultTuple &res, Args&&...args) {
      Ignore( InvokeAndStore<num-1, FunctorTuple, ResultTuple, Args...>
               (funcs, res, std::forward<Args>(args)...) );
      std::get<num-1>(res)
         = std::get<num-1>(funcs)(std::forward<Args>(args)...);
   }
};
template<class FunctorTuple, class ResultTuple, class...Args>
struct InvokeAndStore<1, FunctorTuple, ResultTuple, Args...> {
   InvokeAndStore(FunctorTuple &funcs, ResultTuple &res, Args&&...args) {
      std::get<0>(res) = std::get<0>(funcs)(std::forward<Args>(args)...);
   }
};
template<class FunctorTuple, class ResultTuple, class...Args>
struct InvokeAndStore<0, FunctorTuple, ResultTuple, Args...> {
   InvokeAndStore(FunctorTuple &funcs, ResultTuple &res, Args&&...args) {}
};
template<class...Functors, class...Args> inline auto Invoke
(std::tuple<Functors...> funcs, Args&&...args) noexcept {
   std::tuple<std::invoke_result_t<Functors, Args...>...> res;
   Ignore( InvokeAndStore<sizeof...(Functors),
                          decltype(funcs), decltype(res), Args...>
            (funcs, res, std::forward<Args>(args)...) );
   return res;
}

// Invokes every function in an std::vector with the given arguments
// and returns the vector of results of type Res:
template<class Res, class...Args> inline std::vector<Res> Invoke
(std::vector<std::function<Res(Args...)>> funcs, Args...args) {
   size_t const size = funcs.size();
   std::vector<Res> res(size);
   for(size_t k = 0; k != size; ++k)
      res[k] = funcs[k](std::forward<Args>(args)...);
   return res;
}

// Invokes every function in an std::array with the given arguments
// and returns the array of results of type Res:
template<class Res, size_t N, class...Args> inline std::array<Res, N> Invoke
(std::array<std::function<Res(Args...)>, N> funcs, Args...args) {
   std::array<Res, N> res;
   for(size_t k = 0; k != N; ++k)
      res[k] = funcs[k](std::forward<Args>(args)...);
   return res;
}

// Returns a linearly spaced std::vector of nPoints of type Real from in to fin,
// ends included:
template<class Real> inline std::vector<Real> LinearRange
(Real const in, Real const fin, size_t const n_points) {
   std::vector<Real> ret(n_points);
   if(n_points > 1) {
      Real const step = (fin - in)/(n_points - 1);
      for(size_t k = 0; k != n_points; ++k)
         ret[k] = in + k*step;
   }
   // Midpoint if n_points == 1, empty if n_points == 0:
   else if(n_points == 1) ret[0] = 0.5*(in + fin);
   return ret;
}

// Writes an std::array on an std::ostream:
template<class T, size_t N> inline std::ostream& operator <<
(std::ostream &os, std::array<T,N> const &a) {
   auto it = a.begin();
   os << '[';
   if(it != a.end()) {
      for(; it != a.end()-1; ++it) os << *it << ", ";
      os << *it;
   }
   return os << ']';
}

// Writes an std::vector on an std::ostream:
template<class T> inline std::ostream& operator <<
(std::ostream &os, std::vector<T> const &v) {
   auto it = v.begin();
   os << '[';
   if(it != v.end()) {
      for(; it != v.end()-1; ++it) os << *it << ", ";
      os << *it;
   }
   return os << ']';
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

// Set of template classes and the operator which writes
// an std::tuple on an std::ostream:
template<size_t num, class...Types> struct PrintTupleElements {
   PrintTupleElements(std::tuple<Types...> const &tup, std::ostream &os) {
      Ignore( PrintTupleElements<num-1, Types...>(tup, os) );
      os << std::get<num-1>(tup);
      if(num < sizeof...(Types)) os << ", ";
   }
};
template<class...Types> struct PrintTupleElements<1, Types...> {
   PrintTupleElements(std::tuple<Types...> const &tup, std::ostream &os) {
      os << std::get<0>(tup);
      if(1 < sizeof...(Types)) os << ", ";
   }
};
template<class...Types> struct PrintTupleElements<0, Types...> {
   PrintTupleElements(std::tuple<Types...> const&, std::ostream&) {}
};
template<class...Types> inline std::ostream& operator <<
(std::ostream &os, std::tuple<Types...> const &tup) {
   os << '(';
   Ignore( PrintTupleElements<sizeof...(Types), Types...>(tup, os) );
   return os << ')';
}

// Returns a dynamically-allocated pointer to nPoints linearly-spaced values
// of type Real from in to fin, ends included:
template<class Real> inline Real* NewPointerToLinearRange
(size_t const n_points, Real const in, Real const fin) {
   Real *ret = new Real[n_points];
   if(n_points > 1) {
      Real const step = (fin - in)/(n_points - 1);
      for(size_t k = 0; k != n_points; ++k)
         ret[k] = in + k*step;
   }
   // Midpoint if n_points == 1, empty if n_points == 0:
   else if(n_points == 1) ret[0] = 0.5*(in + fin);
   return ret;
}

// Returns the vector of proper divisors of the argument.
template<class Integer> inline std::vector<Integer> ProperDivisors
(Integer const n) {
   if(n <= 0) throw std::invalid_argument("non-positive integer passed "
                                          "to function ProperDivisors");
   std::vector<Integer> divs;
   for(Integer d = 2; d <= n/2; ++d)
      if(n % d == 0) divs.push_back(d);
   return divs;
}

#endif // UTILITY_FUNCTIONS_HPP




