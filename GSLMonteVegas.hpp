#ifndef GSLMONTEVEGAS_HPP
#define GSLMONTEVEGAS_HPP

#include <vector>
#include <utility>
#include <random> // std::random_device
#include <gsl/gsl_monte_vegas.h>
#include <functional>
#include <string>
#include <stdexcept> // std::runtime_error


// The struct GSLMonteVegas is a C++ wrapper for the gsl_monte_vegas routine,
// which can integrate a functor until the desired error goal is achieved.
// Notice however that the functor must be convertible to an
// std::function<double(double*)>.
// Interface:
   // GSLMonteVegas(std::function<double(double*)> const &f,
   //               std::vector<std::pair<double, double>> const &bounds,
   //               gsl_rng_type const *gen_type = gsl_rng_mt19937)
   // GSLMonteVegas& operator() (size_t const n_calls,
   //                            double const error_goal,
   //                            size_t const max_iter = 1000,
   //                            int const stage = 0)
// The operator() iteratively calls gsl_monte_vegas_integrate with the user-
// provided n_calls, until the desired error_goal is achieved. The chisquare is
// also ensured to be in between 0.5 and 1.5. The optional parameter sets the
// initial integration stage (see the GSL documentation), the default being a
// new clean integration. The return value is *this. If the integration routine
// raises an error, the integer error code is thrown. If it is not able to reach
// the desired precision after max_iter iterations, an appropriate exception is
// thrown.
// The result, absolute error, chisquare and pointer to integrator state can be
// accessed via public member functions Result(), AbsError(), ChiSquare(),
// State(). Params of the integrator can be got and set via the
// homonymous method.

double GSLMonteVegas_Integrand(double*, size_t, void*);

struct GSLMonteVegas_Exception : std::runtime_error {
   using std::runtime_error::runtime_error;
};

struct GSLMonteVegas {
   GSLMonteVegas(std::function<double(double*)> const &f,
                 std::vector<std::pair<double, double>> const &bounds,
                 gsl_rng_type const *gen_type = gsl_rng_mt19937)
   : _dim(bounds.size()), _func(f),
     _gen(gsl_rng_alloc(gen_type)),
     _state(gsl_monte_vegas_alloc(_dim)) {
      for(auto const &b : bounds) {
         _xl.push_back( std::get<0>(b) ),
         _xu.push_back( std::get<1>(b) );
      }
   }
   GSLMonteVegas(const GSLMonteVegas&) = delete;
   GSLMonteVegas& operator= (const GSLMonteVegas&) = delete;
   ~GSLMonteVegas() {
      gsl_rng_free(_gen);
      gsl_monte_vegas_free(_state);
   }
   friend double GSLMonteVegas_Integrand
   (double* x, size_t, void* this_) {
      return (((GSLMonteVegas*)this_)->_func)(x);
   }
   GSLMonteVegas& operator() (size_t const n_calls,
                              double const error_goal,
                              size_t const max_iter = 1000,
                              int const stage = 0) {
      gsl_monte_function integrand
         = {GSLMonteVegas_Integrand, _dim, (void*)this};
      gsl_rng_set(_gen, std::random_device{}());
      auto par = Params();
      par.stage = stage, Params(par);
      if(stage == 0) {
         // Warmup integration, to set the grid:
         gsl_monte_vegas_integrate (&integrand, _xl.data(), _xu.data(), _dim,
                                    n_calls/10, _gen, _state,
                                    &_result, &_abserr);
         // Keep the grid, discard the estimate:
         par = Params(), par.stage = 1, Params(par);
      }
      double err = 0.;
      size_t iter = 0;
      do {
         if(iter >= max_iter) _OutOfTime(error_goal, iter);
         int err_code = gsl_monte_vegas_integrate
                           (&integrand, _xl.data(), _xu.data(), _dim,
                            n_calls, _gen, _state,
                            &_result, &_abserr);
         if(err_code != 0) throw(err_code);
         _chisquare = gsl_monte_vegas_chisq(_state);
         err = std::abs(_abserr/_result);
         par = Params(), par.stage = 3, Params(par);
         ++iter;
      }
      while(std::abs(_chisquare - 1.0) > 0.5 or err > error_goal);
      return *this;
   }
   double AbsError() const { return _abserr; }
   double ChiSquare() const { return _chisquare; }
   double Result() const { return _result; }
   gsl_monte_vegas_params Params() const {
      gsl_monte_vegas_params params;
      gsl_monte_vegas_params_get(_state, &params);
      return params;
   }
   void Params(gsl_monte_vegas_params const &params) {
      gsl_monte_vegas_params_set(_state, &params);
   }
   gsl_monte_vegas_state* State() const { return _state; }
   private:
   void _OutOfTime
   (double const error_goal, size_t const iter) const {
      std::string why = "could not achieve the precision of ";
      why += std::to_string(error_goal),
      why += " after ",
      why += std::to_string(iter),
      why += " iterations,\ncurrent estimate is ",
      why += std::to_string(_result),
      why += ", absolute error is ",
      why += std::to_string(_abserr),
      why += ", and chisquare is ",
      why += std::to_string(_chisquare);
      throw GSLMonteVegas_Exception(why);
   }
   double _abserr, _chisquare, _result;
   size_t _dim;
   std::function<double(double*)> _func;
   gsl_rng *const _gen;
   gsl_monte_vegas_state *const _state;
   std::vector<double> _xl, _xu;
};

#endif   /* GSLMONTEVEGAS_HPP */




