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
// which can integrate a functor convertible to an
// std::function<double(double*)> until the desired error goal is achieved.

double GSLMonteVegas_Integrand(double*, size_t, void*);

struct GSLMonteVegas_Exception : std::runtime_error {
   using std::runtime_error::runtime_error;
};

template<class Functor> struct GSLMonteVegas {
   GSLMonteVegas(Functor *f,
                 std::vector<std::pair<double, double>> const &bounds,
                 gsl_rng_type const *gen_type = gsl_rng_mt19937)
   :  _func(f), _gen(gsl_rng_alloc(gen_type)) {
         _dim = bounds.size(),
         _xl = new double[_dim],
         _xu = new double[_dim],
         _state = gsl_monte_vegas_alloc(_dim);
         std::pair<double, double> bound;
         for(size_t d = 0; d != _dim; ++d)
            bound = bounds[d],
            _xl[d] = std::get<0>(bound),
            _xu[d] = std::get<1>(bound);
   }
   GSLMonteVegas(const GSLMonteVegas&) = delete;
   GSLMonteVegas& operator= (const GSLMonteVegas&) = delete;
   ~GSLMonteVegas() {
      delete[] _xl,
      delete[] _xu,
      gsl_monte_vegas_free(_state),
      gsl_rng_free(_gen);
   }
   friend double GSLMonteVegas_Integrand(double* x, size_t, void* this_) {
      return (((GSLMonteVegas*)this_)->_integrand)(x);
   }
   GSLMonteVegas& operator() (size_t const n_calls,
                              double const error_goal,
                              size_t const max_iter = 1000,
                              int const stage = 0) {
      _integrand = *_func;
      gsl_monte_function integrand
         = {GSLMonteVegas_Integrand, _dim, (void*)this};
      gsl_rng_set(_gen, std::random_device{}());
      auto par = Params();
      par.stage = stage, Params(par);
      if(stage == 0) {
         // Warmup integration, to set the grid:
         gsl_monte_vegas_integrate (&integrand, _xl, _xu, _dim,
                                    n_calls/10, _gen, _state,
                                    &_result, &_abserr);
         // Keep the grid, discard the estimate:
         par = Params(), par.stage = 1, Params(par);
      }
      double err = 0.;
      int err_code = 0;
      size_t iter = 0;
      do {
         if(iter >= max_iter) _OutOfTime(error_goal, iter);
         err_code = gsl_monte_vegas_integrate
                    (&integrand, _xl, _xu, _dim,
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
   void _OutOfTime (double const error_goal, size_t const iter) const {
      std::string why = "could not achieve the precision of "
                        + std::to_string(error_goal)
                        + " after "
                        + std::to_string(iter)
                        + " iterations,\ncurrent estimate is "
                        + std::to_string(_result)
                        + ", absolute error is "
                        + std::to_string(_abserr)
                        + ", and chisquare is "
                        + std::to_string(_chisquare);
      throw GSLMonteVegas_Exception(why);
   }
   double _abserr, _chisquare, _result;
   size_t _dim;
   Functor * const _func;
   std::function<double(double*)> _integrand;
   gsl_rng *const _gen;
   gsl_monte_vegas_state *_state;
   double *_xl, *_xu;
};

#endif   /* GSLMONTEVEGAS_HPP */




