#ifndef GSLMONTEVEGAS_HPP
#define GSLMONTEVEGAS_HPP

#include <vector>
#include <utility>
#include <random> // std::random_device
#include <gsl/gsl_monte_vegas.h>

// The struct GSLMonteVegas is a C++ wrapper for the gsl_monte_vegas routine,
// which can integrate a functor until the desired error goal is achieved.
// Notice however that the functor must accept a double*.
// Interface:
   // GSLMonteVegas(Functor const &f,
   //               std::vector<std::pair<double, double>> const &bounds,
   //               double const error_goal,
   //               gsl_rng_type const *gen_type = gsl_rng_mt19937);
   // GSLMonteVegas& Integrate(size_t const n_calls, int const stage = 0);
// The Integrate method iteratively calls gsl_monte_vegas_integrate with the
// user provided n_calls, until the error_goal is achieved. The chisquare is
// also ensured to be in between 0.5 and 1.5. The optional parameter sets the
// initial integration stage (see the GSL documentation), the default being a
// new clean integration. The return value is *this.
// The result, absolute error, chisquare and pointer to integrator state can be
// accessed via public member functions Result(), AbsErr(), ChiSquare(),
// State(). Params of the integrator can be got and set via the
// homonymous method.

double _GSLIntegrand(double*, size_t, void*);

template<class Functor> struct GSLMonteVegas {
   GSLMonteVegas(Functor const &f,
                 std::vector<std::pair<double, double>> const &bounds,
                 double const error_goal,
                 gsl_rng_type const *gen_type = gsl_rng_mt19937)
   : _dim(bounds.size()), _error_goal(error_goal), _func(f),
     _gen(gsl_rng_alloc(gen_type)),
     _state(gsl_monte_vegas_alloc(_dim)) {
      for(auto const &b : bounds) {
         _xl.push_back( std::get<0>(b) ),
         _xu.push_back( std::get<1>(b) );
      }
   }
   GSLMonteVegas(const GSLMonteVegas&) = delete;
   GSLMonteVegas& operator=(const GSLMonteVegas&) = delete;
   ~GSLMonteVegas() {
      gsl_rng_free(_gen);
      gsl_monte_vegas_free(_state);
   }
   friend double _GSLIntegrand(double *x, size_t, void *this_) {
      return ((GSLMonteVegas*)this_)->_func(x);
   }
   GSLMonteVegas& Integrate(size_t const n_calls, int const stage = 0) {
      gsl_monte_function I = {_GSLIntegrand, _dim, (void*)this};
      gsl_rng_set(_gen, std::random_device{}());
      auto par = Params();
      par.stage = stage, Params(par);
      if(stage == 0) {
         // warmup integration, to set the grid:
         gsl_monte_vegas_integrate (&I, _xl.data(), _xu.data(), _dim,
                                    n_calls/10, _gen, _state,
                                    &_result, &_abserr);
         // keep the grid, discard the estimate:
         par = Params(), par.stage = 1, Params(par);
      }
      double err{};
      do {
         int err_code = gsl_monte_vegas_integrate
                           (&I, _xl.data(), _xu.data(), _dim,
                            n_calls, _gen, _state,
                            &_result, &_abserr);
         if(err_code != 0) throw(err_code);
         _chisquare = gsl_monte_vegas_chisq(_state);
         err = std::abs(_abserr/_result);
         par = Params(), par.stage = 3, Params(par);
      }
      while( std::abs(_chisquare - 1.0) > 0.5 or err > _error_goal);
      return *this;
   }
   double AbsErr() const { return _abserr; }
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
   double _abserr, _chisquare, _result;
   size_t _dim;
   double _error_goal;
   Functor _func;
   gsl_rng *_gen;
   gsl_monte_vegas_state *_state;
   std::vector<double> _xl, _xu;
};

#endif   /* GSLMONTEVEGAS_HPP */




