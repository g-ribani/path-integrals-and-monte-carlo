#include <boost/math/quadrature/naive_monte_carlo.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <UtilityFunctions.hpp>
#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif
#ifndef N_THREADS
   #define N_THREADS 1
#endif

void PrintProgress (const double &progress,
                    const double &error_estimate,
                    const double &current_estimate,
                    const std::chrono::duration<double>
                        &estimated_time_to_completion,
                    std::ofstream &out) {
   std::cout << "["
           << std::setw(3) << int(progress * 100.)
           << "%], estimate = "
           << std::setprecision(5)
           << current_estimate << " +/- "
           << error_estimate*current_estimate
           << ", error = "
           << std::setprecision(3)
           << error_estimate
           << ", time to completion = "
           << estimated_time_to_completion.count()
           << " seconds\r" << std::flush;
   #ifdef PRINT_ON_FILE
      out << std::setprecision(5) << current_estimate << "\t\t"
          << error_estimate << '\n' << std::flush;
   #endif
}

struct EuclidHarmonicOscillator1D {
   EuclidHarmonicOscillator1D(double const &m, double const &f)
    : frequency(f), mass(m) {}
   std::array<double, 2> initial_pos;
   std::array<double, 2> final_pos;
   double const frequency;
   double const mass;
   EuclidHarmonicOscillator1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      step = (final_pos[1] - initial_pos[1])/N;
      return *this;
   }
   protected:
   double n_steps;
   double step;
};

template<class T> struct PathIntegrand;

template<> struct PathIntegrand<EuclidHarmonicOscillator1D>
 : EuclidHarmonicOscillator1D {
   double operator() (std::vector<double> const &x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = x.size();
      std::vector<double> deltaX(dim + 1);
      deltaX[0] = x[0] - initial_pos[1];
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = final_pos[1] - x[dim-1];
      double ret = 1.;
      for(size_type k = 0; k != dim + 1; ++k) {
         double deltaS = mass/2.*deltaX[k]*deltaX[k]/step
                         + step*mass*frequency*frequency/4.
                           *( x[k]*x[k] + x[k-1]*x[k-1] );
         ret *= std::sqrt(mass/(2.*M_PI*step)) * std::exp(-deltaS);
      }
      return ret;
   }
   PathIntegrand(const EuclidHarmonicOscillator1D &osci)
      : EuclidHarmonicOscillator1D(osci) {}
};

double NaiveMCAmplitude (EuclidHarmonicOscillator1D& osci,
                         const std::size_t n_steps,
                         const std::pair<double, double> range,
                         const double error_goal,
                         std::ofstream &out) {
   osci.SetNSteps(n_steps);
   std::vector<std::pair<double, double>> bounds(n_steps-1, range);
   auto F = PathIntegrand<EuclidHarmonicOscillator1D>(osci);
   boost::math::quadrature::naive_monte_carlo<double, decltype(F)>
    mc(F, bounds, error_goal, /*singular = */ false, N_THREADS);
   std::future<double> task = mc.integrate();
   while
    (task.wait_for(std::chrono::seconds(1)) != std::future_status::ready)
      PrintProgress(mc.progress(),
                    mc.current_error_estimate(),
                    mc.current_estimate(),
                    mc.estimated_time_to_completion(),
                    out);
   return task.get();
}


int main() {
   std::cout << std::scientific;

   EuclidHarmonicOscillator1D osci(1., 1.);
   osci.initial_pos = {0., 0.};
   osci.final_pos = {1., 1.};

   std::ifstream input("naive.in");
   double xmin, xmax, error_goal;
   size_t n_steps;
   input >> n_steps >> xmin >> xmax >> error_goal;
   auto min_max = std::make_pair(xmin, xmax);
   input.clear(), input.close();

   std::ofstream output("naive.out", std::ios::app);
   #ifdef PRINT_ON_FILE
      output << "\n\n\t************\t\n\n"
            << "(xmin, xmax) = " << min_max << '\n'
            << "nsteps = " << n_steps << "\n\n"
            << "estimate\t\terror estimate\n";
   #endif
   double amp =
    NaiveMCAmplitude(osci, n_steps, min_max, error_goal, output);
   std::cout << "(xmin, xmax) = " << min_max << '\n'
             << "nsteps = " << n_steps << std::endl
             << "Naive MC Amplitude: " << amp << " +/- " << amp*error_goal
             << std::endl;
   #ifdef PRINT_ON_FILE
      output << "Naive MC Amplitude: " << amp << " +/- " << amp*error_goal
             << std::endl;
   #endif
   output.clear(), output.close();
}




