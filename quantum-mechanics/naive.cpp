#include <boost/math/quadrature/naive_monte_carlo.hpp>
#include <iostream>
#include <fstream>
#include "QMPathIntegrals.hpp"
#include <UtilityFunctions.hpp>
#ifndef N_THREADS
   #define N_THREADS 1
#endif

// PrintProgress function taken from a Boost example
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
   #ifdef PRINT_TO_FILE
      out << std::setprecision(5) << current_estimate << "\t\t"
          << error_estimate << '\n' << std::flush;
   #endif
}

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
   while    // wait and print cycle taken from a Boost example
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
   double exact_amp = osci.ExactAmplitude();

   std::ifstream input("naive.in");
   double xmin, xmax, error_goal;
   size_t n_steps;
   input >> n_steps >> xmin >> xmax >> error_goal;
   auto min_max = std::make_pair(xmin, xmax);
   input.clear(), input.close();

   std::ofstream output("naive.out", std::ios::app);
   #ifdef PRINT_TO_FILE
      output << "\n\n\t************\t\n\n"
            << "(xmin, xmax) = " << min_max << '\n'
            << "nsteps = " << n_steps << '\n'
            << "exact amplitude = " << exact_amp << "\n\n"
            << "estimate\t\terror estimate\n";
   #endif
   double mc_amp = NaiveMCAmplitude(osci, n_steps, min_max, error_goal, output);
   std::cout << "(xmin, xmax) = " << min_max << '\n'
             << "nsteps = " << n_steps << '\n'
             << "exact amplitude = " << exact_amp << std::endl
             << "Naive MC Amplitude: " << mc_amp << " +/- " << mc_amp*error_goal
             << std::endl;
   #ifdef PRINT_TO_FILE
      output << "Naive MC Amplitude: " << mc_amp << " +/- " << mc_amp*error_goal
             << std::endl;
   #endif
   output.clear(), output.close();
}




