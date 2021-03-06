#include <cstdlib>   // std::system
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include "headers/GSLMonteVegas.hpp"
#include "headers/EuclidParticle1d.hpp"
#include <libIntegrate/Integrate.hpp>

// TASK:
//    Compute the energy and wavefunction of the ground state for the
//    anharmonic oscillator using the Vegas algortihm.

// USAGE:
//    For plotting, provide a log file name on the command line.


int main(int narg, char const **args) {

   using std::cout, std::clog, std::cerr;
   using namespace boost::math::double_constants;
   cout << std::setprecision(6) << std::fixed;

   // Setup:
   double mass = 1.,
          frequency = 1.,
          time_step =  0.5;
   auto potential = [&mass, &frequency](double x) {
      return mass*POW_2(frequency*x)/2. + POW_4(x)/2.;
   };
   size_t n_steps = 8;
   double propa_time = n_steps*time_step,
          delta_x = 5.;
          //^ Paths based at x can span the interval ]x-delta_x, x+delta_x[.
   double err = 1e-2;
   size_t n_calls = 100'000;
   EuclidParticle1d osci(mass, potential);

   // Set output stream:
   std::ostream *log_stream;
   std::ofstream file;
   std::string file_name;
   if(narg < 2) log_stream = &cout;
   else {
      file_name = args[1];
      file.open(file_name, std::ios_base::app);
      if(!file) {
      cerr << "Could not open file " << file_name << "... "
           << "using stdout for logging instead. No plotting.\n\n"
           << std::flush;
      log_stream = &cout;
      }
      else log_stream = &file;
   }
   cout << "mass = " << mass << '\n'
        << "frequency = " << frequency << '\n'
        << "potential = harmonic + x^4/2\n"
        << "time_step = " << time_step << '\n'
        << "n_steps = " << n_steps << "\n\n" << std::flush;
   if(file) {
      file  << "mass = " << mass << '\n'
            << "frequency = " << frequency << '\n'
            << "potential = harmonic + x^4/2\n"
            << "time_step = " << time_step << '\n'
            << "n_steps = " << n_steps << "\n\n" << std::flush;
   }

   // Compute the amplitude on several closed paths:
   double const x_begin = 0.,
                x_end = 2.;
   size_t const n_points = 11;
   std::vector<double> xs = LinearRange(x_begin, x_end, n_points),
                       amps(n_points);
   for(size_t k = 0; k != n_points; ++k) {
      double x = xs[k];
      osci.Initial({0., x}),
      osci.Final({propa_time, x});   // closed path
      std::vector<std::pair<double, double>> bounds
         (n_steps-1, {x-delta_x, x+delta_x});
      osci.NSteps(n_steps);
      PathIntegrand<EuclidParticle1d> path_integrand(&osci);
      GSLMonteVegas vegas(&path_integrand, bounds);
      try { vegas(n_calls, err); }
      catch(int err) {
         cout << "Vegas integration routine returned with error code "
              << err << "\n\n" << std::flush;
         continue;
      }
      double result = vegas.Result(),
             abserr = vegas.AbsError(),
             chisquare = vegas.ChiSquare();
      amps[k] = result;
      cout << "Vegas(x = " << x << ") = " << result << " +/- " << 3.*abserr
           << ", err = " << abserr/result
           << ", chi square per dof = " << chisquare << "\n\n" << std::flush;
   }


   // Extract the ground state energy:
   _1D::SimpsonRule<double> simpson;
   double integral = 2.*simpson(xs, amps);
   double vacuum_energy = -std::log(integral)/propa_time;
   cout << "Vacuum energy = " << vacuum_energy << "\n\n" << std::flush;
   if(file) file << "Vacuum energy = " << vacuum_energy
                  << "\n\n" << std::flush;


   // Extract the ground state wavefunction:
   std::vector<double> ground_wf(n_points);
   *log_stream << "Ground state wavefunction:\n";
   for(size_t k = 0; k != n_points; ++k) {
      ground_wf[k] = std::sqrt(amps[k]/integral);
      *log_stream << xs[k] << '\t' << ground_wf[k] << '\n';
   }
   *log_stream << std::endl;

   // Plot via python script, if log file was successfully used:
   if(file) {
      file.clear(), file.close();
      std::string program("python anharmonic.py ");
      program.append(file_name).append(" &");
      Ignore( std::system(program.data()) );   // return value ignored
   }
   else file.clear(), file.close();
}




