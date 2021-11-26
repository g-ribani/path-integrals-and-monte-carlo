#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "GSLMonteVegas.hpp"
#include "PathIntegrals.hpp"
#include <libIntegrate/Integrate.hpp>
#include "UtilityFunctions.hpp"

// TASK:
//    Compute the energy and wavefunction of the ground state for the
//    harmonic oscillator using the Vegas algortihm.

// USAGE:
//    For plotting, provide a log file name on the command line
//    (it will be automatically put in the log/ subdirectory).
//    To get less verbose printings on the screen, redirect stderr.

int main(int narg, char const **args) {
   using std::cout, std::clog, std::cerr;
   using namespace boost::math::double_constants;
   cout << std::setprecision(6) << std::fixed;

   // Setup:
   double mass = 1.,
          frequency = 1.,
          time_step =  0.5;
   size_t n_steps = 8;
   double propa_time = n_steps*time_step,
          delta_x = 5.;
          //^ Paths based at x can span the interval ]x-delta_x, x+delta_x[.
   double err = 1e-2;
   size_t n_calls = 100'000;
   EuclidHarmonicOscillator1D osci(mass, frequency);

   // Set output stream:
   std::ostream *os;
   std::ofstream file;
   std::string file_name;
   if(narg < 2) os = &cout;
   else {
      file_name.append("log/").append(args[1]);
      file.open(file_name, std::ios_base::app);
      if(!file) {
      cerr << "Could not open file " << file_name << "... "
           << "using stdout for logging instead. No plotting.\n\n"
           << std::flush;
      os = &cout;
      }
      else os = &file;
   }
   cout << "mass = " << mass << '\n'
        << "frequency = " << frequency << '\n'
        << "time_step = " << time_step << '\n'
        << "n_steps = " << n_steps << "\n\n" << std::flush;
   if(file) {
      (*os) << "mass = " << mass << '\n'
            << "frequency = " << frequency << '\n'
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
      osci.Initial(0., x),
      osci.Final(propa_time, x);   // closed path
      clog << "Exact(x = " << x << ") = "
           << osci.ExactAmplitude() << '\n';
      std::vector<std::pair<double, double>> bounds
         (n_steps-1, {x-delta_x, x+delta_x});
      GSLMonteVegas vegas_ho
         (PathIntegrand<EuclidHarmonicOscillator1D>(osci, n_steps),
          bounds, err);
      try { vegas_ho.Integrate(n_calls); }
      catch(int err) {
         clog << "Vegas integration routine returned with error code "
              << err << "\n\n" << std::flush;
         continue;
      }
      double result = vegas_ho.Result(),
             abserr = vegas_ho.AbsErr(),
             chisquare = vegas_ho.ChiSquare();
      amps[k] = result;
      clog << "Vegas = " << result << " +/- " << 3.*abserr
           << ", err = " << abserr/result
           << ", chi square per dof = " << chisquare << "\n\n" << std::flush;
   }


   // Extract the ground state energy:
   _1D::SimpsonRule<double> simpson;
   double integral = 2.*simpson(xs, amps);
   double vacuum_energy = -std::log(integral)/propa_time;
   cout << "Vacuum energy = " << vacuum_energy << "\n\n" << std::flush;
   if(file) (*os) << "Vacuum energy = " << vacuum_energy
                  << "\n\n" << std::flush;


   // Extract the ground state wavefunction:
   std::vector<double> ground_wf(n_points);
   (*os) << "Ground state wavefunction:\n";
   for(size_t k = 0; k != n_points; ++k) {
      ground_wf[k] = std::sqrt(amps[k]/integral);
      (*os) << xs[k] << '\t' << ground_wf[k] << '\n';
   }
   (*os) << std::endl;
         //^ this is essential for correct file reading from the python script

   // Plot via python scripts, if log file was successfully used:
   if(file) {
      file.clear(), file.close();
      std::string program("python harmonic.py ");
      program.append(file_name).append(" &");
      Unused( std::system(program.data()) );   // return value ignored
   }
   else file.clear(), file.close();
}




