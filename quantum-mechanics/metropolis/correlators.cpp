#include <cstdlib>   // std::system
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>
#include "PathIntegrals.hpp"
#include "UtilityFunctions.hpp"

// TASK:
//    Compute the first excitation energy for the
//    harmonic oscillator using the Metropolis algortihm.

// USAGE:
//    For plotting, provide a log file name on the command line.

struct Kernel1 {
   Kernel1(size_t const n, size_t const i) : interval(i), n_steps(n) {}
   template<class Vector> double operator() (Vector const &x) {
      double ker = 0.;
      for(size_t site = 0; site != n_steps; ++site)
         ker += x[site]*x[(site+interval) % n_steps];
      return ker/n_steps;
   }
   size_t interval;
   size_t n_steps;
};

struct Kernel3 {
   Kernel3(size_t const n, size_t const i) : interval(i), n_steps(n) {}
   template<class Vector> double operator() (Vector&& x) {
      double ker = 0.;
      for(size_t site = 0; site != n_steps; ++site)
         ker += POW_3(x[site]*x[(site+interval) % n_steps]);
      return ker/n_steps;
   }
   size_t interval;
   size_t n_steps;
};

template<class Vector> bool IsTrue(Vector&& x) {
   for(auto const b : x) if((bool)b == false) return false;
   return true;
}

int main(int narg, char const **args) {

   using std::cin, std::cout, std::clog, std::cerr, std::flush, std::endl;
   cout << std::setprecision(6) << std::fixed;
   clog << std::setprecision(6) << std::fixed;

   // Input params:
   double mass, frequency, time_step;
   size_t n_steps;
   double epsi;
   size_t n_conf, n_corr;

   // Read inputs:
   std::string word;
   std::vector<bool> read_params(7, false);
   do {
      cin >> word;
      if(!cin) {
         cin.clear();
         break;
      }
      else if(word == "mass") {
         cin >> word,
         cin >> mass;
         if(cin) read_params[0] = true;
      }
      else if(word == "frequency") {
         cin >> word,
         cin >> frequency;
         if(cin) read_params[1] = true;
      }
      else if(word == "time_step") {
         cin >> word,
         cin >> time_step;
         if(cin) read_params[2] = true;
      }
      else if(word == "n_steps") {
         cin >> word,
         cin >> n_steps;
         if(cin) read_params[3] = true;
      }
      else if(word == "epsi") {
         cin >> word,
         cin >> epsi;
         if(cin) read_params[4] = true;
      }
      else if(word == "n_conf") {
         cin >> word,
         cin >> n_conf;
         if(cin) read_params[5] = true;
      }
      else if(word == "n_corr") {
         cin >> word,
         cin >> n_corr;
         if(cin) read_params[6] = true;
      }
      word.erase();
   }
   while(!IsTrue(read_params));

   // Set log stream:
   std::ostream *out_stream;
   std::ofstream file;
   std::string file_name;
   if(narg < 2) {
      out_stream  = &cout;
      file.setstate(std::ios_base::badbit);
   }
   else {
      file_name = args[1];
      file.open(file_name, std::ios_base::app);
      if(!file) {
      cerr << "Could not open file " << file_name << "... "
           << "using stdout for logging instead. No plotting.\n\n"
           << flush;
      out_stream = &cout;
      }
      else out_stream  = &file;
   }

   // Print infos:
   cout  << "mass = " << mass << '\n'
         << "frequency = " << frequency << '\n'
         << "time_step = " << time_step << '\n'
         << "n_steps = " << n_steps << '\n'
         << "n_conf = " << n_conf << '\n'
         << "n_corr = " << n_corr << '\n'
         << "epsi = " << epsi << "\n\n" << flush;
   if(file) {
      file  << "mass = " << mass << '\n'
            << "frequency = " << frequency << '\n'
            << "time_step = " << time_step << '\n'
            << "n_steps = " << n_steps << '\n'
            << "n_conf = " << n_conf << '\n'
            << "n_corr = " << n_corr << '\n'
            << "epsi = " << epsi << "\n\n" << flush;
   }

   // Setup:
   double propa_time = n_steps*time_step;
   EuclidHarmonicOscillator1d osci(mass, frequency);
   // auto potential
   // = [&](double x) { return 0.5*mass*POW_2(frequency*x) + 0.5*POW_4(x); };
   // EuclidParticle1d osci(mass, potential);
   //^ could also try with some parity invariant anharmonic potential
   osci.Initial({0., 0.}),
   osci.Final({propa_time, 0.});
   Metropolis<EuclidParticle1d> metropolis(&osci, n_steps);

   // Compute the propagator for the first time intervals on the lattice:
   size_t max_interval = /*n_steps/4*/ 8;
   std::vector<double> propas, errs;
   std::vector<std::function<double(std::vector<double>)>> kernels;
   for(size_t i = 0; i <= max_interval; ++i)
      kernels.push_back( Kernel1(n_steps, i) );
   metropolis(kernels,
              n_conf, n_corr, 100*n_corr, epsi,
              BC_Type::periodic);
   propas = metropolis.Results(),
   errs = metropolis.AbsErrors();
   cout << '\n';
   cout << "last accept rate: " << metropolis.AcceptRate() << "\n\n";

   // Elaborate:
   std::vector<double> energies(max_interval), errors(max_interval);
   for(size_t i = 0; i <= max_interval-1; ++i) {
      energies[i]
         = std::log(propas[i]/propas[i+1])/time_step,
      errors[i]
         = std::sqrt(POW_2(errs[i]/propas[i])
                     + POW_2(errs[i+1]/propas[i+1]))/time_step;
   }
   std::vector<double> slopes(max_interval-1);
   for(size_t i = 0; i != max_interval-1; ++i)
      slopes[i] = std::abs(energies[i+1] - energies[i]);
   size_t const best_interval = FindMinIndex(slopes);
   double const best_time = (best_interval+0.5)*time_step,
                best_energy = 0.5*(energies[best_interval]
                                    + energies[best_interval+1]);

   // Print out:
   *out_stream << "excitation energy estimates:\n";
   for(size_t i = 0; i != max_interval-1; ++i)
      *out_stream << i << " / " << i+1 << " -> "
           << energies[i] << " +/- " << errors[i] << '\n';
   *out_stream << '\n' << flush;
   cout << "best_estimate = " << best_energy << '\n'
        << "around t = " << best_time << endl;
   if(file) file << "best_estimate = " << best_energy << '\n'
                 << "around t = " << best_time << endl;

   // Plot via python script, if log file was successfully used:
   if(file) {
      std::string program("python correlators.py ");
      program.append(file_name).append(" &");
      Ignore( std::system(program.data()) );   // return value ignored
   }
   file.clear(), file.close();
}




