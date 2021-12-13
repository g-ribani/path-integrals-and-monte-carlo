#include <cstdlib>   // std::system
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>
#include <boost/math/statistics/univariate_statistics.hpp>
#include "headers/EuclidHarmonicOscillator1d.hpp"

// TASK:
//    Compute the uncertainties on the first excitation energy for the
//    harmonic oscillator obtained via several different binnings
//    of the Metropolis results.

// USAGE:
//    For plotting, provide a log file name on the command line.


struct Kernel1 {
   Kernel1(size_t const n, size_t const i) : interval(i), n_steps(n) {}
   template<class Vector> double operator() (Vector &&x) {
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
   template<class Vector> double operator() (Vector &&x) {
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
   osci.Initial({0., 0.}),
   osci.Final({propa_time, 0.});
   osci.NSteps(n_steps);
   Metropolis<EuclidParticle1d> metropolis(&osci);

   // Compute the propagator for the first time intervals on the lattice:
   size_t const interval = 2;
   std::vector<std::function<double(double*)>> kernels
      { Kernel1(n_steps, interval),
         Kernel1(n_steps, interval+1) };
      //^ try also with Kernel3
   metropolis.BoundaryConditions(BC_Type::periodic);
   auto configs = metropolis(kernels, n_conf, n_corr, 100*n_corr, epsi);

   // Compute energy estimate:
   using namespace boost::math::statistics;
   double propa0 = mean(configs[0]),
          propa1 = mean(configs[1]),
          estimate = std::log(propa0/propa1)/time_step;
   cout << "estimate = " << estimate << '\n'
        << "accept rate = " << metropolis.AcceptRate() << "\n\n";

   // Incremental binning:
   std::vector<size_t> bin_size, n_bins;
   for(size_t s = 1; s <= n_conf/10; ++s) if(n_conf % s == 0)
      bin_size.push_back(s),
      n_bins.push_back(n_conf/s);
   size_t const n_binnings = n_bins.size();
   std::vector<std::vector<double>>
      binned_propas(2, std::vector<double>(n_binnings)),
      binned_errs(2, std::vector<double>(n_binnings));
   size_t N, s;
   for(size_t b = 0; b != n_binnings; ++b) {
      N = n_bins[b], s = bin_size[b];
      std::vector<double> bin_averages(N);
      for(size_t i = 0; i != 2; ++i) {
         auto configs_begin = configs[i].data();
         for(size_t n = 0; n != N; ++n)
            bin_averages[n] = mean(configs_begin+n*s, configs_begin + (n+1)*s);
         binned_propas[i][b] = mean(bin_averages);
         binned_errs[i][b] = std::sqrt(sample_variance(bin_averages)/N);
      }
   }

   // Elaborate:
   std::vector<double> binned_rel_errors(n_binnings);
      for(size_t b = 0; b != n_binnings; ++b)
         binned_rel_errors[b]
            = std::sqrt(POW_2(binned_errs[0][b]/binned_propas[0][b])
                  + POW_2(binned_errs[1][b]/binned_propas[1][b]))
               /(time_step*estimate);

   // Print out:
   *out_stream << "excitation energy errors ( bin_size ):\n";
   for(size_t b = 0; b != n_binnings; ++b)
      *out_stream << binned_rel_errors[b]
         << " ( " << bin_size[b] << " )\n";
   *out_stream << '\n' << flush;

   // Plot via python script, if log file was successfully used:
   if(file) {
      std::string program("python binning.py ");
      program.append(file_name).append(" &");
      Ignore( std::system(program.data()) );   // return value ignored
   }
   file.clear(), file.close();
}




