#include <iostream>
#include <fstream>
#include "MonteCarlo.hpp"

int main(int nargs, char const **args) {
   if(nargs<3) {
      std::cerr << "tell me the name of two output files, the first "
                << "for Crude and the second for Vegas\n";
      exit(1);
   }
   std::ofstream outputs[2];
   outputs[0].open(args[1], std::ios::out | std::ios::app),
   outputs[1].open(args[2], std::ios::out | std::ios::app);

   EuclidHarmonicOscillator1D osci(1., 1.);
   osci.AddBoundaryCondition(0.,0),
   osci.AddBoundaryCondition(1.,1);

   double xmin = -5.,
          xmax = 5.;
   std::size_t nEvals = 10'000'000,
                nSteps = 6;
   MCResult amp;
   std::mt19937_64 gen;

   double exact_amp = osci.ExactAmplitude();
   std::cout << std::scientific
             << "Exact amplitude = " << exact_amp << "\n\n" << std::flush;
   for(auto &out : outputs) {
      out << std::scientific
          << "\n\n\t***********\t\n\n"
          << "Exact amplitude = " << exact_amp << '\n'
          << "nEvals = " << (float)nEvals << '\n'
          << "xmin xmax = " << xmin << ' ' << xmax << "\n\n" << std::flush;
   }

   while(true) {
      #ifdef BE_CRUDE
         amp = CrudeMCAmplitude(&osci, nSteps, xmin, xmax, gen, nEvals);
         outputs[0] << "Crude MonteCarlo amplitude = "
                    << amp.res << " +/- " << amp.err
                    << " // " << nSteps << " steps"
                    << "\n\n" << std::flush;
         std::cout << "crude = " << amp.res << '\n' << std::flush;
      #endif

      #ifdef BE_VEGAS
         try { amp = VegasMCAmplitude(&osci, nSteps, xmin, xmax, nEvals); }
         catch(int errcode) {
            std::cout << "Vegas error code " << errcode << '\n' << std::flush;
            continue;
         }
         outputs[1] << "VEGAS MonteCarlo amplitude = "
                    << amp.res << " +/- " << amp.err
                    << " // " << nSteps << " steps"
                    << "\n\n" << std::flush;
         std::cout << "vegas = " << amp.res << '\n' << std::flush;
      #endif
   }

   std::cout << '\n';
   for(auto &out : outputs) out.clear(), out.close();
}




