#include <iostream>
#include "MonteCarlo.hpp"

int main() {
   std::mt19937_64 gen;

   EuclidHarmonicOscillator1D osci(1., 1.);
   osci.AddBoundaryCondition(0.,0),
   osci.AddBoundaryCondition(1.,1);

   std::cout << std::scientific;
   std::cout << "Exact amplitude = " << osci.ExactAmplitude() << std::endl;
   std::clog << std::scientific;
   std::clog << "\n\n\t***********\t\n\n"
             << "Exact amplitude = " << osci.ExactAmplitude()
             << "\n\n" << std::flush;

   double xmin = -5.,
          xmax = 5.;
   MCResult amp;

   std::size_t nEvals = 10'000'000,
                nSteps = 10;

   std::clog << "nEvals = " << (float)nEvals << '\n'
               << "xmin xmax = " << xmin << ' ' << xmax << "\n\n";

   for(; true;) { for(int repeat = 0; true; ++repeat) {
      std::cout << "\nrepeat = " << repeat << '\n' << std::flush;

   #ifdef BE_CRUDE
      #pragma omp parallel
      amp = CrudeMCAmplitude(&osci, nSteps, xmin, xmax, gen, nEvals);
      std::clog << "Crude MonteCarlo amplitude = "
               << amp.res << " +/- " << amp.err << " // " << nSteps << " steps"
               << "\n\n" << std::flush;
      std::cout << "crude = " << amp.res << '\n' << std::flush;
   #endif

   #ifdef BE_VEGAS
      try{
         #pragma omp parallel
         amp = VegasMCAmplitude(osci, nSteps, xmin, xmax, nEvals);
      }
      catch(int errcode) {
         std::cout << "Vegas error code " << errcode << '\n' << std::flush;
         continue;
      }
      std::clog << "VEGAS MonteCarlo amplitude = "
            << amp.res << " +/- " << amp.err << " // " << nSteps << " steps"
            << "\n\n" << std::flush;
      std::cout << "vegas = " << amp.res << '\n' << std::flush;
   #endif
      std::clog << '\n';
   }
} }




