#include <iostream>
#include <iomanip>
#include <OutputManager/OutputManager.h>
#include "DynamicalSystems.hpp"
using boost::core::demangle;

template<class C, class D> void CheckSolvability(DynamicalSystem<C,D>* s) {
   if( s->IsExactlySolvable() and s->IsDeterministic() )
      std::cout << "The " << demangle(typeid(*s).name()) <<
                     " is exactly solvable with these BCs!\n";
}

int main() {
   EuclideanHarmonicOscillator osci;
   std::map<const char*,double> vars =
      {std::pair("frequency",osci.Frequency()),
       std::pair("mass",osci.Mass())};
   std::cout << "Properties of the " << demangle(typeid(osci).name()) << ": "
      << vars << '\n';
   osci.SetCondition(0.,0.);
   try { osci.SetDiscretisedPath(LinearRange(0.,1.,11)); }
   catch(...) { osci.SetCondition(1.,1.); }
   CheckSolvability(&osci);
   osci.SetDiscretisedPath(LinearRange(0.,1.,11));
   osci.PrintDiscretisedPath("Discretized path is:\n");
}





