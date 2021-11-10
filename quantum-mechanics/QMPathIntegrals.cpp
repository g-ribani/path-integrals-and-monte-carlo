#include <array>
#include <cmath>
#include <vector>
#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

template<class T> struct PathIntegrand;


//// one-dimensional particle in a potential in euclidean time

template<class PotentialFunc> struct EuclidParticle1D {
   EuclidParticle1D(double const &m,
                     PotentialFunc const &V = [](double const&) { return 0.; })
    : mass(m), potential(V) {}
   std::array<double, 2> initial_pos;
   std::array<double, 2> final_pos;
   double const mass;
   EuclidParticle1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      step = (final_pos[1] - initial_pos[1])/N;
      return *this;
   }
   PotentialFunc potential;
   protected:
   double n_steps;
   double step;
};

template<class PotentialFunc>
 struct PathIntegrand<EuclidParticle1D<PotentialFunc>>
 : EuclidParticle1D<PotentialFunc> {
   double operator() (std::vector<double> const &x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = x.size();
      std::vector<double> deltaX(dim + 1);
      deltaX[0] = x[0] - this->initial_pos[1];
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = this->final_pos[1] - x[dim-1];
      double const &mass = this->mass;
      double &step = this->step;
      PotentialFunc &potential = this->potential;
      double ret = 1.;
      for(size_type k = 0; k != dim + 1; ++k) {
         double deltaS = mass/2.*deltaX[k]*deltaX[k]/step
                         + step*( potential(x[k]) + potential(x[k-1]) );
         ret *= std::sqrt(mass/(2.*M_PI*step)) * std::exp(-deltaS);
      }
      return ret;
   }
   PathIntegrand(const EuclidParticle1D<PotentialFunc> &part)
      : EuclidParticle1D<PotentialFunc>(part) {}
};


//// one-dimensional free particle in euclidean time

struct EuclidFreeParticle1D {
   EuclidFreeParticle1D(double const &m) : mass(m) {}
   std::array<double, 2> initial_pos;
   std::array<double, 2> final_pos;
   double const mass;
   EuclidFreeParticle1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      step = (final_pos[1] - initial_pos[1])/N;
      return *this;
   }
   double ExactAmplitude() const {
      double const deltaT = final_pos[0] - initial_pos[0],
                   deltaX = final_pos[1] - initial_pos[1];
      return std::sqrt( mass/(2.*M_PI*deltaT) )
         *std::exp( -mass*deltaX*deltaX/(2.*deltaT) );
   }
   protected:
   double n_steps;
   double step;
};

template<> struct PathIntegrand<EuclidFreeParticle1D>
 : EuclidFreeParticle1D {
   double operator() (std::vector<double> const &x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = x.size();
      std::vector<double> deltaX(dim + 1);
      deltaX[0] = x[0] - initial_pos[1];
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = final_pos[1] - x[dim-1];
      double ret = 1.;
      for(size_type k = 0; k != dim + 1; ++k) {
         double deltaS = mass/2.*deltaX[k]*deltaX[k]/step;
         ret *= std::sqrt(mass/(2.*M_PI*step)) * std::exp(-deltaS);
      }
      return ret;
   }
   PathIntegrand(const EuclidFreeParticle1D &part)
      : EuclidFreeParticle1D(part) {}
};


//// one-dimensional harmonic oscillator in euclidean time

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
   double ExactAmplitude() const {
      double const deltaT = final_pos[0] - initial_pos[0],
                   x0 = initial_pos[1],
                   x1 = final_pos[1];
      return std::sqrt( mass/(2.*M_PI*deltaT) )
         *std::exp( -mass*frequency/2.
                     *( (x0*x0 + x1*x1)/std::tanh(frequency*deltaT)
                        - 2.*x0*x1/std::sinh(frequency*deltaT) ) );
   }
   protected:
   double n_steps;
   double step;
};

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




