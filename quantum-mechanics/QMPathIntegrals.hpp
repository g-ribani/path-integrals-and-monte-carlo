#include <cmath>
#include <vector>
#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

template<class T> struct PathIntegrand;


//// one-dimensional particle in a potential in euclidean time

template<class PotentialFunc> struct EuclidParticle1D {
   struct Point {
      double t;
      double x;
   };
   EuclidParticle1D(double const &m,
                     PotentialFunc const &V = [](double const&) { return 0.; })
    : mass(m), potential(V) {}
   Point initial;
   Point final;
   double const mass;
   EuclidParticle1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      time_step = (final.t - initial.t)/N;
      return *this;
   }
   PotentialFunc potential;
   protected:
   double n_steps;
   double time_step;
};

template<class PotentialFunc>
 struct PathIntegrand<EuclidParticle1D<PotentialFunc>>
 : EuclidParticle1D<PotentialFunc> {
   double operator() (std::vector<double> const &x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = x.size();
      std::vector<double> deltaX(dim + 1);
      deltaX[0] = x[0] - this->initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = this->final.x - x[dim-1];
      double const &mass = this->mass;
      double &time_step = this->time_step;
      PotentialFunc &potential = this->potential;
      double ret = 1., deltaS;
      deltaS = mass/2.*deltaX[0]*deltaX[0]/time_step
               + time_step*( potential(x[0]) + potential(this->initial.x) )/2.;
      ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      for(size_type k = 1; k != dim; ++k) {
         double deltaS = mass/2.*deltaX[k]*deltaX[k]/time_step
                         + time_step*( potential(x[k]) + potential(x[k-1]) )/2.;
         ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      }
      deltaS = mass/2.*deltaX[dim]*deltaX[dim]/time_step
              + time_step*( potential(this->final.x) + potential(x[dim-1]) )/2.;
      ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      return ret;
   }
   PathIntegrand(const EuclidParticle1D<PotentialFunc> &part)
      : EuclidParticle1D<PotentialFunc>(part) {}
};


//// one-dimensional free particle in euclidean time

struct EuclidFreeParticle1D {
   struct Point {
      double t;
      double x;
   };
   EuclidFreeParticle1D(double const &m) : mass(m) {}
   Point initial;
   Point final;
   double const mass;
   EuclidFreeParticle1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      time_step = (final.t - initial.t)/N;
      return *this;
   }
   double ExactAmplitude() const {
      double const deltaT = final.t - initial.t,
                   deltaX = final.x - initial.x;
      return std::sqrt( mass/(2.*M_PI*deltaT) )
         *std::exp( -mass*deltaX*deltaX/(2.*deltaT) );
   }
   protected:
   double n_steps;
   double time_step;
};

template<> struct PathIntegrand<EuclidFreeParticle1D>
 : EuclidFreeParticle1D {
   double operator() (std::vector<double> const &x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = x.size();
      std::vector<double> deltaX(dim + 1);
      deltaX[0] = x[0] - initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = final.x - x[dim-1];
      double ret = 1.;
      for(size_type k = 0; k != dim + 1; ++k) {
         double deltaS = mass/2.*deltaX[k]*deltaX[k]/time_step;
         ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      }
      return ret;
   }
   PathIntegrand(const EuclidFreeParticle1D &part)
      : EuclidFreeParticle1D(part) {}
};


//// one-dimensional harmonic oscillator in euclidean time

struct EuclidHarmonicOscillator1D {
   struct Point {
      double t;
      double x;
   };
   EuclidHarmonicOscillator1D(double const &m, double const &f)
    : frequency(f), mass(m) {}
   Point initial;
   Point final;
   double const frequency;
   double const mass;
   EuclidHarmonicOscillator1D& SetNSteps(std::size_t const &N) {
      n_steps = N;
      time_step = (final.t - initial.t)/N;
      return *this;
   }
   double ExactAmplitude() const {
      double const deltaT = final.t - initial.t,
                   x0 = initial.x,
                   x1 = final.x;
      return std::sqrt( mass/(2.*M_PI*deltaT) )
         *std::exp( -mass*frequency/2.
                     *( (x0*x0 + x1*x1)/std::tanh(frequency*deltaT)
                        - 2.*x0*x1/std::sinh(frequency*deltaT) ) );
   }
   double n_steps;
   double time_step;
};

template<> struct PathIntegrand<EuclidHarmonicOscillator1D>
 : EuclidHarmonicOscillator1D {
   double operator() (std::vector<double> const &x) const {
      using size_type = std::vector<double>::size_type;
      auto const dim = x.size();
      std::vector<double> deltaX(dim + 1);
      deltaX[0] = x[0] - initial.x;
      for(size_type k = 1; k != dim; ++k) deltaX[k] = x[k] - x[k-1];
      deltaX[dim] = final.x - x[dim-1];
      double ret = 1., deltaS;
      deltaS = mass/2.*deltaX[0]*deltaX[0]/time_step
               + time_step*mass*frequency*frequency/4.
                  *( x[0]*x[0] + (initial.x)*(initial.x) );
      ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      for(size_type k = 1; k != dim; ++k) {
         deltaS = mass/2.*deltaX[k]*deltaX[k]/time_step
                  + time_step*mass*frequency*frequency/4.
                     *( x[k]*x[k] + x[k-1]*x[k-1] );
         ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      }
      deltaS = mass/2.*deltaX[dim]*deltaX[dim]/time_step
               + time_step*mass*frequency*frequency/4.
                  *( (final.x)*(final.x) + x[dim-1]*x[dim-1] );
      ret *= std::sqrt(mass/(2.*M_PI*time_step)) * std::exp(-deltaS);
      return ret;
   }
   PathIntegrand(const EuclidHarmonicOscillator1D &osci)
      : EuclidHarmonicOscillator1D(osci) {}
};




