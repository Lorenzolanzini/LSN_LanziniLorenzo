/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__
#include <vector>

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
  
  
  double square_modulus(double x, double mean, double sigma);
  double function(double x, double mean, double sigma);
  
  //progressive data blocking
  std::vector<double> progressive( std::vector <double> &Numbers, int N, int blocchi);
  std::vector<double> progressive2( std::vector <double> &Numbers, int N, int blocchi);
  
  //calculaation of the integral
  double integral(std::vector<double> sampling, double mean, double sigma);
  double integral2(std::vector<double> sampling, double mean, double sigma);
  
  
  std::vector<double> simulated_annealing(std::vector<double> sampling, double Tin, double Tfinal, double step_m, double step_s, double mean, double sigma);
  //S_A algorithm takes samplig, initial and final Ts, mean, sigma and steps as parameters. It returns
  //a vector with the value of mean and sigma
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
