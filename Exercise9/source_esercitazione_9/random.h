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

using namespace std;

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

//TSP METHODS

//void initial_population(int n);
void Check(std::vector<int> &v );
void Check_c(std::vector<int> &v );

double cost(std::vector<vector<int>> pop,std::vector<double> &coordinates_x, std::vector<double> &coordinates_y, int n);
void sort(std::vector <int> &v);
int choose(int M, int p);

//Genetic Mutation Methods
std::vector<int> permutation (std::vector<int> &v);
std::vector<int> shift (std::vector<int> &v, int m);
std::vector<int> inversion (std::vector<int> &v);
std::vector<int> permutation_m (std::vector<int> &v , unsigned int m);

//Read from file
std::vector<double> ReadAll(const char *filename);
std::vector<vector<int>> ReadAllint(const char *filename, int jtot);

//Ordering
std::vector<vector<int>> order(std::vector<vector<int>> &pop, int n, const char*filenamex, const char *filenamey);
void order2(std::vector<vector<int>> &pop, int n, const char*filenamex, const char *filenamey);

//CROSSOVER
std::vector<int> crossover1(std::vector<int> &gen1, std::vector<int> &gen2);

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
