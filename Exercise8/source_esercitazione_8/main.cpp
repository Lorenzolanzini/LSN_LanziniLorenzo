/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
//methods and functions are in class random. See random.h and random.cpp
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   double step = 2;
   double y;
   double x_old = 0.;
   double x_new;
   double mean;
   mean = 0.815302; //these values for sigma and mean are the best values found. 
   double sigma = 0.638483 ; //At the beginning I started with mean=1, sigma=1
   int counter = 0;
   int tot=100000;

   std:: vector<double> sampling (tot);
   //start sampling 100000 points
   for(int i=0; i<100000; i++){
   
      x_new = rnd.Rannyu(-1, 1) * step + x_old;

      if(rnd.square_modulus(x_new, mean, sigma) > rnd.square_modulus(x_old, mean, sigma)){
         x_old = x_new;
         counter ++;
      }
      else{
         y=rnd.Rannyu();
         if(y<rnd.square_modulus(x_new, mean, sigma)/rnd.square_modulus(x_old, mean, sigma)){
            x_old = x_new;
            counter ++;
         } 
      }
      sampling[i] = x_old;
}
   

   std::ofstream fout ("Gaussian_sampling_final.txt");
   for(int i=0; i<100000; i++){
      fout<< sampling[i] <<std::endl;
   }
   //std::cout<<rnd.integral(sampling, mean, sigma)<<endl;
   
   //double val = rnd.integral(sampling, mean, sigma);
   //std::cout<<val<<endl;

   std::vector<double> v(2);
   //simulated_annealing //Tin=1, Tf=0.005, step_m=1, step_s=0.25, mean=2.5, sigma=0.75
   v=rnd.simulated_annealing(sampling, 1.0, 0.005, 1, 0.25, 2.5, 0.75);
   //print results for parameters
   std::cout<<v[0]<<" "<<v[1]<<endl;
   //print integral value
   std::cout<<rnd.integral(sampling, v[0], v[1]);
   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
