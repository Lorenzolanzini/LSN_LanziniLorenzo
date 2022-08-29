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
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>


using namespace std;
 
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

   //data
   double sigma = 0.25;
   double T = 1;
   double K = 100;
   double r = 0.1;
   int N = 1000;
   double S_0 = 100;
   double Zi;
   double Si_T;
   std::vector<double> CN(10000);
   int blocks = 100;
   double profit;

   std::vector<double> C(N);

   //2 methods for both call and put options are here;
   // comment the ones that you don't want to study

   //call 1
   for(int k =0; k < 10000; k++){
      for(int i=0; i<N; i++){
         Zi=rnd.Gauss(0, 1);
         Si_T=S_0*pow(M_E, (r-sigma*sigma*0.5)*T +sigma*Zi*sqrt(T));
         if(Si_T < K){
            profit = 0;
         }
         else{
            profit = Si_T - K;
         }
         C[i]=pow(M_E, -r*T)*profit;
      }
      for(int j=0; j<N; j++){
         CN[k] += C[j];
      }
      CN[k]/=(double)N;
   }

   //results with progrssive averages
   std::vector<double> prog_call1(blocks);
   std::vector<double> prog_call1_inc(blocks);
   for(int l=0; l<blocks; l++){
      prog_call1_inc[l] = rnd.progressive2(CN, 10000, blocks)[l];
      prog_call1[l] = rnd.progressive(CN, 10000, blocks)[l];
   }

   ofstream call1 ("call1.txt");
   for(int l = 0; l < blocks; l++){
      call1 << l*100<<" "<<prog_call1[l]<<" "<<prog_call1_inc[l]<<endl;
   }   
   CN.clear();
   C.clear();
   std::vector<double> S (100);
   
   //call 2
   for(int k =0; k < 10000; k++){
      for(int i=0; i<N; i++){
         S[0]=S_0;
         for(int m=1; m<100; m++){
            Zi=rnd.Gauss(0, 1);
            S[m]=S[m-1]*pow(M_E, (r-sigma*sigma*0.5)*0.01 + sigma*Zi*sqrt(0.01));
            Si_T = S[m];
         }

         if(Si_T < K){
            profit = 0;
         }
         else{
         profit = Si_T - K;
         }
         C[i]=pow(M_E, -r*T)*profit;
      }
      for(int j=0; j<N; j++){
         CN[k] += C[j];
      }
      CN[k]/=(double)N;
   }

   //progressive averages

   std::vector<double> prog_call2(blocks);
   std::vector<double> prog_call2_inc(blocks);
   for(int l=0; l<blocks; l++){
      prog_call2_inc[l] = rnd.progressive2(CN, 10000, blocks)[l];
      prog_call2[l] = rnd.progressive(CN, 10000, blocks)[l];
   }

   ofstream call2 ("call2.txt");
   for(int l = 0; l < blocks; l++){
      call2 << l*100<<" "<<prog_call2[l]<<" "<<prog_call2_inc[l]<<endl;
   }
 
   //put 1
   CN.clear();
   C.clear();
   for(int k =0; k < 10000; k++){
      for(int i=0; i<N; i++){
         Zi=rnd.Gauss(0, 1);
         Si_T=S_0*pow(M_E, (r-sigma*sigma*0.5)*T +sigma*Zi*sqrt(T));
         if(Si_T > K){
            profit = 0;
         }
         else{
            profit = -Si_T + K;
         }
         C[i]=pow(M_E, -r*T)*profit;
      }
      for(int j=0; j<N; j++){
         CN[k] += C[j];
      }
      CN[k]/=(double)N;}

   std::vector<double> prog_put1(blocks);
   std::vector<double> prog_put1_inc(blocks);
   for(int l=0; l<blocks; l++){
   prog_put1_inc[l] = rnd.progressive2(CN, 10000, blocks)[l];
   prog_put1[l] = rnd.progressive(CN, 10000, blocks)[l];
   }

   ofstream put1 ("put1.txt");
   for(int l = 0; l < blocks; l++){
      put1 << l*100<<" "<<prog_put1[l]<<" "<<prog_put1_inc[l]<<endl;
   }  

   //put 2
   CN.clear();
   C.clear();
   std::vector<double> S2(100);
   for(int k =0; k < 10000; k++){
      for(int i=0; i<N; i++){
         S2[0]=S_0;
         for(int m=1; m<100; m++){
            Zi=rnd.Gauss(0, 1);
            S2[m]=S2[m-1]*pow(M_E, (r-sigma*sigma*0.5)*0.01 + sigma*Zi*sqrt(0.01));
            Si_T = S2[m];
         }
         if(Si_T > K){
            profit = 0;
         }  
         else{
            profit = -Si_T + K;
         }
         C[i]=pow(M_E, -r*T)*profit;
      }
      for(int j=0; j<N; j++){
         CN[k] += C[j];
      }
   CN[k]/=(double)N;}

   std::vector<double> prog_put2(blocks);
   std::vector<double> prog_put2_inc(blocks);
   for(int l=0; l<blocks; l++){
      prog_put2_inc[l] = rnd.progressive2(CN, 10000, blocks)[l];
      prog_put2[l] = rnd.progressive(CN, 10000, blocks)[l];
   }

   ofstream put2 ("call2.txt");
   for(int l = 0; l < blocks; l++){
      put2 << l*100<<" "<<prog_put2[l]<<" "<<prog_put2_inc[l]<<endl;
   }

   rnd.SaveSeed();
   prog_put2.clear();
   prog_put2_inc.clear();
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
