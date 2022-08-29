/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "funzioni.h"
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>

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

  

   int M=1000000;
   int N=100;
   int L=M/N;

   std::vector<double> numeri(M);
   for(int k=0; k<=M; k++){
      numeri[k]=rnd.Rannyu();  //filling the vector with random numbers
   }
   //creazione di blocchi
   std::vector<double> medie(N);
   std::vector<double> medie2(N);
   double sum;
   int k=0;
   for(int i=0; i<N; i++){
      sum=0;
      for(int j=0; j<=L; j++){
         k=j+i*L;
         sum+=numeri[k];
      }
      medie[i]=sum/L;
      medie2[i]= medie[i]*medie[i];
   }
   //cumulative
   std::vector<double> prog(N);
   std::vector<double> prog2(N);
   std::vector<double>error(N);
   
   for(int i=0; i<N; i++){
      prog[i]=0;
      prog2[i]=0;
      for(int j=0; j<i+1; j++){
         prog[i] += medie[j];
         prog2[i] += medie2[j];
      }
      prog[i]/=((i+1));  //progressive averages
      prog2[i]/=((i+1));
      if(i==0){
         error[i]=0;
      }
      else{
      error[i]=sqrt((-(prog[i]*prog[i])+(prog2[i]))/(i)); //progressive errors
      }
   }
    ofstream fout ("data.txt");
   for(int i=0; i<N; i++){
      fout<<i*L<<" "<<prog[i]<<" " <<error[i]<<endl;
   }
   fout.close();


   //varianza
   std::vector<double> varianze(N);
   std::vector<double> varianze2(N);
   for(int i=0; i<N; i++){
      sum=0;
      for(int j=0; j<=L; j++){
         k=j+i*L;
         sum+=(numeri[k]-0.5)*(numeri[k]-0.5);
      }
      varianze[i]=sum/L;
      varianze2[i]= varianze[i]*varianze[i];
   }
   //cumulative
   std::vector<double> prog_var(N);
   std::vector<double> prog_var2(N);
   std::vector<double>error_var(N);
   
   for(int i=0; i<N; i++){
      prog_var[i]=0;
      prog_var2[i]=0;
      for(int j=0; j<i+1; j++){
         prog_var[i] += varianze[j];
         prog_var2[i] += varianze2[j];
      }
      prog_var[i]/=((i+1));
      prog_var2[i]/=((i+1));
      if(i==0){
         error_var[i]=0;
      }
      else{
      error_var[i]=sqrt((-(prog_var[i]*prog_var[i])+(prog_var2[i]))/(i));
      }
   }
   ofstream fout2 ("Varianza.txt");
   for(int i=0; i<N; i++){
      fout2<<i*L<<" "<<prog_var[i]<<" " <<error_var[i]<<endl;
   }
   fout.close();
   
   //test chi-quadro
   int blocchi=100;
   int n;

   std::vector<int> occupazioni(blocchi);
   std::vector<int> chi2(100);
   for(int i=0; i<100; i++){
      chi2[i]=0;
         for(int j=0; j<blocchi; j++){
            n=0;
            for(int l=0; l<10000; l++){
               k=l+(i*10000);
               if(numeri[k]<=0.01*(j+1) && numeri[k]>0.01*j){
                  n++;
               }
            }
            occupazioni[j]=(n-100)*(n-100)/(100);
            chi2[i]+=occupazioni[j];
      }
   }
   ofstream chi ("chi_quadro.txt");

   for(int i=0; i<N; i++){
      chi<<chi2[i]<<" " <<endl;
      }

   //Exponential distribution   
   double lambda=1;
   ofstream exponential("exp.txt");
   for(int i=0; i<10000; i++) {
      exponential<<rnd.Exponential(lambda)<<endl;
   }

   //Cauchy Distribution
   double gamma=1;
   double mean=0;
   ofstream Cauchy("Cauchy.txt");
   for(int i=0; i<10000; i++){
      Cauchy<<rnd.Cauchy(mean, gamma)<<endl;
   }
   rnd.SaveSeed();
   std::vector <double> m(10000);
   
   //Exponential (2-average)
   N=2;
   for(int i=0; i<10000; i++){
      m[i]=0;
      for(int j=0; j<N; j++){
         m[i]+=rnd.Exponential(lambda)/N;
      }
   }
   std::ofstream exp2("exp2.txt");
   for(int i=0; i<10000; i++){
      exp2<<m[i]<<endl;
   }
   m.clear();
   N=10;
   for(int i=0; i<10000; i++){
      m[i]=0;
      for(int j=0; j<N; j++){
         m[i]+=rnd.Exponential(lambda)/N;
      }
   }
   std::ofstream exp10("exp10.txt");
   for(int i=0; i<10000; i++){
      exp10<<m[i]<<endl;
   }
   m.clear();

   N=100;
   for(int i=0; i<10000; i++){
      m[i]=0;
      for(int j=0; j<N; j++){
         m[i]+=rnd.Exponential(lambda)/N;
      }
   }
    std::ofstream exp100("exp100.txt");
   for(int i=0; i<10000; i++){
      exp100<<m[i]<<endl;
   }
   m.clear();

   //Cauchy
   N=2;
   for(int i=0; i<10000; i++){
      m[i]=0;
      for(int j=0; j<N; j++){
         m[i]+=rnd.Cauchy(mean, gamma)/N;
      }
   }
   std::ofstream cauchy2 ("Cauchy2.txt");
   for(int i=0; i<10000; i++){
      cauchy2<<m[i]<<endl;
   }
   m.clear();

   N=10;
   for(int i=0; i<10000; i++){
      m[i]=0;
      for(int j=0; j<N; j++){
         m[i]+=rnd.Cauchy(mean, gamma)/N;
      }
   }
   std::ofstream cauchy10 ("Cauchy10.txt");
   for(int i=0; i<10000; i++){
      cauchy10<<m[i]<<endl;
   }
   m.clear();

   N=100;
   for(int i=0; i<10000; i++){
      m[i]=0;
      for(int j=0; j<N; j++){
         m[i]+=rnd.Cauchy(mean, gamma)/N;
      }
   }
   std::ofstream cauchy100 ("Cauchy100.txt");
   for(int i=0; i<10000; i++){
      cauchy100<<m[i]<<endl;
   }
   
   m.clear();


   ///Buffon's Experiment
   double l = 0.9; //needle
   double d = 1; //distance of the grating
   double angle;
   double x;
   double y;
   //number between 0 to angle and starting point
   std::vector<double> starting_x(100000);
   std::vector<double> theta(100000);
   vector<double>pigreco(10000);
   int counter;
   for(int j=0; j<10000; j++){
   counter=0;
   for(int i=0; i<100000; i++){
      starting_x[i] = rnd.Rannyu();
      do{
         x = rnd.Rannyu()*2 -1;
      y = rnd.Rannyu()*2 - 1;}while(x*x +y*y > 1);
      
      if(y>=0){
         angle = acos(x/(sqrt(x*x + y*y)));
      }
      else{
         angle = - acos(x/(sqrt(x*x + y*y)));
      }
      theta[i] = angle; //theta extracted without using pi
      if(l*cos(theta[i])+starting_x[i]>d || l*cos(theta[i])+starting_x[i]<0){
         counter++;
      }
   }
   pigreco[j]=2*100000*l/(counter*d);
   }
   
   N = 10000;
   blocchi=100;
   std::vector<double> pigrecoav(blocchi);
   std::vector<double> pigrecoRMSQ(blocchi);

   for(int i=0; i<blocchi; i++){
      pigrecoav[i]=rnd.progressive(pigreco, N, blocchi)[i];
      pigrecoRMSQ[i]= rnd.progressive2(pigreco, N, blocchi)[i];
   }

   ofstream buffon("Buffon.txt");
   for(int j=0; j<blocchi; j++){
      buffon<<j*N/blocchi<<" "<<pigrecoav[j]<<" "<<pigrecoRMSQ[j]<<endl;
   }
   ofstream buffon2("theta.txt");
   for(int j=0; j<100000; j++){
      buffon2<<theta[j]<<endl;
   }
   pigrecoav.clear();
   pigrecoRMSQ.clear();

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
