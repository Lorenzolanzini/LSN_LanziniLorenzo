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

   
   rnd.SaveSeed();
   
   //metodo della media
   int N=100;
   
   std::vector<double> casuali(1000000);
   for(int i=0; i<1000000; i++){
      casuali[i]=0.5*M_PI*cos(M_PI*rnd.Rannyu()*0.5);
   }
   ofstream fout ("Metodo_Media.txt");
   std::vector<double> devstd_Media(N);
   devstd_Media=rnd.progressive2(casuali, 1000000, N);
   for(int i=0; i<N; i++){
      fout<<i*10000<<" "<<rnd.progressive(casuali, 1000000, N)[i]<<" "<<devstd_Media[i]<<endl;
   }
   fout.close();
   //devstd
   ofstream std_Media("Devstd_Media.txt");
   for(int i=0; i<N; i++){
      std_Media<<i*10000<<" "<<devstd_Media[i]<<" "<<rnd.progressive2(devstd_Media, 100, 100)[i]<<endl;
   }

   //importance sampling
   std::vector<double> numeri(1000000);

   int count=0;
   double x;
   double y;
   while(count<1000000){
   x = rnd.Rannyu();
   y = rnd.Rannyu()*(M_PI/2);
   if(y<=abs(M_PI*0.5-M_PI*M_PI*M_PI*x*x/(16))/((4*sqrt(2)/3)-M_PI*0.5+M_PI*M_PI*M_PI/(48))){
      numeri[count]=x;
      count++;
   }
   }
   ofstream fout2("Sampling.txt");
   for(int i=0; i<1000000;i++){
      fout2<<numeri[i]<<endl;
   }



   int Ntot=1000000;
   int blocks=100;
   int n=Ntot/(int)blocks;
   
   std::vector<double> values(Ntot);
   for(int i=0; i<Ntot; i++){
      values[i] = abs(M_PI*0.5*cos(M_PI*numeri[i]*0.5)*((4*sqrt(2)/3)-M_PI*0.5 +M_PI*M_PI*M_PI/(48))/(M_PI*0.5-((M_PI*M_PI*M_PI*numeri[i]*numeri[i])/16)));
   }
   std::vector<double> Devstd_IS(100);
   Devstd_IS=rnd.progressive2(values, Ntot, blocks);

   ofstream fout3 ("Sampling_Method.txt");
   for(int i=0; i<blocks; i++){
      fout3<<i*n<<" "<<rnd.progressive(values, Ntot, blocks)[i]<<" "<<Devstd_IS[i]<<endl;
   }
   fout3.close();

   ofstream Devstd_sampling("Devstd_I_S.txt");
   for(int i=0; i<blocks; i++){
      Devstd_sampling<<i*n<<" "<<Devstd_IS[i]<<" "<<rnd.progressive2(Devstd_IS, 100, 100)[i]<<endl;
   }


   //exercise 02.2 discrete random walk
   int k;
   double sign;
   double choice; 
   blocks=100;
   int Nu=10000;
   int steps =100;
   int nb=Nu/(int)blocks;
   std::cout<<nb<<endl;
   std::vector<int>riga{0, 0, 0, 0, 0}; //1st: step; 2nd, 3rd, 4th:coordinates; 5th: radius
   // std::vector<vector<int>> random_walk;
   std::vector<vector<double>> radius(steps, vector<double>(steps));

   
   for(int j=0; j<blocks; j++){
   for(int l=0; l<nb; l++){
   for(int i=0; i<steps; i++){
      riga[0] += 1;
      choice=rnd.Rannyu(0.0, 3.0);
      if(choice<=1){
         k=1;
      }
      else{
         if(choice>1 && choice<=2){
         k=2;
         }
         else{
         k=3;
         }
      }
      sign = rnd.Rannyu(0, 2);
      if(sign>=1){
         riga[k]+=1;
      }
      else{
         riga[k] -= 1;
      }
      riga[4] = (riga[1]*riga[1] + riga[2]*riga[2] + riga[3]*riga[3]);
      //random_walk.push_back(riga);
      radius[i][j] += riga[4]/(double)nb;
         }
      
      //random_walk.clear();
      riga.clear();
      riga={0, 0, 0, 0, 0};   
      }
   }
   for(int i=0; i<steps; i++){
   for(int j=0; j<blocks; j++){
      radius[i][j]=sqrt(radius[i][j]);
   }
   }
   std::vector<double> averages(steps);
   std::vector<double> std(steps);
   for(int i=0; i<steps; i++){
      //for(int j=0; j<blocks; j++){
        // averages[i] += radius[i][j];
      //}
      // averages[i] /= (double)blocks;
      // averages[i]=sqrt(averages[i]);
      averages[i] = (rnd.progressive(radius[i], nb, blocks )[blocks-1]);
      std[i] = rnd.progressive2(radius[i], nb, blocks)[blocks-1];
        }

   std::ofstream RMSQ("RMSQ.txt");
   for(int i=0; i<steps; i++){
      RMSQ<<(i+1)<<" "<<averages[i]<<" "<<std[i]<<endl;
   }

   //continous random walk
   double theta;
   double phi; 
   std::vector<double> riga_c{0, 0, 0, 0, 0};
   // std::vector<vector<double>> C_RandomWalk;
   std::vector<vector<double>> radius_c(steps, vector<double>(steps));
   
   for(int j=0; j<blocks; j++){
   for(int l=0; l<nb; l++){
   for(int i=0; i<steps; i++){
      theta=rnd.Rannyu(0, M_PI);
      phi=rnd.Rannyu(0, 2*M_PI);

      riga_c[0] += 1;
      riga_c[1] += cos(phi)*sin(theta);
      riga_c[2] += sin(phi)*sin(theta);
      riga_c[3] += cos(theta);
      riga_c[4] = (riga_c[1]*riga_c[1]+riga_c[2]*riga_c[2]+riga_c[3]*riga_c[3]);
      // C_RandomWalk.push_back(riga_c);
      radius_c[i][j] += riga_c[4]/(double)nb;
         }
      //random_walk.clear();
      riga_c.clear();
      riga_c={0, 0, 0, 0, 0};   
      }
   }
   for(int i=0; i<steps; i++){
   for(int j=0; j<blocks; j++){
      radius_c[i][j]=sqrt(radius_c[i][j]);
   }
   }
   std::vector<double> averages_c(steps);
   std::vector<double> std_c(steps);

   for(int i=0; i<steps; i++){
     
      averages_c[i] = (rnd.progressive(radius_c[i], nb, blocks )[blocks-1]);
      std_c[i] = rnd.progressive2(radius_c[i], nb, blocks)[blocks-1];
   }
   


   std::ofstream RMSQC("RMSQC.txt");
   for(int i=0; i<steps; i++){
      RMSQC<<(i+1)<<" "<<averages_c[i]<<" "<<std_c[i]<<endl;
   }

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
