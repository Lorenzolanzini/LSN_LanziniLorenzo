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
#include <cmath>
#include <cstdlib>
#include "random.h"
#include <fstream>
#include <vector>

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}
//double gaussian, probability distribution
double Random :: square_modulus(double x, double mean, double sigma){
   return pow(pow(M_E, -pow(x-mean, 2)/(2*sigma*sigma)) + pow(M_E, -pow(x + mean, 2)/(2*sigma*sigma)), 2);
}


double Random :: function(double x, double mean, double sigma){
   return -0.5*1/(pow(sigma, 4))*(pow(M_E, -pow((x+mean)/sigma, 2)))*(4*(1+pow(M_E, 4*mean*x/(sigma*sigma)))*mean*mean-8*(-1+pow(M_E, 4*mean*x/(sigma*sigma)))*mean*x-2*pow((1+pow(M_E, 2*mean*x/(sigma*sigma))),2)*(sigma*sigma-2*x*x))/(pow(M_E, -0.5*pow((x-mean)/sigma, 2))+pow(M_E, -0.5*pow((x+mean)/sigma, 2)))+pow(x,4)-5*0.5*x*x;}

//progressive average, average
std:: vector<double> Random :: progressive( vector <double> &Numbers, int N, int blocchi){
   int n= Numbers.size()/(int)blocchi;
   int k;
   std::vector<double> blocks(blocchi);
   std::vector<double> progressive(blocchi);
   for(int i=0; i<blocchi; i++){
      for(int j=0; j<n; j++){
         k= i*n + j;
         blocks[i] += Numbers[k]; 
      }
      blocks[i] /= n;
   }
   for(int i = 0; i<blocchi; i++){
      progressive[i]=0;
      for(int j=0; j<i+1; j++){
         progressive[i] += blocks[j];
      }
      progressive[i] /= (i+1);
   }
   return progressive;
   }
   

//progressive averages--for uncertainty
vector<double> Random :: progressive2(vector <double> &Numbers, int N, int blocchi){
   int n = N/(int)blocchi;
   int k;
   std::vector<double> blocks(blocchi);
   std::vector<double> blocks2(blocchi);
   std::vector<double> progressive(blocchi);
   std::vector<double> progressive2(blocchi);
   
   for(int i=0; i<blocchi; i++){
      blocks[i]=0;
      for(int j=0; j<n; j++){
         k= i*n + j;
         blocks[i] += Numbers[k]; 
      }
      blocks[i] /= (double)n;
      blocks2[i]= blocks[i]*blocks[i];
   }
   for(int i=0; i<blocchi; i++){
      progressive[i]=0;
      for(int j=0; j<i+1; j++){
         progressive[i] += blocks[j];
         progressive2[i] += blocks2[j];
      }
      progressive[i] /= (double)(i+1);
      progressive2[i] /=(double) (i+1);
      if(i==0){
         progressive2[i]=0;
      }
      else{
      progressive2[i] = sqrt((-progressive[i]*progressive[i]+progressive2[i])/(i));
   }
   }
   
   return progressive2; 
}

//Integral calculation, average with data blocking
double Random :: integral(std::vector<double> sampling, double mean, double sigma){
   std::vector<double>functions(100000);
   double x_old=0;
   double y;
   double step = 2*mean;
   for(int i=0; i<100000; i++){
      double x_new = Rannyu(-1, 1) * step + x_old;
      if(square_modulus(x_new, mean, sigma) > square_modulus(x_old, mean, sigma)){
         x_old = x_new;
      }
      else{
         y=Rannyu();
         if(y<square_modulus(x_new, mean, sigma)/square_modulus(x_old, mean, sigma)){
            x_old = x_new;
         } 
      }
      sampling[i] = x_old;
   }
   for(int i=0; i<100000; i++){
      functions[i]=function(sampling[i], mean, sigma);
   }
   std::vector<double>values = progressive(functions, 100000, 200);
   std::vector<double>values2 = progressive2(functions, 100000, 200);
   std::ofstream integral("Integral_final.txt");
   for(int j=0; j<200; j++){
      integral<<j*500<<" "<<values[j]<<" "<<values2[j]<<endl;
   }
   functions.clear();
   return values[199];
   }

   //integral calculation, uncertainty
   double Random :: integral2(std::vector<double> sampling, double mean, double sigma){
      std::vector<double>functions(100000);
      double x_old=0;
      double y;
      double step = 2*mean;
      for(int i=0; i<100000; i++){
         double x_new = Rannyu(-1, 1) * step + x_old;
         if(square_modulus(x_new, mean, sigma) > square_modulus(x_old, mean, sigma)){
            x_old = x_new;
         }
         else{
            y=Rannyu();
            if(y<square_modulus(x_new, mean, sigma)/square_modulus(x_old, mean, sigma)){
            x_old = x_new;
            } 
         }
         sampling[i] = x_old;
      }
      for(int i=0; i<100000; i++){
         functions[i]=function(sampling[i], mean, sigma);
      }
      std::vector<double>values = progressive(functions, 100000, 200);
      std::vector<double>values2 = progressive2(functions, 100000, 200);
      /*std::ofstream integral("Integral.txt");
      for(int j=0; j<200; j++){
         integral<<j*100<<" "<<values[j]<<" "<<values2[j]<<endl;
      }*/
      functions.clear();
      return values2[199];
   }

//SIMULATED ANNEALING

std::vector<double>Random::simulated_annealing(std::vector<double> sampling, double Tin, double Tf, double step_m, double step_s, double mean, double sigma){
      std::vector<double> v(2);
      double T=Tin;   
      double number;
      double sigma_n;
      double mean_n;
      double p;
      double accepted;
      double attempted;
      double accepted_s;
      double attempted_s;
      ofstream H ("H.txt");
      for(int i=0; i<45; i++){
         accepted=0;
         attempted=0;
         accepted_s=0;
         attempted_s=0;
         if(T>0.1){
            T = T-0.05;
         } //decidendo di quanto far scendere T
         else{
            if(T>0.01){
               T=T-0.005;
            }
            else{
               T=T-0.001;
            }
         }
         for(int j=0; j<100; j++){ //100 steps each temperature
            
            number=Rannyu();
            if(number<0.5){// sigma
               attempted_s+=1;
               sigma_n = sigma + Rannyu(-1, 1)*step_s;
               if(integral(sampling, abs(mean), abs(sigma_n)) < integral(sampling, abs(mean), abs(sigma))){
                  sigma = sigma_n;
                  accepted_s+=1;//accepted move for sigma
               }
               else{
                  p=Rannyu();
                  if(p < pow(M_E, -1/T *((integral(sampling, abs(mean), abs(sigma_n)))-integral(sampling, abs(mean), abs(sigma))))){
                     sigma=sigma_n;
                     accepted_s+=1;//accepted move for sigma
                  }
                  else{ //not accepted move for sigma
                  sigma=sigma;
                  }
               
               }
            }
            else{//mean
               attempted+=1;
               mean_n = mean + Rannyu(-1, 1)*step_m;
               if(integral(sampling, abs(mean_n), abs(sigma)) < integral(sampling, abs(mean), abs(sigma))){
                  mean = mean_n;
                  accepted+=1; //accepted for mean
               }
               else{
                  p=Rannyu();
                  if(p < pow(M_E, -1/T *((integral(sampling, abs(mean_n), abs(sigma)))-integral(sampling, abs(mean), abs(sigma))))){
                     mean=mean_n;
                     accepted+=1; //accepted for mean
                  }
                  else{//not accepted for mean
                  }
               }
            }
            H << j+i*100<<" "<<mean<<" "<<sigma<<" "<<integral(sampling, abs(mean), abs(sigma))<<" "<<integral2(sampling, abs(mean), abs(sigma))<<endl;
            }
            std::cout<<T<<" acceptance ratio mean = "<<accepted/attempted<<" Acceptance Ratio sigma = "<<accepted_s/attempted_s<<std::endl;
      }
         v[0]=mean;
         v[1]=sigma;
         return v;
      }



double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
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
