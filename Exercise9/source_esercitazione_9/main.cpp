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
#include <algorithm>
#include <random>
#include <cmath>


using namespace std;

int niterations = 500;  //number of iterations of the GA
//TSP problem
//Methods for genetic algorithm (crossover, order, mutations,..) are
//as methods in class random
//See random.h and random.cpp

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

   //creation of 34 random cities coordinates, commented


   /*ofstream fout("Squarex.dat");
   for(int i=0; i<34; i++){
      fout<<rnd.Rannyu()<<endl;
   }
   ofstream fouty("Squarey.dat");
   for(int i=0; i<34; i++){
      fouty<<rnd.Rannyu()<<endl;
   }
   
   double angolo;
    ofstream fout3("Circumferencey.dat");
    ofstream fout2("Circumferencex.dat");
    for(int i=0; i<34; i++){
      angolo = rnd.Rannyu()*2*M_PI;
      fout2<<0.5*cos(angolo)+1<<endl;
      fout3<<0.5*sin(angolo)+1<<endl;
   }*/
  

   /*//Population
   std::vector<int> in(35);
   for(int i=0; i<34; i++){
      in[i]=i+1;
   }
   in[34]= 1;
   Check(in);
   std::vector<int> population(35);
   ofstream fout3("population.txt", ios::app);
   population = shift(in, 34);
   for(int i=0; i<10; i++){
   fout3<<endl;
   population = shift(population, 34);
   Check(population);
   for(int i=0; i<35; i++){
      fout3<<population[i]<<" ";
   }
   }
   fout3<<endl;
   population = permutation_m(in, 12);
   for(int i=0; i<35; i++){
      fout3<<population[i]<<" ";
   }*/
    

   int counter=0;

   //reading coordinates from file 
   std::vector<double> coordinates_x;
   std::vector<double> coordinates_y;
   const char *filenamex = "Squarex.dat"; //or Circumferencex.dat
   const char *filenamey = "Squarey.dat"; ////or Circumferencey.dat

   //reading data from file
   coordinates_x = rnd.ReadAll(filenamex); 
   coordinates_y = rnd.ReadAll(filenamey);

   //creating population
   std::vector<std::vector<int>> popolazione;
   std::vector<int> popolazione_aus(35);

   //initial population
   for(int k =0; k<34; k++){
      popolazione_aus[k] = k + 1; //first sequence 1, 2,3, 4.. , 34, 1
   }
   popolazione_aus[34] = 1;
   popolazione.push_back(popolazione_aus);
   for(int j=1; j<=15; j++){
      counter=0;
      random_shuffle(popolazione_aus.begin(), popolazione_aus.end()); 
      popolazione.push_back(popolazione_aus);
      for(int l=0; l<34;l++){
         if(counter==0){
            if(l!=0 && popolazione[j][l]==1){
               popolazione[j][l] = popolazione[j][0];
               popolazione[j][0]=1;
               counter+=1;
            }
         }
         else{
            if(l!=34 && popolazione[j][l]==1){
               popolazione[j][l]=popolazione[j][34];
               popolazione[j][34]=1;
            }
         }
      }
      rnd.Check(popolazione[j]); //Checking if the sequence of cities respects costrains
   }


   
  
   int n;
   n=16; //size of the initial population
   //std::vector<vector<int>> sequence;

   //parents
   std::vector<int>gen1; 
   std::vector<int>gen2;
   //sons
   std::vector<int>fig1; 
   std::vector<int>fig2;

   int f=0;
   cout<<endl;
   int j=0;
   int k=0;
   int niterations = 20000;
   //sequence = rnd.order(popolazione, n, filenamex, filenamey);
   double prob;
   for(int i=0; i<niterations; i++){
      rnd.order2(popolazione, n, filenamex, filenamey);
      prob=rnd.Rannyu();
      if(prob<0.6){ //60% probability for crossover
         do{  
            j = (int)n*pow(rnd.Rannyu(), 5) +1;
            k = (int)n*pow(rnd.Rannyu(), 5) +1;
         }while(k==j);
         
         for(int s=0; s<35; s++){
            gen1.push_back(popolazione[j][s]);
            gen2.push_back(popolazione[k][s]);
         }

         fig1=rnd.crossover1(gen1, gen2);
         fig2=rnd.crossover1(gen2, gen1);
         popolazione.push_back(fig1);
         popolazione.push_back(fig2);
         //sequence.clear();
         n+=2;
      }
      else{
      f++;
      if(prob< 0.7){ //10% probability for permutation among m contiguos cities
         j = (int)n*pow(rnd.Rannyu(), 1) +1;
         for(int s=0; s<35; s++){
            gen1.push_back(popolazione[j][s]);
         }
         fig1=rnd.permutation_m(gen1, (int)rnd.Rannyu(1, 33));
         popolazione.push_back(fig1);
         //sequence.clear();
         n+=1;
      }
      else{
         if(prob<0.8){//10% probability of permutation(5 cities)
            j = (int)n*pow(rnd.Rannyu(), 1) +1;
            for(int s=0; s<35; s++){
               gen1.push_back(popolazione[j][s]);
            }
            fig1=rnd.permutation(gen1);
            popolazione.push_back(fig1);
            //sequence.clear();
            n+=1;
         }
         else{
            if(prob<0.9){ //10% probabilityu of shift
               j = (int)n*pow(rnd.Rannyu(), 1) +1;
               for(int s=0; s<35; s++){
                  gen1.push_back(popolazione[j][s]);
               }
               fig1=rnd.shift(gen1, (int)rnd.Rannyu(2, 10));
               popolazione.push_back(fig1);
               //sequence.clear();
               n+=1;
            }
            else{ //10% probability of inversion
               j = (int)n*pow(rnd.Rannyu(), 1) +1;
               for(int s=0; s<35; s++){
                  gen1.push_back(popolazione[j][s]);
               }
               fig1=rnd.inversion(gen1);
               popolazione.push_back(fig1);
               //sequence.clear();
               n+=1;

            }

         }
      }
   }
      
      if(n>101){ //erasing "expensive" members of population
         rnd.order2(popolazione, n, filenamex, filenamey);
         std::vector<vector<int>>::iterator it = popolazione.end();
      do{
         popolazione.erase(it);
         it--;
      }while(popolazione.size()>101);
      n=100; 
      }

   
         gen1.clear();
         gen2.clear();
         fig1.clear();
         fig2.clear();
 
   }
   //ordering
   rnd.order2(popolazione, n, filenamex, filenamey);

   //saving final order pop
   ofstream final ("final_ordered_pop.txt");
   for(int i=0; i<n; i++){
      for(int k=0; k<35; k++){
         final<<popolazione[i][k]<< " ";}
         final<<endl;
      }

   std::cout<<rnd.cost(popolazione, coordinates_x, coordinates_y, 0)<<endl;
   std::cout<<f<<endl;

   /*ofstream cordinatefinali("final-cord-sq.txt");
   for(int i=0; i<35; i++){
      cordinatefinali<<i<<" "<<coordinates_x[popolazione[0][i]]<<" "<<coordinates_y[popolazione[0][i]]<<endl;
   }*/
   //rnd.order2(popolazione, n, filenamex, filenamey);
   //std::cout<<rnd.cost(popolazione, coordinates_x, coordinates_y, 0)<<endl;
   double average=0;
   for(int i=0; i<100; i++){
      average+=rnd.cost(popolazione, coordinates_x, coordinates_y, i);
   }
   average/=(double)100;
   ofstream square("Square-costs.txt", ios::app);
   square<<niterations<<" "<<rnd.cost(popolazione, coordinates_x, coordinates_y, 0)<<" "<<average<<endl;
   /*
   ofstream circ("Circ-costs.txt", ios::app);
   circ<<niterations<<" "<<rnd.cost(popolazione, coordinates_x, coordinates_y, 0)<<" "<<average<<endl;
   */

   popolazione.clear();
   

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
