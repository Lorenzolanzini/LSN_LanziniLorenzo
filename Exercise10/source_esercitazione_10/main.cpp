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
#include "mpi.h"
#include <fstream>
#include <string>
#include "random.h"
#include <algorithm>
#include <random>
#include <cmath>

//make to compile
// mpiexec -np "NContinents" ./main.exe
//TSP methods in random.h, random.cpp

using namespace std;

int niterations = 15000; // # iterations of the program
int n_migr = 125; //The number of migrations will be niterations/n_migr 
std::vector<int> migrazione1;
std::vector<int> migrazione2;
int mig1[51];
int mig2[51];
int world_best[51];
std::vector<int> best_path;
int best_cont;
std::vector<int> aus2;
std::vector<int> aus1;
std::vector<vector<int>>::iterator it2;
std::vector<vector<int>>::iterator it;
std::vector<int>gen1;
std::vector<int>gen2;
std::vector<int>fig1;
std::vector<int>fig2;


int main (int argc, char *argv[]){
   int size, rank;

   double min_length; //World best path's length  
   

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   double lengths[size]; //Here I will save the lengths of best paths for each continent
   double tstart= MPI_Wtime();

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("primes32001.in");
   if (Primes.is_open()){
      for(int i=0;i<=rank;i++) Primes >> p1 >> p2 ; //each thread has its primes
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

   std::vector<double> coordinates_x;
   std::vector<double> coordinates_y;
   const char *filenamex = "capitali americane x.txt";
   const char *filenamey = "capitali americane y.txt";
   coordinates_x = rnd.ReadAll(filenamex);
   coordinates_y = rnd.ReadAll(filenamey);
   std::vector<std::vector<int>> popolazione;
   std::vector<int> popolazione_aus(51);

   int counter=0;
   int n;
   n=16; // # To create the initial population with random shuffle

   for(int k =0; k<50; k++){
      popolazione_aus[k] = k + 1;
   }
   popolazione_aus[50] = 1;
   popolazione.push_back(popolazione_aus);
   for(int j=1; j<n; j++){
      counter=0;
      random_shuffle(popolazione_aus.begin(), popolazione_aus.end());
      popolazione.push_back(popolazione_aus);
      for(int l=0; l<50;l++){
         if(counter==0){
            if(l!=0 && popolazione[j][l]==1){
               popolazione[j][l] = popolazione[j][0];
               popolazione[j][0]=1;
               counter+=1;
            }
         }
         else{
            if(l!=50 && popolazione[j][l]==1){
               popolazione[j][l]=popolazione[j][50];
               popolazione[j][50]=1;
            }
         }
      }
      rnd.Check(popolazione[j]);}
      //std::cout<<"rank "<<rank<<" popolazione fatta "<<" "<<popolazione[1][10]<<endl;
   

      int f=0;
      int j=0;
      int k=0;
   
      double prob;

      //TSP algortithm 
      for(int i=2; i<niterations; i++){
         //cout<<rank<<" start running program"<<i<<endl;
         n=popolazione.size();
         rnd.order2(popolazione, n, coordinates_x, coordinates_y);
      
         prob=rnd.Rannyu();
      
         if(prob<0.70){ // 70% probabilty for crossover
            do{  
               j = (int)n*pow(rnd.Rannyu(), 5) ;
               k = (int)n*pow(rnd.Rannyu(), 5) ;
            }while(k==j);
            for(int s=0; s<51; s++){
               gen1.push_back(popolazione[j][s]);
               gen2.push_back(popolazione[k][s]);
            }

            fig1=rnd.crossover1(gen1, gen2); //two sons are added to the population
            fig2=rnd.crossover1(gen2, gen1);
            rnd.Check_c(fig1);
            rnd.Check_c(fig2);
            popolazione.push_back(fig1);
            popolazione.push_back(fig2);
            n+=2;
         }
         else{
            f++;
            if(prob< 0.80){ // 10 % probability of Permutation among  m  contiguous cities 
               j = (int)n*pow(rnd.Rannyu(), 2) ;
               for(int s=0; s<51; s++){
                  gen1.push_back(popolazione[j][s]);
               }
               fig1=rnd.permutation_m(gen1);
               rnd.Check(fig1);
               popolazione.push_back(fig1);
               n+=1;
            }
            else{
               if(prob<0.90){ // 10% probability of permutation (5 cities)
                  j = (int)n*pow(rnd.Rannyu(), 2) ;
                  for(int s=0; s<51; s++){
                     gen1.push_back(popolazione[j][s]);
                  }
                  fig1=rnd.permutation(gen1);
                  rnd.Check(fig1);
                  popolazione.push_back(fig1);
                  n+=1;
               }
               else{
                  if(prob<0.92){ // 2% probability of shift
                     j = (int)n*pow(rnd.Rannyu(), 2) ;
                     for(int s=0; s<51; s++){
                        gen1.push_back(popolazione[j][s]);
                     }
                     fig1=rnd.shift(gen1, (int)rnd.Rannyu(2, 10));
                     //rnd.Check(fig1);
                     popolazione.push_back(fig1);
                     //sequence.clear();
                     n+=1;
                  }
                  else{ // 8% probability of inversion of the order
                     j = (int)n*pow(rnd.Rannyu(), 2) ;
                     for(int s=0; s<51; s++){
                        gen1.push_back(popolazione[j][s]);
                     }
                     fig1=rnd.inversion(gen1);
                     //rnd.Check(fig1);
                     popolazione.push_back(fig1);
            
                     n+=1;
                  }
               }
            }
         }
         //std::cout<<"rank "<<rank<<"it. done"<<endl;
      
         // Erase "Expensive" members of population
         if(n>101){
            rnd.order2(popolazione, n, coordinates_x, coordinates_y);
            it = popolazione.end();
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
   
   
      //migration
      if(i % n_migr == 0 && size>1){
			
         cout<<"starting migration"<<endl;
         int cont1,cont2;
			if(not rank){//Choose two random continents that will echange their best path
			   cont1=floor(rnd.Rannyu()*size);
			   do  cont2=floor(rnd.Rannyu()*size); while(cont2==cont1);
			}
			
			MPI_Bcast(&cont1,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(&cont2,1,MPI_INT,0,MPI_COMM_WORLD);
			
			if(rank==cont1){
            //migrazione1=rnd.best_p(popolazione);
				for(int k=0;k<51;k++){
					mig1[k] = popolazione[0][k];
				}
            //migrazione1.clear();
			}
			if(rank==cont2){
				//migrazione2=rnd.best_p(popolazione);
				for(int k=0;k<51;k++){
					mig2[k] = popolazione[0][k];
				}
            //migrazione2.clear();
			}
			MPI_Bcast(mig1,51,MPI_INT,cont1,MPI_COMM_WORLD); //send genes to all the other continents
			MPI_Bcast(mig2,51,MPI_INT,cont2,MPI_COMM_WORLD);
         cout<<"sending members"<<endl;
			if(rank==cont1){//on cont1, add the best of cont2 and erase the worst element
				cout<<"--"<<mig2[0]<<"--"<<endl;
            for(int k=0;k<51;k++){
               aus1.push_back(mig2[k]);
				}
            cout<<"--"<<aus1[0]<<"--"<<endl;
            popolazione.push_back(aus1);
            aus1.clear();
            rnd.order2(popolazione, popolazione.size(), coordinates_x, coordinates_y);
            it2 = popolazione.end();
            popolazione.erase(it2);
            n=popolazione.size();        
            cout<<"gene added"<<endl;
			}
			if(rank==cont2){//on cont2, add the best of cont1 and erase the worst element
				
            for(int k=0;k<51;k++){
               aus2.push_back(mig1[k]);
				}
            popolazione.push_back(aus2);
            aus2.clear();
            rnd.order2(popolazione, popolazione.size(), coordinates_x, coordinates_y);
            it2 = popolazione.end();
            popolazione.erase(it2);
            n=popolazione.size();
            cout<<"gene added"<<endl;
			}
         cout<<"------"<<i<<" "<<rank<<" migration done"<<endl;
         //migration done
		}

   }   
      double single_length;
      single_length=rnd.cost(popolazione, coordinates_x, coordinates_y, 0);
      //builing lengths[size]
	   MPI_Gather(&single_length, 1 ,MPI_DOUBLE, lengths ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	   
	   cout<<rank<<"gather 1"<<endl;
      if (not rank){//core 0

		   //finding the minimum, and find the corresponding index
		
		   min_length=lengths[0];
         for(int k=1; k<size; k++){
            if(lengths[k]<min_length){
               min_length=lengths[k];
            }   
            cout<<"lengths done"<<endl;  
         }

	   for(int l=0;l<size;l++) if(lengths[l]==min_length){ best_cont=l;}//Writing which continent found best path
	}
      //comunicating the best length and continent to the world
	MPI_Bcast(&min_length,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&best_cont,1,MPI_INT,0,MPI_COMM_WORLD);
	   
   //cout<<rank<<" Bcast both done"<<endl;
   //if(not rank) cout<<best_cont<<endl;
	   
   if (rank==best_cont){
      //best_path=rnd.best_p(popolazione); 
      for(int k=0; k<51; k++){
            world_best[k]=popolazione[0][k];
      }
   }      
   cout<<rank<<"wbest done"<<endl;
    
   cout<<"Time taken : "<<MPI_Wtime() - tstart<<endl;

   if(rank==best_cont){
      ofstream fout("best_path_cord_15000_5.txt");
      for(int i=0; i<51; i++){
         fout<<coordinates_x[popolazione[0][i]]<<" "<<coordinates_y[popolazione[0][i]]<<endl;
      }
   }

   rnd.SaveSeed();
   if(not rank){
      ofstream fout2 ("Min_lengths_15000.txt", ios::app);
      fout2<<min_length<<endl;
   }
   MPI_Finalize();   
   
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
