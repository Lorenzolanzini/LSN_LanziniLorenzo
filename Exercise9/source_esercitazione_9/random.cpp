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
#include<string.h>
#include <vector>
#include  <bits/stdc++.h>

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

//Check member of population
void Random :: Check(std::vector<int> &v){
    int counter = 0;
    if(v[0]==1){}
    else{counter++;}

    if(v[34]==v[0]){}
    else{counter++;}

    for(int i=0; i<34; i++){
        for(int j=i+1; j<34; j++){
            if(v[i]==v[j]){
                counter ++;
            }
        }
    }
    if(counter!=0){
        cout<<"Problem, check this member of the population"<<endl;
    }
}

//Check to see if crossover works
void Random :: Check_c(std::vector<int> &v){
    int counter = 0;
    if(v[0]==1){}
    else{counter++;
    cout<<"1"<<endl;}

    if(v[34]==v[0]){}
    else{counter++;
    cout<<"34"<<endl;}

    for(int i=0; i<34; i++){
        for(int j=i+1; j<34; j++){
            if(v[i]==v[j]){
                counter ++;
                cout<<"el ug"<<endl;
            }
        }
    }
    if(counter!=0){
        cout<<"Problem, check this member of the population crossover"<<endl;
    }
}

//Method for total Length
double Random::cost(std::vector<vector<int>> pop , std::vector<double> &x,std::vector<double> &y, int n){
    double length2=0;
    for(int i=0; i<=34; i++){
        length2 += sqrt((x[pop[n][i]]-x[pop[n][i-1]])*(x[pop[n][i]]-x[pop[n][i-1]]) + (y[pop[n][i]]-y[pop[n][i-1]])*(y[pop[n][i]]-y[pop[n][i-1]]));
    }
    return length2;
}
//function to choose gene
int Random::choose(int M, int p){
    double c;
    c = Rannyu();
    cout<<c<<endl;
    return int(M*pow(c, p) + 1);
}

//permutation
std::vector<int> Random::permutation (std::vector<int> &v){
    unsigned int c;
    unsigned int d;
    int e;
    for(int i =0; i<5; i++ ){
        c = (unsigned int)Rannyu(1, 33);
        d = (unsigned int)Rannyu(1, 33);
        e = v[c];
        //std::cout<<c<<" "<<d<<endl;
        v[c] = v[d];
        v[d] = e;   }

    v[0] = 1;
    v[34] = 1;
    Check(v);
    return v; 
}
//shift genetic mutation
std::vector<int> Random:: shift (std::vector<int> &v, int m){
    std::vector<int> v2(35);
    int h;
    do{ h = (int)Rannyu(1, 16);}while(m+h>33);
    for(int i=0; i<34; i++){
        v2[i]=v[i];
    }
    int e;
    
    for(int i = 0; i<m; i++){
        if(i+h+m<34){
        e=v2[h+m+i];
        v[h+m+i]=v2[i+h];
        v[i+h]=e;}
        else{
            e=v2[i + h + m -33];
            v[i + h + m -33]=v2[i+h];
            v[i+h]=e;
        }
    }
    v[0]=1;
    v[34]=1;
    //cout<<"s"<<endl;
    Check(v);
    return v;
}

//permutation among m contiguos cities
std::vector<int> Random::permutation_m (std::vector<int> &v , unsigned int m){
    std::vector<int> v2(35);
    for(int i=0; i<=34; i++){
        v2[i] = v[i];
    }
    unsigned int c;
    unsigned int d;
    int e;
        c = Rannyu(1, 34);
        do{
        d = Rannyu(1, 34);}while(d-c<m);
        for(unsigned int i=0; i<m; i++){
        if(c+i < 34){
            e = v2[c+i];
            if(d+i<34){
                v2[c+i] = v2[d+i];
                v2[d+i] = e;   }
            else{
                v2[c+i] = v2[d+i - 33];
                v2[d+i-33] = e;
            }
        }
        else{
            e=v2[c+i-33];
            if(d+i<34){
                v2[c+i-33]=v2[d+i];
                v2[d+i]=e;
            }
            else{
                v2[c+i-33]=v2[d+i-33];
                v2[d+i-33]=e;
            }
        }
    }
    v2[0] = 1;
    v2[34] = 1;
    //cout<<"pm"<<endl;
    Check(v2);
    return v2;
    
}

// inversion of the order genetic mutation
std::vector<int> Random:: inversion (std::vector<int> &v){
    std::vector<int> inversion(35);
    int m;
    int j;
    j=(int)Rannyu(1, 33);
    m=(int)Rannyu(1, 33-j);
    for(int i=0; i<35; i++){
        if(i<j || i>j+m){
            inversion[i] = v[i];
        }
        else{
        inversion[i] = v[2*j + m -i];}
    }
    //cout<<"inv"<<endl;
    Check(inversion);
    return inversion;
}


std::vector<double> Random:: ReadAll(const char *filename){
vector<double> v;
ifstream fin(filename);
if(!fin){cout<<"Cannot open file"<<endl;
		exit(11);}
else{
	while(!fin.eof()){
		double val;
		fin>> val;
		v.push_back(val);
	}
}
return v;
}

std::vector<vector<int>> Random:: ReadAllint(const char *filename, int jtot){
vector<vector<int>> v;
ifstream fin(filename);

if(!fin){cout<<"Cannot open file"<<endl;
		exit(11);}
else{
	for(int j=0;j<jtot;j++){
        for(int i=0; i<34; i++){
		int val;
		fin>> val;
        v[j][i] = val;
        }
	}
}
return v;
}

//order function that returns the population (Not useful, use void order2)
std::vector<vector<int>> Random:: order(std::vector<vector<int>> & pop, int n, const char *filenamex, const char *filenamey){
   std::vector<double> x = ReadAll(filenamex);
   std::vector<double> y = ReadAll(filenamey);
   std::vector<vector<int>> pop_aus;
   std::vector<int> temp;
   double costo;
   double costo2;
   for(int j=0; j<n-1; j++){
    for(int l=0; l<n; l++){
    pop_aus.push_back(pop[l]);
   }
         costo = cost(pop, x, y, j );
      for(int i =j+1 ; i<n; i++){
        costo2=cost(pop,x, y, i);
         if(costo2 < costo){
           pop.clear();
           for(int k=0; k<n; k++){
            if(k==i){
                pop.push_back(pop_aus[j]);
            }
            else{
                if(k==j){
                    pop.push_back(pop_aus[i]);
                }
            else{
            pop.push_back(pop_aus[k]);}}
           }
         }
      }
   pop_aus.clear();}
   return pop;
   }
 
//CROSSOVER 
std::vector<int> Random:: crossover1(std::vector<int> &gen1, std::vector<int> &gen2){
    int l = int(Rannyu(1, 33));
    std::vector<int> crossover1;
    std::vector<int> missing_p;
    for(int y=0; y<34-l; y++){
        missing_p.push_back(gen1[l+y]);
    }
    for(int i=0; i<l; i++){
        crossover1.push_back(gen1[i]);
    }
    for(int i=1; i<34; i++){
        for(int y=0; y<34-l; y++){
            if(missing_p[y]==gen2[i]){
                crossover1.push_back(missing_p[y]);
            }
        }
    }
    missing_p.clear();
    crossover1.push_back(1);
    Check_c(crossover1);

    return crossover1;
}  

//Void method for order (BETTER)
void Random::order2(std::vector<vector<int>> &pop, int n, const char*filenamex, const char *filenamey){
    std::vector<double> x = ReadAll(filenamex);
    std::vector<double> y = ReadAll(filenamey);
    std::vector<vector<int>> pop_aus;
    std::vector<int> temp;
    double costo;
    double costo2;
    for(int j=0; j<n-1; j++){
        for(int l=0; l<n; l++){
            pop_aus.push_back(pop[l]);
        }
        costo = cost(pop, x, y, j );
        for(int i =j+1 ; i<n; i++){
            costo2=cost(pop,x, y, i);
            if(costo2 < costo){
            pop.clear();
            for(int k=0; k<n; k++){
                if(k==i){
                    pop.push_back(pop_aus[j]);
                }
                else{
                    if(k==j){
                    pop.push_back(pop_aus[i]);
                }
                    else{
                    pop.push_back(pop_aus[k]);
                    }
                }
            }
        }
    }
    pop_aus.clear();
    }

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
