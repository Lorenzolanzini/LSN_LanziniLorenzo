#include "funzioni.h"

using namespace std;

double CalcolaMedia(const double *data, int size){
	double media=0;
	for(int k=0;k<size;k++){
	media+= data[k];
	}
	media=media/double(size);
	return media;
}
double CalcolaVarianza(const double *data, int size){
	double varianza=0;
	for(int k=0; k<size; k++){
		varianza+= (data[k]-CalcolaMedia(data, size))*(data[k]-CalcolaMedia(data, size));
	}
	varianza=sqrt(varianza/double(size-1));
	return varianza;
}
double *ReadDataFromFile(const char *Filename, int ndata){
	double*data= new double[ndata];
	ifstream fin(Filename);
	if(!fin) {
		cout <<"Cannot open file"<<endl;
		  exit(0);}
        else{
	for(int k=0; k<ndata ; k++){
		fin>>data[k];
		if(fin.eof()){cout<<"End of file reached Exiting"<<endl;
		exit(0);}
		}
	}return data;}

void ScambiaByValue(double a, double b){
	double temp=a;
	a=b;
	b=temp;
}
void ScambiaByReference(double &a, double &b){
	double temp=a;
	a=b;
	b=temp;
}
void ScambiaByPointer(double*a, double*b){
	double temp=*a;
	*a=*b;
	*b=temp;
}

void selection_sort(double *vec, int size){
	int imin=0;
	double min=0;
for(int j=0; j<size-1; j++){
	imin=j;
	min=vec[imin];
	for(int i=j+1; i<size; i++){
		if(vec[i]<min){
			min=vec[i];
			imin=i;
		}
	}
	double c= vec[j];
	vec[j]=vec[imin];
	vec[imin]=c;}
}


double CalcolaMediana(const double *vec, int size){
	double *vcopy= new double[size];
	for(int k=0; k<size; k++){
		vcopy[k]=vec[k];
	}
	double mediana=0;
	selection_sort (vcopy, size);
	if(size%2==0){
		mediana=(vcopy[(size/2) -1] + vcopy[size/2])/2;
	}else{
		mediana=vcopy[size/2];
	}
	delete[] vcopy;
	return mediana;
}

void Print(const char *Filename, double *data, int ndata){
	ofstream fout(Filename);
	for(int k=0; k<ndata; k++) fout <<data[k] <<endl;
}

void Print(double *data, int ndata){
for(int k=0; k<ndata; k++){
	cout<<data[k]<<endl;}}

bool is_close(double a, double b, double epsilon=1e-7){
	return fabs(a-b)<epsilon;}

void test_media(){
double x[]={1, 2, 3, 4};
assert(is_close(CalcolaMedia(x, 4), 2.5));
}

void test_mediana(){
	
double x[]={1, 2, 3, 4};
assert(is_close(CalcolaMediana(x, 4), 2.5));		
}

void test_varianza(){

double x[]={1, 2, 3, 4};
assert(is_close(CalcolaVarianza(x, 4), 1.29099445));
}