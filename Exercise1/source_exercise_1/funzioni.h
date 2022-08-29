#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

double CalcolaMedia(const double *, int);
double CalcolaVarianza(const double *, int);

double CalcolaMediana(const double *, int);

double *ReadDataFromFile(const char *, int);

void scambiaByValue(double, double);
void scambiaByRef(double &, double &);
void scambiaByPointer(double *, double *);

void selection_sort(double *, int);
void Print(double *data, int ndata);
void Print(const char *Filename, double *data, int ndata);
void test_media();
void test_varianza();
void test_mediana();
