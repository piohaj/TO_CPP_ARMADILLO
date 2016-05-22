#ifndef data_model_h
#define data_model_h

#define ARMA_64BIT_WORD
#include <armadillo>
#include <iostream>
#include <string>
#include<tbb/tbb.h>
#include<mkl_service.h>
#include <vector>

using namespace std;
using namespace arma;
using namespace tbb;

/// ========= VF =========

// struktura wyjsciowa z algorytmu VF
struct SER
{
    cx_mat res;
    cx_mat poles;
    mat h;
    double err;
};

struct input_data
{
    cx_mat f;
    cx_vec s;
};

//============= NETWORK MODEL ==============

// struktura z wartosciami R i L dla bieguna rzeczywistego
struct real_pole_net
{
    double R, L;
};

struct imag_pole_net
{
    double R, L, C, G; 
};


struct Y_network_data
{
   double R, C;
   vector <real_pole_net> real_pole_nets;
   vector <imag_pole_net> imag_pole_nets;
};

#endif
