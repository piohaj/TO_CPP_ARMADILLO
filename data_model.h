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

#define NON_SPLITING 1
#define ALL_SPLITING 2
#define COLUMN_SPLITING 3
#define NO_ELEM 5.0


struct vf_opts
{
    string out_file_name;
    double tol;
    int start_row;
    int end_row;
    int max_iters;
    double R_max;
    double C_min;
    int split_strat;
    int pasivity_check;
};

// global configuration structure
extern vf_opts global_conf;
/// ========= VF =========

// struktura wyjsciowa z algorytmu VF
struct SER
{
    cx_mat res;
    cx_mat poles;
    mat h;
    mat d;
    double err;
    mat err_table;
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
