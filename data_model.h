#ifndef data_model_h
#define data_model_h

#define ARMA_64BIT_WORD
#include <armadillo>
#include <iostream>
#include <string>
#include<tbb/tbb.h>
#include<mkl_service.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <sys/stat.h>
#include <ctime>

using namespace std;
using namespace arma;
using namespace tbb;

#define NON_SPLITING 1
#define ALL_SPLITING 2
#define COLUMN_SPLITING 3
#define NO_ELEM datum::inf

#define TOUCHSTONE_FILE 0
#define RAW_FILE 1
#define PI 3.1416


struct vf_opts
{
    string out_file_name;
    string in_file_name;
    double tol;
    int optim_model;
    int min_row;
    int max_row;
    int max_iters;
    int calc_parallel_RC;
    double R_max;
    double C_min;
    int split_strat;
    int pasivity_check;
    int spice_simulation;
    string spice_program_loc;
    int gnuplot_generation;
};

/// ========= VF =========

// struktura wyjsciowa z algorytmu VF
struct SER
{
    cx_mat res;
    cx_mat poles;
    mat h;
    mat d;
    double err;
    vec err_table;
};

struct input_data
{
    cx_mat f;
    cx_vec s;
    vec freq;
};

struct raw_params
{
    int Nc_ports;
    int Ns;
    int last_line;
};

//============= NETWORK MODEL ==============
// struktura z wartosciami R i L dla bieguna rzeczywistego
struct real_pole_net
{
    double R, L;
    double res, pole;
};

struct imag_pole_net
{
    double R, L, C, G; 
    cx_double res, pole;
};


struct Y_network_data
{
   double R, C;
   vector <real_pole_net> real_pole_nets;
   vector <imag_pole_net> imag_pole_nets;
};


// ===================== KONFIGURACJA TOUCHSTONE =========

struct touchstone_conf
{
    bool is_touchstone;
    string freq_unit;
    string data_type;
    string data_type2;
    double R0;
    int Ns;
};

// ================= Dane do gnuplot ====================
struct gnuplot_data
{
    cx_mat input_data;
    cx_mat simulation_data;
    vec freq;
};


#endif
