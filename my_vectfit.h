#define ARMA_64BIT_WORD
#include <armadillo>
#include <iostream>
#include <string>
#include<mkl_service.h>

using namespace std;
using namespace arma;

struct opts
{
    int relax;
    int stable;
    int asymp;
    int spy1;
    int spy2;
    int logx;
    int logy;
    int errplot;
    int phaseplot;
    int skip_pole;
    int skip_res;
    int complx_ss;
    int legend;
};

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

SER my_vectorfit(const cx_mat& f, const cx_vec& s, cx_vec poles);

SER my_vf_all_splitting(const cx_mat *f, const cx_vec *s, cx_mat *poles);

cx_mat logspace(double a, double b, int n);
int sign( double x );

//funkcje do generowania lub ladowaniai z pliku danych testowych
input_data prepare_sample_data();
input_data load_data_from_file( int N, int Nc, int Ns );

void calculate_poles_res_set(const cx_mat *f, const cx_vec *s, cx_mat *poles1, SER *wynik, int rr);
