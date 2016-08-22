#define ARMA_64BIT_WORD
#include <armadillo>
#include <iostream>
#include <string>

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
    double qr_time;
};

struct input_data
{
    cx_mat f;
    cx_mat s;
};

SER my_vectorfit3(const cx_mat& f, const cx_mat& s, cx_vec poles, cx_mat weight);
cx_mat logspace(double a, double b, int n);
int sign( double x );

//funkcje do generowania lub ladowaniai z pliku danych testowych
input_data prepare_sample_data();
input_data load_data_from_file( int N, int Nc, int Ns );
