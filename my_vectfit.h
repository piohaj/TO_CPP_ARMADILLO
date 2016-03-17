#define ARMA_64BIT_WORD
#include <armadillo>
#include <iostream>
#include <string>
#include<tbb/tbb.h>
#include<mkl_service.h>
#include<cstdlib>

using namespace std;
using namespace arma;
using namespace tbb;

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

SER my_vf_column_splitting(const cx_mat *f, const cx_vec *s, cx_mat *poles);

cx_mat logspace(double a, double b, int n);
int sign( double x );

//funkcje do generowania lub ladowaniai z pliku danych testowych
input_data prepare_sample_data();
input_data load_data_from_file( int N, int Nc, int Ns );

class vf_column
{
    const cx_mat *f;
    const cx_vec *s;
    cx_mat *poles1;
    SER *wynik;

public:
    vf_column( const cx_mat *f_in, const cx_vec *s_in, cx_mat *poles_in, SER *wynik_in )
         : f(f_in),
           s(s_in),
           poles1(poles_in),
           wynik(wynik_in)
         {}

    void operator() ( const blocked_range<int>& r ) const;
};
