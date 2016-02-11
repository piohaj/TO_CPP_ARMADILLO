#define ARMA_64BIT_WORD
#include <armadillo>
#include <iostream>
#include <string>
#include<tbb/tbb.h>

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
    cx_mat s;
};

SER my_vectorfit3(const cx_mat& f, const cx_mat& s, cx_vec poles, cx_mat weight);
cx_mat logspace(double a, double b, int n);
int sign( double x );

//funkcje do generowania lub ladowaniai z pliku danych testowych
input_data prepare_sample_data();
input_data load_data_from_file( int N, int Nc, int Ns );

// klasa potrzebna do zrownoleglenia dekompozycji QR w TBB
class QR_calculation
{
    const cx_mat *A;
    const cx_mat *f;
    int N;
    int Ns;
    mat *AA_poles;
    mat *bb_poles;
    
public:
    QR_calculation( const cx_mat *A_in, const cx_mat *f_in, int N_in, int Ns_in, mat *AA_poles_out, mat *bb_poles_out)
                  : A(A_in), 
                    f(f_in), 
                    N(N_in), 
                    Ns(Ns_in), 
                    AA_poles(AA_poles_out), 
                    bb_poles(bb_poles_out)
                  {}

    // zadanie dla jednego watku - w TBB definiowany poprzez przeciazenie ()
    void operator() ( const blocked_range<int>& r ) const;
};
