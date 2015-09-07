#include <armadillo>
#include <iostream>

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
    double h;
    double err;
};

cx_mat my_vectorfit3(cx_mat f, cx_mat s, cx_vec poles, cx_mat weight);
cx_mat logspace(double a, double b, int n);
int sign( double x );
