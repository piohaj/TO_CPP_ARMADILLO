#define ARMA_64BIT_WORD
#include <armadillo>
#include <iostream>
#include <string>
#include<tbb/tbb.h>
#include<mkl_service.h>
#include<vector>


using namespace std;
using namespace arma;
using namespace tbb;

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


// obliczanie wartosci modelu i wstawienie go do 
void parse_SER(SER *input_SER, Y_network_data *output_network_data);
real_pole_net parse_real_pole( cx_mat res, cx_mat poles );
imag_pole_net parse_imag_pole( cx_mat res, cx_mat poles );

