#include "my_vectfit3.h"

int main()
{

    // przygotowanie danych testowych
    cx_mat f, s, poles;
    cx_mat weight;

    
    int Ns = 101;
    f = zeros<cx_mat>(1, Ns);
    s = 2 * 3.14 * 1i * logspace(0, 4, Ns);


    for ( int k = 0; k < Ns ; k++ )
    {
        cx_double sk = s(k);
        
        f(k) = cx_double(2,0)/(sk+cx_double(5,0)) + cx_double(30, 40)/(sk - cx_double(-100,500)) + cx_double(30,-40)/(sk-cx_double(-100, -500)) + cx_double(0.5, 0);
    } 

    weight = ones<cx_mat>(1, Ns);
    poles = -2 * 3.14 * logspace(0,4,3);


// wlaczenie algorytmu
    vectorfit3(f, s, poles, weight); 

    return 0;
}


