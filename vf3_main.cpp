#include "my_vectfit.h"
#include <time.h>

int main()
{

    // przygotowanie danych testowych
    mat f_real, f_imag, s_real, s_imag;
    cx_mat f, s;
    cx_vec poles;
    cx_mat weight;
    SER wynik;

    
    int Ns = 1010;
    int N = 3;
    f = zeros<cx_mat>(1, Ns);
    s = 2 * 3.14 * 1.0I * linspace<cx_mat>(0, 350, Ns);

    for ( int k = 0; k < Ns ; k++ )
    {
        cx_double sk = s(k);
        
        f(k) = cx_double(2,0)/(sk+cx_double(5,0)) + cx_double(30, 40)/(sk - cx_double(-100,500)) + cx_double(30,-40)/(sk-cx_double(-100, -500)) + cx_double(0.5, 0);
    } 
/*
    f_real.load( "f_real.dat", raw_ascii );
    f_imag.load( "f_imag.dat", raw_ascii );
    f = cx_mat(f_real, f_imag);
    f_imag.reset();
    f_real.reset();

	s_real.load( "s_real.dat", raw_ascii );
	s_imag.load( "s_imag.dat", raw_ascii );
	s = cx_mat(s_real, s_imag).st();
    s_imag.reset();
    s_real.reset();
*/

    poles = -2 * 3.14 * logspace(0,4,N);

// wlaczenie algorytmu
    clock_t tStart = clock();
    wynik = my_vectorfit3(f, s, poles, weight); 
    double executionTime = (double)(clock() - tStart)/CLOCKS_PER_SEC;

    printf("Czas wykonania algorytmu: %.6fs \n", executionTime); 

    cout << "Poles: \n" <<  wynik.poles << endl;
    cout << "Res: \n" <<  wynik.res << endl;
    cout << "h: " <<  wynik.h << endl;
    cout << "err: " <<  wynik.err << endl;
    return 0;
}


