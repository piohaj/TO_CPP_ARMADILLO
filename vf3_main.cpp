#include "my_vectfit.h"

int main()
{

    // przygotowanie danych testowych
    mat f_real, f_imag, s_real, s_imag;
    cx_mat f, s;
    cx_vec poles;
    cx_mat weight;
    SER wynik;

    
    int Ns = 1010000;
    int N = 3;
    f = zeros<cx_mat>(2, Ns);
    s = 2 * 3.14 * 1.0I * linspace<cx_mat>(0, 550, Ns);

    for ( int k = 0; k < Ns ; k++ )
    {
        cx_double sk = s(k);
        
        f(0, k) = cx_double(2,0)/(sk+cx_double(5,0)) + cx_double(30, 40)/(sk - cx_double(-100,500)) + cx_double(30,-40)/(sk-cx_double(-100, -500)) + cx_double(0.5, 0);
    } 

    for ( int kk = 0; kk < Ns ; kk++ )
    {
        cx_double sk = s(kk);
        
        f(1, kk) = cx_double(3,0)/(sk+cx_double(5,0)) + cx_double(300, 40)/(sk - cx_double(-100,500)) + cx_double(300,-40)/(sk-cx_double(-100, -500)) + cx_double(0.9, 0);
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

    wall_clock timer;
    // wlaczenie algorytmu
    cout << "Vector fitting" << endl;
    timer.tic();
    int iter = 1;
    for ( iter = 1; i < 11; i++ )
    {
        wynik = my_vectorfit3(f, s, poles, weight); 
        
        if ( wynik.err.max() < 1e-5 )
        {
            break;
        }
    }
    double executionTime = timer.toc();

    printf("Czas wykonania algorytmu: %.6fs \n", executionTime); 

    wynik.poles.print("poles=");
    wynik.res.print("residues=");
    wynik.h.print("h=");
    wynik.err.print("RMS-err=");
    cout << "Iter: " << iter << endl;

    return 0;
}


