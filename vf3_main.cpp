#include "my_vectfit.h"


// program na wejsciu przyjmuje 3 dane (w celu wczytania odpowiedniego benczmarka):
// $1 - N rzad przyblizenia
// $2 - Nc liczba portow badanego ukladu
// $3 - Ns liczba probek pomiarowych
// np.
// vf3 4 2 100 - bedzie probowalo czytac pliki: f_real_N4_Nc2_Ns100.dat itd.
// w przydapku niepodania wartosci wejsciowych zostanie uruchomiony przebieg testowy 
int main(int argc, char* argv[])
{

    // przygotowanie danych testowych
    input_data data;
    cx_vec poles;
    cx_mat weight;
    SER wynik;
    int N = 0,
        Ns = 0,
        Nc = 0;

    if ( argc == 1 )
    {
        cout << "===Test for sample data===\n";
        data = prepare_sample_data();
        N = 3;
        Ns = 101;
    }
    else if ( argc == 4 )
    {
        N = atoi( argv[1] );
        Nc = atoi( argv[2] );
        Ns = atoi( argv[3] );
        cout << N << " " << Nc << " " << Ns << endl;
        return 2;
    }
    else
    {
        cout << "TODO\n";
        return 1;
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
//    f_real.load("test_real.dat");
//    f_imag.load("test_imag.dat");
//    f = cx_mat(f_real, f_imag);

//    f.print("f=");


    //poles = -2 * 3.14 * logspace(0,4,N);
    // complex starting poles
    poles = zeros<cx_vec>(N);
    mat bet = linspace<mat>(1, imag(data.s(Ns-1)), N/2);
    int m = 0;
    for ( int n = 0; n < N-1; n=n+2 )
    {
        double alf = -bet(m)*1e-2;
        poles(n) = cx_double(alf, bet(m));
        poles(n+1) = conj(poles(n));
        m++;
    }


    wall_clock timer;
    // wlaczenie algorytmu
    cout << "Vector fitting" << endl;
    timer.tic();
    int iter = 1;
    for ( iter = 1; iter < 11; iter++ )
    {
        poles.print("Input poles: ");
        wynik = my_vectorfit3(data.f, data.s, poles, weight); 
        poles = wynik.poles;
        
        cout << "Iter: " << iter << endl;
        cout << "Err: " << wynik.err.max() << endl;
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


