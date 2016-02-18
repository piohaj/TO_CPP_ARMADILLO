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
    cx_mat poles;
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
        Nc = 2;
    }
    else if ( argc == 4 )
    {
        cout << "=== Dane z pliku ===\n";
        N = atoi( argv[1] );
        Nc = atoi( argv[2] );
        Ns = atoi( argv[3] );

        try 
        {
            data = load_data_from_file(N, Nc, Ns);
        }
        catch(string err)
        {
           cout << err << endl;
           return 1;
        }
    }
    else
    {
        cout << "Nieprawidlowa liczba argumentow wejsciowych\n";
        return 1;
    }

    poles = zeros<cx_mat>(Nc, N);
    mat bet = linspace<mat>(imag(data.s(0)), imag(data.s(Ns-1)), N/2);
    
    for ( int mm = 0; mm < Nc; mm++ )
    {
        int m = 0;
        for ( int n = 0; n < N-1; n=n+2 )
        {
            double alf = -bet(m)*1e-2;
            poles(mm, n) = cx_double(alf, bet(m));
            poles(mm, n+1) = conj(poles(n));
            m++;
        }
    }

    wall_clock timer;
    // wlaczenie algorytmu
    cout << "Vector fitting" << endl;
    timer.tic();
    int iter = 1;
    for ( iter = 1; iter < 11; iter++ )
    {
//        poles.print("Input poles: ");
        wynik = my_vf_all_splitting(data.f, data.s, poles); 
        poles = wynik.poles;
        
        cout << "Iter: " << iter << endl;
        cout << "Err: " << wynik.err << endl;
        if ( wynik.err < 1e-5 )
        {
            break;
        }
    }
    double executionTime = timer.toc();

    printf("Czas wykonania algorytmu: %.6fs \n", executionTime); 

    wynik.poles.print("poles=");
    wynik.res.print("residues=");
    wynik.h.print("h=");
    cout << "RMS-err= " << wynik.err << endl;
    cout << "Iter: " << iter << endl;

    //zapis statystyk do pliku
    fstream plik;
    plik.open("stats_cpp_all_split.txt", ios::out | ios::app);
    plik << N << ";" << Nc << ";" << Ns << ";" << iter << ";" << wynik.err << ";" << executionTime << endl;
    plik.flush();

    plik.close();

    return 0;
}


