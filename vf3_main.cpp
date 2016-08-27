#include "my_vectfit.h"
#define VF_REPEAT 3


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
    MKL_Set_Num_Threads(1);
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

    // wlaczenie algorytmu
    cout << "Vector fitting " << VF_REPEAT << " times" << endl;
    wall_clock timer;
    double exec_time=0;
    int iter = 0;

    for ( int k=0; k < VF_REPEAT ; k++ )
    {
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
    
        timer.tic();
        iter = 1;
        for ( iter = 1; iter < 4; iter++ )
        {
    //        poles.print("Input poles: ");
            wynik = my_vf_all_splitting( &data.f, &data.s, &poles); 
            poles = wynik.poles;
            
            cout << "Iter: " << iter << endl;
            cout << "Err: " << wynik.err << endl;
//            if ( wynik.err < -40 )
//            {
//                break;
//            }
        }
        double executionTime = timer.toc();
        exec_time = exec_time + executionTime;
        cout << "Sample exec time: "<< executionTime << endl;
    }

    exec_time = exec_time / VF_REPEAT;
    printf("Sredni czas wykonania algorytmu po %d wywolaniach: %.6fs \n", VF_REPEAT, exec_time); 

//    wynik.poles.print("poles=");
//    wynik.res.print("residues=");
//    wynik.h.print("h=");
    cout << "RMS-err= " << wynik.err << endl;
    cout << "Iter: " << iter << endl;

    //zapis statystyk do pliku
    fstream plik;
    plik.open("stats_cpp_all_split_parallel.txt", ios::out | ios::app);
    plik << N << ";" << Nc << ";" << Ns << ";" << iter << ";" << wynik.err << ";" << exec_time << endl;
    plik.flush();

    plik.close();

    return 0;
}


