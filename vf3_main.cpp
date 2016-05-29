#include "data_model.h"
#include "network_model.h"
#include "additional_features.h"
#include "my_vectfit_non.h"
#define VF_REPEAT 1

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
        Ns = 10;
        Nc = 4;
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


    double exec_time=0;
    int iter;
    wall_clock timer;
    // wlaczenie algorytmu
    cout << "Vector fitting " << VF_REPEAT << " times" << endl;
    for ( int k = 0; k < VF_REPEAT ; k++ )
    {
        timer.tic();
        wynik = vf_high_level( data.f, data.s, ALL_SPLITING, N-2, N+2 ); 
        double executionTime = timer.toc();
        cout<< "Exec one: "<< executionTime << endl;
        exec_time = exec_time + executionTime;
    }

    exec_time = exec_time / VF_REPEAT;

    printf("Sredni czas wykonania algorytmu po %d wywolaniach: %.6fs \n", VF_REPEAT, exec_time); 

    cout << "\n\n\n";
    wynik.poles.print("poles=");
    wynik.res.print("residues=");
    wynik.h.print("h=");
    wynik.d.print("d=");
    cout << "RMS-err(wybrany)= " << wynik.err <<endl;

    // utworzenie modelu cir i zapis do pliku
    //ofstream myfile;
    //myfile.open("test_mgr.cir");
    //create_model_netlist( &wynik, Nc, myfile);

    //myfile.close();

    //zapis statystyk do pliku
/*
    fstream plik;
    plik.open("stats_cpp_parallel.txt", ios::out | ios::app);
    plik << N << ";" << Nc << ";" << Ns << ";" << iter << ";" << wynik.err << ";" << exec_time << endl;
    plik.flush();

    plik.close();
*/
    return 0;
}


