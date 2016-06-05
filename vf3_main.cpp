#include "data_model.h"
#include "network_model.h"
#include "additional_features.h"
#include "my_vectfit_non.h"
#define VF_REPEAT 1

vf_opts global_conf;
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
    int split_strat = 1;

    read_conf();

    if ( argc == 1 )
    {
        cout << "===Test for sample data===\n";
        data = prepare_sample_data();
        N = 3;
        Ns = 10;
        Nc = 4;
    }
    else if ( argc == 2)
    {
        split_strat = atoi( argv[1] );
    }
    else if ( argc == 5 )
    {
        cout << "=== Dane z pliku ===\n";
        N = atoi( argv[1] );
        Nc = atoi( argv[2] );
        Ns = atoi( argv[3] );
        split_strat = atoi( argv[4] );

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
        wynik = vf_high_level( data.f, data.s,
                               global_conf.split_strat, global_conf.start_row,
                               global_conf.end_row, global_conf.max_iters ); 
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
    ofstream myfile;
    string file_name = "test_mgr.cir";
    if ( split_strat == NON_SPLITING )
    {
        file_name = "test_mgr_non.cir";
    }
    else if ( split_strat == ALL_SPLITING )
    {
        file_name = "test_mgr_all.cir";
    }
    else if ( split_strat == COLUMN_SPLITING )
    {
        file_name = "test_mgr_column.cir";
    }

    myfile.open( file_name.c_str() );
    create_model_netlist( &wynik, Nc, myfile);
    myfile.close();

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


