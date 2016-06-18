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
    vf_opts global_conf;
    int N = 0,
        Ns = 0,
        Nc = 0;

    // read config
    int err = read_conf( global_conf );
    // something went wrong, exit
    if ( err )
    {
        return err;
    }

    int split_strat = global_conf.split_strat;

    if ( argc == 1 )
    {
        if ( global_conf.in_file_name == "test" )
        {
            cout << "===Test for sample data===\n";
            data = prepare_sample_data();
            N = 3;
            Ns = 10;
            Nc = 4;
        }
        else
        {
            try
            {
                data = load_vf_data( global_conf.in_file_name );
                Nc = data.f.n_rows;
            }
            catch (string err)
            {
                cout << err << endl;
                return 1;
            }
        }
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

    cout << "Nc = " << Nc << endl;
    double exec_time=0;
    int iter;
    wall_clock timer;
    // wlaczenie algorytmu
    cout << "Vector fitting " << VF_REPEAT << " times" << endl;
    for ( int k = 0; k < VF_REPEAT ; k++ )
    {
        timer.tic();
        wynik = vf_high_level( data.f, data.s, global_conf ); 
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

    myfile.open( global_conf.out_file_name.c_str() );
    create_model_netlist( &wynik, Nc, data.freq, myfile, global_conf );
    myfile.close();

    if ( global_conf.ngspice_simulation )
    {
        string system_cmd = "ngspice -b " + global_conf.out_file_name;
        cout << "\n\n##### Symulacja wygenerowanego modelu w ng_spice: #####\n"
             << "##### " << system_cmd << " #####\n";
        if ( system( system_cmd.c_str() ) )
        {
            check_model_simulation_results( data.f, global_conf );
        }
    }
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


