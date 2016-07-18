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
    MKL_Set_Num_Threads(1);
    input_data data;
    cx_mat poles;
    SER wynik;
    vf_opts global_conf;
    string conf_file_name;
    int split_strat = 1;
    int N = 0,
        Ns = 0,
        Nc = 0;

    // jeden argument wywolania - plik z konfiguracja
    if ( argc == 2 )
    {
        conf_file_name = string( argv[1] );

        // read config
        int err = read_conf( global_conf, conf_file_name );
        // something went wrong, exit
        if ( err )
        {
            cout << err << endl;
            return err;
        }
        split_strat = global_conf.split_strat;

        try
        {
            data = read_raw_file( global_conf.in_file_name );
            Nc = data.f.n_rows;
        }
        catch (const int& err)
        {
            switch (err)
            {
                case 1:
                    cout << "Nie mozna otworzyc pliku .raw" << endl;
                break;

                case 2:
                    cout << "Nieprawidlowa zawartosc pliku .raw\n"
                         << "Prawidlowa zawartosc: analiza AC, prady z kazdego ze zrodel" << endl;
                break;

                case 3:
                    cout << "Blad podczas operacji na plikach, przy wczytywaniu wejsciowego .raw" << endl;
                break;
            }
            return err;
        }
    }
    else
    {
        cout << "Nieprawidlowa liczba argumentow wejsciowych\nJedynym argumentem jest sciezka do pliku z konfiguracja.\n";
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

    if ( global_conf.spice_simulation )
    {
        string system_cmd = global_conf.spice_program_loc + " -b " + global_conf.out_file_name;
        cout << "\n\n##### Symulacja wygenerowanego modelu w ng_spice: #####\n"
             << "##### " << system_cmd << " #####\n";
        system( system_cmd.c_str() );

        if ( check_spice_log( global_conf.out_file_name ) )
        {
            check_model_simulation_results( data.f, global_conf );
        }
    }

    return 0;
}


