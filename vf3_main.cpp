#include "data_model.h"
#include "network_model.h"
#include "additional_features.h"
#include "my_vectfit_non.h"

// program jako jedyny argument przyjmuje plik z konfiguracja
// brak tego pliku spowoduje bledne zakonczenie dzialania aplikacji
int main(int argc, char* argv[])
{
    MKL_Set_Num_Threads(0); // ustawienie najwiekszej mozliwej liczby watkow dla MKL

    input_data data;
    cx_mat poles;
    SER wynik;
    vf_opts global_conf;
    string conf_file_name;
    int split_strat = 1;
    int N = 0,
        Ns = 0,
        Nc = 0;
    gnuplot_data gp_data;
    time_t startup_time;

    time( &startup_time );
    cout << "## Czas rozpoczecia działania aplikacji: " << ctime( &startup_time ) << endl;

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
            if ( recognize_file( global_conf.in_file_name ) )
            {
                data = read_raw_file( global_conf.in_file_name );
            }
            else
            {
                read_touchstone( global_conf.in_file_name, data );
            }
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

                case 10:
                    cout << "To nie jest plik touchstone, albo nieprawidlowy format pliku" << endl;
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
    cout << "## Przyblizanie funkcji ukladowych algorytmem VF\n";
    timer.tic();
    wynik = vf_high_level( data.f, data.s, global_conf );
    double executionTime = timer.toc();
    cout<< "## Czas przyblizania funkcji ukladowych: "<< executionTime <<  "s" << endl;
    exec_time = exec_time + executionTime;

    cout << "\n\n";
    cout << "## Otrzymane przyblizenie (VF):\n";
    wynik.poles.print("## Bieguny=");
    wynik.res.print("## Residua=");

    if ( global_conf.calc_parallel_RC )
    {
        wynik.h.print("## h=");
        wynik.d.print("## d=");
    }

    cout << "## Blad RRMS uzyskany po aproksymacji VF i wyborze optymalnego rzedu:\n" 
         << "## RRMS_VF = " << wynik.err << " dB\n" << endl;

    // zapis otrzymanych wynikow do pliku - do wczytania w matlabie
    cout << "## Zapis otrzymanego modelu (VF) do plikow:\n";
    save_results_mats( wynik, global_conf.out_file_name );

    cout << "\n## Generowanie modelu ukladowego .cir" << endl;
    // utworzenie modelu cir i zapis do pliku
    ofstream myfile;

    myfile.open( global_conf.out_file_name.c_str() );
    create_model_netlist( &wynik, Nc, data.freq, myfile, global_conf );
    myfile.close();

    if ( global_conf.spice_simulation )
    {
        string system_cmd = global_conf.spice_program_loc + " -b " + global_conf.out_file_name;
        cout << "\n\n##### Symulacja wygenerowanego modelu w LTspice: #####\n"
             << "##### " << system_cmd << " #####\n";
        system( system_cmd.c_str() );

        if ( check_spice_log( global_conf.out_file_name ) )
        {
            check_model_simulation_results( data.f, global_conf, gp_data );
            // rysowanie wykresow
            if ( global_conf.gnuplot_generation )
            {
                gp_data.freq = data.freq;
                try
                {
                    prepare_gnuplot_script( gp_data, global_conf.out_file_name );
                }
                catch( int & err )
                {
                    cout<< "Nie mozna utowrzyc katalogu gnuplot\n";
                    return err;
                }
            }
        }
    }

    return 0;
}


