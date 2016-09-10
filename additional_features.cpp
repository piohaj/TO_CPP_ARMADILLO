#include "additional_features.h"

SER vf_high_level( cx_mat& f, const cx_vec& s, vf_opts conf )
{
    SER wynik;
    int row_iterations_num = conf.max_row - conf.min_row + 1;
    SER *wynik_iter = new SER[row_iterations_num];
    cx_mat poles;
    int Nc = f.n_rows;
    int Ns = f.n_cols;
    mat result_idx;
    int split_strat = conf.split_strat;
    
    if ( conf.max_row < conf.min_row )
    {
        cout << "++ Podano nieprawidlowy przedzial rzedow przyblizenia" << endl;
        delete[] wynik_iter;
        return wynik;
    }

    // sprawdzenie i ewentualne wymuszenie pasywnosci
    if ( conf.pasivity_check )
    {
        if ( check_and_make_passive( f ) != 0 )
        { 
            cout << "++ Problem podczas wymuszania pasywnosci" << endl;
            delete[] wynik_iter;
            return wynik;
        }
    }

    cout << "\n";
    cout << "++ Aplikacja szuka optymalnego przyblizenia wsrod rzedow z przedzialu:\n"
         << "++ <" << conf.min_row << ";" << conf.max_row << ">\n" << endl;
    
    if ( split_strat == NON_SPLITING )
    {
        int high_iter = 0;
        for ( int row = conf.min_row; row <= conf.max_row; row++ )
        {
            //biguny poczatkowe
            poles = prepare_input_poles(s, split_strat, Nc, row, Ns);

            int iter = 0;
            // wywolanie algorytmu
            for ( iter = 1; iter <= conf.max_iters; iter++ )
            {
                wynik_iter[high_iter] = my_vf_non_splitting(f, s, poles, conf); 
    	        poles = wynik_iter[high_iter].poles;
    		
    	        if ( wynik_iter[high_iter].err < conf.tol )
                {
    	            break;
    	        }
            }
            cout << "++ Rzad przyblizenia: " << row << ", RRMS: "
                 << wynik_iter[high_iter].err << " dB" << endl;
            high_iter++;
        }
        result_idx = choose_best_aprox( wynik_iter, row_iterations_num, Nc, split_strat, conf.rms_diff );
        cout << "\n++ Optymalne rozwiazanie (zgodnie z konfiguracja):\n"
             << "++ Rzad: " << conf.min_row + int(result_idx(0)) << "\n" << endl;
        wynik = wynik_iter[int(result_idx(0))];
    }
    else if ( split_strat == ALL_SPLITING )
    {
        int high_iter = 0;
        for ( int row = conf.min_row; row <= conf.max_row; row++ )
        {
            //bieguny poczatkowe
            poles = prepare_input_poles(s, split_strat, Nc, row, Ns);

            int iter = 0;
            // wywolanie algorytmu
            for ( iter = 1; iter <= conf.max_iters; iter++ )
            {
                wynik_iter[high_iter] = my_vf_all_splitting(&f, &s, &poles, conf); 
    	        poles = wynik_iter[high_iter].poles;
    		
    	        if ( wynik_iter[high_iter].err < conf.tol )
                {
    	            break;
    	        }
            }
            cout << "++ Rzad przyblizenia: " << row << ", RRMS: "
                 << wynik_iter[high_iter].err << " dB" << endl;
            wynik_iter[high_iter].err_table.print("++ Bledy RRMS dla kolejnych elementow Y (dB)=");

            high_iter++;
        }
       
       result_idx = choose_best_aprox( wynik_iter, row_iterations_num, Nc, split_strat, conf.rms_diff );
       cout << "\n++ Optymalne rozwiazanie (zgodnie z konfiguracja):\n";
       (result_idx + conf.min_row).print("++ Rzedy (dla kolejnych elementow Y (kolumnowo):");

       wynik = cumulate_model( split_strat, result_idx, wynik_iter, Nc, conf.max_row);
    }
    else if ( split_strat == COLUMN_SPLITING )
    {
        int high_iter = 0;
        for ( int row = conf.min_row; row <= conf.max_row; row++ )
        {
            //bieguny poczatkowe
            poles = prepare_input_poles(s, split_strat, Nc, row, Ns);

            int iter = 0;
            // wywolanie algorytmu
            for ( iter = 1; iter <= conf.max_iters; iter++ )
            {
                wynik_iter[high_iter] = my_vf_column_splitting(&f, &s, &poles, conf); 
    	        poles = wynik_iter[high_iter].poles;
    		
    	        if ( wynik_iter[high_iter].err < conf.tol )
                {
    	            break;
    	        }
            }
            cout << "++ Rzad przyblizenia: " << row << ", RRMS: "
                 << wynik_iter[high_iter].err << " dB" << endl;
            wynik_iter[high_iter].err_table.print("++ Bledy RRMS dla kolejnych kolumn Y=");

            high_iter++;
        }
        result_idx = choose_best_aprox( wynik_iter, row_iterations_num, Nc, split_strat, conf.rms_diff );

        cout << "\n++ Optymalne rozwiazanie (zgodnie z konfiguracja):\n";
        (result_idx + conf.min_row).print("++ Rzedy (dla kolejnych kolumn Y (pionowo):");

        wynik = cumulate_model( split_strat, result_idx, wynik_iter, Nc, conf.max_row);
    }
    else
    {
        cout << "Nieznana strategia podzialu Y" << endl;
    }


    delete[] wynik_iter;
    return wynik;
} 


cx_mat prepare_input_poles( const cx_vec& s, int split_strat, int Nc, int N, int Ns )
{
    cx_mat poles;
    if ( split_strat == NON_SPLITING )
    {
        poles = zeros<cx_mat>(1, N);
        mat bet = linspace<mat>(imag(s(0))+10, imag(s(Ns-1)), N/2);
    
        int m = 0;
        int n = 0;
        for ( n = 0; n < N-1; n=n+2 )
        {
            double alf = -bet(m)*1e-2;
            poles(0, n) = cx_double(alf, bet(m));
            poles(0, n+1) = conj(poles(n));
            m++;
        }

        // zabezpieczenie przed zerowymi biegunami, ktore sa niebezpieczne gdy dane wejsciowe maja probke freq=0
        for ( int nn = 0; nn < N; nn++ )
        {
            if ( poles(0, nn) == cx_double(0, 0) )
            {
                poles(0, nn) = cx_double(10, 0);
            }
        }
    }
    else if ( split_strat == ALL_SPLITING )
    {
        poles = zeros<cx_mat>(Nc, N);
        mat bet = linspace<mat>(imag(s(0))+10, imag(s(Ns-1)), N/2);
        
        for ( int mm = 0; mm < Nc; mm++ )
        {
            int m = 0;
            int n = 0;
            for ( n = 0; n < N-1; n=n+2 )
            {
                double alf = -bet(m)*1e-2;
                poles(mm, n) = cx_double(alf, bet(m));
                poles(mm, n+1) = conj(poles(mm, n));
                m++;
            }
        }

        // zabezpieczenie przed zerowymi biegunami, ktore sa niebezpieczne gdy dane wejsciowe maja probke freq=0
        for ( int mm = 0; mm < Nc; mm++ )
        {
            for ( int nn = 0; nn < N; nn++ )
            {
                if ( poles(mm, nn) == cx_double(0, 0) )
                {
                    poles(mm, nn) = cx_double(10, 0);
                }
            }
        }
    }
    else if ( split_strat == COLUMN_SPLITING )
    {
        int column_num = sqrt(Nc);
        poles = zeros<cx_mat>(column_num, N);
        mat bet = linspace<mat>(imag(s(0))+10, imag(s(Ns-1)), N/2);
        
        for ( int mm = 0; mm < column_num; mm++ )
        {
            int m = 0;
            int n = 0;
            for ( n = 0; n < N-1; n=n+2 )
            {
                double alf = -bet(m)*1e-2;
                poles(mm, n) = cx_double(alf, bet(m));
                poles(mm, n+1) = conj(poles(n));
                m++;
            }
        }

        // zabezpieczenie przed zerowymi biegunami, ktore sa niebezpieczne gdy dane wejsciowe maja probke freq=0
        for ( int mm = 0; mm < column_num; mm++ )
        {
            for ( int nn = 0; nn < N; nn++ )
            {
                if ( poles(mm, nn) == cx_double(0, 0) )
                {
                    poles(mm, nn) = cx_double(10, 0);
                }
            }
        }
    }
 
    return poles;
}

mat choose_best_aprox( SER *input, int size, int Nc, int split_strat, double rms_diff )
{
    mat result;
    if ( split_strat == NON_SPLITING )
    {
        result = zeros<mat>(1,1);
        double err = 0;
        for ( int i = 0; i < size; i++ )
        {
            double err_temp = input[i].err;
            if ( ( err - err_temp ) >= rms_diff )
            {
                err = input[i].err;
                result(0) = i;
            }
        }
    }
    else if ( split_strat == ALL_SPLITING )
    {
        result = zeros<mat>(Nc, 1);
        for ( int i = 0; i < Nc; i++ )
        {
            double err = 0;
            for ( int j = 0; j < size; j++ )
            {
                double err_temp = input[j].err_table[i];

                if ( ( err - err_temp ) >= rms_diff )
                {
                    err = input[j].err_table[i];
                    result(i) = j;
                }
            }
        }
    }
    else if ( split_strat == COLUMN_SPLITING )
    {
        result = zeros<mat>( sqrt(Nc), 1 );
        for ( int i = 0; i < sqrt(Nc); i++ )
        {
            double err = 0;
            for ( int j = 0; j < size; j++ )
            {
                double err_temp = input[j].err_table[i];

                if ( (err - err_temp) >= rms_diff )
                {
                    err = input[j].err_table[i];
                    result(i) = j;
                }
            }
        }
    
    }
    return result;
}  


SER cumulate_model( int split_strat, mat& indexes, SER *iter_models, int Nc, int max_row )
{
    SER wynik;

    if ( split_strat == ALL_SPLITING )
    {
        wynik.poles = zeros<cx_mat>(Nc, max_row);
        wynik.res = zeros<cx_mat>(Nc, max_row);
        wynik.d = zeros<mat>(Nc, 1);
        wynik.h = zeros<mat>(Nc, 1);
        double err_temp = 0.0;
    
        for ( int n = 0; n < Nc ; n++ )
        {
            cx_mat poles_temp = iter_models[int(indexes(n))].poles.row(n);
            cx_mat res_temp = iter_models[int(indexes(n))].res.row(n);
            mat d_temp = iter_models[int(indexes(n))].d.row(n);
            mat h_temp = iter_models[int(indexes(n))].h.row(n);

            wynik.poles( n, span(0,poles_temp.n_cols-1) ) = poles_temp;
            wynik.res( n, span(0, res_temp.n_cols-1) ) = res_temp;
            wynik.d.row(n) = d_temp;
            wynik.h.row(n) = h_temp;
            err_temp += iter_models[int(indexes(n))].err_table(n); 
        }
        wynik.err = err_temp/Nc;
    }
    else if ( split_strat == COLUMN_SPLITING )
    {
        wynik.poles = zeros<cx_mat>(sqrt(Nc), max_row);
        wynik.res = zeros<cx_mat>(Nc, max_row);
        wynik.d = zeros<mat>(Nc, 1);
        wynik.h = zeros<mat>(Nc, 1);
        double err_temp = 0.0;
    
        for ( int n = 0; n < Nc ; n++ )
        {
            cx_mat poles_temp = iter_models[int(indexes(n/sqrt(Nc)))].poles.row(n/sqrt(Nc));
            cx_mat res_temp = iter_models[int(indexes(n/sqrt(Nc)))].res.row(n);
            mat d_temp = iter_models[int(indexes(n/sqrt(Nc)))].d.row(n);
            mat h_temp = iter_models[int(indexes(n/sqrt(Nc)))].h.row(n);

            wynik.poles( n/sqrt(Nc), span(0, poles_temp.n_cols-1) ) = poles_temp;
            wynik.res( n, span(0, res_temp.n_cols-1) ) = res_temp;
            wynik.d.row(n) = d_temp;
            wynik.h.row(n) = h_temp;
            err_temp += iter_models[int(indexes(n/sqrt(Nc)))].err_table(n/sqrt(Nc)); 
        }
        wynik.err = err_temp / Nc;
    }

    return wynik;
}



bool ispassive_s( cx_cube& s_params )
{
    bool result = true;
    int Ns = s_params.n_slices;

    int ii = 0;
    while ( result && ii < Ns )
    {
        result = ispassive_slice( s_params.slice(ii) );
        ii++;
    }
  
    return result; 
} 

bool ispassive_slice(cx_mat& s_slice)
{
    return !( norm(s_slice, 2) > (1+sqrt(1e-16)) );
}

int makepassive( cx_cube& s_params )
{
    double threshold = 1 - sqrt(1e-16);
    int Ns = s_params.n_slices;
    bool result = true;

    for ( int j = 0; j < Ns; j++ )
    {
        result = ispassive_slice( s_params.slice(j) );    
        int idx = 0;
        
        while ( result == false && idx < 10 )
        {
            cout << "Brak pasywnosci w plastrze " << j << endl;
            cx_mat U, V;
            vec zigm;
            svd( U, zigm, V, s_params.slice(j) );
            int nzigma = zigm.n_elem;
            mat zigma = diagmat( zigm );
            cx_mat upsilon = eye<cx_mat>(nzigma, nzigma);
            cx_mat psi = eye<cx_mat>(nzigma, nzigma);
           
            for ( int ii = 0; ii < nzigma; ii++ )
            {
                if ( zigma(ii, ii) <= threshold )
                {
                    upsilon(ii, ii) = 0;
                    psi(ii, ii) = 0;
                }
                else
                {
                    upsilon(ii, ii) = 1;
                    psi(ii, ii) = threshold;
                }
             }
             cx_mat zigma_viol = zigma*upsilon - psi;
             cx_mat s_viol = U*zigma_viol*V.t();

             // compensate data
             s_params.slice(j) = s_params.slice(j) - s_viol;
             idx = idx + 1;
             result = ispassive_slice( s_params.slice(j) );
         }

         if ( result == false )
         {
             cout << "Nie udalo sie wymusic pasywnosci" << endl;
             return 1;
         }
    }

    return 0;
}

int check_and_make_passive( cx_mat& y ) 
{
    int result = 0;
    cx_cube s_params = y2s( y );

    if ( !ispassive_s( s_params ) )
    {
        cout << "Dane nie sa pasywne, wymuszanie pasywnosci..." << endl;
        result = makepassive( s_params );
        y = s2y( s_params );
    }
    else
    {
        cout << "Dane dla modelu pasywnego, wymuszanie pasywnosci nie jest potrzebne" << endl;
    }

    return result;
}


cx_mat s2y(cx_cube s_params, int z0)
{
    int Ns = s_params.n_slices;
    int Nc_ports = s_params.n_cols;
    
    cx_cube yy = zeros<cx_cube>(Nc_ports, Nc_ports, Ns);

    cx_mat I = eye<cx_mat>(Nc_ports, Nc_ports);

    for ( int k = 0; k < Ns; k++ )
    {
        yy.slice(k) = ( I - s_params.slice(k) ) * inv( z0 * (I + s_params.slice(k) ) );
    }

    return cube2mat( yy );
}

cx_mat cube2mat( cx_cube& cube )
{
    int Ns = cube.n_slices;
    int Nc_ports = cube.n_cols;
    int Nc = pow( Nc_ports, 2 );
    cx_mat y = zeros<cx_mat>(Nc, Ns);

    for ( int i = 0; i < Ns; i++ )
    {
        // kolumna z danymi dla danej probki czestotliwosciowej
        cx_mat temp = zeros<cx_mat>(Nc, 1);

        for ( int k = 0; k < Nc_ports; k++ )
        {
            temp.rows(Nc_ports * k, Nc_ports*k + Nc_ports - 1) = cube.slice(i).col(k);
        }

        y.col(i) = temp;
    }

    return y;
}

cx_cube y2s(cx_mat y, int z0)
{
    int Ns = y.n_cols;
    int Nc = y.n_rows;
    int Nc_ports = sqrt(Nc);
    cx_cube s_params = zeros<cx_cube>(Nc_ports, Nc_ports, Ns);
    cx_cube yy = make_cube( y );

    cx_mat I = eye<cx_mat>(Nc_ports, Nc_ports);

    for ( int k = 0; k < Ns; k++ )
    {
        s_params.slice(k) = solve( I + z0 * yy.slice(k), I - z0 * yy.slice(k) );
    }

    return s_params;
}

cx_cube make_cube( cx_mat& y )
{
    int Ns = y.n_cols;
    int Nc = y.n_rows;

    cx_cube yy = zeros<cx_cube>(sqrt(Nc), sqrt(Nc), Ns); 

    for ( int i = 0; i < Ns; i++ )
    {
        cx_mat temp = zeros<cx_mat>(sqrt(Nc), sqrt(Nc));
        for ( int j = 0 ; j < sqrt(Nc) ; j++ )
        {
           temp.col(j) = y( span(sqrt(Nc)*j, sqrt(Nc) + sqrt(Nc)*j - 1), i );   
        }
        yy.slice(i) = temp;
    }

    return yy;
}


int check_model_simulation_results( const cx_mat& f, const vf_opts& conf, gnuplot_data & gp_data )
{
    string simulation_data_file = conf.out_file_name;
    simulation_data_file.replace(simulation_data_file.end()-3,simulation_data_file.end(), "raw");
    int Ns = f.n_cols;
    input_data spice_simulation_data;

    try
    {
        spice_simulation_data = read_raw_file( simulation_data_file );
    }
    catch (const int& err)
    {
        cout << err << endl;
        return err;
    }

    gp_data.input_data = f;
    gp_data.simulation_data = spice_simulation_data.f;

    // obliczanie rms miedzy danymi wejsciowymi a modelem
    mat diff_real = real(f - spice_simulation_data.f);
    mat diff_imag = imag(f - spice_simulation_data.f);

    //double rms_err = sqrt( ( accu( pow(diff_real, 2) + pow(diff_imag, 2) ) ) / Ns );
    double rms_err = sqrt( accu( pow( abs( f - spice_simulation_data.f ), 2 ) ) /
                     accu ( pow ( abs(f), 2 ) ) );
    double rms_err_db = 20 * log10( rms_err );
    cout << "## Blad RRMS (pomiedzy danymi wejsciowymi a danymi uzyskanymi z symulacji LT Spice) = " << rms_err << endl;
    cout << "## Blad RRMS (pomiedzy danymi wejsciowymi a danymi uzyskanymi z symulacji LT Spice) = " << rms_err_db << " db" << endl;

    return 0;
}

void force_stable_poles( cx_mat& poles )
{
    int N = poles.n_elem;

    for ( int i = 0; i < N; i++ )
    {
        double pole_real = real(poles(i));

        if ( pole_real > 0 )
        {
            poles(i) = poles(i) - cx_double(2*pole_real, 0);
        }
    }
}


// wczytywanie konfiguracji z pliku ciala funkcji
int read_conf( vf_opts& global_conf, string file_name )
{
    map <string,string> conf_map = read_conf_file(file_name);

    try
    {
        global_conf.out_file_name = read_param("out_file_name", conf_map); 
        global_conf.in_file_name = read_param("in_file_name", conf_map); 
        global_conf.tol = atof( read_param("min_rms", conf_map).c_str() ); 
        global_conf.rms_diff = atof( read_param("rms_diff", conf_map).c_str() ); 
        global_conf.min_row = atoi( read_param("min_row", conf_map).c_str() ); 
        global_conf.max_row = atoi( read_param("max_row", conf_map).c_str() ); 
        global_conf.max_iters = atoi( read_param("max_iters", conf_map).c_str() ); 
        global_conf.R_max = atof( read_param("R_max", conf_map).c_str() ); 
        global_conf.C_min = atof( read_param("C_min", conf_map).c_str() ); 
        global_conf.split_strat = atoi( read_param("spliting_strategy", conf_map).c_str() ); 
        global_conf.pasivity_check = atoi( read_param("pasivity_check", conf_map).c_str() ); 
        global_conf.spice_simulation = atoi( read_param("spice_simulation", conf_map).c_str() );
        global_conf.spice_program_loc = read_param("spice_program_loc", conf_map);
        global_conf.calc_parallel_RC = atoi( read_param("calc_parallel_RC", conf_map).c_str() );
        global_conf.gnuplot_generation = atoi( read_param("gnuplot_generation", conf_map).c_str() );
    }
    catch( const string& err )
    {
        cout << err;
        return 1;
    }

    return 0;
}

vector<string> my_split(string str, const char delim)
{
    vector<string> v;
    string tmp;
    string::const_iterator i = str.begin();

    for(i ; i <= str.end(); ++i) 
    {
        if ( *i != delim && !isspace(*i) && i != str.end())
        {
            tmp += *i;
        } 
        else if ( *i == delim || i == str.end() )
        {
            v.push_back(tmp);
            tmp = "";
        }
    }
    return v;
}

map <string, string>read_conf_file(string file_name)
{
    fstream plik(file_name.c_str(), std::ios::in);
    string single_line;
    map <string, string>conf_map;
    vector<string> vec;
    
    if ( plik.good() == false )
    {
        cout << "Nie udalo sie otworzyc pliku\n";
    }

    while ( getline( plik, single_line ) )
    {
        if ( strncmp(single_line.c_str(), "//", 2) != 0 && single_line.find('=') != string::npos )
        {
            vec = my_split(single_line, '=');
            conf_map[vec[0]] = vec[1];
            vec.clear();
        }
    }
    
    plik.close();
    return conf_map;
}


string read_param( string param, map <string,string> &conf_map)
{
    string temp = conf_map[param];
    if ( temp.empty() )
    {
        throw "Brak parametru " + param + " w pliku konfiguracyjnym!!!\n";
    }
    return temp;
}


raw_params read_input_data_params( string file_name )
{
    fstream plik( file_name.c_str(), std::ios::in );
    string single_line;
    vector<string> vec;
    raw_params params;
    params.last_line = 0;
    params.Nc_ports = 0;
    params.Ns = 0;
    bool correct_raw_file = false;
    bool check_variables = false;
    int check_var_num = 0;

    if ( plik.good() == false )
    {
        throw 1;
    }

    while ( getline(plik, single_line ) )
    {
        params.last_line++;
        if ( single_line.find("Plotname:") != string::npos )
        {
            if ( single_line.compare("Plotname: AC Analysis") != 0 )
            {
                throw 2; 
            }
            correct_raw_file = true;
        }
        else if ( single_line.find("No. Variables:") != string::npos )
        {
            vec = my_split (single_line, ':');
            params.Nc_ports = atoi( vec[1].c_str() );
        }
        else if ( single_line.find("No. Points:") != string::npos )
        {
            vec = my_split (single_line, ':');
            params.Ns = atoi( vec[1].c_str() );
        }
        else if ( single_line.find("Variables") != string::npos )
        {
            check_variables = true;
        }
        else if ( single_line.find("Binary:") != string::npos )
        {
            break;
        }
        else if ( check_variables && single_line.find("device_current") != string::npos )
        {
            check_var_num++;
        }
    }
    
    // input file is not spice raw file
    if ( correct_raw_file == false || check_var_num != params.Nc_ports-1 )
    {
        throw 2;
    }
 
    plik.close();

    return params;
}

string create_temp_bin_file( int last_no_bin_data_line, string file_name )
{
    string temp_file_name = file_name + "_temp";
    string single_line;
    int line_number = 0;
    
    fstream plik( file_name.c_str(), std::ios::in );
    fstream temp_file( temp_file_name.c_str(), std::ios::out );

    if ( plik.good() == false || temp_file.good() == false )
    {
        throw 1;
    }

    while ( getline(plik, single_line) )
    {
        line_number++;
        if ( line_number > last_no_bin_data_line )
        {
            single_line = single_line + "\n";
            temp_file.write ( &single_line[0], single_line.length() );
        }
    }

    plik.close();
    temp_file.close();
    
    return temp_file_name;
}

input_data read_raw_file( string file_name )
{
    string temp_file;
    cx_mat temp_data;
    raw_params input_raw_params;

    try
    {
        // reading params from raw file
        input_raw_params = read_input_data_params( file_name );
        // create temp file to be loaded by armadillo
        temp_file = create_temp_bin_file( input_raw_params.last_line, file_name ); 
        temp_data.load(temp_file);
        if ( std::remove( temp_file.c_str() ) != 0 )
        {
            throw 3;
        }
    }
    catch( const int& err )
    {
        throw err;
    }

    return parse_raw_data( temp_data, input_raw_params );
}

input_data parse_raw_data( const cx_mat& temp_data, raw_params params )
{
    input_data data;
    int temp_data_size = temp_data.n_rows;
    int Ns = params.Ns / ( params.Nc_ports - 1 );
    cx_mat f_temp;
    data.freq = zeros<vec>(Ns);
    data.f = zeros<cx_mat>(pow(params.Nc_ports - 1, 2), Ns);

    int j = 0;
    for ( int i = 0; i < temp_data_size/(params.Nc_ports - 1); i = i+params.Nc_ports )
    {
        cx_double freq_temp = temp_data(i);
        data.freq(j) = real(freq_temp);
        j++;
    }
    data.s = 2*PI*cx_double(0,1) * data.freq;

    for ( int i = 0; i < temp_data_size; i = i+params.Nc_ports )
    {
         cx_mat temp = temp_data.rows(i+1, i+params.Nc_ports-1);
         f_temp = join_horiz(f_temp, temp);
    }
    f_temp = f_temp.st();

    j = 0;
    for ( int x = params.Nc_ports-2; x >= 0; x-- )
    {
        for ( int z = 0; z < params.Ns; z=z+Ns )
        {
            cx_mat temp = f_temp(span(z, z+Ns-1), x).st();
            data.f.row(j) = temp;
            j++;
        }
    } 

    return data;
}


bool check_spice_log( string file_name )
{
     file_name.replace(file_name.end()-3, file_name.end(), "log");
     string single_line;
     int line_count = 0;

     fstream file( file_name.c_str(), ios::in );
     if ( file.good() == false )
     {
         return false;
     }

     while ( getline(file, single_line) )
     {
         cout << single_line << endl;
         if ( single_line.find("Error") != string::npos )
         {
             return false;
         }
         line_count++;
     }
     file.close();

     if ( line_count <= 1 )
     {
         return false;
     }
     
     return true;
}


int get_ports_num_touchstone( string file_name )
{
    int Nc_ports = 0;

    string::const_iterator i = file_name.begin();
    string tmp;
    bool is_dot = false;

    for(i ; i <= file_name.end(); ++i)
    {
        if ( isdigit(*i) && is_dot )
        {
            tmp += *i;
        }

        if ( *i == '.' ) is_dot = true;
    }

    Nc_ports = atoi(tmp.c_str());

    if ( Nc_ports == 0 )
    {
        throw 1; // nieprawidlowa nazwa pliku
    }

    return Nc_ports;
}

int count_spaces_in_header( string header )
{
    int spaces_count = 0;
    string::const_iterator i = header.begin();
    
    for(i ; i <= header.end(); ++i)
    {
        if ( isspace(*i) )
        {   
            spaces_count++;
        } 
    }

    return spaces_count;
}


touchstone_conf check_header_touchstone( string file_name )
{
     touchstone_conf conf;
     conf.is_touchstone = false;
     string single_line;
     conf.Ns = 0;
     vector<string> conf_line;
     fstream file( file_name.c_str(), ios::in );
     if ( file.good() == false )
     {
         return conf;
     }

     while ( getline(file, single_line) )
     {
         if ( single_line.find("#") != string::npos )
         {
             if ( count_spaces_in_header( single_line ) < 7 )
             {
                 conf_line = my_split(single_line, ' ');
                 conf.is_touchstone = true;
             }
             else
             {
                 return conf;
             }
         }
         else if ( single_line.find("!") == string::npos )
         {
             conf.Ns++;
         }
     }

     conf.freq_unit = str_toupper( conf_line[1] );
     conf.data_type = str_toupper( conf_line[2] );
     conf.data_type2 = str_toupper( conf_line[3] );
     conf.R0 = atof( conf_line[5].c_str() );

     return conf;
}

string str_toupper( string & str )
{
    string out;
    string::const_iterator i = str.begin();

    for(i ; i <= str.end()-1; ++i)
    {
        out += toupper(*i);
    }
   
    return out;
}


void read_touchstone( string file_name, input_data & data )
{
    map<string, double> touchstone_freq_unit;
    touchstone_freq_unit["GHZ"] = 10e9;
    touchstone_freq_unit.insert(std::pair<string, double>("MHZ",10e6));
    touchstone_freq_unit.insert(std::pair<string, double>("KHZ",10e3));
    touchstone_freq_unit.insert(std::pair<string, double>("HZ",1));

    touchstone_conf conf = check_header_touchstone( file_name );

    if ( conf.is_touchstone == false )
    {
        throw 10;
    }

    int Nc_ports = get_ports_num_touchstone( file_name );
    int Nc = pow( Nc_ports, 2 );
    int Ns = conf.Ns;

    std::ifstream ifile( file_name.c_str(), std::ios::in);
    int line_it = 0;
    std::string single_line;
    double n = 0;
    arma::mat slice;

    // ===== dla dwuwrotnikow
    if ( Nc_ports == 2 )
    {
        for ( int i = 0; getline(ifile,single_line); i++ )
        {
           std::stringstream stream(single_line);
           arma::mat temp = zeros<mat>(1, Nc*2+1);
    
           if ( single_line[0] != '!' && single_line[0] != '#' )
           {
               for ( line_it = 0; stream >> n; line_it++ )
               {
                   temp(line_it) = n;
               }
               slice = join_vert( slice, temp );
           }
        }
        ifile.close();
    
        data.freq = slice.col(0) * touchstone_freq_unit[conf.freq_unit];
    
        cx_mat ri_touchstone;
    
        for ( int i = 1; i <= Nc*2; i=i+2 )
        {
            mat data_col1 = slice.col(i);
            mat data_col2 = slice.col(i+1);
            cx_mat temp;
     
            if ( conf.data_type2 == "DB" ) { cout<<"Dane w postaci db\n"; db2magnitude( data_col1 ); }
            if ( conf.data_type2 == "DB" || conf.data_type2 == "MA")
            {
                cout<< "Dane w postaci wykladniczej\n";
                temp = angle2canonic( data_col1, data_col2 );
            }
            else
            {
                temp = cx_mat(data_col1, data_col2);
            }
    
            ri_touchstone = join_vert( ri_touchstone, temp.st() );
        }
     
        if ( conf.data_type == "S" )
        {
            cx_cube ri_cube = make_cube( ri_touchstone );
            data.f = s2y( ri_cube, conf.R0 );
        }
        else
        {
            data.f = ri_touchstone;
        }
    }
    else // dla wielowrotnikow
    {
        data.freq = zeros<vec>(Ns/Nc_ports);
        int slice_iter = 0;
        int freq_iter = 0;
        cx_cube f_data_cube;

        for ( int i = 0; getline(ifile,single_line); i++ )
        {
            std::stringstream stream(single_line);
            arma::mat temp = zeros<mat>(1, Nc_ports*2+1);
    
            if ( single_line[0] != '!' && single_line[0] != '#' )
                {
                    for ( line_it = 0; stream >> n; line_it++ )
                    {
                        temp(line_it) = n;
                    }
                
                if ( slice_iter == 0 ) // jesli linia z czestotliwoscia
                {
                    data.freq(freq_iter++) = temp(0) * touchstone_freq_unit[conf.freq_unit];
                    slice = join_vert( slice, temp.cols(1, temp.n_cols - 1) );
                }
                else
                {
                    slice = join_vert(slice, temp.cols(0, temp.n_cols - 2) );
                }
                slice_iter++;
    
                if ( slice_iter == Nc_ports )
                {
                    slice_iter = 0;
                    cx_mat ri_touchstone;
                    for ( int i = 0; i < Nc_ports*2; i=i+2 )
                    {
                        mat data_col1 = slice.col(i);
                        mat data_col2 = slice.col(i+1);
                        cx_mat temp;
            
                        if ( conf.data_type2 == "DB" ) { cout<<"Dane w postaci db\n"; db2magnitude( data_col1 ); }
                        if ( conf.data_type2 == "DB" || conf.data_type2 == "MA")
                        {   
                            cout<< "Dane w postaci wykladniczej\n";
                            temp = angle2canonic( data_col1, data_col2 );
                        }
                        else
                        {   
                            temp = cx_mat(data_col1, data_col2);
                        }
                        
                        ri_touchstone = join_horiz( ri_touchstone, temp );
                   }
                   f_data_cube = join_slices(f_data_cube, ri_touchstone);
                   slice.reset();
                }                
            }
        }
        ifile.close();

        if ( conf.data_type == "S" )
        {
            data.f = s2y( f_data_cube, conf.R0 );
        }
    }
    data.s = data.freq * PI * 2 * cx_double(0,1);
}

cx_mat angle2canonic( const mat& mag, const mat& angle )
{
    mat real =  mag % cos( angle * PI/180 );
    mat imag =  mag % sin( angle * PI/180 );

    return cx_mat(real, imag);
}


void db2magnitude( mat& db )
{
    for ( int i = 0; i < db.n_elem; i++ )
    {
        double x = db(i);
        db(i) = pow(10, x/20);
    }
}

int recognize_file ( string file_name )
{
    vector<string> vec;
    vec = my_split(file_name, '.');
    string file_type = vec[1];

    if ( file_type == "raw" )
    {
        return RAW_FILE;
    }
    else
    {
        return TOUCHSTONE_FILE;
    }
}

void save_results_mats( SER & results, string file_name )
{
    vector<string> vec = my_split(file_name, '.');
    string file_prefix = vec[0];
    string poles_file_name = file_prefix + "_poles.mat";
    string res_file_name = file_prefix + "_res.mat";
    string d_file_name = file_prefix + "_d.mat";
    string h_file_name = file_prefix + "_h.mat";

    int Nc_poles = results.poles.n_rows;
    int Nc = results.res.n_rows;
    int N = results.res.n_cols;

    mat poles_real = real(results.poles);
    mat poles_imag = imag(results.poles);
    mat poles_mat = join_horiz( poles_real, poles_imag);

    mat res_real = real(results.res);
    mat res_imag = imag(results.res);
    mat res_mat = join_horiz( res_real, res_imag );

    cout << poles_file_name << endl;
    poles_mat.save(poles_file_name, raw_ascii);

    cout << res_file_name << endl;
    res_mat.save(res_file_name, raw_ascii);

    cout << d_file_name << endl;
    results.d.save(d_file_name, raw_ascii);

    cout << h_file_name << endl;
    results.h.save(h_file_name, raw_ascii);
}

void prepare_gnuplot_script( gnuplot_data & data, string name )
{
    cout<<"\n========== Wykresy porownawcze ========= " << endl;
    cout << "Przygotowanie danych do wykresow oraz skryptu gnuplot" << endl;
    vector<string> vec = my_split(name, '.');
    string dir_name = vec[0] + "_gnuplot";

    int dir_err = mkdir(dir_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
    if ( dir_err && errno != EEXIST )
    {
        throw dir_err;
    }
    cout << "Folder roboczy wykresow (gnuplot): " << dir_name << endl;

    string gp_script_name = dir_name + "/" + vec[0] + "_script.gp";
    cout << "Nazwa skryptu do uruchomienia w gnuplot " << gp_script_name << endl;

    ofstream gpfile;
    gpfile.open( gp_script_name.c_str() );

    int Nc = data.simulation_data.n_rows;
    int Nc_ports = sqrt(Nc);
    int Ns = data.simulation_data.n_cols;
    mat file_mat = zeros<mat>(Ns, 4);
    file_mat.col(0) = data.freq;
    string file_name = dir_name + "/";

    gpfile << "set xlabel \"freq\"" << endl;
    gpfile << "set linestyle 1 lt 2 lw 2 pt 12 ps 2" << endl;
    gpfile << "set linestyle 2 lt 1 lw 2 pt 2 ps 1\n" << endl;
    int gp_it = 0;
    int data_it = 0;
    for ( int i = 0; i < Nc_ports ; i++ )
    {
       for ( int j = 0; j < Nc_ports; j++ )
       {
           ostringstream ss_abs, ss_abs_file, ss_index;
           cx_mat simulation_temp = data.simulation_data.row(data_it).st();
           cx_mat input_temp = data.input_data.row(data_it).st();
           data_it++;
    
           // przygowanie odpowiedniego indeksu Y, tj. 11, 21, itd.
           ss_index << j+1 << i+1;

           ss_abs_file << "Y" << ss_index.str() << "_abs";
           ss_abs << file_name << ss_abs_file.str();
           file_mat.col(1) = abs(input_temp);
           file_mat.col(2) = abs(simulation_temp);
           file_mat.col(3) = abs(input_temp - simulation_temp);
            
           file_mat.save(ss_abs.str(), raw_ascii);
           gpfile << "set title \"abs(Y" << ss_index.str() <<")\"\n";
           gpfile << "set logscale y\n";
           gpfile << "set term x11 " << gp_it++ << endl;
           gpfile << "plot \'" << ss_abs_file.str() << "\' u 1:2 title \'Input data\' with linespoints ls 1, \'"
                  << ss_abs_file.str()
                  << "\' u 1:3 title \'Simulation data\' with linespoints ls 2, \'"
                  << ss_abs_file.str()
                  << "\' u 1:4 title \'Blad bezwzgledny\' with linespoints\n" << endl;
    
           ostringstream ss_angle, ss_angle_file;
           ss_angle_file << "Y" << ss_index.str() << "_angle";
           ss_angle << file_name << ss_angle_file.str();
           file_mat.col(1) = gp_angle(input_temp);
           file_mat.col(2) = gp_angle(simulation_temp);
    
           file_mat.save(ss_angle.str(), raw_ascii);
           gpfile << "set title \"angle(Y" << ss_index.str() <<")\"\n";
           gpfile << "unset logscale y\n";
           gpfile << "set term x11 " << gp_it++ << endl;
           gpfile << "plot \'" << ss_angle_file.str() << "\' u 1:2 title \'Input data\' with linespoints ls 1, \'"
                  << ss_angle_file.str()
                  << "\' u 1:3 title \'Simulation data\' with linespoints ls 2\n" << endl;
       }
    }
    gpfile << "pause -1\n";

    gpfile.close();

}

mat gp_angle( cx_mat & data )
{
    mat data_real = real(data);
    mat data_imag = imag(data);
    int Ns = data.n_rows;
    mat result = zeros<mat>(Ns, 1);

    for ( int i = 0; i < Ns; i++ )
    {
        result(i) = std::atan2(data_imag(i), data_real(i));
    }

    return result;
}
