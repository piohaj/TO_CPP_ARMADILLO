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
        cout << "Podano nieprawidlowy przedzial rzedow przyblizenia" << endl;
        delete[] wynik_iter;
        return wynik;
    }

    // sprawdzenie i ewentualne wymuszenie pasywnosci
    if ( conf.pasivity_check )
    {
        if ( check_and_make_passive( f ) != 0 )
        { 
            cout << "Problem podczas wymuszania pasywnosci" << endl;
            delete[] wynik_iter;
            return wynik;
        }
    }

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
                wynik_iter[high_iter] = my_vf_non_splitting(f, s, poles); 
    	        poles = wynik_iter[high_iter].poles;
    		
    	        if ( wynik_iter[high_iter].err < conf.tol )
                {
    	            break;
    	        }
            }
    	    cout << "Row: " << row << endl;
    	    cout << "Err: " << wynik_iter[high_iter].err << endl;
            high_iter++;
        }
        result_idx = choose_best_aprox( wynik_iter, row_iterations_num, Nc, split_strat, conf.rms_diff );
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
                wynik_iter[high_iter] = my_vf_all_splitting(&f, &s, &poles); 
    	        poles = wynik_iter[high_iter].poles;
    		
    	        if ( wynik_iter[high_iter].err < conf.tol )
                {
    	            break;
    	        }
            }
    	    cout << "Row: " << row << endl;
    	    cout << "Err: " << wynik_iter[high_iter].err << endl;
            wynik_iter[high_iter].poles.print("poles=");
            wynik_iter[high_iter].res.print("residues=");
            wynik_iter[high_iter].h.print("h=");
            wynik_iter[high_iter].d.print("d=");

            high_iter++;
        }
       result_idx = choose_best_aprox( wynik_iter, row_iterations_num, Nc, split_strat, conf.rms_diff );
       result_idx.print("result_idx");
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
                wynik_iter[high_iter] = my_vf_column_splitting(&f, &s, &poles); 
    	        poles = wynik_iter[high_iter].poles;
    		
    	        if ( wynik_iter[high_iter].err < conf.tol )
                {
    	            break;
    	        }
            }
    	    cout << "Row: " << row << endl;
    	    cout << "Err: " << wynik_iter[high_iter].err << endl;
            wynik_iter[high_iter].poles.print("poles=");
            wynik_iter[high_iter].res.print("residues=");
            wynik_iter[high_iter].h.print("h=");
            wynik_iter[high_iter].d.print("d=");

            high_iter++;
        }
        result_idx = choose_best_aprox( wynik_iter, row_iterations_num, Nc, split_strat, conf.rms_diff );
        result_idx.print("result_idx");
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
        mat bet = linspace<mat>(imag(s(0)), imag(s(Ns-1)), N/2);
    
        int m = 0;
        for ( int n = 0; n < N-1; n=n+2 )
        {
            double alf = -bet(m)*1e-2;
            poles(0, n) = cx_double(alf, bet(m));
            poles(0, n+1) = conj(poles(n));
            m++;
        }
    }
    else if ( split_strat == ALL_SPLITING )
    {
        poles = zeros<cx_mat>(Nc, N);
        mat bet = linspace<mat>(imag(s(0)), imag(s(Ns-1)), N/2);
        
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
    }
    else if ( split_strat == COLUMN_SPLITING )
    {
        int column_num = sqrt(Nc);
        poles = zeros<cx_mat>(column_num, N);
        mat bet = linspace<mat>(imag(s(0)), imag(s(Ns-1)), N/2);
        
        for ( int mm = 0; mm < column_num; mm++ )
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
    }
 
    return poles;
}

mat choose_best_aprox( SER *input, int size, int Nc, int split_strat, double rms_diff )
{
    mat result;
    if ( split_strat == NON_SPLITING )
    {
        result = zeros<mat>(1,1);
        double err = 10e6;
        for ( int i = 0; i < size; i++ )
        {
            if ( input[i].err < err )
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
            double err = 10e6;
            for ( int j = 0; j < size; j++ )
            {
                double err_temp = input[j].err_table[i];

                if ( abs( err - err_temp ) > rms_diff )
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
            double err = 10e6;
            for ( int j = 0; j < size; j++ )
            {
                if ( input[j].err_table[i] < err )
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
    
        for ( int n = 0; n < Nc ; n++ )
        {
            cx_mat poles_temp = iter_models[int(indexes(n))].poles.row(n);
            cx_mat res_temp = iter_models[int(indexes(n))].res.row(n);
            mat d_temp = iter_models[int(indexes(n))].d.row(n);
            mat h_temp = iter_models[int(indexes(n))].h.row(n);

            wynik.poles.row(n) = poles_temp;
            wynik.res.row(n) = res_temp;
            wynik.d.row(n) = d_temp;
            wynik.h.row(n) = h_temp;
        }
    }
    else if ( split_strat == COLUMN_SPLITING )
    {
        wynik.poles = zeros<cx_mat>(sqrt(Nc), max_row);
        wynik.res = zeros<cx_mat>(Nc, max_row);
        wynik.d = zeros<mat>(Nc, 1);
        wynik.h = zeros<mat>(Nc, 1);
    
        for ( int n = 0; n < Nc ; n++ )
        {
            cx_mat poles_temp = iter_models[int(indexes(n/sqrt(Nc)))].poles.row(n/sqrt(Nc));
            cx_mat res_temp = iter_models[int(indexes(n/sqrt(Nc)))].res.row(n);
            mat d_temp = iter_models[int(indexes(n/sqrt(Nc)))].d.row(n);
            mat h_temp = iter_models[int(indexes(n/sqrt(Nc)))].h.row(n);

            wynik.poles.row(n/sqrt(Nc)) = poles_temp;
            wynik.res.row(n) = res_temp;
            wynik.d.row(n) = d_temp;
            wynik.h.row(n) = h_temp;
        }
    }

    return wynik;
}


int read_conf( vf_opts& global_conf )
{
    Config cfg;

    // Read the file. If there is an error, report it and exit.
    try
    {
        cfg.readFile("file.conf");
    }
    catch(const FileIOException &fioex)
    {
        std::cerr << "I/O error while reading file." << std::endl;
        return 1;
    }
    catch(const ParseException &pex)
    {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
       return 1;
    }

    try
    {
        const Setting &root = cfg.getRoot();
        root.lookupValue("out_file_name", global_conf.out_file_name);
        root.lookupValue("in_file_name", global_conf.in_file_name);
        root.lookupValue("min_rms", global_conf.tol);
        root.lookupValue("rms_diff", global_conf.rms_diff);
        root.lookupValue("min_row", global_conf.min_row);
        root.lookupValue("max_row", global_conf.max_row);
        root.lookupValue("max_iters", global_conf.max_iters);
        root.lookupValue("R_max", global_conf.R_max);
        root.lookupValue("C_min", global_conf.C_min);
        root.lookupValue("spliting_strategy", global_conf.split_strat);
        root.lookupValue("pasivity_check", global_conf.pasivity_check);
    }
    catch (const SettingNotFoundException &nfex)
    {
        cerr << "No 'name' setting in configuration file." << endl;
        return 1;
    }

    return 0;

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
        cout << "Brak pasywnosci w plastrze " << j << endl;
        int idx = 0;
        
        while ( result == false && idx < 10 )
        {
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
