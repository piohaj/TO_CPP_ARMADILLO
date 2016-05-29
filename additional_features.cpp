#include "additional_features.h"

SER vf_high_level( const cx_mat& f, const cx_vec& s, int split_strat, int min_row, int max_row, int max_iters )
{
    SER wynik;
    int row_iterations_num = max_row - min_row + 1;
    SER *wynik_iter = new SER[row_iterations_num];
    cx_mat poles;
    int Nc = f.n_rows;
    int Ns = f.n_cols;
    mat result_idx;
    
    if ( max_row < min_row )
    {
        cout << "Podano nieprawidlowy przedzial rzedow przyblizenia" << endl;
        delete[] wynik_iter;
        return wynik;
    }

    if ( split_strat == NON_SPLITING )
    {
        int high_iter = 0;
        for ( int row = min_row; row <= max_row; row++ )
        {
            //biguny poczatkowe
            poles = prepare_input_poles(s, split_strat, Nc, row, Ns);

            int iter = 0;
            // wywolanie algorytmu
            for ( iter = 1; iter < max_iters; iter++ )
            {
                wynik_iter[high_iter] = my_vf_non_splitting(f, s, poles); 
    	        poles = wynik_iter[high_iter].poles;
    		
    	        if ( wynik_iter[high_iter].err < 1e-10 )
                {
    	            break;
    	        }
            }
    	    cout << "Row: " << row << endl;
    	    cout << "Err: " << wynik_iter[high_iter].err << endl;
            high_iter++;
        }
        result_idx = choose_best_aprox( wynik_iter, row_iterations_num, Nc, split_strat );
        wynik = wynik_iter[int(result_idx(0))];
    }
    else if ( split_strat == ALL_SPLITING )
    {
        int high_iter = 0;
        for ( int row = min_row; row <= max_row; row++ )
        {
            //bieguny poczatkowe
            poles = prepare_input_poles(s, split_strat, Nc, row, Ns);

            int iter = 0;
            // wywolanie algorytmu
            for ( iter = 1; iter < max_iters; iter++ )
            {
                wynik_iter[high_iter] = my_vf_all_splitting(&f, &s, &poles); 
    	        poles = wynik_iter[high_iter].poles;
    		
    	        if ( wynik_iter[high_iter].err < 1e-10 )
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
       result_idx = choose_best_aprox( wynik_iter, row_iterations_num, Nc, split_strat );
       // wynik = wynik_iter[result_idx];dd
       result_idx.print("result_idx");
       wynik = cumulate_model( result_idx, wynik_iter, Nc, max_row);
    }
    else
    {
        cout << "Nieznana strategia podzialu Y" << endl;
        delete[] wynik_iter;
        return wynik;
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
 
    return poles;
}

mat choose_best_aprox( SER *input, int size, int Nc, int split_strat )
{
    mat result = zeros<mat>(Nc,1);
    if ( split_strat == NON_SPLITING )
    {
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
        for ( int i = 0; i < Nc; i++ )
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


SER cumulate_model( mat& indexes, SER *iter_models, int Nc, int max_row )
{
    SER wynik;
    wynik.poles = zeros<cx_mat>(Nc, max_row);
    wynik.res = zeros<cx_mat>(Nc, max_row);
    wynik.d = zeros<mat>(Nc, 1);
    wynik.h = zeros<mat>(Nc, 1);

    for ( int n = 0; n < Nc ; n++ )
    {
        cout << "index="<< indexes(n) << endl;
        wynik.poles.row(n) = iter_models[int(indexes(n))].poles.row(n);
        wynik.res.row(n) = iter_models[int(indexes(n))].res.row(n);
        wynik.d.row(n) = iter_models[int(indexes(n))].d.row(n);
        wynik.h.row(n) = iter_models[int(indexes(n))].h.row(n);
    }

    return wynik;
}
