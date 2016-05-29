#include "additional_features.h"
#include "my_vectfit_non.h"

SER vf_high_level( const cx_mat& f, const cx_vec& s, int split_strat, int min_row, int max_row, int max_iters )
{
    SER wynik;
    int row_iterations_num = max_row - min_row + 1;
    SER *wynik_iter = new SER[row_iterations_num];
    cx_mat poles;
    int Nc = f.n_rows;
    int Ns = f.n_cols;
    
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
    }
    else
    {
        cout << "Nieznana strategia podzialu Y" << endl;
        delete[] wynik_iter;
        return wynik;
    }

    int result_idx = choose_best_aprox( wynik_iter, row_iterations_num );
    wynik = wynik_iter[result_idx];

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
 
    return poles;
}

int choose_best_aprox( SER *input, int size )
{
    int result = 0;
    double err = 10e6;
    for ( int i = 0; i < size; i++ )
    {
        if ( input[i].err < err )
        {
            err = input[i].err;
            result = i;
        }
    }
    return result;
}  
