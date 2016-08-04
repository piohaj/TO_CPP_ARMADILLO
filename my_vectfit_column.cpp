#include "my_vectfit_column.h"

void vf_column::operator() ( const blocked_range<int>& r ) const
{
    for ( int rr = r.begin(); rr < r.end() ; rr++ )
    {
        int N = poles1->n_cols, //rzad rozwiazania
            Ns = f->n_cols, // liczba probek pomiarowych
            column_elems = sqrt(f->n_rows); // liczba kolumn macierzy Y
        mat imag_check = zeros<mat>(1, N); //wektor informujacy czy dany biegun jest zespolony
        cx_mat poles = poles1->row(rr);
    
        for ( int i = 0; i < N; i++ )
        {
            if ( imag(poles(i)) != 0 ) 
            {
                if ( i == 0 )
                {
                    imag_check(i) = 1;
                }
                else
                {
                    if ( imag_check(i-1) == 0 || imag_check(i-1) == 2 )
                    {
                        imag_check(i) = 1; imag_check(i+1) = 2;
                    }
                    else
                    {
                        imag_check(i) = 2;
                    }
                }
            }
        }

        cx_mat A = zeros<cx_mat>(Ns, 2*N+RC_offset);
    
        // wypelnienie lewej strony macierzy A
        for ( int m = 0; m < N ; m++ )
        {
            if ( imag_check(m) == 0 )
            {
                A.col(m) = cx_double(1,0) / ( *s - poles(m));
            }
            else if ( imag_check(m) == 1 )
            {
                A.col(m) = cx_double(1,0) / ( *s - poles(m)) + cx_double(1,0) / ( *s - conj(poles(m)) );
                A.col(m+1) = cx_double(0,1) / ( *s - poles(m) ) - cx_double(0,1) / ( *s - conj(poles(m)) );
            }
        }
    
        if ( RC_offset )
        {
            A.col(N) = ones<cx_mat>(1,Ns).st();
            A.col(N+1) = *s;
        }

        // przygotowanie macierzy pod wszystkie elementy Y dla danej kolumny
        mat AA_poles = zeros<mat>(column_elems*N, N);
        mat bb_poles = zeros<mat>(column_elems*N, 1);
    
        int n = 0;
        //zmienna m - po zakresie danej kolumny Y w f
        for ( int m = rr*column_elems; m < column_elems*(rr+1); m++ ) 
        {
            // wypelnienie prawej strony macierzy A
            for ( int i = 0; i < N; i++ )
            {
                A.col(i+N+RC_offset) = -strans(f->operator()(m, span(0, Ns-1))) % A( span(0, Ns-1), i);
            }
        
            mat A_real = join_vert( real(A), imag(A) );

            cx_mat f_lsp = f->row(m).st();
            mat f_lsp_real = join_vert( real(f_lsp), imag(f_lsp) );
        
            // dokompozycja QR macierzy A
            mat Q, R;
            qr_econ(Q, R, A_real);
        
            mat bb = Q.st() * f_lsp_real;
        
            bb_poles.rows(n*N, (n+1)*N-1) = bb.rows(N+RC_offset, 2*N+RC_offset-1);
            AA_poles( span(n*N, (n+1)*N-1 ), span( 0, N-1 ) ) = R( span(N+RC_offset, 2*N+RC_offset-1), span(N+RC_offset, 2*N+RC_offset-1) );

            n++;
        } 
        // rozwiazanie ukladu rownan
        mat x = solve(AA_poles, bb_poles);

        // przy pomocy metody wartosci wlasnych macierzy obliczanie zer funkcji sigma - szukane bieguny
        cx_mat poles_diag = diagmat(poles);
        mat b_ones = ones<mat>(N,1);
        mat x_trans = x.st();
        //int x_trans_end = x_trans.n_elem - 1;
        //x_trans = x_trans( 1, span( x_trans_end-N, x_trans_end ) );
        int m = 0;
        mat poles_diag_real = zeros<mat>(N, N);
        for ( int n = 0 ; n < N ; n++ )
        {
            if ( m < N )
            {
                if ( abs(poles_diag(m,m)) > abs(real(poles_diag(m,m))) ) // sprawdzenie czy otrzymany biegun jest zespolony
                {
    			poles_diag_real(m+1,m)=-imag(poles_diag(m,m));
    			poles_diag_real(m,m+1)=imag(poles_diag(m,m));
    			poles_diag_real(m,m)=real(poles_diag(m,m));
    			poles_diag_real(m+1,m+1)=poles_diag_real(m,m);
    			b_ones(m,0)=2; b_ones(m+1,0)=0;
    			m=m+1;
                }
                else
                {
                    poles_diag_real(m,m) = real(poles_diag(m, m));
                }
            }
            m++;
        }
    
        //obliczanie wartosci wlasnych macierzy
        mat H = poles_diag_real - b_ones * x_trans;
        //cx_mat Hi = cx_mat(H, zeros<mat>(N,N));
    	poles = eig_gen(H);

        // force stable poles
        force_stable_poles( poles );

       //=============================================
       // obliczanie residuów szukanej funkcji
       // ============================================
        // sparwdzenie ktore bieguny sa zespolone
        imag_check = zeros<mat>(1,N);
        for ( int i = 0; i < N; i++ )
        {
            if ( imag(poles(i)) != 0 ) 
            {
                if ( i == 0 )
                {
                    imag_check(i) = 1;
                }
                else
                {
                    if ( imag_check(i-1) == 0 || imag_check(i-1) == 2 )
                    {
                        imag_check(i) = 1; imag_check(i+1) = 2;
                    }
                    else
                    {
                        imag_check(i) = 2;
                    }
                }
            }
        }
    
        cx_mat AA_res = zeros<cx_mat>(Ns, N+RC_offset);
    
        // wypelnienie lewej strony macierzy AA_res
        for ( int m = 0; m < N; m++ )
        {
            if ( imag_check(m) == 0 )
            {
                AA_res.col(m) = 1 / ( *s - poles(m));
            }
            else if ( imag_check(m) == 1 )
            {
                AA_res.col(m) = 1 / ( *s - poles(m)) + 1 / ( *s - conj(poles(m)));
                AA_res.col(m+1) = cx_double(0,1) / ( *s - poles(m)) - cx_double(0,1) / ( *s - conj(poles(m)));
            }
        }
    
        if ( RC_offset )
        {
            AA_res.col(N) = ones<cx_mat>(1, Ns).st();
            AA_res.col(N+1) = *s;
        }
    
        mat AA_res_real = join_vert( real(AA_res), imag(AA_res) );
    
        // obliczenie metoda najmniejszych kwadratow residuów szukanej funkcji
        // dla kazdego z elementow danej kolumny macierzy Y
        mat f_lsp_res = zeros<mat>(2*Ns, column_elems);

        int nn = 0;
        for ( int m = rr*column_elems; m < column_elems*(rr+1); m++ )
        {
            f_lsp_res.col(nn) = join_vert( real(f->row(m)).st(), imag(f->row(m)).st() ); 
            nn++;
        }

        x = solve(AA_res_real, f_lsp_res);

        if ( RC_offset )
        {
            for ( int m = 0; m < column_elems; m++ )
            {
                wynik->d(rr*column_elems+m,0) = x(N, m);
            }
    
            for ( int m = 0; m < column_elems; m++ )
            {
                wynik->h(rr*column_elems+m,0) = x(N+1, m);
            }
        }
        wynik->poles.row(rr) = poles.st();

        // zapis residuów w postaci zespolonej (jesli takie istnieja) i wpisanie ich do obiektu
        // z wynikami
        for ( int m_port = 0; m_port < column_elems; m_port++ )
        {
            m = 0;
            for ( int i = 0; i < imag_check.n_elem; i++ )
            {
                if ( imag_check(i) == 0 )
                {
                    wynik->res(rr*column_elems+m_port, m) = x(i, m_port);
                }
                else if ( imag_check(i) == 1 )
                {
                    wynik->res(rr*column_elems+m_port, m) = x(i, m_port) + cx_double(0,1) * x(i+1, m_port);
                    wynik->res(rr*column_elems+m_port, m+1) = conj(wynik->res(rr*column_elems+m_port, m));
                }
                m++;
            }
        }
    }
} 


SER my_vf_column_splitting(const cx_mat *f, const cx_vec *s, cx_mat *poles, vf_opts & conf )
{
    SER wynik;
    int Nc = f->n_rows;
    int N = poles->n_cols;
    int Ns = s->n_elem;
    int column_num = sqrt(Nc);
    int RC_offset = 0;
    if ( conf.calc_parallel_RC ) RC_offset = 2;

    if ( pow(column_num, 2) != Nc )
    {
        cout << "Bledne dane wejsciowe, macierz Y nie jest kwadratowa" << endl;
        return wynik;
    }

    // prealokacja
    wynik.res = zeros<cx_mat>(Nc, N);
    wynik.poles = zeros<cx_mat>(column_num, N);
    wynik.h = zeros<mat>(Nc,1);
    wynik.d = zeros<mat>(Nc,1);
    wynik.err = 0.0;
    wynik.err_table = zeros<mat>(column_num, 1);

    // wielowatkowe uruchomienie algorytmu VF
    task_scheduler_init init();
    parallel_for( blocked_range<int>(0, column_num),
              vf_column(f, s, poles, &wynik, RC_offset) );

    // obliczanie bledu metody najmniejszych kwadratow do obliczonego modelu
    rms_err_calculation(&wynik, f, s, N);

    return wynik;
}

void rms_err_calculation(SER *wynik, const cx_mat *f, const cx_vec *s, int N)
{
    int Nc = f->n_rows;
    int Ns = s->n_elem;
    int column_num = sqrt(Nc);
    cx_mat f_check = zeros<cx_mat>(Nc, Ns);

    for ( int m = 0; m < Nc; m++ )
    {
        for ( int i = 0; i < Ns; i++ )
        {
            cx_double sk = s->operator()(i);
            for ( int j = 0; j < N; j++ )
            {
                f_check(m, i) = f_check(m, i) + wynik->res(m, j) / ( sk - wynik->poles(m/column_num, j) );
            }
            f_check(m, i) = f_check(m, i) + sk * wynik->h(m, 0) + wynik->d(m,0);
        } 
    }

    // mat diff_real = real(*f - f_check);
    // mat diff_imag = imag(*f - f_check);

    //wynik->err = sqrt( ( accu( pow(diff_real, 2) + pow(diff_imag, 2) ) ) / Ns );

    double rms_err_db = sqrt( accu( pow( abs( *f - f_check ), 2 ) ) /
                     accu ( pow ( abs(*f), 2 ) ) );
    wynik->err = 20 * log10( rms_err_db );

    // wypelnienie macierzy z RMS dla kazdego z elementow Y
    for ( int j = 0; j < Nc; j=j+2 )
    {
//        wynik->err_table.row(j/sqrt(Nc)) += sqrt( accu( pow(diff_real.row(j), 2)
//                                            + pow(diff_imag.row(j), 2) ) / Ns ) ;
        double rms_err_row_db = sqrt( accu( pow( abs( f->rows(j,j+1) - f_check.rows(j,j+1) ), 2) ) /
                     accu ( pow ( abs(f->rows(j,j+1)), 2) ) );
        wynik->err_table.row(j/sqrt(Nc)) += 20 * log10( rms_err_row_db );
    }
}
