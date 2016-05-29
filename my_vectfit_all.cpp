#include "my_vectfit_all.h"

void vf_all::operator() ( const blocked_range<int>& r ) const
{
    for ( int rr = r.begin(); rr != r.end() ; rr++ )
    {
        int N = poles1->n_cols, //rzad rozwiazania
            Ns = f->n_cols; // liczba probek pomiarowych
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

        cx_mat A = zeros<cx_mat>(Ns, 2*N+2);
    
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
                A.col(m+1) = 1.0i / ( *s - poles(m) ) - 1.0i / ( *s - conj(poles(m)) );
            }
        }
    
        A.col(N) = ones<cx_mat>(1,Ns).st();
        A.col(N+1) = s->st();
    
        // wypelnienie prawej strony macierzy A
        for ( int i = 0; i < N; i++ )
        {
            A.col(i+N+2) = -strans(f->row(rr)) % A( span(0, Ns-1), i);
        }
    
        // obliczanie x metoda najmniejszych kwadratow Ax=b
        mat A_real = join_vert( real(A), imag(A) );
    
        cx_mat f_lsp = f->row(rr).st();
        mat f_lsp_real = join_vert( real(f_lsp), imag(f_lsp) );
    
        // dokompozycja QR macierzy A
        mat Q, R;
        qr_econ(Q, R, A_real);
    
        mat bb = Q.st() * f_lsp_real;
        bb = bb.rows(N+2, 2*N+1);
    
        mat AA = R( span(N+2, 2*N+1), span(N+2, 2*N+1) );

        // rozwiazanie ukladu rownan
        mat x = solve(AA, bb, solve_opts::fast);

        // przy pomocy metody wartosci wlasnych macierzy
        // obliczanie zer funkcji sigma - szukane bieguny
        cx_mat poles_diag = diagmat(poles);
        mat b_ones = ones<mat>(N,1);
        mat x_trans = x.st();

        int m = 0;
        mat poles_diag_real = zeros<mat>(N, N);
        for ( int n = 0 ; n < N ; n++ )
        {
            if ( m < N )
            {
                // sprawdzenie czy otrzymany biegun jest zespolony
                if ( abs(poles_diag(m,m)) > abs(real(poles_diag(m,m))) )
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
    
        cx_mat AA_res = zeros<cx_mat>(Ns, N+2);
    
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
                AA_res.col(m+1) = 1.0i / ( *s - poles(m)) - 1.0i / ( *s - conj(poles(m)));
            }
        }
    
        AA_res.col(N) = ones<cx_mat>(1, Ns).st();
        AA_res.col(N+1) = s->st();
    
        mat AA_res_real = join_vert( real(AA_res), imag(AA_res) );
    
        // obliczenie metoda najmniejszych kwadratow residuów szukanej funkcji
        x = solve(AA_res_real, f_lsp_real);
    
        // zapis residuów w postaci zespolonej (jesli takie istnieja) do obiektu z wynikami
        m = 0;
        for ( int i = 0; i < N ; i++ )
        {
            if ( imag_check(i) == 0 )
            {
                wynik->res(rr,m) = x(i);
            }
            else if ( imag_check(i) == 1 )
            {
                wynik->res(rr, m) = x(i) + 1.0i * x(i+1);
                wynik->res(rr, m+1) = conj(wynik->res(rr, m));
            }
            m++;
        }
    
        wynik->d(rr,0) = x(x.n_elem - 2);
        wynik->h(rr,0) = x(x.n_elem - 1);
        wynik->poles.row(rr) = poles.st();
    }
} 


SER my_vf_all_splitting(const cx_mat *f, const cx_vec *s, cx_mat *poles)
{
    SER wynik;
    int Nc = f->n_rows;
    int N = poles->n_cols;
    int Ns = s->n_elem;

    // prealokacja
    wynik.res = zeros<cx_mat>(Nc, N);
    wynik.poles = zeros<cx_mat>(Nc, N);
    wynik.h = zeros<mat>(Nc,1);
    wynik.d = zeros<mat>(Nc,1);
    wynik.err = 0.0;
    wynik.err_table = zeros<mat>(Nc, 1);

    // wielowatkowe uruchomienie algorytmu VF
    task_scheduler_init init();
    parallel_for( blocked_range<int>(0, Nc),
              vf_all(f, s, poles, &wynik) );

    // obliczanie bledu metody najmniejszych kwadratow dla kazdego z portow
    cx_mat f_check = zeros<cx_mat>(Nc, Ns);
    for ( int m = 0; m < Nc; m++ )
    {
        for ( int i = 0; i < Ns; i++ )
        {
            cx_double sk = s->operator()(i);
            for ( int j = 0; j < N; j++ )
            {
                f_check(m, i) = f_check(m, i) + wynik.res(m, j)
                                / ( sk - wynik.poles(m,j));
            }
            f_check(m, i) = f_check(m, i) + wynik.d(m, 0) + sk * wynik.h(m,0);
        } 
    }
     
    cx_mat diff = *f - f_check;

    wynik.err = sqrt( accu ( accu( pow(abs(diff), 2) ) ) );

    // wypelnienie macierzy z RMS dla kazdego z elementow Y
    for ( int j = 0; j < Nc; j++ )
    {
        wynik.err_table.row(j) = sqrt( accu ( pow(abs(diff.row(j)), 2)) ) ;
    }

    return wynik;
}
