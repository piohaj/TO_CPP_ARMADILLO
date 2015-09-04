#include "my_vectfit.h"

cx_mat logspace(double a, double b, int n)
{
    cx_mat return_mat = linspace<cx_mat>(a, b, n);
    double wykladnik;
    cx_double simple_result;

    for (int i = 0; i < n; i++)
    {
        wykladnik = real(return_mat(i));
        simple_result = pow(10, wykladnik);
        return_mat(i) = simple_result;
    }

    return return_mat;
}

int sign( double x )
{
    if ( x == 0 )
    {
        return 0;
    }
    else if ( x < 0 )
    {
        return -1;
    }
    else 
    {
        return 1;
    }
}

// funkcja oblicza wspolczynniki modelu 
// mozliwe zespolone dane
// do rozwiazania rownania metoda najmniejszych kwadratow uzyto dekompozycji QR
// algorymt zaklada istnienie stalej h
cx_mat my_vectorfit3(cx_mat f, cx_mat s, cx_mat poles, cx_mat weight)
{
    int N = poles.n_elem, //rzad rozwiazania
        Ns = f.n_elem; // liczba probek pomiarowych
    mat imag_check = zeros<mat>(1, N); //wektor informujacy czy dany biegun jest zespolony

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

    cx_mat A = zeros<cx_mat>(Ns, 2*N+1);

    // wypelnienie lewej strony macierzy A
    for ( int m = 0; m < N ; m++ )
    {
        if ( imag_check(m) == 0 )
        {
            A.col(m) = 1 / (s - poles(m));
        }
        else if ( imag_check(m) = 1 )
        {
            A.col(m) = 1 / (s - poles(m)) + 1 / (s - conj(poles(m)) );
            A.col(m+1) = 1i / (s-poles(m)) - 1i / (s - conj(poles(m)) );
        }
    }

    A.col(N+1) = ones<cx_mat>(1,Ns).st();
    
    // wypelnienie prawej strony macierzy A

    for ( int i = 0; i < N; i++ )
    {
        A.col(i+N+1) = - strans(f(0, span(0, Ns-1))) % A( span(0, Ns-1), i);
    }

    // obliczanie x metoda najmniejszych kwadratow Ax=b
    mat A_real = join_vert( real(A), imag(A) );
    cx_mat f_lsp = f.st();
    mat f_lsp_real = join_vert( real(f_lsp), imag(f_lsp) );

    // dokompozycja QR macierzy A
    mat Q, R;
    qr(Q, R, A_real);

    mat AA = Q.st() * f_lsp_real;
    AA = AA.rows(N+1, 2*N);

    mat bb = R( span(N+1, 2*N), span(N+1, 2*N) );
    mat x = solve(AA, bb);
    
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
				b_ones(m,1)=2; b_ones(m+1,1)=0;
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
    mat H = poles_diag_real - b_ones * x;
    cx_mat Hi = cx_mat(H, zeros<mat>(N,N));
    vec eigval2;
    eig_sym(eigval2, poles, Hi);

    cout << eigval2;

    return poles;
}
