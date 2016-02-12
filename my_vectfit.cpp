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

void QR_calculation::operator() ( const blocked_range<int>& r ) const
{
    for ( int m = r.begin(); m != r.end(); m++ )
    {
        cx_mat AA_port = *A;
        // wypelnienie prawej strony macierzy A
        cx_mat part;
        for ( int i = 0; i < N; i++ )
        {
            part = -strans(f->operator()(m, span(0, Ns-1))) % A->operator()( span(0, Ns-1), i);
            AA_port.col(i+N+1) = part;
        }

        // obliczanie x metoda najmniejszych kwadratow Ax=b
        mat A_real = join_vert( real(AA_port), imag(AA_port) );
        //AA_port.reset();
  
        cx_mat f_lsp = f->row(m).st();
        mat f_lsp_real = join_vert( real(f_lsp), imag(f_lsp) );

        // dokompozycja QR macierzy A
        //cout<< "QR " << m << endl; 
        mat Q, R;
        qr_econ(Q, R, A_real);
        //A_real.reset();

        mat bb = Q.st() * f_lsp_real;
        //Q.reset();

        bb_poles->rows(m*N, (m+1)*N-1) = bb.rows(N+1, 2*N); 
        AA_poles->operator()( span(m*N, (m+1)*N-1), span( 0, N-1 ) ) = R( span(N+1, 2*N), span(N+1, 2*N) );

        //R.reset();
        //bb.reset();
    }
}


// funkcja oblicza wspolczynniki modelu 
// mozliwe zespolone dane
// do rozwiazania rownania metoda najmniejszych kwadratow uzyto dekompozycji QR
// algorymt zaklada istnienie stalej h
// algorytm przystosowany do obliczenia modelu wieloportowego
SER my_vectorfit(const cx_mat& f, const cx_vec& s, cx_vec poles)
{
    int N = poles.n_elem, //rzad rozwiazania
        Ns = s.n_elem, // liczba probek pomiarowych
        Nc = f.n_rows; // liczba portow do symulacji
    SER wynik; // struktura z wynikiem dzialania algorytmu

    if ( N > Ns )
    {
       cout << "Podane zle dane do algorymtu" << endl;
       return wynik;
    }

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
            A.col(m) = cx_double(1,0) / (s - poles(m));
        }
        else if ( imag_check(m) == 1 )
        {
            A.col(m) = cx_double(1,0) / (s - poles(m)) + cx_double(1,0) / (s - conj(poles(m)) );
            A.col(m+1) = 1.0i / (s-poles(m)) - 1.0i / (s - conj(poles(m)) );
        }
    }

    A.col(N) = ones<cx_mat>(1,Ns).st();

    // przygotowanie macierzy pod wiele portow
    mat AA_poles = zeros<mat>(Nc*N, N);
    mat bb_poles = zeros<mat>(Nc*N, 1);

    // wielowatkowe (TTB) oblicznie wspolczynnikow AA_poles - QR rownolegle
    task_scheduler_init init();
    parallel_for(blocked_range<int>(0, Nc),
           QR_calculation( &A, &f, N, Ns, &AA_poles, &bb_poles) );

    // obliczenie x dla wszystkich portow badanego ukladu
    mat x = solve(AA_poles, bb_poles);

    //x.print("x=");

    
    // przy pomocy metody wartosci wlasnych macierzy obliczanie zer funkcji sigma - szukane bieguny
    cx_mat poles_diag = diagmat(poles); // tu jeszcze bieguny wejsciowe
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
      //  cout << "m="<<m<<endl;
     //   cout << "n="<<n<<endl;
    }

    //obliczanie wartosci wlasnych macierzy
    mat H = poles_diag_real - b_ones * x_trans;
    //poles_diag_real.print("poles_diag_real=");
    //b_ones.print("b_ones=");
    //x_trans.print("x_trans=");

    poles = eig_gen(H);
    //cout << "poles: "<< poles << endl;
  
    H.reset();

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
     //  cout << "imag_check("<<i<<")="<<imag_check(i)<<endl;
    }


    //cout << "imag_check " << imag_check << endl;
    cx_mat AA_res = zeros<cx_mat>(Ns, N+1);

    // wypelnienie lewej strony macierzy AA_res
    for ( int m = 0; m < N; m++ )
    {
        if ( imag_check(m) == 0 )
        {
            AA_res.col(m) = 1 / (s - poles(m));
        }
        else if ( imag_check(m) == 1 )
        {
            AA_res.col(m) = 1 / (s - poles(m)) + 1 / (s - conj(poles(m)));
            AA_res.col(m+1) = 1.0i / (s - poles(m)) - 1.0i / (s - conj(poles(m)));
        }
    }

    AA_res.col(N) = ones<cx_mat>(1, Ns).st();

    mat AA_res_real = join_vert( real(AA_res), imag(AA_res) );

    // obliczenie metoda najmniejszych kwadratow residuów szukanej funkcji
    // dla kazdego portu ukladu
    mat f_lsp_res = zeros<mat>(2*Ns, Nc);

    for (int m = 0; m < Nc; m++)
    {
        f_lsp_res.col(m) = join_vert( real(f.row(m)).st(), imag(f.row(m)).st() );
    }

    //x = inv( (AA_res_real.st() * AA_res_real ) ) * AA_res_real.st() * f_lsp_real;
    x = solve(AA_res_real, f_lsp_res);

    // zapis residuów w postaci zespolonej (jesli takie istnieja)
    wynik.res = zeros<cx_mat>(Nc, N);
    
    for ( int m_port = 0; m_port < Nc; m_port++ )
    {
        m = 0;
        for ( int i = 0; i < imag_check.n_elem ; i++ )
        {
            if ( imag_check(i) == 0 )
            {
                wynik.res(m_port, m) = x(i, m_port);
            }
            else if ( imag_check(i) == 1 )
            {
                wynik.res(m_port, m) = x(i, m_port) + 1.0i * x(i+1, m_port);
                wynik.res(m_port, m+1) = conj(wynik.res(m_port, m));
            }
            m++;
        }
    }

    // rzeczywiste wspolczynniki h
    wynik.h = zeros<mat>(Nc, 1);
    for ( int m = 0; m < Nc; m++ )
    {
        wynik.h(m, 0) = x(N, m);
    }


//    cout << "Wynik.h: " << wynik.h << endl;
    // wstawienie biegunow do struktury wynikow
    wynik.poles = poles;

     // obliczanie bledu metody najmniejszych kwadratow dla kazdego z portow
     cx_mat f_check = zeros<cx_mat>(Nc, Ns);
     for ( int m = 0; m < Nc; m++ )
     {
         for ( int i = 0; i < Ns; i++ )
         {
             for ( int j = 0; j < N; j++ )
             {
                 f_check(m, i) = f_check(m, i) + wynik.res(m, j) / ( s(i) - wynik.poles(j));
             }
             f_check(m, i) = f_check(m, i) + wynik.h(m, 0);
         } 
     }
     
     cx_mat diff = f - f_check;

     wynik.err = sqrt( accu ( accu( pow(abs(diff), 2) ) ) );

     return wynik;
}


input_data prepare_sample_data()
{
    input_data data;
    int Ns = 101;
    int N = 3;
    data.f = zeros<cx_mat>(2, Ns);
    data.s = 2 * 3.14 * 1.0I * linspace<cx_mat>(1, 550, Ns);

    for ( int k = 0; k < Ns ; k++ )
    {
        cx_double sk = data.s(k);
        
        data.f(0, k) = cx_double(2,0)/(sk+cx_double(5,0)) + cx_double(30, 40)/(sk - cx_double(-100,500)) + cx_double(30,-40)/(sk-cx_double(-100, -500)) + cx_double(0.5, 0);
    } 

    for ( int kk = 0; kk < Ns ; kk++ )
    {
        cx_double sk = data.s(kk);
        
        data.f(1, kk) = cx_double(3,0)/(sk+cx_double(5,0)) + cx_double(300, 40)/(sk - cx_double(-100,500)) + cx_double(300,-40)/(sk-cx_double(-100, -500)) + cx_double(0.9, 0);
    } 

    return data;
}
 

input_data load_data_from_file( int N, int Nc, int Ns )
{
    input_data data;
    mat f_real, f_imag, s_real, s_imag;
    ostringstream oss;
    oss << "N" << N << "_Nc" << Nc << "_Ns" << Ns << ".dat";
    string file_sufix = oss.str();

    string f_real_file = "./benczmarki/f_real_" + file_sufix;
    string f_imag_file = "./benczmarki/f_imag_" + file_sufix;
    string s_real_file = "./benczmarki/s_real_" + file_sufix;
    string s_imag_file = "./benczmarki/s_imag_" + file_sufix;
   
    if ( f_real.load(f_real_file) == false )
    {
        throw "Brak pliku " + f_real_file;
    }
    if ( f_imag.load(f_imag_file) == false )
    {
        throw "Brak pliku " + f_imag_file;
    }
    if ( s_real.load(s_real_file) == false )
    {
        throw "Brak pliku " + s_real_file;
    }
    if ( s_imag.load(s_imag_file) == false )
    {
        throw "Brak pliku " + s_imag_file;
    }
     
    data.f = cx_mat(f_real, f_imag);
    data.s = cx_mat(s_real, s_imag).st();

    return data;
}

