#include "vectfit3.h"

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

void vectorfit3(cx_mat f, cx_mat s, cx_mat poles, cx_mat weight)
{
    //struktura z opcjami algorytmu
    opts opt;
 
    // pierwsza wersja optaultowe ustawienia
    opt.relax=1;      //Use vector fitting with relaxed non-triviality constraint
    opt.stable=1;     //Enforce stable poles
    opt.asymp=3;      //Include only D in fitting (not E)   
    opt.skip_pole=0;  //Do NOT skip pole identification
    opt.skip_res=0;   //Do NOT skip identification of residues (C,D,E) 
    opt.complx_ss=1;   //Create complex state space model
    opt.spy1=0;       //No plotting for first stage of vector fitting
    opt.spy2=1;       //Create magnitude plot for fitting of f(s) 
    opt.logx=1;       //Use logarithmic abscissa axis
    opt.logy=1;       //Use logarithmic ordinate axis 
    opt.errplot=1;    //Include deviation in magnitude plot
    opt.phaseplot=1;  //exclude plot of phase angle (in addition to magnitiude)
    opt.legend=1;     //Do include legends in plots

    // zmienne zakresu
    double TOLlow = 1e-18,
           TOLhigh = 1e18;

    // rozmiar macierzy biegunow
    int a, b;

    a = poles.n_rows;
    b = poles.n_cols;

    if ( s(0) == cx_double(0,0) && a == 1 )
    {
        cout << "Wszedlem w ifa" << endl;
        if ( poles(0) == cx_double(0, 0) && poles(1) != cx_double(0, 0) )
        {
			cout << "Wszedlem w pierwszego podifa" << endl;
			poles(0)=-1;
		}
		else if ( poles(1) == cx_double(0, 0) && poles(0) != cx_double(0, 0) )
		{
			cout << "Wszedlem w drugiego podifa" << endl;
			poles(1) = -1;
		}
		else if ( poles(0) == cx_double(0, 0) && poles(1) == cx_double(0, 0) )
		{
			cout << "Wszedlem w trzeciego podifa" << endl;
			poles(1) = cx_double(-1, 10);
			poles(2) = cx_double(-1, -10);
		}
	}

    a = s.n_rows;
    b = s.n_cols;

    if ( a < b )
    {
        s = inv(s);
    }

    cx_mat LAMBD = diagmat(poles);
    int Ns = s.n_rows,
        N = LAMBD.n_rows,
        Nc = f.col(1).n_rows;
    cx_mat B = ones<cx_mat>(N, 1);
    cx_mat SERA = poles,
           SERC = zeros<cx_mat>(Nc, N),
           SERD = zeros<cx_mat>(Nc, 1),
           SERE = zeros<cx_mat>(Nc, 1);

    int common_weight = 0;
    weight = weight.t();
    if ( weight.row(0).n_cols == 1 )
    {
        common_weight = 1;
    }
    else if ( weight.row(0).n_cols == Nc )
    {
        common_weight = 0;
    }
    else
    {
        cout << "ERROR in vectfit3: Invalid size of array weight" << endl;
    }

    int offs = 0;
    if ( opt.asymp == 1 )
    {
        offs = 0;
    }
    else if ( opt.asymp == 2 )
    {
        offs = 1;
    }
    else
    {
        offs = 2;
    }

// ==================================================================================
// POLE INDENTIFICATION
// ==================================================================================
   cx_mat Escale, Dk;
      mat cindex;

   if ( opt.skip_pole != 1 )
   {
       Escale = zeros<cx_mat>(1, Nc+1 );

       //==========================================================================
       // Finding out which starting poles are complex
       // ========================================================================
       cindex = zeros<mat>(1,N);

       for ( int m = 0; m < N; m++ )
       {
           if ( imag( LAMBD(m,m)) != 0 )
           {
               if ( m == 0 )
               {
                   cindex(m)==1;
               }
               else
               {
                   if ( cindex(m-1) == 0 || cindex(m-1) == 2 )
                   {
                       cindex(m) = 1;
                       cindex(m+1) = 2;
                   }
                   else
                   {
                       cindex(m) = 2;
                   }
               }
           }
       }

       cout << "LAMBD:" << LAMBD << endl;

       // ===============================================================================
       // Building system - matrix :
       // ===============================================================================
       Dk = zeros<cx_mat>(Ns, N);

       for ( int m = 0; m < N; m++ )
       {
           if ( cindex(m) == 0 ) //real pole
           {
               Dk.col(m) = cx_double(1,0) / (s-LAMBD(m,m));
           }
           else if ( cindex(m) == 1 ) //complex pole 1st part
           {
              Dk.col(m) = 1 / (s - LAMBD(m,m)) + 1 / (s - conj(LAMBD(m,m)));
              Dk.col(m+1) = 1i / (s-LAMBD(m,m)) - 1i / (s - conj(LAMBD(m,m)));
           }    
       }


        if ( opt.asymp == 1 || opt.asymp == 2 )
        {
            Dk = join_horiz(Dk, ones<cx_mat>(1,Ns).t());
            cout << Dk << endl;
        }
        else if ( opt.asymp == 3 )
        {
            Dk = join_horiz(Dk, ones<cx_mat>(1,Ns).t());
            Dk = join_horiz(Dk, s);
            cout << Dk << endl;
        }

        // Scaling for las row of LS - problem (pole identyfication)

        double scale = 0;

        for ( int m = 0; m < Nc; m++ )
        {
            if ( weight.row(0).n_cols == 1 )
            {
                scale = scale + ( norm(weight.col(0) % f.row(m).t()) ) * ( norm(weight.col(0) % f.row(m).t()) );
            }
            else
            {
                scale = scale + ( norm(weight.col(m) % f.row(m).t()) ) * ( norm(weight.col(m) % f.row(m).t()) );
            }    
        }

        scale = sqrt(scale)/Ns;

        cx_mat weig;

        if ( opt.relax == 1 )
        {
            cx_mat AA = zeros<cx_mat>(Nc*(N+1), N+1);
            cx_mat bb = zeros<cx_mat>(Nc*(N+1), 1);

            for ( int n = 0; n < Nc ; n++ )
            {
                cx_mat A = zeros<cx_mat>(Ns, N+offs+N+1);

                if ( common_weight == 1 )
                {
                    weig = weight;    
                }
                else
                {
                    weig = weight.col(n);
                }

                for ( int m = 0; m < N +offs ; m++ )
                {
                    
                }
            }

        }











    }           
}
