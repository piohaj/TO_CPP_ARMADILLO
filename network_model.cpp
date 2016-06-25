#include "network_model.h"


void parse_SER(SER *input_SER, Y_network_data *output_network_data )
{
    int N = input_SER->poles.n_cols;
    int Nc = input_SER->res.n_rows;
    cx_mat poles = input_SER->poles;

    int N_poles = poles.n_rows;

    mat imag_check = NO_ELEM*ones<mat>(N_poles, N); //wektor informujacy czy dany biegun jest zespolony

    for ( int j = 0 ; j < N_poles ; j++ )
    {
        for ( int i = 0; i < N; i++ )
        {
            if ( imag(poles(j, i)) != 0 && abs(poles(j,i)) != 0.0 ) 
            {
                if ( i == 0 )
                {
                    imag_check(j, i) = 1;
                }
                else
                {
                    if ( imag_check(j, i-1) == 0 || imag_check(j, i-1) == 2 )
                    {
                        imag_check(j, i) = 1; imag_check(j, i+1) = 2;
                    }
                    else
                    {
                        imag_check(j, i) = 2;
                    }
                }
            }
            else if ( abs(poles(j,i)) != 0.0 )
            {
                imag_check(j, i) = 0;
            }
        }
    }

    imag_check.print("imag_check=");

    // sprawdzenie strategi podzialu macierzy Y
    // SPLLITTING_STRT
    // 1 - no-splitting
    // 2 - all-splitting
    // 3 - column-splitting
    int SPLITING_STRT = NON_SPLITING; // no-splitting
    if ( N_poles == Nc ) // all_splitting
    {
        SPLITING_STRT = ALL_SPLITING;
    }
    else if ( N_poles == sqrt(Nc) ) // column splliting    
    {
        SPLITING_STRT = COLUMN_SPLITING;
    }
    
    // zmienna do prawidlowego wyboru biegunow dla danego elementu macierzy Y
    int nn = 0;

    for ( int i = 0; i < Nc ; i++ )
    {
        if ( SPLITING_STRT == NON_SPLITING ) { nn = 0; }
        else if ( SPLITING_STRT == ALL_SPLITING ) { nn = i; }
        else if ( SPLITING_STRT == COLUMN_SPLITING ) { nn = i/sqrt(Nc); } 

        Y_network_data network_data_sample; // struktura z danymi sieci dla pojedynczego elementu macierz Y
        // wyczyszczenie struktury pomocniczej
        network_data_sample.real_pole_nets.clear();
        network_data_sample.imag_pole_nets.clear();
        network_data_sample.R = 0;
        network_data_sample.C = 0;

        real_pole_net real_pole_net_sample;
        imag_pole_net imag_pole_net_sample;

        // czy element z przekatnej macierzy Y
        int is_diag = -1; // nie jest
        if ( ( i - floor(i/(sqrt(Nc)+1))*(sqrt(Nc)+1) ) == 0 ) is_diag = 1; // jest

        // rozwiazanie problemu z zerowym R rownoleglym
        if ( input_SER->d(i) < 1e-16 ) 
        {
            network_data_sample.R = 1e16; 
        }
        else
        {
            network_data_sample.R = is_diag * 1/input_SER->d(i); // rownolegle R
        }

        // rownolegle C
        if ( input_SER->d(i) < 1e-16 )
        {
            network_data_sample.C = 1e-16;
        }
        else
        {
            network_data_sample.C = is_diag * input_SER->h(i);
        }

        for ( int j = 0; j < N ; j++ ) // po wszystkich residuach
        {

            if ( imag_check(nn, j) == 0 && imag_check(nn, j) != NO_ELEM ) //biegun rzeczywisty
            {
                 // obliczanie parametrow dla galezi od bieguna real
                 real_pole_net_sample = parse_real_pole( input_SER->res(i, j), poles(nn, j), is_diag );
                 network_data_sample.real_pole_nets.push_back(real_pole_net_sample);
            }
            else if ( imag_check(nn, j) == 1 && imag_check(nn, j) != NO_ELEM ) // biegun zespolony
            {
                 // obliczanie parametrow dla galezi od bieguna imag
                 imag_pole_net_sample = parse_imag_pole( input_SER->res(i, j), poles(nn, j), is_diag );
                 network_data_sample.imag_pole_nets.push_back(imag_pole_net_sample);
            }
        }
        output_network_data[i] = network_data_sample;
    }

}


real_pole_net parse_real_pole( cx_double res, cx_double poles, int is_diag )
{
    real_pole_net net;

    cout << "\nBiegun rzeczywisty" << endl;
    cout << "Pole " << poles << endl;
    cout << "Res " << res << endl;

    net.R = -real(poles)/real(res) * is_diag;
    net.L = 1/real(res) * is_diag;

    net.res = real(res);
    net.pole = real(poles);

    cout << "R " << net.R << endl;
    cout << "L " << net.L << endl;

    return net;
}


imag_pole_net parse_imag_pole( cx_double res, cx_double poles, int is_diag )
{
    imag_pole_net net;
    double res_real = real(res) * is_diag;
    double res_imag = imag(res) * is_diag;
    double poles_real = real(poles);
    double poles_imag = imag(poles);
    cout << "\nBiegun zespolony" << endl;
    cout << "Pole " << poles << endl;
    cout << "Res " << res << endl;

    net.L = 1/(2 * res_real);
    net.R = 2*net.L*(net.L*(res_real*poles_real + res_imag*poles_imag) - poles_real); 
    net.C = 1 / ( net.L*( pow(poles_real,2) + pow(poles_imag,2) + 2*net.R*(res_real*poles_real + res_imag*poles_imag)) ); 
    net.G = -2*( res_real*poles_real + res_imag*poles_imag) * net.C * net.L;

    net.res = cx_double(res_real, res_imag);
    net.pole = cx_double(poles_real, poles_imag);

    cout << "\nR " << net.R << endl;
    cout << "L " << net.L << endl;
    cout << "C " << net.C << endl;
    cout << "G " << net.G << endl;
    return net;
}


void print_network_data( Y_network_data *Y, int i )
{
     cout << "\n\nNetwork data for " << i << " Y element" << endl;

     cout << "R0: " << Y[i].R << "om" << endl;
     cout << "C0: " << Y[i].C << "F" << endl;

     cout << "Real pole branches: " << endl;
     for ( int j = 0 ; j < Y[i].real_pole_nets.size() ; j++ )
     {
         cout << "R: " << Y[i].real_pole_nets[j].R << " om" <<endl; 
         cout << "L: " << Y[i].real_pole_nets[j].L << " henr" <<endl; 
     }
     
     cout << "Imag pole branches:" <<endl;
     for ( int j = 0 ; j < Y[i].imag_pole_nets.size() ; j++ )
     {
         cout << "R: " << Y[i].imag_pole_nets[j].R << " om" <<endl; 
         cout << "L: " << Y[i].imag_pole_nets[j].L << " henr" <<endl; 
         cout << "C: " << Y[i].imag_pole_nets[j].C << " F" <<endl; 
         cout << "G: " << Y[i].imag_pole_nets[j].G << " Si" <<endl; 
     }
}


void create_model_netlist( SER *input_SER, int Nc, const vec& freq, ofstream &cir_file, vf_opts& conf )
{
    int Nc_port = sqrt(Nc);

    if ( pow(Nc_port, 2) != Nc )
    {
        cout << "Podana macierz nie jest kwadratowa - model nie bedzie syntezowany" <<endl;
        return;
    }

    // macierz ze strukturami z wartosciami elementow modelu
    Y_network_data *data;
    data = new Y_network_data[Nc];

    // przygotowanie wartosci elementow
    parse_SER( input_SER, data );

    // wyswietlenie wartosci elementow
    for ( int i = 0; i < Nc ; i++ )
    {
        print_network_data( data, i );
    }

    //pierwsza linika cira
    cir_file << "Generated netlist: " << conf.out_file_name << "\n" << endl;

    if ( conf.ngspice_simulation )
    {
        cir_file << "*Initial values for AC ports\n";
        for ( int i = 1 ; i <= Nc_port; i++ )
        {
            cir_file << ".param Vg" << i << "=0\n";
        }
    }

    // przygotowanie subckt dla kazdego elementu macierzy Y
    for ( int j = 1; j <= Nc_port; j++ )
    {
        for ( int i = 1; i <= Nc_port; i++ )
        {
            Y_network_data Y_temp = get_Y( data, i, j, Nc_port);

            // wygenerowanie odpowiedniego indeksu elementu macierzy Y
            ostringstream ss;
            ss << i << j;
            string y_inx = ss.str();
            create_subckt( Y_temp, y_inx, cir_file, conf );
        }
    }       
    
    cir_file << endl;
    cir_file << "*** Complete cir ***" <<endl;

    int current_index = 1; // indeks zrodel pradowych sterowanych pradem 
    int voltage_index = 1; // indeks zrodel napieciowych sterowanych napieciem

    int node_volt = 0; // nr wezla od zrodla napieciowego sterowanego napieciem

    // zlozenie calego modelu obwodoweg do cira (subckt i zrodla sterowane)
    for ( int i = 0; i < Nc_port; i++ )
    {
        node_volt = (i+1)*1000;

        for ( int j = 0; j < Nc_port; j++ )
        {
            if ( j == i ) // element z przekatnej Y
            {
                int node = i+1;
                cir_file << "V" << node << " " << node << " 0 DC 0 AC {Vg" << node << "}" << endl; // port
                cir_file << "X_Y" << node << node << " " << node << " 0 " << "Y" << node << node << endl; // element z przekatnej macierzy Y 
              
                // dodanie zrodel pradowych sterowanych pradem
                for ( int k = 1; k <= Nc_port; k++ )
                {
                    if ( k != node )
                    {
                        cir_file << "F" << current_index << " " << node << " 0 " << "E" << k << "_" << node << " 1" << endl;
                        current_index++;
                    }
                }
            }
            else // element spoza przekatnej
            {
                int e_idx = node_volt / 1000;
                cir_file << "X_Y" << i+1 << j+1 << " " << node_volt << " 0 " << "Y" << i+1 << j+1 << endl;
                // zrodlo napieciowe sterowane napieciem
                cir_file << "E" << j+1 << "_" << e_idx << " 0 " << node_volt << " " << j+1 << " 0 1" <<endl;

                node_volt++;
            }
        }
    }

    // przygotowanie komendy symulacji
    double freq_start = freq(0);
    double freq_end = freq( freq.n_elem - 1 );
    cir_file << "\n.ac lin " << freq.n_elem << " " << freq_start << " " << freq_end << endl;

    // przygotowanie .control dla ngspice
    if ( conf.ngspice_simulation )
    {
        cir_file << "\n.control\nset filetype=ascii\n\nlet start=1\nlet step=1\n"
                 << "let end=" << sqrt(Nc) << "\nlet iter=start\n\nwhile iter le end\n\n";

        // ustawienie wszystkich wrot AC na 0
        for ( int i = 1; i <= sqrt(Nc); i++ )
        {
            cir_file << "    alter V" << i << " AC 0" << endl;
        }

        cir_file << "\n";

        // ustawienie AC -1 na odpowiednim wrocie
        // ujemna wartosc po to aby uzyskane prady na zrodlach byly parametrami Y
        for ( int i = 1; i <= sqrt(Nc); i++ )
        {
            cir_file << "    if iter = " << i 
                     << "\n        alter V" << i << " AC -1\n    end\n";
        }

        cir_file << "\n    run\n"
		 << "\n    wrdata " << conf.out_file_name;
        for ( int i = 1; i <= sqrt(Nc); i++ )
        {
            cir_file << " I(V" << i << ")";
        }
        cir_file << "\n    set appendwrite"
                 << "\n    let iter = iter + step"
                 << "\nend";
        cir_file << "\n\n.endc" << endl;
    }
    else
    {
        // dane do symulacji lt spice
        cir_file << "\n.step param run 1 " << Nc_port << " 1" << endl;
        for ( int i = 1; i <= Nc_port; i++ )
        {
            cir_file << ".param Vg" << i << "=if(run==" << i << ",-1,0)" << endl;
        }
    }

    cir_file << "\n.end";

    delete[] data;
}


void create_subckt( Y_network_data data, string index, ofstream &cir_file, vf_opts& conf )
{
     int node = 3;
     int R_index = 1;
     int L_index = 0;
     int C_index = 1;

     cir_file << endl;
     cir_file << "*** Subcircuit for Y" << index << endl;
     cir_file << ".subckt Y" << index << " 1 2" <<endl;

     if ( data.R > conf.R_max )
     {
         cout << "Rownolegle R dla Y" << index << " pominiete (zgodnie z konfiguracja)" << endl; 
     }
     else
     {
         cir_file << "R0 1 2 " << data.R << endl;
     }
  
     if ( data.C < conf.C_min )
     {
         cout << "Rownolegle C dla Y" << index << " pominiete (zgodnie z konfiguracja)" << endl; 
     }
     else
     {
         cir_file << "C0 1 2 " << data.C << endl;
     }

     // real poles
     for ( int i = 0; i < data.real_pole_nets.size(); i++ )
     {
         cir_file << "*** Real pole: (res=" << data.real_pole_nets[i].res 
                  << " pole=" << data.real_pole_nets[i].pole << ") ***" << endl;
         cir_file << "R" << R_index << " 1 " << node << " " << data.real_pole_nets[i].R << endl; 
         cir_file << "L" << L_index << " " << node << " 2 " << data.real_pole_nets[i].L << endl; 
         R_index++;
         L_index++;
         node++;
     }

     // imag poles
     for ( int i = 0; i < data.imag_pole_nets.size(); i++ )
     {
         cir_file << "*** Imag pole: para (res=" << data.imag_pole_nets[i].res 
                  << " pole=" << data.imag_pole_nets[i].pole << ") ***" << endl;
         cir_file << "R" << R_index++ << " 1 " << node << " " << data.imag_pole_nets[i].R << endl; 
         cir_file << "L" << L_index << " " << node << " " << ++node << " " << data.imag_pole_nets[i].L << endl; 
         cir_file << "C" << C_index << " " << node << " 2 " << data.imag_pole_nets[i].C << endl; 
         cir_file << "R" << R_index << " " << node << " 2 " << 1 / data.imag_pole_nets[i].G << endl; // G modelu przedstawione jako R 

         R_index++;
         L_index++;
         C_index++;
         node++;
     }
     cir_file << ".ends" << endl;
}


Y_network_data get_Y( Y_network_data *input, int i, int j, int Nc_port )
{
    i=i-1;
    j=j-1;
    return input[i+j*Nc_port];
}







