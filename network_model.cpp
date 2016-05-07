#include "network_model.h"


void parse_SER(SER *input_SER, Y_network_data *output_network_data)
{
    int N = input_SER->poles.n_rows;
    int Nc = input_SER->res.n_rows;
    cx_mat poles = input_SER->poles;
    //output_network_data = new Y_network_data[Nc];

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


    for ( int i = 0 ; i < Nc ; i++ )
    {
        Y_network_data network_data_sample; // struktura z danymi sieci dla pojedynczego elementu macierz Y
        // wyczyszczenie struktury pomocniczej
        network_data_sample.real_pole_nets.clear();
        network_data_sample.imag_pole_nets.clear();
        network_data_sample.R = 0;

        real_pole_net real_pole_net_sample;
        imag_pole_net imag_pole_net_sample;

        network_data_sample.R = 1/input_SER->h(i); // rownolegle R

        for ( int j = 0; j < N ; j++ ) // po wszystkich residuach
        {
            if ( imag_check(j) == 0 ) //biegun rzeczywisty
            {
                 // obliczanie parametrow dla galezi od bieguna real
                 real_pole_net_sample = parse_real_pole( input_SER->res(i, j), poles(j) );
                 network_data_sample.real_pole_nets.push_back(real_pole_net_sample);
            }
            else if ( imag_check(j) == 1 ) // biegun zespolony
            {
                 // obliczanie parametrow dla galezi od bieguna imag
                 imag_pole_net_sample = parse_imag_pole( input_SER->res(i, j), poles(j) );
                 network_data_sample.imag_pole_nets.push_back(imag_pole_net_sample);
            }
        }
        output_network_data[i] = network_data_sample;
    }

}


real_pole_net parse_real_pole( cx_double res, cx_double poles )
{
    real_pole_net net;

    cout << "\nBiegun rzeczywisty" << endl;
    cout << "Pole " << poles << endl;
    cout << "Res " << res << endl;

    net.R = -real(poles)/real(res);
    net.L = 1/real(res);

    cout << "R " << net.R << endl;
    cout << "L " << net.L << endl;

    return net;
}


imag_pole_net parse_imag_pole( cx_double res, cx_double poles )
{
    imag_pole_net net;
    double res_real = real(res);
    double res_imag = imag(res);
    double poles_real = real(poles);
    double poles_imag = imag(poles);
    cout << "\nBiegun zespolony" << endl;
    cout << "Pole " << poles << endl;
    cout << "Res " << res << endl;

    net.L = 1/(2 * res_real);
    net.R = 2*net.L*(net.L*(res_real*poles_real + res_imag*poles_imag) - poles_real); 
    net.C = 1 / ( net.L*( pow(poles_real,2) + pow(poles_imag,2) + 2*net.R*(res_real*poles_real + res_imag*poles_imag)) ); 
    net.G = -2*( res_real*poles_real + res_imag*poles_imag) * net.C * net.L;

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


void create_cir( Y_network_data *data, int N, int Nc)
{
    int Nc_port = sqrt(Nc);

    if ( pow(Nc_port, 2) != Nc )
    {
        cout << "Podana macierz nie jest kwadratowa - model nie bedzie syntezowany" <<endl;
        return;
    }

    //pierwsza linika cira
    cout << "Generated model" << endl;

    // przygotowanie subckt dla kazdego elementu macierzy Y
    for ( int j = 1; j <= Nc_port; j++ )
    {
        for ( int i = 1; i <= Nc_port; i++ )
        {
            Y_network_data Y_temp = get_Y( data, i, j, Nc_port);
            int y_inx = i*10+j;
            create_subckt( Y_temp, y_inx );
        }
    }       
    
    cout << endl;
    cout << "*** Complete cir ***" <<endl;

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
                cout << "V" << node << " " << node << " 0 AC 1" << endl; // port
                cout << "X_Y" << node << node << " " << node << " 0 " << "Y" << node << node << endl; // element z przekatnej macierzy Y 
              
                // dodanie zrodel pradowych sterowanych pradem
                for ( int k = 1; k <= Nc_port; k++ )
                {
                    if ( k != node )
                    {
                        cout << "F" << current_index << " " << node << " 0 " << "E" << k << " 1" << endl;
                        current_index++;
                    }
                }
            }
            else // element spoza przekatnej
            {
                cout << "X_Y" << i+1 << j+1 << " " << node_volt << " 0 " << "Y" << i+1 << j+1 << endl;
                // zrodlo napieciowe sterowane napieciem
                cout << "E" << j+1 << " 0 " << node_volt << " " << j+1 << " 0 1" <<endl;

                node_volt++;
            }
        }
    }
    cout << ".end";
}


void create_subckt( Y_network_data data, int index )
{
     int node = 3;
     int R_index = 1;
     int L_index = 0;
     int C_index = 1;

     cout << endl;
     cout << "*** Subcircuit for Y" << index << endl;
     cout << ".subckt Y" << index << " 1 2" <<endl;
     cout << "R0 1 2 " << data.R << endl;

     // real poles
     for ( int i = 0; i < data.real_pole_nets.size(); i++ )
     {
         cout << "*** Real pole ***" <<endl;
         cout << "R" << R_index << " 1 " << node << " " << data.real_pole_nets[i].R << endl; 
         cout << "L" << L_index << " " << node << " 2 " << data.real_pole_nets[i].L << endl; 
         R_index++;
         L_index++;
         node++;
     }

     // imag poles
     for ( int i = 0; i < data.imag_pole_nets.size(); i++ )
     {
         cout << "*** Imag pole ***" <<endl;
         cout << "R" << R_index++ << " 1 " << node << " " << data.imag_pole_nets[i].R << endl; 
         cout << "L" << L_index << " " << node << " " << ++node << " " << data.imag_pole_nets[i].L << endl; 
         cout << "C" << C_index << " " << node << " 2 " << data.imag_pole_nets[i].C << endl; 
         cout << "R" << R_index << " " << node << " 2 " << 1 / data.imag_pole_nets[i].G << endl; // G modelu przedstawione jako R 

         R_index++;
         L_index++;
         C_index++;
         node++;
     }
     cout << ".ends" << endl;
}


Y_network_data get_Y( Y_network_data *input, int i, int j, int Nc_port )
{
    i=i-1;
    j=j-1;
    return input[i+j*Nc_port];
}







