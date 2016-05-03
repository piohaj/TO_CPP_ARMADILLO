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

    cout << Nc << endl;
    cout << N << endl;

    for ( int i = 0; i < Nc_port; i++ )
    {
        int index_help=1;
        for ( int j = i*Nc_port; j < i*Nc_port+Nc_port; j++ )
        {
            cout << "Y" << index_help*10 + i + 1 <<endl;
            index_help++;
        }
    }
}
