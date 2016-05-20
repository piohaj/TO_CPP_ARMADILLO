#include "my_vectfit.h"

// obliczanie wartosci modelu i wstawienie go do 
void parse_SER(SER *input_SER, Y_network_data *output_network_data);
real_pole_net parse_real_pole( cx_double res, cx_double poles, int is_diag );
imag_pole_net parse_imag_pole( cx_double res, cx_double poles, int is_diag );

void print_network_data( Y_network_data *Y, int i);

void create_subckt( Y_network_data data, string index );
void create_cir( Y_network_data *data, int N, int Nc);

// funkcja do wyciagania danych z jednowymiarowej tablicy z danymi do cira
// wyciagniecie Y11 to get_Y( data, 1, 1, 2 ); w przypadku dwuwrotnika
Y_network_data get_Y( Y_network_data *input, int i, int j, int Nc_port );
