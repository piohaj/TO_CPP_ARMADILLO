#include "my_vectfit.h"

// obliczanie wartosci modelu i wstawienie go do 
void parse_SER(SER *input_SER, Y_network_data *output_network_data);
real_pole_net parse_real_pole( cx_double res, cx_double poles );
imag_pole_net parse_imag_pole( cx_double res, cx_double poles );

void print_network_data( Y_network_data *Y, int i);

void create_subckt( Y_network_data data, int index );
void create_cir( Y_network_data *data, int N, int Nc);
