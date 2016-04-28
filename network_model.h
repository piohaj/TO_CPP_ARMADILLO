#include "my_vectfit.h"

// obliczanie wartosci modelu i wstawienie go do 
void parse_SER(SER *input_SER, Y_network_data *output_network_data);
real_pole_net parse_real_pole( cx_mat res, cx_mat poles );
imag_pole_net parse_imag_pole( cx_mat res, cx_mat poles );

