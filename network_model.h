#ifndef network_model_h
#define network_model_h

#include "data_model.h"
// obliczanie wartosci modelu i wstawienie go do 
void parse_SER(SER *input_SER, Y_network_data *output_network_data );
real_pole_net parse_real_pole( cx_double res, cx_double poles, int is_diag );
imag_pole_net parse_imag_pole( cx_double res, cx_double poles, int is_diag );

void print_network_data( Y_network_data *Y, int i);

void create_subckt( Y_network_data data, string index, ofstream &cir_file, vf_opts& conf );

// create_model_netlist - funkcja do tworzenia pliku cir do symulacji w LTspice na podstawie danych uzyskanych z algorytmu vector fitting
// input_SER - dane uzyskane z vector fitting
// Nc - liczba elementow macierzy Y 
void create_model_netlist( SER *input_SER, int Nc, ofstream &cir_file, vf_opts& conf );

// funkcja do wyciagania danych z jednowymiarowej tablicy z danymi do cira
// wyciagniecie Y11 to get_Y( data, 1, 1, 2 ); w przypadku dwuwrotnika
Y_network_data get_Y( Y_network_data *input, int i, int j, int Nc_port );

#endif
