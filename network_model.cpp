#include "network_model.h"


void parse_SER(SER *input_SER, Y_network_data *output_network_data)
{
    int N = input_SER->poles.n_rows;
    int Nc = input_SER->res.n_rows;

    cout << "Ilosc elementow Y przyblizenia "<< Nc << endl;
}


