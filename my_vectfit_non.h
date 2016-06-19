#ifndef my_vectfit_non_h
#define my_vectfit_non_h

#include "data_model.h"
#include "additional_features.h"

// Funckja my_vectorfit
// Funkcja znajduje model ukladu w postaci funkcji wymiernej
// na podstawie danych pomiarowych.
// Dane wejsciowe:
// f - macierz zawierajaca wartosci pomiarow dla kazdego elementu macierzy Y
// s - wektor zawierajacy probki wartosci czestotliwosci dla wartosci pomiaru z f
// poles - wektor z biegunami poczatkowymi
SER my_vf_non_splitting(const cx_mat& f, const cx_vec& s, cx_mat poles);


cx_mat logspace(double a, double b, int n);
int sign( double x );

//funkcje do generowania lub ladowaniai z pliku danych testowych
input_data prepare_sample_data();
input_data load_data_from_file( int N, int Nc, int Ns );

input_data load_vf_data( string file_name );
// klasa potrzebna do zrownoleglenia dekompozycji QR w TBB
class QR_calculation
{
    const cx_mat *A;
    const cx_mat *f;
    int N;
    int Ns;
    mat *AA_poles;
    mat *bb_poles;
    
public:
    QR_calculation( const cx_mat *A_in, const cx_mat *f_in, int N_in, int Ns_in,
                    mat *AA_poles_out, mat *bb_poles_out)
                  : A(A_in),
                    f(f_in),
                    N(N_in),
                    Ns(Ns_in),
                    AA_poles(AA_poles_out),
                    bb_poles(bb_poles_out)
                  {}

    // zadanie dla jednego watku - w TBB definiowany poprzez przeciazenie ()
    void operator() ( const blocked_range<int>& r ) const;
};

#endif
