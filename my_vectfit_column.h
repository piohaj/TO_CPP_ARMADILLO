#ifndef my_vectfit_column_h
#define my_vectfit_column_h

#include "data_model.h"
#include "additional_features.h"

SER my_vf_column_splitting(const cx_mat *f, const cx_vec *s, cx_mat *poles);
void rms_err_calculation(SER *wynik, const cx_mat *f, const cx_vec *s, int N);

class vf_column
{
    const cx_mat *f;
    const cx_vec *s;
    cx_mat *poles1;
    SER *wynik;

public:
    vf_column( const cx_mat *f_in, const cx_vec *s_in, cx_mat *poles_in,
               SER *wynik_in )
         : f(f_in),
           s(s_in),
           poles1(poles_in),
           wynik(wynik_in)
         {}

    void operator() ( const blocked_range<int>& r ) const;
};

#endif
