#ifndef my_vectfit_all_h
#define my_vectfit_all_h

#include "data_model.h"
#include "additional_features.h"

SER my_vf_all_splitting(const cx_mat *f, const cx_vec *s, cx_mat *poles, vf_opts& conf);

class vf_all
{
    const cx_mat *f;
    const cx_vec *s;
    cx_mat *poles1;
    SER *wynik;
    int RC_offset;

public:
    vf_all( const cx_mat *f_in, const cx_vec *s_in, cx_mat *poles_in,
            SER *wynik_in, int RC_offset_in )
         : f(f_in),
           s(s_in),
           poles1(poles_in),
           wynik(wynik_in),
           RC_offset(RC_offset_in)
         {}

    void operator() ( const blocked_range<int>& r ) const;
};

#endif
