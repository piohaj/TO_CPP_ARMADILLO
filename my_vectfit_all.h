#include "data_model.h"

SER my_vf_all_splitting(const cx_mat *f, const cx_vec *s, cx_mat *poles);

class vf_all
{
    const cx_mat *f;
    const cx_vec *s;
    cx_mat *poles1;
    SER *wynik;

public:
    vf_all( const cx_mat *f_in, const cx_vec *s_in, cx_mat *poles_in,
            SER *wynik_in )
         : f(f_in),
           s(s_in),
           poles1(poles_in),
           wynik(wynik_in)
         {}

    void operator() ( const blocked_range<int>& r ) const;
};
