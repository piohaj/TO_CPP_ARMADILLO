#ifndef additional_features_h
#define additional_features_h

#include "data_model.h"
#include "my_vectfit_non.h"
#include "my_vectfit_all.h"
#include "my_vectfit_column.h"

SER vf_high_level( const cx_mat& f, const cx_vec& s, int split_strat, int min_row, int max_row, int max_iters = 10 );
cx_mat prepare_input_poles( const cx_vec& s, int split_strat, int Nc, int N, int Ns );
mat choose_best_aprox( SER *input, int size, int Nc, int split_strat );
SER cumulate_model( int split_strat, mat& indexes, SER *iter_models, int Nc, int max_row );

void read_conf();

cx_cube make_cube( cx_mat& y );
cx_cube y2s(cx_mat y, int z0 = 50);
cx_mat s2y(cx_mat s_params, int z0 =50);
bool ispassive(cx_mat& y);

#endif
