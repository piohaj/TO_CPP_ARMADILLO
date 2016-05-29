#ifndef additional_features_h
#define additional_features_h

#include "data_model.h"

SER vf_high_level( const cx_mat& f, const cx_vec& s, int split_strat, int min_row, int max_row, int max_iters = 10 );
cx_mat prepare_input_poles( const cx_vec& s, int split_strat, int Nc, int N, int Ns );
int choose_best_aprox( SER *input, int size );

#endif
