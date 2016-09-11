#ifndef additional_features_h
#define additional_features_h

#include "data_model.h"
#include "my_vectfit_non.h"
#include "my_vectfit_all.h"
#include "my_vectfit_column.h"

SER vf_high_level( cx_mat& f, const cx_vec& s, vf_opts conf );
cx_mat prepare_input_poles( const cx_vec& s, int split_strat, int Nc, int N, int Ns );
mat choose_best_aprox( SER *input, int size, int Nc, int split_strat, double rms_diff );
SER cumulate_model( int split_strat, mat& indexes, SER *iter_models, int Nc, int max_row );


cx_cube make_cube( cx_mat& y );
cx_cube y2s(cx_mat y, int z0 = 50);
cx_mat s2y(cx_cube s_params, int z0 = 50);
cx_mat cube2mat( cx_cube& cube );
bool ispassive_slice(cx_mat& s_slice);
bool ispassive_s( cx_cube& s_params );

int makepassive( cx_cube& s_params );

int check_and_make_passive ( cx_mat& y );

int check_model_simulation_results( const cx_mat& f, const vf_opts& conf, gnuplot_data & data );

void force_stable_poles( cx_mat& poles ); 
mat gp_angle( cx_mat & data );

// wczytywanie konfiguracji
vector<string> my_split(string str, const char delim);
map <string, string>read_conf_file(string file_name);
string read_param( string param, map <string,string> &conf_map);
int read_conf( vf_opts& global_conf, string file_name );

// wczytywanie plikow
int recognize_file ( string file_name );

// wczytywanie plikow raw
raw_params read_input_data_params( string file_name );
string create_temp_bin_file( int last_no_bin_data_line, string file_name );
input_data read_raw_file( string file_name );
input_data parse_raw_data( const cx_mat& temp_data, raw_params params );

// wczytywanie touchstone
string str_toupper( string & str );
int get_ports_num_touchstone( string file_name );
int count_spaces_in_header( string header );
touchstone_conf check_header_touchstone( string file_name );
void read_touchstone( string file_name, input_data & data );
cx_mat angle2canonic( const mat& mag, const mat& angle);
void db2magnitude( mat& db );

// funkcje spice
bool check_spice_log( string file_name );

// zapis wynikow do plikow
void save_results_mats( SER & results, string file_name );

// wykresy gnuplot
void prepare_gnuplot_script( gnuplot_data & data, string name );

// odwracalnosc ukladu
void reciprocal_make_mat( cx_mat &f );
void reciprocal_fix_results( SER &wynik, cx_mat &f, int Nc_ports );
cx_mat get_Y_elem(cx_mat &f, int i, int j, int Nc_ports);
#endif
