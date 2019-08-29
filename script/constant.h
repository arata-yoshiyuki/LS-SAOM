#ifndef CONSTANT_H_
#define CONSTANT_H_

#include <iostream>
#include <random>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <tuple>
#include <ctime>
#include <unordered_map>
#include <set>

using namespace std;
using mat = vector<vector<double>>;
using net = unordered_map<int, unordered_map<int, int>>;
using cov = unordered_map<int, unordered_map<int, double>>;
using h_set = unordered_map<int, set<int>>;

extern const int n_proc;
extern const int net_size;
extern const int n_item;
extern const int pp;

extern const int n_sub;
extern const int n_phase_1;
extern const int n_phase_2;
extern const int n_phase_3;
extern const double eps;
extern const double a_gain;

extern const int kj_ind;
extern const int log_age;
extern const int log_sale;
extern const int d_lab_pro;
extern const int g_rate;
extern const int pf_rate;
extern const int loc_X;
extern const int loc_Y;

extern const string f_name_kj, f_name_sk_1, f_name_sk_2;

extern vector<string> theta_name;
extern vector<double> pre_theta;
extern const int sim_flag, ph_1_flag, ph_2_flag, ph_3_flag, gof_flag, gof_2_flag;
extern vector<double> sim_theta_A, sim_theta_B, sim_theta_C;

extern const mat pre_D_hat_temp;

#endif /* CONSTANT_H_ */
