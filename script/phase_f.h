#ifndef PHASE_F_H_
#define PHASE_F_H_

extern "C"{
  void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
};

random_device rd;

void display(const mat& sample, string name);
mat Inverse_matrix(const mat& A);
mat Product_matrix(mat A, mat B);
tuple <h_set, h_set> initialize_f(const net& net_y, const net& net_x, const net& net_x_b, const net& cov_sec);
void simulation(const vector<double>& theta, net& net_x, net& net_x_b,
  const cov& covar, const net& cov_sec,
  const h_set& f_subset_ind, const h_set& ind_subset_f,
  const unsigned& rand);
tuple <mat, vector<double>, mat, mat, double> phase3(
  const vector<double>& theta,
  const net& net_y, const net& net_y_b,
  const net& net_x, const net& net_x_b,
  const cov& covar, const net& cov_sec,
  const h_set& f_subset_ind, const h_set& ind_subset_f,
  // const net& two_dis_1_t1, const net& two_dis_1_t3,
  int n_phase, int my_rank, int p);
void phase2(vector<double>& theta,
  const net& net_y, const net& net_y_b,
  const net& net_x, const net& net_x_b,
  const cov& covar, const net& cov_sec,
  const h_set& f_subset_ind, const h_set& ind_subset_f,
  // const net& two_dis_1_t1, const net& two_dis_1_t3,
  const vector<vector<double> >& D_hat_temp,
  int my_rank, int p);

void display(const mat& sample, string name){
	int number = sample.size();
	cout << name << endl;
	for(int i = 0; i < number; ++i) {
		for(int j = 0; j < number; ++j) cout << sample[i][j] << ", ";
		cout << endl;
	}
}
mat Inverse_matrix(const mat& A){
  double temp_A[pp*pp];
  for(int j=0; j < pp; ++j){
    for(int i=0; i < pp; ++i){
      temp_A[pp*j + i] = A[i][j];
    }
  }
  int ipiv[pp];

	mat B(pp, vector<double>(pp, 0) );
	for(int i=0; i < pp; ++i) B[i][i] = 1.0;

  double temp_B[pp*pp];
  for(int j=0; j < pp; ++j){
    for(int i=0; i < pp; ++i){
      temp_B[pp*j + i] = B[i][j];
    }
  }
  int info;

  dgesv_(&pp,&pp,temp_A,&pp,ipiv,temp_B,&pp,&info);

	mat C(pp, vector<double>(pp, 0));
  for(int j=0; j < pp; ++j){
    for(int i=0; i < pp; ++i){
      C[i][j] = temp_B[pp*j + i];
    }
  }

	return C;
}
mat Product_matrix(mat A, mat B){
	int n, m, l;
	n = A.size();
	m = A[0].size();
	l = B[0].size();
	mat C(n, vector<double>(l, 0) );
    for(int i=0; i < n; ++i) {
            for(int j=0; j < l; ++j) {
                    for(int k=0; k < m; ++k) {
                            C[i][j] += A[i][k] * B[k][j];
                    }
            }
    }
    return C;
}
tuple <h_set, h_set> initialize_f(const net& net_y, const net& net_x, const net& net_x_b, const net& cov_sec){

  h_set f_subset_ind;
  for(int i = 1; i <= net_size; ++i){
    f_subset_ind[i];
    for(auto itr = net_x.at(i).begin(); itr != net_x.at(i).end(); ++itr) {
      int sup_i_ind = cov_sec.at(itr->first).at(kj_ind);
      f_subset_ind.at(i).insert(sup_i_ind);
    }
    for(auto itr = net_y.at(i).begin(); itr != net_y.at(i).end(); ++itr){
      int sup_i_ind = cov_sec.at(itr->first).at(kj_ind);
      f_subset_ind.at(i).insert(sup_i_ind);
    }
  }

  h_set ind_subset_f;
  for(int i = 1; i <= net_size; ++i){
    int i_ind = cov_sec.at(i).at(kj_ind);
    ind_subset_f[i_ind].insert(i);
  }

  /*
  net two_dis_t1, two_dis_t3;

  // indirect network type 1
  for(int i = 1; i <= net_size; ++i){
    two_dis_t1[i];
    for(auto itr_i = net_x.at(i).begin(); itr_i != net_x.at(i).end(); ++itr_i){
      int temp_j = itr_i->first;
      for(auto itr_temp_j = net_x.at(temp_j).begin(); itr_temp_j != net_x.at(temp_j).end(); ++itr_temp_j) {
        int indir_k = itr_temp_j->first;
        two_dis_t1.at(i)[indir_k] = 1;
      }
    }
  }

  // indirect network type 3
  for(int i = 1; i <= net_size; ++i){
    two_dis_t3[i];
    for(auto itr_i = net_x_b.at(i).begin(); itr_i != net_x_b.at(i).end(); ++itr_i){
      int temp_j = itr_i->first;
      for(auto itr_temp_j = net_x_b.at(temp_j).begin(); itr_temp_j != net_x_b.at(temp_j).end(); ++itr_temp_j) {
        int indir_k = itr_temp_j->first;
        two_dis_t3.at(i)[indir_k] = 1;
      }
    }
  }

  */
  // return make_tuple(f_subset_ind, ind_subset_f, two_dis_t1, two_dis_t3);
  return make_tuple(f_subset_ind, ind_subset_f);
}
void simulation(const vector<double>& theta,
  net& net_x, net& net_x_b,
  const cov& covar, const net& cov_sec,
  const h_set& f_subset_ind, const h_set& ind_subset_f,
  const unsigned& rand){

    double rate_p = theta[0];
    double rate_outdeg_p = theta[1];
    double rate_age_p = theta[2];

    double outde_p = theta[3];
    double ego_sale_p = theta[4];

    double alt_lab_pro_p = theta[5];
    double alt_grow_p = theta[6];
    double alt_pf_p = theta[7];
    double alt_sale_p = theta[8];

    double inPop_p = theta[9];

    double sim_size_p = theta[10];
    double sim_tech_p = theta[11];
    double geo_dis_p = theta[12];

    double recipro_p = theta[13];
    double two_dis_t1_p = theta[14];
    double two_dis_t3_p = theta[15];

    if(rate_p < 0.0001) cout << "negative value of rate_para " << endl;


    mt19937 engine_2(rand);
    double time_var = 0.0;
    double time_1 = 1.0;

    // initial rate weight
    vector<double> rate_w((net_size + 1), 1.0); // Note that firm ID starts from 1
    rate_w[0] = 0.0; //there is no probability that "0" is chosen
    for(int i = 1; i <= net_size; ++i) rate_w[i] = rate_p*exp(rate_outdeg_p*log(net_x.at(i).size() + 1) + rate_age_p*covar.at(i).at(log_age)); //to avoid log(0), 1 is added.

    // the following parts of weight does not change over time.
    unordered_map<int, double> fixed_var;
    for(int sup_i = 1; sup_i <= net_size; ++sup_i) fixed_var[sup_i] = outde_p + alt_lab_pro_p*covar.at(sup_i).at(d_lab_pro) + alt_grow_p*covar.at(sup_i).at(g_rate) + alt_pf_p*covar.at(sup_i).at(pf_rate) + alt_sale_p*covar.at(sup_i).at(log_sale);

    while(time_var < time_1){
      double rate_sum = accumulate(rate_w.begin(), rate_w.end(), 0.0); // Note 0.0 rather than 0.

      exponential_distribution<> dist_net(rate_sum);
  		time_var += dist_net(engine_2);
      discrete_distribution<int> dist_2(rate_w.begin(), rate_w.end());
      int actor_id = dist_2(engine_2);

      //Alternative
      if(f_subset_ind.at(actor_id).size() != 0){
        vector<int> vec_ind;//the industry subset of acotr_id
        for(auto itr = f_subset_ind.at(actor_id).begin(); itr != f_subset_ind.at(actor_id).end(); ++itr) vec_ind.push_back(*itr);
        vector<double> eq_w(vec_ind.size(), 1.0);
        discrete_distribution<> vec_ind_dist(eq_w.begin(), eq_w.end());//equal weight
        int c_ind = vec_ind[vec_ind_dist(engine_2)];

        vector<int> vec_f;//the industry subset of acotr_id
        for(auto itr = ind_subset_f.at(c_ind).begin(); itr != ind_subset_f.at(c_ind).end(); ++itr) vec_f.push_back(*itr);

        vector<double> u_weight;
        for(auto itr = vec_f.begin(); itr != vec_f.end(); ++itr){
          int sup_id = *itr;

          double tie_gen = 0.0;
          if(net_x.at(actor_id).count(sup_id) == 0) tie_gen = 1.0;
          double d_tie_gen = 2*tie_gen - 1; // 1 or -1

          double size_simil = abs(covar.at(actor_id).at(log_sale) - covar.at(sup_id).at(log_sale));
          double tech_simil = 0.0;
          if(cov_sec.at(actor_id).at(kj_ind) == cov_sec.at(sup_id).at(kj_ind)) tech_simil = 1.0;
          double dis_X = covar.at(sup_id).at(loc_X) - covar.at(actor_id).at(loc_X);
          double dis_Y = covar.at(sup_id).at(loc_Y) - covar.at(actor_id).at(loc_Y);

          double recipro = 0.0;
          if(net_x.at(sup_id).count(actor_id) != 0) recipro = 1.0;

          double two_dis_var1 = 0.0;
          for(auto itr = net_x.at(actor_id).begin(); itr != net_x.at(actor_id).end(); ++itr){
    				int temp_id = itr->first;
    				if(net_x.at(temp_id).count(sup_id)){
    					two_dis_var1 = 1.0;
    					break;
    				}
    			}
          double two_dis_var3 = 0.0;
          for(auto itr = net_x_b.at(actor_id).begin(); itr != net_x_b.at(actor_id).end(); ++itr){
            int temp_id = itr->first;
            if(net_x_b.at(temp_id).count(sup_id)){
              two_dis_var3 = 1.0;
              break;
            }
          }

          if(actor_id == sup_id){
            u_weight.push_back(0.0);
          }else{
            u_weight.push_back(exp(d_tie_gen*(fixed_var.at(sup_id)
              + ego_sale_p*covar.at(actor_id).at(log_sale)
              + inPop_p*log(net_x_b.at(sup_id).size() + 1)
              + sim_size_p*size_simil + sim_tech_p*tech_simil + geo_dis_p*sqrt(pow(dis_X, 2) + pow(dis_Y, 2))
              + recipro_p*recipro
              )
              + tie_gen*(two_dis_t1_p*two_dis_var1 + two_dis_t3_p*two_dis_var3))
            );
          }

        }
  			u_weight.push_back(1.0);
        int last_index = u_weight.size() - 1;//last_index means that actor does nothing

  			//#choice
        discrete_distribution<int> dist_3(u_weight.begin(), u_weight.end());
        int alter_index = dist_3(engine_2);
        int alter_id;
        if(alter_index == last_index ){
          alter_id = actor_id;// no change
        }else{
          alter_id = vec_f[alter_index];
        }

  		  if(actor_id != alter_id){
          if(net_x.at(actor_id).count(alter_id) == 0){ // a new tie is generated
            net_x.at(actor_id)[alter_id] = 1;
            net_x_b.at(alter_id)[actor_id] = 1;
            // cout << "Creation " << endl;

            // for actor_id
            rate_w.at(actor_id) = exp(rate_outdeg_p*(log(net_x.at(actor_id).size() + 1) - log(net_x.at(actor_id).size())))*rate_w.at(actor_id);
          }else{
            // cout << "DELETION " << endl;
            net_x.at(actor_id).erase(alter_id);
  				  net_x_b.at(alter_id).erase(actor_id);

            // for actor_id
            rate_w.at(actor_id) = exp(rate_outdeg_p*(log(net_x.at(actor_id).size() + 1) - log(net_x.at(actor_id).size() + 2) ) )*rate_w.at(actor_id);
          }
  		  }
      }
  	}
  }
tuple <mat, vector<double>, mat, mat, double> phase3(
  const vector<double>& theta,
  const net& net_y, const net& net_y_b,
  const net& net_x, const net& net_x_b,
  const cov& covar, const net& cov_sec,
  const h_set& f_subset_ind, const h_set& ind_subset_f,
  // const net& two_dis_1_t1, const net& two_dis_1_t3,
  int n_phase, int my_rank, int p){

  	unsigned rand_seed[n_phase];
  	for(int seed_i = 0; seed_i < n_phase; ++seed_i) rand_seed[seed_i]= rd();

  	mat Sum_Sq(pp, vector<double>(pp, 0));
    mat Sum_d(pp, vector<double>(pp, 0));
  	vector<double> Sum_S(pp, 0);
  	mat S_sim(pp, vector<double>(pp, 0));

    int ib = n_phase/n_proc;
    int istat, iend;

    istat = my_rank*ib;
    iend = (my_rank + 1)*ib;

  	for(int i = istat; i < iend; ++i){
  		unsigned rand = rand_seed[i];
      vector<double> S_org(pp, 0);

      #pragma omp parallel for schedule(dynamic, 1)
      for(int k = 0; k <= pp; ++k){
        net sim_net(net_x);
        net sim_net_b(net_x_b);
        // net sim_two_dis_t1(two_dis_1_t1);
        // net sim_two_dis_t3(two_dis_1_t3);

        if(k != pp){
          vector<double> theta_sim(theta);
          theta_sim[k] += eps;
          simulation(theta_sim, sim_net, sim_net_b, covar, cov_sec, f_subset_ind, ind_subset_f, rand);
          #pragma omp critical
          S_sim[k] = stat_gen(sim_net, sim_net_b, net_x, net_x_b, covar, cov_sec);
        }else{
          simulation(theta, sim_net, sim_net_b, covar, cov_sec, f_subset_ind, ind_subset_f, rand);
          #pragma omp critical
          S_org = stat_gen(sim_net, sim_net_b, net_x, net_x_b, covar, cov_sec);
        }
      }

  		mat d(pp, vector<double>(pp, 0));
  		for(int k = 0; k < pp; ++k){
  			for(int j = 0; j < pp; ++j) d[j][k] = (S_sim[k][j] - S_org[j])/eps;
  		}

  		for(int j = 0; j < pp; ++j) Sum_S[j] += S_org[j];
  		for(int j = 0; j < pp; ++j){
  			for(int k = 0; k < pp; ++k)	Sum_d[j][k] += d[j][k];
  		}

  		mat S_org_col(pp, vector<double>(1, 0));
  		for(int j = 0; j < pp; ++j) S_org_col[j][0] = S_org[j];
  		mat S_org_row(1, vector<double>(pp, 0));
  		for(int j = 0; j < pp; ++j) S_org_row[0][j] = S_org[j];
  		mat Sq = Product_matrix(S_org_col, S_org_row);
  		for(int j = 0; j < pp; ++j){
  			for(int k = 0; k < pp; ++k){
  				Sum_Sq[j][k] += Sq[j][k];
  			}
  		}
  	}

  	double comm_Sum_Sq[pp*pp];
    double comm_Sum_d[pp*pp];
  	double comm_all_Sum_Sq[pp*pp*n_proc];
  	double comm_all_Sum_d[pp*pp*n_proc];
    double comm_all_Sum_S[pp*n_proc];

  	for(int i = 0; i < pp; ++i){
  		for(int j = 0; j < pp; ++j){
  			comm_Sum_Sq[i*pp + j] = Sum_Sq[i][j];
  			comm_Sum_d[i*pp + j] = Sum_d[i][j];
  		}
  	}

  	MPI_Allgather(&comm_Sum_Sq[0], pp*pp, MPI_DOUBLE, &comm_all_Sum_Sq[0], pp*pp, MPI_DOUBLE, MPI_COMM_WORLD);
  	MPI_Allgather(&comm_Sum_d[0], pp*pp, MPI_DOUBLE, &comm_all_Sum_d[0], pp*pp, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&Sum_S[0], pp, MPI_DOUBLE, &comm_all_Sum_S[0], pp, MPI_DOUBLE, MPI_COMM_WORLD);

  	for(int i = 0; i < pp; ++i){
  		for(int j = 0; j < pp; ++j){
  			Sum_Sq[i][j] = 0.0;
  			Sum_d[i][j] = 0.0;
  			for(int k = 0; k < n_proc; ++k){
  				Sum_Sq[i][j] += comm_all_Sum_Sq[i*pp + j + k*pp*pp];
  			 	Sum_d[i][j] += comm_all_Sum_d[i*pp + j + k*pp*pp];
  			}
  		}
  	}

    for(int i = 0; i < pp; ++i){
      Sum_S[i] = 0;
      for(int k = 0; k < n_proc; ++k) Sum_S[i] += comm_all_Sum_S[i + pp*k];
    }

  	mat D_hat(pp, vector<double>(pp, 0) );
  	for(int i = 0; i < pp; ++i){
  		for(int j = 0; j < pp; ++j){
  			D_hat[i][j] = Sum_d[i][j]/n_phase;
  		}
  	}

    mat s_bar(pp, vector<double>(1, 0));
  	for(int i = 0; i < pp; ++i) s_bar[i][0] = Sum_S[i]/n_phase;

    mat Sum_Sq_temp(pp, vector<double>(pp, 0));
  	for(int i = 0; i < pp; ++i){
  		for(int j = 0; j < pp; ++j){
  			Sum_Sq_temp[i][j] = Sum_Sq[i][j]/n_phase;
  		}
  	}
  	mat s_bar_trans(1, vector<double>(pp, 0));
  	for(int i = 0; i < pp; ++i) s_bar_trans[0][i] = s_bar[i][0];
  	mat s_bar_temp = Product_matrix(s_bar, s_bar_trans);

  	mat s_sigma(pp, vector<double>(pp, 0));
  	for(int i = 0; i < pp; ++i){
  		for(int j = 0; j < pp; ++j){
  			s_sigma[i][j] = Sum_Sq_temp[i][j] - s_bar_temp[i][j];
  		}
  	}

  	mat inverse_s_sigma = Inverse_matrix(s_sigma);
    vector<double> s_obs = stat_gen(net_y, net_y_b, net_x, net_x_b, covar, cov_sec);
  	mat s_dev(pp, vector<double>(1, 0) );
  	for(int i = 0; i < pp; ++i) s_dev[i][0] = s_bar[i][0] - s_obs[i];

  	vector<double> tstat(pp);
  	for(int i = 0; i < pp; ++i) tstat[i] = (s_bar[i][0] - s_obs[i])/sqrt(s_sigma[i][i]);

  	mat s_dev_trans(1, vector<double>(pp, 0));
  	for(int i = 0; i < pp; ++i) s_dev_trans[0][i] = s_dev[i][0];
  	mat s_dev_temp = Product_matrix(inverse_s_sigma, s_dev);
  	mat max_test_temp = Product_matrix(s_dev_trans, s_dev_temp);

  	double max_test = sqrt(max_test_temp[0][0]);
  	mat inverse_D_hat = Inverse_matrix(D_hat);

  	mat inverse_D_hat_trans(pp, vector<double>(pp, 0) );
  	for(int i = 0; i < pp; ++i){
  		for(int j = 0; j < pp; ++j) inverse_D_hat_trans[j][i] = inverse_D_hat[i][j];
  	}
  	mat inverse_D_hat_temp = Product_matrix(s_sigma, inverse_D_hat_trans);
  	mat Cov_theta = Product_matrix(inverse_D_hat, inverse_D_hat_temp);

  	return make_tuple(D_hat, tstat, Cov_theta, s_sigma, max_test);
}
void phase2(vector<double>& theta,
  const net& net_y, const net& net_y_b,
  const net& net_x, const net& net_x_b,
  const cov& covar, const net& cov_sec,
  const h_set& f_subset_ind, const h_set& ind_subset_f,
  // const net& two_dis_1_t1, const net& two_dis_1_t3,
  const vector<vector<double> >& D_hat_temp,
  int my_rank, int p){

  	unsigned rand_seed[n_sub*n_phase_2];
  	for(int seed_i = 0; seed_i < n_sub*n_phase_2; ++seed_i) rand_seed[seed_i]= rd();

  	for(int sub = 0; sub < n_sub; ++sub){
  		if(my_rank == 0) cout << "current sub is " << sub << endl;
  		if(theta[0] < 0.0001){
  			if(my_rank == 0) cout << "negative value of theta[0] " << endl;
  			break;
  		}

      mat Sum_Sq(pp, vector<double>(pp, 0));
      vector<mat> store_Sq(n_phase_2, mat(pp, vector<double>(pp, 0) ) );
    	vector<double> Sum_S(pp, 0);

      int ib = n_phase_2/n_proc;
      int istat, iend;

      istat = my_rank*ib;
      iend = (my_rank + 1)*ib;

      #pragma omp parallel for schedule(dynamic, 1)
    	for(int sim_i = istat; sim_i < iend; ++sim_i){
    		unsigned rand = rand_seed[sub*n_phase_2 + sim_i];

    		//S_org
    		vector<double> S_org(pp, 0);
        net sim_net(net_x);
        net sim_net_b(net_x_b);
        simulation(theta, sim_net, sim_net_b, covar, cov_sec, f_subset_ind, ind_subset_f, rand);
        S_org = stat_gen(sim_net, sim_net_b, net_x, net_x_b, covar, cov_sec);

        #pragma omp critical
    		for(int j = 0; j < pp; ++j) Sum_S[j] += S_org[j];

        mat S_org_col(pp, vector<double>(1, 0));
    		for(int j = 0; j < pp; ++j) S_org_col[j][0] = S_org[j];

        mat S_org_row(1, vector<double>(pp, 0));
    		for(int j = 0; j < pp; ++j) S_org_row[0][j] = S_org[j];
    		mat Sq = Product_matrix(S_org_col, S_org_row);

        #pragma omp critical
        store_Sq[sim_i] = Sq;
    	}


      for(int j = 0; j < pp; ++j){
        for(int k = 0; k < pp; ++k){
          for(int sim_i = istat; sim_i < iend; ++sim_i) Sum_Sq[j][k] += store_Sq[sim_i][j][k];
        }
      }

    	double comm_Sum_Sq[pp*pp];
    	double comm_all_Sum_Sq[pp*pp*n_proc];
      double comm_all_Sum_S[pp*n_proc];

    	for(int i = 0; i < pp; ++i){
    		for(int j = 0; j < pp; ++j) comm_Sum_Sq[i*pp + j] = Sum_Sq[i][j];
    	}

    	MPI_Allgather(&comm_Sum_Sq[0], pp*pp, MPI_DOUBLE, &comm_all_Sum_Sq[0], pp*pp, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgather(&Sum_S[0], pp, MPI_DOUBLE, &comm_all_Sum_S[0], pp, MPI_DOUBLE, MPI_COMM_WORLD);

    	for(int i = 0; i < pp; ++i){
    		for(int j = 0; j < pp; ++j){
    			Sum_Sq[i][j] = 0.0;
    			for(int k = 0; k < n_proc; ++k) Sum_Sq[i][j] += comm_all_Sum_Sq[i*pp + j + k*pp*pp];
    		}
    	}

      for(int i = 0; i < pp; ++i){
        Sum_S[i] = 0;
        for(int k = 0; k < n_proc; ++k) Sum_S[i] += comm_all_Sum_S[i + pp*k];
      }

    	mat s_bar(pp, vector<double>(1, 0));
    	for(int i = 0; i < pp; ++i) s_bar[i][0] = Sum_S[i]/n_phase_2;

    	mat Sum_Sq_temp(pp, vector<double>(pp, 0));
    	for(int i = 0; i < pp; ++i){
    		for(int j = 0; j < pp; ++j) Sum_Sq_temp[i][j] = Sum_Sq[i][j]/n_phase_2;
    	}
    	mat s_bar_trans(1, vector<double>(pp, 0));
    	for(int i = 0; i < pp; ++i) s_bar_trans[0][i] = s_bar[i][0];
    	mat s_bar_temp = Product_matrix(s_bar, s_bar_trans);

    	mat s_sigma(pp, vector<double>(pp, 0));
    	for(int i = 0; i < pp; ++i){
    		for(int j = 0; j < pp; ++j)	s_sigma[i][j] = Sum_Sq_temp[i][j] - s_bar_temp[i][j];
    	}
    	mat inverse_s_sigma = Inverse_matrix(s_sigma);

    	vector<double> s_obs = stat_gen(net_y, net_y_b, net_x, net_x_b, covar, cov_sec);
    	mat s_dev(pp, vector<double>(1, 0));
    	for(int i = 0; i < pp; ++i) s_dev[i][0] = s_bar[i][0] - s_obs[i];

    	vector<double> tstat(pp);
    	for(int i = 0; i < pp; ++i) tstat[i] = (s_bar[i][0] - s_obs[i])/sqrt(s_sigma[i][i]);

    	mat s_dev_trans(1, vector<double>(pp, 0));
    	for(int i = 0; i < pp; ++i) s_dev_trans[0][i] = s_dev[i][0];
    	mat s_dev_temp = Product_matrix(inverse_s_sigma, s_dev);
    	mat max_test_temp = Product_matrix(s_dev_trans, s_dev_temp);

    	double max_test = sqrt(max_test_temp[0][0]);


  		if(my_rank == 0){
  			cout << "new paramters are ";
  			for(int i = 0; i < pp; ++i) cout << theta[i] << ", ";
  			cout << endl;

  			cout << "tstat are ";
  			for(int i = 0; i < pp; ++i) cout << tstat[i] << ", ";
  			cout << endl;

  			cout << "max_test is " << max_test << endl;

  			// display(D_hat, "D_hat_phase1");
  		}

  		mat inv_D_hat = Inverse_matrix(D_hat_temp);
  		// double a_gain = 0.3;

  		//parameter update
  		for(int i =0; i < pp; ++i){
  			for(int j=0; j < pp; ++j){
  				theta[i] -= a_gain * inv_D_hat[i][j] * s_dev[j][0];
  			}
  		}

    }
}


#endif /* PHASE_F_H_ */
