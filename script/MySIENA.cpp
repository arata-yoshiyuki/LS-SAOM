#include "constant.h"
#include "mpi.h"
#include "omp.h"
#include "Data.h"
#include "Statistics.h"
#include "phase_f.h"

int main(int argc, char** argv) {

  cov cov_1;
  net cov_sec_1;
  net net_1, net_1_b, net_2, net_2_b;

  data_cov(f_name_kj, cov_1);
  data_sec(f_name_kj, cov_sec_1);

  data(f_name_sk_1, net_1);
  data_b(f_name_sk_1, net_1_b);
  data(f_name_sk_2, net_2);
  data_b(f_name_sk_2, net_2_b);

  //	Parameter settings
  vector<double> theta(pre_theta);
  auto initial = initialize_f(net_2, net_1, net_1_b, cov_sec_1);
  auto f_subset_ind = get<0>(initial);
  auto ind_subset_f = get<1>(initial);

  // for gof_flag is true
  vector<unsigned int> sim_seed;
  random_device sim_rd;
  if(gof_flag || gof_2_flag){
    int n_sim = n_proc;
    for(int sim_i = 0; sim_i < n_sim; ++sim_i) sim_seed.push_back(sim_rd());
  }

  //Parallel
  double t_start, t_finish;
	int my_rank, p;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

  if(my_rank == 0) cout << "files: " << f_name_kj << ", " << f_name_sk_1 << ", " << f_name_sk_2 << "n_node " << n_proc << ", net_size " << net_size << ", n_phase_1 " << n_phase_1 << ", n_phase_2 " << n_phase_2 << ", n_phase_3 " << n_phase_3 << endl;

  //test one simulation
  if(sim_flag){
    net sim_net(net_1);
    net sim_net_b(net_1_b);

    t_start = MPI_Wtime();
    if(my_rank == 0) simulation(theta, sim_net, sim_net_b, cov_1, cov_sec_1, f_subset_ind, ind_subset_f, 1234567);
    t_finish = MPI_Wtime();
    if(my_rank == 0) cout << "One simulation " << t_finish - t_start << endl;
    /*
    if(my_rank == 0){
      vector<double> test_S(stat_gen(sim_net, sim_net_b, net_1, net_1_b, cov_1, cov_sec_1));
      for(int i = 0; i < pp; ++i) cout << "after test_S[" << i << "]: " << test_S[i] << endl;
    }
    */
  }

  //Phase 1. Get D_hat. d_ij is the derivative of statistics of i w.r.t. parameter j.
  mat D_hat_temp;
  if(ph_1_flag){
    t_start = MPI_Wtime();
    auto ini_phase = phase3(theta, net_2, net_2_b, net_1, net_1_b, cov_1, cov_sec_1, f_subset_ind, ind_subset_f, n_phase_1, my_rank, p);
    D_hat_temp = get<0>(ini_phase);
    if(my_rank == 0) display(D_hat_temp, "D_hat_temp");
    t_finish = MPI_Wtime();
    if(my_rank == 0) cout << "Phase 1 time is " << t_finish - t_start << endl;
  }else{
    D_hat_temp = pre_D_hat_temp;
  }

  //Phase 2. Convergence of parameters.
  if(ph_2_flag){
    t_start = MPI_Wtime();
    phase2(theta, net_2, net_2_b, net_1, net_1_b, cov_1, cov_sec_1, f_subset_ind, ind_subset_f, D_hat_temp, my_rank, p);
    if(my_rank == 0){
      cout << "Last paramters are " << endl;
      for(int i = 0; i < pp; ++i) cout << theta_name[i] << ": " << theta[i] << endl;
      cout << endl;
    }
    t_finish = MPI_Wtime();
    if(my_rank == 0) cout << "Convergence time " << t_finish - t_start << endl;
  }

  // Phase 3. Covariance of theta.
  if(ph_3_flag){
    t_start = MPI_Wtime();
    auto phase3_tuple = phase3(theta, net_2, net_2_b, net_1, net_1_b, cov_1, cov_sec_1, f_subset_ind, ind_subset_f, n_phase_3, my_rank, p);
    if(my_rank == 0){
      auto D_hat_phase3 = get<0>(phase3_tuple);
      auto tstat = get<1>(phase3_tuple);
      auto Cov_theta = get<2>(phase3_tuple);
      auto max_test = get<4>(phase3_tuple);

      display(D_hat_phase3, "D_hat_phase3");
      cout << "Last tstat are ";
      for(int i = 0; i < pp; ++i) cout << tstat[i] << ", ";
      cout << endl;
      cout << "Last max_test is " << max_test << endl;
      display(Cov_theta, "Cov_theta");

      cout << "Last paramters are " << endl;
      for(int i = 0; i < pp; ++i) cout << theta_name[i] << ", " << theta[i] << ", " << sqrt(Cov_theta[i][i]) << ", " << endl;
      cout << endl;
    }
    t_finish = MPI_Wtime();
    if(my_rank == 0) cout << "Phase 3 time is " << t_finish - t_start << endl;
  }

  // Simulations for Goodness of Fit
  if(gof_flag){
    net sim_net(net_1);
    net sim_net_b(net_1_b);
    simulation(theta, sim_net, sim_net_b, cov_1, cov_sec_1, f_subset_ind, ind_subset_f, sim_seed[my_rank]);

    ofstream fout("sim_test.txt");
    if(!fout) cout << "False FOUT" << endl;

    for(int i = 1; i <= net_size; ++i){
      for(auto itr = sim_net.at(i).begin(); itr != sim_net.at(i).end(); ++itr){
        fout << i << ", " << itr->first << endl;
      }
    }
    fout.close();
  }

  // Simulations for Goodness of Fit 2 for the entire period
  if(gof_2_flag){
    mt19937 gof_2_mt(sim_seed[my_rank]);
    net sim_net(net_1);
    net sim_net_b(net_1_b);
    simulation(sim_theta_A, sim_net, sim_net_b, cov_1, cov_sec_1, f_subset_ind, ind_subset_f, gof_2_mt());
    simulation(sim_theta_B, sim_net, sim_net_b, cov_1, cov_sec_1, f_subset_ind, ind_subset_f, gof_2_mt());
    simulation(sim_theta_C, sim_net, sim_net_b, cov_1, cov_sec_1, f_subset_ind, ind_subset_f, gof_2_mt());

    ofstream fout("sim_entire.txt");
    if(!fout) cout << "False FOUT" << endl;

    for(int i = 1; i <= net_size; ++i){
      for(auto itr = sim_net.at(i).begin(); itr != sim_net.at(i).end(); ++itr){
        fout << i << ", " << itr->first << endl;
      }
    }
    fout.close();
  }

  MPI_Finalize();
}
