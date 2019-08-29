#ifndef STATISTICS_H_
#define STATISTICS_H_

double change(const net& net_y, const net& net_x);
double RateX_network(const net& net_y, const net& net_x);
double RateX(const net& net_y, const net& net_x, const cov& covar, const int& x_var);
double outdegree(const net& net_y, const int& actor_j);
double egoX(const net& net_y, const cov& covar, const int& x_var, const int& actor_j);
double altX(const net& net_y, const cov& covar, const int& x_var, const int& actor_j);
double inPop(const net& net_y, const net& net_x, const net& net_x_b, const int& actor_j);
double cont_similar(const net& net_y, const net& net_x, const cov& covar, const int& actor_j);
double disc_similar(const net& net_y, const net& net_x, const net& cov_sec, const int& actor_j);
double geo_dis(const net& net_y, const net& net_x, const cov& covar, const int& actor_j);
double recipro(const net& net_y, const net& net_x, const int& actor_j);
double two_dis_t1(const net& net_y, const net& net_x, const int& actor_j);
double two_dis_t3(const net& net_y, const net& net_x, const net& net_x_b, const int& actor_j);

vector<double> stat_gen(const net& net_y, const net& net_y_b,
		const net& net_x, const net& net_x_b,
		const cov& covar, const cov& cov_sec);


double change(const net& net_y, const net& net_x){
	double change = 0.0;
	for(int i = 1; i <= net_size; ++i){
	    for(auto itr = net_x.at(i).begin(); itr != net_x.at(i).end(); ++itr) {
	    	int j = itr->first;
				if(net_y.at(i).count(j) == 0) change += 1.0;
	    }
			for(auto itr = net_y.at(i).begin(); itr != net_y.at(i).end(); ++itr) {
	    	int j = itr->first;
				if(net_x.at(i).count(j) == 0) change += 1.0;
	    }
	}

	return change;
}
double RateX_network(const net& net_y, const net& net_x){
	double Rate_stat = 0.0;
	for(int i = 1; i <= net_size; ++i){
		for(auto itr = net_x.at(i).begin(); itr != net_x.at(i).end(); ++itr) {
			int j = itr->first;
			if(net_y.at(i).count(j) == 0) Rate_stat += log(net_x.at(i).size() + 1);
		}
		for(auto itr = net_y.at(i).begin(); itr != net_y.at(i).end(); ++itr) {
			int j = itr->first;
			if(net_x.at(i).count(j) == 0) Rate_stat += log(net_x.at(i).size() + 1);
		}
	}
	return Rate_stat;
}
double RateX(const net& net_y, const net& net_x, const cov& covar, const int& x_var){
	double Rate_stat = 0.0;
	for(int i = 1; i <= net_size; ++i){
		for(auto itr = net_x.at(i).begin(); itr != net_x.at(i).end(); ++itr) {
			int j = itr->first;
			if(net_y.at(i).count(j) == 0) Rate_stat += covar.at(i).at(x_var);
		}
		for(auto itr = net_y.at(i).begin(); itr != net_y.at(i).end(); ++itr) {
			int j = itr->first;
			if(net_x.at(i).count(j) == 0) Rate_stat += covar.at(i).at(x_var);
	  }
	}
	return Rate_stat;
}
double outdegree(const net& net_y, const int& actor_j){
	double out_deg = (double)net_y.at(actor_j).size();
	return out_deg;
}
double egoX(const net& net_y, const cov& covar, const int& x_var, const int& actor_j){
	double new_egoX = 0.0;
	new_egoX = net_y.at(actor_j).size()*covar.at(actor_j).at(x_var);
	return new_egoX;
}
double altX(const net& net_y, const cov& covar, const int& x_var, const int& actor_j){
	double new_altX = 0.0;
	for(auto itr = net_y.at(actor_j).begin(); itr != net_y.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		new_altX += covar.at(sup_id).at(x_var);
	}
	return new_altX;
}
double inPop(const net& net_y, const net& net_x, const net& net_x_b, const int& actor_j){
	double inPop = 0.0;
	for(auto itr = net_y.at(actor_j).begin(); itr != net_y.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_x.at(actor_j).count(sup_id) == 0){ // newly created
			inPop += log(net_x_b.at(sup_id).size() + 1);
		}
	}
	for(auto itr = net_x.at(actor_j).begin(); itr != net_x.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_y.at(actor_j).count(sup_id) == 0){ // ties deleted
			inPop -= log(net_x_b.at(sup_id).size() + 1);
		}
	}
	return inPop;
}
double cont_similar(const net& net_y, const net& net_x, const cov& covar, const int& actor_j){
	double new_cont_similar = 0.0;
	for(auto itr = net_y.at(actor_j).begin(); itr != net_y.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_x.at(actor_j).count(sup_id) == 0 )	new_cont_similar += abs(covar.at(sup_id).at(log_sale) - covar.at(actor_j).at(log_sale));
	}
	for(auto itr = net_x.at(actor_j).begin(); itr != net_x.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_y.at(actor_j).count(sup_id) == 0 ) new_cont_similar -= abs(covar.at(sup_id).at(log_sale) - covar.at(actor_j).at(log_sale));
	}
	return new_cont_similar;
}
double disc_similar(const net& net_y, const net& net_x, const net& cov_sec, const int& actor_j){
	double new_disc_similar = 0.0;
	for(auto itr = net_y.at(actor_j).begin(); itr != net_y.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_x.at(actor_j).count(sup_id) == 0){ // newly created
			if(cov_sec.at(actor_j).at(kj_ind) == cov_sec.at(sup_id).at(kj_ind)) new_disc_similar += 1.0;
		}
	}

	for(auto itr = net_x.at(actor_j).begin(); itr != net_x.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_y.at(actor_j).find(sup_id) == net_y.at(actor_j).end()){ // ties deleted
			if(cov_sec.at(actor_j).at(kj_ind) == cov_sec.at(sup_id).at(kj_ind)) new_disc_similar -= 1.0;
		}
	}
	return new_disc_similar;
}
double geo_dis(const net& net_y, const net& net_x, const cov& covar, const int& actor_j){
	double new_geo_dis = 0.0;
	for(auto itr = net_y.at(actor_j).begin(); itr != net_y.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_x.at(actor_j).count(sup_id) == 0 ){ // newly created
			double dis_X = covar.at(sup_id).at(loc_X) - covar.at(actor_j).at(loc_X);
			double dis_Y = covar.at(sup_id).at(loc_Y) - covar.at(actor_j).at(loc_Y);
			new_geo_dis += sqrt(pow(dis_X, 2) + pow(dis_Y, 2));
		}
	}
	for(auto itr = net_x.at(actor_j).begin(); itr != net_x.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_y.at(actor_j).count(sup_id) == 0 ){ // ties deleted
			double dis_X = covar.at(sup_id).at(loc_X) - covar.at(actor_j).at(loc_X);
			double dis_Y = covar.at(sup_id).at(loc_Y) - covar.at(actor_j).at(loc_Y);
			new_geo_dis -= sqrt(pow(dis_X, 2) + pow(dis_Y, 2));
		}
	}
	return new_geo_dis;
}
double recipro(const net& net_y, const net& net_x, const int& actor_j){
	double new_recipro = 0.0;
	for(auto itr = net_y.at(actor_j).begin(); itr != net_y.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_x.at(actor_j).count(sup_id) == 0){ // newly created
			if(net_x.at(sup_id).count(actor_j) == 1) new_recipro += 1.0;
		}
	}
	for(auto itr = net_x.at(actor_j).begin(); itr != net_x.at(actor_j).end(); ++itr){
		int sup_id = itr->first;
		if(net_y.at(actor_j).count(sup_id) == 0){ // a tie deleted
			if(net_x.at(sup_id).count(actor_j) == 1) new_recipro -= 1.0;
		}
	}
	return new_recipro;
}
double two_dis_t1(const net& net_y, const net& net_x, const int& actor_j){
	double new_two_dis = 0.0;
	for(auto itr_y = net_y.at(actor_j).begin(); itr_y != net_y.at(actor_j).end(); ++itr_y){
		int indir_k = itr_y->first;
		if(net_x.at(actor_j).count(indir_k) == 0){	//a tie to indir_id is newly created
			for(auto itr_x = net_x.at(actor_j).begin(); itr_x != net_x.at(actor_j).end(); ++itr_x){
				int temp_id = itr_x->first;
				if(net_x.at(temp_id).count(indir_k)){
					new_two_dis += 1.0;
					break;
				}
			}
		}
	}
	return new_two_dis;
}
double two_dis_t3(const net& net_y, const net& net_x, const net& net_x_b, const int& actor_j){
	double new_two_dis = 0.0;
	for(auto itr_y = net_y.at(actor_j).begin(); itr_y != net_y.at(actor_j).end(); ++itr_y){
		int indir_k = itr_y->first;
		if(net_x.at(actor_j).count(indir_k) == 0){	//a tie to indir_id is newly created
			for(auto itr_x = net_x_b.at(actor_j).begin(); itr_x != net_x_b.at(actor_j).end(); ++itr_x){
				int temp_id = itr_x->first;
				if(net_x_b.at(temp_id).count(indir_k)){
					new_two_dis += 1.0;
					break;
				}
			}
		}
	}
	return new_two_dis;
}

vector<double> stat_gen(const net& net_y, const net& net_y_b,
		const net& net_x, const net& net_x_b,
		const cov& covar, const net& cov_sec){

	vector<double> S(pp, 0);
	S[0] = change(net_y, net_x);
	S[1] = RateX_network(net_y, net_x);
	S[2] = RateX(net_y, net_x, covar, log_age);

	for(int i = 1; i < (net_size + 1); ++i) S[3] += outdegree(net_y, i);
	for(int i = 1; i < (net_size + 1); ++i) S[4] += egoX(net_y, covar, log_sale, i);
	for(int i = 1; i < (net_size + 1); ++i) S[5] += altX(net_y, covar, d_lab_pro, i);
	for(int i = 1; i < (net_size + 1); ++i) S[6] += altX(net_y, covar, g_rate, i);
	for(int i = 1; i < (net_size + 1); ++i) S[7] += altX(net_y, covar, pf_rate, i);
	for(int i = 1; i < (net_size + 1); ++i) S[8] += altX(net_y, covar, log_sale, i);
	for(int i = 1; i < (net_size + 1); ++i) S[9] += inPop(net_y, net_x, net_x_b, i);

	for(int i = 1; i < (net_size + 1); ++i) S[10] += cont_similar(net_y, net_x, covar, i);
	for(int i = 1; i < (net_size + 1); ++i) S[11] += disc_similar(net_y, net_x, cov_sec, i);
	for(int i = 1; i < (net_size + 1); ++i) S[12] += geo_dis(net_y, net_x, covar, i);

	for(int i = 1; i < (net_size + 1); ++i) S[13] += recipro(net_y, net_x, i);
	for(int i = 1; i < (net_size + 1); ++i) S[14] += two_dis_t1(net_y, net_x, i);
	for(int i = 1; i < (net_size + 1); ++i) S[15] += two_dis_t3(net_y, net_x, net_x_b, i);

	return S;
}

#endif /* STATISTICS_H_ */
