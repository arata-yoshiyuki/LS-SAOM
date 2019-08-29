#ifndef DATA_H_
#define DATA_H_

#include <fstream>
#include <sstream>

void data(string f_name, net& network);
void data_b(string f_name, net& network);
void data_sec(string f_name, net& cov_sec);
void data_cov(string f_name, cov& covariate);
vector<string> split(string &s, char delim);

void data(string f_name, net& network){
	ifstream ifs(f_name);
	string sample;

	for(int i = 1; i <= net_size; ++i) network[i];
	while(getline(ifs, sample) ){
		vector<string> sample_split = split(sample, ',');
		vector<int> test_vec(2);
		for(int j = 0; j < 2; ++j) test_vec[j] = stoi(sample_split[j]);
		network[test_vec[0]][test_vec[1]] = 1;
	}
	ifs.close();
}
void data_b(string f_name, net& network){
	ifstream ifs(f_name);
	string sample;

	for(int i = 1; i <= net_size; ++i) network[i];
	while(getline(ifs, sample) ){
		vector<string> sample_split = split(sample, ',');
		vector<int> test_vec(2);
		for(int j = 0; j < 2; ++j) test_vec[j] = stoi(sample_split[j]);
		network[test_vec[1]][test_vec[0]] = 1;
	}
	ifs.close();
}
void data_sec(string f_name, net& cov_sec){
	ifstream ifs(f_name);
	string sample;
	int i = 1;	// firm ID starts from 1

	while(getline(ifs, sample)){
		vector<string> sample_split = split(sample, ',');
		for(int j = 0; j < 2; ++j) cov_sec[i][j] = stoi(sample_split[j]);
		i += 1;
	}
	ifs.close();
}
void data_cov(string f_name, cov& covariate){
	ifstream ifs(f_name);
	string sample;
	int i = 1;	// firm ID starts from 1

	while(getline(ifs, sample) ){
		vector<string> sample_split = split(sample, ',');
		for(int j = 0; j < n_item; ++j) covariate[i][j] = stod(sample_split[j]);
		i += 1;
	}
	ifs.close();
}

vector<string> split(string &s, char delim){
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
    if (!item.empty()) {
            elems.push_back(item);
        }
    }
    return elems;
}

#endif /* DATA_H_ */
