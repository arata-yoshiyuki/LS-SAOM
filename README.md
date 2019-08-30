# LS-SAOM
Large-scale stochastic actor-oriented model (LS-SAOM) implementation.

* LS-SAOM is a C++ program to implement the SAOM estimation for a large network data.
* I applied this program to Japanese inter-firm transaction network data in [Arata and Mundt (2019)](https://www.rieti.go.jp/jp/publications/dp/19e027.pdf) and uploaded it here for replication purpose.
* This program is still under construction and any comments are welcome.

## 1. Overview
### 1.1. SAOM

* The SAOM is an empirical model to estimate how a network is dynamically generated (see Snijders (1996)).
* Suppose that a network consists of firms and ties between the firms, and the ties change over time.
* The SAOM assumes that the network evolution within a period is the accumulation of many micro-steps, at each of which a randomly selected firm can change its ties to others such that its objective function is maximized.
[](* The objective function is assumed to depend on a firm's characteristics and local network structure.)
* An important feature is that the network evolution is sequential, i.e., a firm changes its ties given the *current* network, which is the consequence of changes that have occurred before.
* Main parameters to be estimated is the ones that characterize the objective function.

### 1.2. Estimation

* The estimation is based on methods of moments, i.e., we search for a set of parameters that fit a set of empirical moments most closely.
* For calculation of moments given a set of proposed parameters, we use a simulation method.
* In practice, the estimation process consists of the following three steps:
  1. Given a set of proposed parameters, we simulate the model and calculate the average of moments across the simulations (i.e., simulated moments).
  2. The set of parameters are updated such that simulated moments with the new parameters would become closer to the empirical ones.
  3. Repeat 1. and 2. until convergence.
* To reduce the computation time, simulations in Step 2 are parallelized.

## 2. C++ codes

* *MySIENA.cpp* is the main file, in which Step 1-3 are implemented.
* *constant.cpp* is a file in which a set of parameters are given, including a set of constants (e.g., the number of firms) and initial parameters used in estimation.
* Compiling and linking the two files together, we obtain the object file. File *s_script.sh* describes the compile options in my environment.
* By executing the object file, we get an estimation output.

## 3. Results

* Output files in my analysis are given under the folder *result*.

## 4. References

* Arata, Y. and Mundt. P., 2019. Topology and Formation of Production Input Interlinkages: Evidence from Japanese microdata, RIETI DP, 19-E-027.
* Snijders, T., 1996. Stochastic actor-oriented models for network change. The Journal of Mathematical Sociology 21, 149-172.

