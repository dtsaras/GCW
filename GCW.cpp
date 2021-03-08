#include <sys/time.h>
#include "option.h"
#include "hypergraph.hpp"
#include "sfmt/SFMT.h"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

double print_time_main(struct timeval &start, struct timeval &end){
	double usec;

	usec = (end.tv_sec*1000 + (end.tv_usec/1000)) - (start.tv_sec*1000 + start.tv_usec/1000);
	return usec;
}

bool calculateInfluence_DSSA(HyperGraph & hg, Graph & g, vector<int> & seeds, int t, double & deg, float epsilon, float delta, int m, long long int maxSamples, int iter, sampler<int> &samp, int p){
	long long  counter = 0;
	int n = g.getSize();
	// cout << "!!!!!!!!!!!!!!!! 1" << endl;
	unsigned k = seeds.size();
	vector<unsigned int> link(n + 1, seeds.size());
	double f = (log(6/delta)+lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))*n/(k*log(6*log2(n)/delta));
	double lambda1 = 1+(1+epsilon)*(2+2*epsilon/3)*log(3*log2(f)/delta)/(epsilon*epsilon);
	double degree=0;
	// cout << "!!!!!!!!!!!!!!!! 2" << endl;
	for (unsigned int i = 0; i < k;++i){
		link[seeds[i]] = i;
	}
	vector<bool> maxSeed(t, false);
	omp_set_num_threads(t);
	// cout << "!!!!!!!!!!!!!!!! 3" << endl;
	#pragma omp parallel
	{
		vector<bool> visit(n+1,false);
		vector<int> visit_mark(n,0);
		int id = omp_get_thread_num();

		if (m == 0){
			while(counter < maxSamples){
				// cout << "iteration: " << counter << endl; // del
				// cout << "!!!!!!!!!!!!!!!! 4" << endl;
				maxSeed[id]=hg.pollingLT2(g,link,k,visit,visit_mark, samp, p);
				// cout << "!!!!!!!!!!!!!!!! 5" << endl;
				#pragma omp critical
				{
					counter += 1;
					if (maxSeed[id]){
						degree++;
					}
				}
				// cout << "!!!!!!!!!!!!!!!! 6" << endl;
				// cin.ignore(numeric_limits <streamsize>::max(), '\n');
			}
		}
		else {
			while(counter < maxSamples){
				maxSeed[id]=hg.pollingIC2(g,link,k,visit,visit_mark, samp, p);
	            #pragma omp critical
	            {
					counter += 1;
	                if (maxSeed[id]){
	                    degree++;
	                }
	            }
        	}
		}	
	}
	// cout << "!!!!!!!!!!!!!!!! 7" << endl;

    if (degree >= lambda1){
		double epsilon_1 = (deg*n/maxSamples)/(degree*n/counter) - 1;
        double epsilon_2 = epsilon*sqrt(n*(1+epsilon)/(degree*n*pow(2,iter-1)/counter));
		double epsilon_3 = epsilon*sqrt(n*(1+epsilon)*(1-1/exp(1)-epsilon)/((1+epsilon/3)*degree*n*pow(2,iter-1)/counter));
	    if ((epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1-1/exp(1)) <= epsilon){
			return true;
        }
    }

	hg.updateDeg();
	// cout << "!!!!!!!!!!!!!!!! 8" << endl;
	return false;
}

bool calculateInfluence_SSA(Graph & g, HyperGraph & hg, vector<int> & seeds, int t, double expected, double epsilon_1, float epsilon_2, float delta, int m, long long int maxSamples, long long int & checkSam, sampler<int> &samp, int p){
	long long  counter = 0;
	int n = g.getSize();
	unsigned int k = seeds.size();
	double f = log2(n*(log(2/delta)+lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/(k*log(2/delta)));
	double upperBound = 1+(2+2*epsilon_2/3)*(1 + epsilon_2)*(log(3/delta)+log(f))/(epsilon_2*epsilon_2);
	vector<unsigned int> link(n+1,k);
	int degree=0;
	vector<float> benefits(n,0);
	for (unsigned int i = 0; i < k;++i){
		link[seeds[i]] = i;
	}
	vector<bool> maxSeed(t, false);
	omp_set_num_threads(t);
	#pragma omp parallel
	{
		vector<bool> visit(n+1,false);
		vector<int> visit_mark(n,0);
		int id = omp_get_thread_num();
		if (m == 0){
		    while (counter <= maxSamples && degree < upperBound){
	            maxSeed[id]=hg.pollingLT(g,link,k,visit,visit_mark, samp, p);
				#pragma omp critical
				{
                	counter++;
					if (maxSeed[id]){
						degree++;
					}
				}
        	}
		} 
		else {
			while (counter <= maxSamples && degree < upperBound){
                maxSeed[id]=hg.pollingIC(g,link,k,visit,visit_mark, samp, p);
                #pragma omp critical
                {
                    counter++;
                    if (maxSeed[id]){
                        degree++;
                    }
                }
            }
		}
	}
	checkSam += counter;

	if (expected <= (1 + epsilon_1)*((double)degree*n/(double)counter)){
		return true;
	}
	return false;
}

void DSSA(Graph & g, int n, int k, int t, vector<int> &seeds, float epsilon, float delta, sampler<int> &samp, int p, int m){
	// cout << "################## 1" << endl;
	HyperGraph hg(n);
	vector<double> degree(k+1,0);

	double f = (log(6/delta)+lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))*n/(k*log(6*log2(n)/delta));
	// cout << "################## 2" << endl;
	double lambda = (2+2*epsilon/3)*log(3*log2(f)/delta)/(epsilon*epsilon);

	long long int totalSamples = (long long int)lambda;

	int iter = 1;

	// cout << "################## 2" << endl;
	addHyperedge(g,hg,t,totalSamples,m, samp, p);
	double nmax = (2+2*epsilon/3)*(lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) + log(6/delta))*n/(epsilon*epsilon*k);

	// cout << "################## 3" << endl;
	while (totalSamples < nmax){
		seeds.clear();
		totalSamples = hg.getNumEdge();
		buildSeedSet(hg,seeds,n,k,degree);
		// cout << "$$$$$$$$$$$$$$ 1" << endl;
		if (calculateInfluence_DSSA(hg,g,seeds,t,degree[k],epsilon,delta,m,totalSamples,iter, samp, p)){
           	break;
        }
        // cout << "$$$$$$$$$$$$$$ 2" << endl;
		iter++;
	}
	// cout << "################## 3" << endl;
}

void SSA(Graph & g, int n, int k, int t, vector<int> &seeds, float epsilon, float delta, sampler<int> &samp, int p, int m){
	HyperGraph hg(n);
	vector<double> degree(k+1,0);

	double epsilon_1= epsilon/(2*(1-1/exp(1)-epsilon));
	double epsilon_2= epsilon/(3*(1-1/exp(1)-epsilon));
	double epsilon_3= (epsilon - (epsilon_1 - epsilon_2 - epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon))/(1-1/exp(1));

	double f = log2(n*(log(2/delta)+lgammal(n+1)-lgammal(k+1)-lgammal(n-k+1))/(k*log(2/delta)));

	double degreeRequired = (2+2*epsilon_3/3)*(1+epsilon_1)*(1+epsilon_2)*log(3*f/delta)/(epsilon_3*epsilon_3);

	long long int curSamples = (long long) degreeRequired;
	long long int totalSamples = 0;
	long long int checkSam = 0;

	while (true){
		seeds.clear();
		addHyperedge(g,hg,t,curSamples,m, samp, p);
		curSamples *= 2;
		totalSamples = hg.getNumEdge();
		buildSeedSet(hg,seeds,n,k,degree);
		if (degree[k] < degreeRequired){
			continue;
		}
		
		if (calculateInfluence_SSA(g, hg, seeds, t, (double)degree[k]*n/(double)hg.getNumEdge(), epsilon_1, epsilon_2, delta, m, totalSamples*2*(1+epsilon_2)/(1-epsilon_2)*(epsilon_3*epsilon_3)/(epsilon_2*epsilon_2), checkSam, samp, p)){
			break;
		}
	}
}

void create_dist(Graph & g, int p, vector<float> &RR_user_Prs){
	float IP;
	float not_AP;
	// cout << "Event: " << p << endl; 
	for(int v = 1; v <= g.getSize(); v++){
		not_AP = 1.0;
		map<double, int>::iterator it = g.A[v].begin();
		while (it->second != p and not_AP != 0){
			not_AP *= 1 - g.AP[it->second][v]; //Probability not influenced by the stronger competitors
			++it;
		}
		// if (v <= g.getSize()){
		// 	cout << v << " " << not_AP << endl;
		// }
		IP = (g.similarities[v][p]/g.tsim_V_p[p]) * not_AP;
		RR_user_Prs[v - 1] = IP;
	}
}

int count_change_strat(Graph & g, Locations * usersLocations, int p){
	int changedSeedsets = 0;
	for (int seed : g.seedSets[p]){
		// Did not find it in the previous seedset
		if (g.seedSets_before[p].find(seed) == g.seedSets_before[p].end()){
			bool not_same_location = true;
			for (int old_seed : g.seedSets_before[p]){
				if (usersLocations->same_location(seed, old_seed)){
					not_same_location = false;
					break;	
				}
			}
			if (not_same_location)
				++changedSeedsets;
		}
	}		
	if(changedSeedsets > 0){
		g.seedSets_before[p].clear();
		for (int seed : g.seedSets[p])
			g.seedSets_before[p].insert(seed);
	}
}

void output_stats_k_and_p(ostream * output, double *** data, string variable, string algorithms, int reps, int * eventsArray, int * budgetsArray){
	string names[3] = {"total similarity", "time", "memory"};
	for (int d = 0; d < 3; ++d){
		*output << names[d] << endl;
		*output << setw(16) << "#";
		for (int i = 0; i < reps; ++i){
			if (variable == "c")
				*output << setw(16) << eventsArray[i];
			else
				*output << setw(16) << budgetsArray[i];
		}
		*output << endl;

		for (int i = 0; i < algorithms.size(); i += 2){
			if (algorithms.substr(i,2) == "EG")
				*output << setw(16) << "GCW";
			else
				*output << setw(16) << algorithms.substr(i,2);
			for (int j = 0; j < reps; ++j){
				*output << setw(16) << data[d][i/2][j];
			}
			*output << endl;
		}
	}
}

void output_stats_q_and_t(ostream * output, string algorithms, vector<double> * round_data_time, vector<double> * round_data_score, int rounds){
	for (int i = 0; i < algorithms.size(); i += 2){
		*output << algorithms.substr(i,2) << endl;
		for (int j = 0; j <= rounds; ++j){
			*output << round_data_time[i/2][j] << "\t" << round_data_score[i/2][j] << endl;
		}
	}
}

void output_basic_utility(ostream * output, vector<double> sim_Sj, vector<double> sim_S_S_Sj, vector<double> DSS, vector<double> final_sim, string variable, int reps, int * eventsArray, int * budgetsArray){
	if (variable == "c")
		*output << setw(35) << "Competitors";
	else
		*output << setw(35) << "Budgets";
	*output << setw(35) << "avg[Sigma(Sj)]" << setw(35) << "avg[Sigma(S) - Sigma(S - {Sj})]" << setw(35) << "D(S,S)" << setw(35) << "Total Similarity" << endl;
	for (int i = 0; i < reps; ++i){
		if (variable == "c")
			*output << setw(35) << eventsArray[i];
		else
			*output << setw(35) << budgetsArray[i];

		*output << setw(35) << sim_Sj[i] << setw(35) << sim_S_S_Sj[i] << setw(35) << DSS[i] << setw(35) << final_sim[i] << endl;
	}
}

void output_mult_rounds(ostream * output, map <string, vector<double>> data_tsim_vec, map <string, vector<double>> data_time_vec, map <string, vector<double>> data_memo_vec, string variable, int reps, int * eventsArray, int * budgetsArray, int rounds){
	string names[3] = {"total similarity", "time", "memory"};

	cout << "###############1" << endl;
	*output << names[0] << endl;
	*output << setw(16) << "#";
	for (int i = 0; i < reps; ++i){
		if (variable == "c")
			*output << setw(16) << eventsArray[i];
		else
			*output << setw(16) << budgetsArray[i];
	}
	*output << endl;

	cout << "###############2" << endl;
	// for (int i = 0; i < algorithms.size(); i += 2){
	for (pair<string, vector<double>> it : data_tsim_vec){
		if (it.first.substr(0,2) == "EG")
			*output << setw(16) << "GCW" + it.first.substr(2,string::npos);
		else
			*output << setw(16) << it.first;
		for (int j = 0; j < reps; ++j){
			*output << setw(16) << it.second[j];
		}
		*output << endl;
	}

	cout << "###############3" << endl;
	*output << names[1] << endl;
	*output << setw(16) << "#";
	for (int i = 0; i < reps; ++i){
		if (variable == "c")
			*output << setw(16) << eventsArray[i];
		else
			*output << setw(16) << budgetsArray[i];
	}
	*output << endl;

	cout << "###############4" << endl;
	// for (int i = 0; i < algorithms.size(); i += 2){
	for (map <string, vector<double>>::iterator it = data_time_vec.begin(); it != data_time_vec.end(); ++it){
		if (it->first.substr(0,2) == "EG")
			*output << setw(16) << "GCW" + it->first.substr(2,string::npos);
		else
			*output << setw(16) << it->first;
		for (int j = 0; j < reps; ++j){
			*output << setw(16) << it->second[j];
		}
		*output << endl;
	}

	cout << "###############5" << endl;
	*output << names[2] << endl;
	*output << setw(16) << "#";
	for (int i = 0; i < reps; ++i){
		if (variable == "c")
			*output << setw(16) << eventsArray[i];
		else
			*output << setw(16) << budgetsArray[i];
	}
	*output << endl;

	cout << "###############6" << endl;

	// for (int i = 0; i < algorithms.size(); i += 2){
	for (map <string, vector<double>>::iterator it = data_memo_vec.begin(); it != data_memo_vec.end(); ++it){
		if (it->first.substr(0,2) == "EG")
			*output << setw(16) << "GCW" + it->first.substr(2,string::npos);
		else
			*output << setw(16) << it->first;
		for (int j = 0; j < reps; ++j){
			*output << setw(16) << it->second[j];
		}
		*output << endl;
	}
	cout << "###############7" << endl;
}

void output_clust(ostream * output, double *** data, string variable, string algorithms){
	for (int i = 0; i < algorithms.size(); i += 2){
		if (algorithms.substr(i,2) == "EG"){
			*output << setw(16) << "GCW";
		}
		else{
			*output << setw(16) << algorithms.substr(i,2);
		}
		*output << setw(16) << data[0][i/2][0] << setw(16) << data[1][i/2][0] << setw(16) << data[2][i/2][0] << endl;
	}
}

int main(int argc, char ** argv){
	struct timeval total_start, total_end, read_graphs_start, read_graphs_end, start_round, end_round, start_init, end_init;
	srand(time(NULL));

	///////////////////////////// Parameters /////////////////////////////
	int m; // Model (Default IC)
	int k; //Budget
	int P; //Number of events
	int t; // num_of_threads (default 1)
	int reps, rounds;
	int num_of_monte_carlo_iterations; //BFS mc1
	int monte_carlo_iterations_tsim; // mc2
	float epsilon, delta, percent, mu, sigma;



	string variable, algorithms, dataset, usersLocationsFile, eventsLocationsFile, cond, fileName, ris_alg, dt_type, sim_file;
	ostream * output;
	bool exp1, exp2, exp3, exp4, exp5, sim, rand_edges;
	// int budgetsArray[] = {1, 2, 4, 8, 16, 32, 64, 128}; // {1, 2, 4, 10, 20, 30, 50, 70, 100} 
	int budgetsArray[] = {1, 5, 10, 20, 50};
	// int budgetsArray[] = {16, 80, 160, 320, 800};
	// int eventsArray[] = {2, 4, 8, 16, 32, 64, 128}; // 
	int eventsArray[] = {4, 8, 16, 32, 64};
	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}
	op.setPara(algorithms, dataset, usersLocationsFile, eventsLocationsFile, sim_file, dt_type, epsilon, delta, k, P, num_of_monte_carlo_iterations, monte_carlo_iterations_tsim, t, percent, reps, variable, cond, rounds, m, exp1, exp2, exp3, exp4, exp5, fileName, ris_alg, sim, rand_edges, mu, sigma);
	
	// Output file
	std::ofstream file;
	char * outFile = op.getPara("-o");
	if (outFile == NULL){
		output = &std::cout;
	}
	else{
		file.open(fileName);
		output = &file;
	}

	cout << "dataset: " << dataset << endl;
	cout << "users Locations file: " << usersLocationsFile << endl;
	cout << "events Locations file: " << eventsLocationsFile << endl;
	cout << "similarities file: " << sim_file << endl;
	cout << "dataset type: " << dt_type << endl;


	//////////////////// Quality Experiments ////////////////////////////
	double ** data_tsim = new double*[algorithms.size()/2];
	double ** data_time = new double*[algorithms.size()/2];
	double ** data_memo = new double*[algorithms.size()/2];

	for (int i = 0; i < algorithms.size()/2; ++i){
		data_tsim[i] = new double[reps];
		data_time[i] = new double[reps];
		data_memo[i] = new double[reps];
	}

	string algorithm;
	map<string, vector<double>> data_tsim_vec;
	map<string, vector<double>> data_time_vec;
	map<string, vector<double>> data_memo_vec;
	for (int a = 0; a < algorithms.size(); a += 2){
		algorithm = algorithms.substr(a,2);
		for (int i = 1; i <= rounds; i++){
			if (!(algorithm == "NA" && i > 1)){
				data_tsim_vec[algorithm + "." + to_string(i)] = vector<double> ();
				data_time_vec[algorithm + "." + to_string(i)] = vector<double> ();
				data_memo_vec[algorithm + "." + to_string(i)] = vector<double> ();
			}
		}
	}

	//////////////////// Rounds Stats ////////////////////////////////////
	vector<double> round_time, round_score;
	vector<double> * round_data_time = new vector<double>[algorithms.size()/2];
	vector<double> * round_data_score = new vector<double>[algorithms.size()/2];
		

	//////////////////// For Basic Ultility Stats ///////////////////////
	double stats[] = {0, 0, 0};
	double avg_tot_sim = 0;
	double avg_similarity_Sj = 0;
	double avg_similarity_S_Sj = 0;
	double total_DSS = 0;
	vector<double> sim_Sj, sim_S_S_Sj, DSS, final_sim; 

	////////////////// Algorithm //////////////////////

	

	Graph g(dataset, m, sim, rand_edges, mu, sigma);
	g.set_events(P);
	Utilities util;
	Locations* usersLocations = new Locations(usersLocationsFile);
	Locations* eventsLocations = new Locations(eventsLocationsFile);
	if (dt_type == "gsn" || dt_type == "geosocial" || dt_type == "GSN" || dt_type == "Geosocial"){
		g.create_AV_sets(usersLocations, eventsLocations, P, util);
	}
	else{
		g.read_similarities(sim_file, P);
	}
		

	cout << "Read graph" << endl;
	gettimeofday(&total_start, NULL);

	for(int event = 0; event < reps; event++){
		if (variable == "c"){ 
			P = eventsArray[event];
			g.change_events(P);
			if (dt_type == "gsn" || dt_type == "geosocial" || dt_type == "GSN" || dt_type == "Geosocial"){
				g.create_AV_sets(usersLocations, eventsLocations, P, util);
			}
			else{
				g.read_similarities(sim_file, P);
			}
		}
		else if (variable == "k"){
			k = budgetsArray[event];
		}

		////////////////////// Event Permutation + AP initialization ///////////////////////
		int* eventsRandomPermutation = new int[P];
		for (int p = 0; p < P; ++p){
			eventsRandomPermutation[p] = p;
			for (int v = 0; v <= g.getSize(); ++v){
		        g.AP[p][v] = 0;
			}
		}
		random_shuffle(eventsRandomPermutation, eventsRandomPermutation + P);


		/////////////////////////////////////////// Run for each algorithm ///////////////////////////////////////
		for (int a = 0; a < algorithms.size(); a += 2){
			algorithm = algorithms.substr(a,2);
			round_time.clear();
			round_score.clear();

			//////////////////////// Change for each algorithm ///////////////////////////////////////////////
			gettimeofday(&start_init, NULL);
			// if (algorithm != "EG" && algorithm != "NA")
			// 	g.initial_seedset(algorithm, util, usersLocations, eventsLocations, P, k);

			// if (algorithm == "CG" || algorithm == "RG"){
			// 	for (int p = 0; p < P; ++p){
			// 		for(int iteration = 0; iteration < num_of_monte_carlo_iterations; iteration++)
			// 			if (m == 0){
			// 				g.rand_LT();
			// 				g.BFS_LT(p, g.AP);
			// 			}
			// 			else
			// 				g.BFS_IC(p, g.AP);
			// 		for (int v = 1; v <= g.getSize(); v++)
			// 			g.AP[p][v] /= num_of_monte_carlo_iterations;
			// 	}
			// }
			gettimeofday(&end_init, NULL);
			float tot_time = print_time_main(start_init, end_init)/1000;

			double old_tsim = 0;
			double new_tsim = 0;
			total_DSS = 0;
			avg_tot_sim = 0;
			avg_similarity_Sj = 0;
			avg_similarity_S_Sj = 0;
			// if (algorithm == "CG" || algorithm == "RG"){
			// 	#pragma omp parallel for
			// 	for(int MC_sim_it = 0; MC_sim_it < monte_carlo_iterations_tsim; MC_sim_it++)
			// 		new_tsim += g.calculate_tsim(P);
			// 	new_tsim = new_tsim/(double)monte_carlo_iterations_tsim;
			// }
			cout << "Round 0. (initialization) (alg: " << algorithm << ", |C| = " << P << ", k = " << k << ") Total Similarity: " << new_tsim << endl;

			round_time.push_back(tot_time);
			round_score.push_back(new_tsim);
			bool condition = true;
			bool sameSeedsets = false;
			int changedSeedsets;
			int round = 0; 

			///////////// AV sets ////////////////////
			// for (int i = 1; i < g.getSize(); i++){
			// 	map<double, int>::iterator it = g.A[i].begin();
			// 	cout << i << ":" << endl;
			// 	while (it != g.A[i].end()){
			// 		cout << " (" << it->first << ", " << it->second << ')';
			// 		++it;
			// 	}
			// 	cout << endl;
			// }

			//////////////////////////////////////////////////////////////////////////////////////////
			// if (algorithm == "CG" || algorithm == "RG" || algorithm == "EG" || algorithm == "NA"){
			if (algorithm == "EG" || algorithm == "NA"){
				while (condition){
					++round;
					changedSeedsets = 0;
					sameSeedsets = true;
					int rand_event;
					gettimeofday(&start_round, NULL);
					//////////////////////////////// Best Response //////////////////////////////////
					for (int best_response = 0; best_response < P; best_response++){
						cout << "\r" << "Best Response " << best_response + 1 << "/" << P << flush;
						// cout << "Best Response " << best_response + 1 << endl;
						rand_event = eventsRandomPermutation[best_response];
						// cout << "*********** 1" << endl;
						///////////////// Create the sampler for sampling nodes ////////////////////
						vector<int> users(g.getSize(), 0);
						for (int v = 1; v <= g.getSize(); v++)
							users[v - 1] = v;
						vector<float> RR_user_Prs(g.getSize(), 1.0);
						if (algorithm != "NA")
							create_dist(g, rand_event, RR_user_Prs);
						sampler<int> samp(users, RR_user_Prs);

						// cout << "*********** 2" << endl;
						vector<int> seeds;
						if (m == 0)
							g.rand_LT();
						// cout << "*********** 2" << endl;
						if (ris_alg == "DSSA")
							DSSA(g, g.getSize(), k, 1, seeds, epsilon, delta, samp, rand_event, m);
						else if (ris_alg == "SSA")
							SSA(g, g.getSize(), k, 1, seeds, epsilon, delta, samp, rand_event, m);
						// cout << "*********** 3" << endl;
						g.seedSets[rand_event].clear();
						for (unsigned int i = 0; i < k; ++i){
							g.seedSets[rand_event].insert(seeds[i]);
						}
						// cout << "*********** 4" << endl;
						//////////////////////////// Update AP //////////////////////////////////////
						if (algorithm != "NA"){
							for (int v = 0; v <= g.getSize(); ++v){
								g.AP[rand_event][v] = 0;
							}
							for(int iteration = 0; iteration < num_of_monte_carlo_iterations; iteration++)
								if (m == 0){
									g.rand_LT();
									g.BFS_LT(rand_event, g.AP);
								}
								else
									g.BFS_IC(rand_event, g.AP);
							for (int v = 1; v <= g.getSize(); v++)
								g.AP[rand_event][v] /= num_of_monte_carlo_iterations;
						}
						// cout << "*********** 5" << endl;
						// for (int i = 0; i < P; i++){
						// 	cout << i << ":";
						// 	for (int j = 1; j <= g.getSize(); j++){
						// 		cout << "  " << g.AP[i][j];
						// 	}
						// 	cout << endl;
						// }
					    

						if (cond == "seeds"){
							changedSeedsets += count_change_strat(g, usersLocations, rand_event);
							if (changedSeedsets != 0)
								sameSeedsets = false;
						}

						////////////////////////////// Basic Utility stats /////////////////////////
						if (exp3){
							double tot_sim = 0;
							double similarity_Sj = 0;
							double similarity_S_Sj = 0;
							for(int MC_sim_it = 0; MC_sim_it < monte_carlo_iterations_tsim; MC_sim_it++){
								g.calculate_sim_stats(P, rand_event, stats);
								tot_sim += stats[0];
								similarity_Sj += stats[1];
								similarity_S_Sj += stats[2];
							}
							tot_sim /= monte_carlo_iterations_tsim;
							similarity_Sj /= monte_carlo_iterations_tsim;
							similarity_S_Sj /= monte_carlo_iterations_tsim;
							total_DSS += similarity_Sj - (tot_sim - similarity_S_Sj);

							avg_tot_sim += tot_sim/P;
							avg_similarity_Sj += similarity_Sj/P;
							avg_similarity_S_Sj += similarity_S_Sj/P;
						}
					}
					gettimeofday(&end_round, NULL);
					tot_time += print_time_main(start_round, end_round)/1000;
		
					///////////////////////////// End of Best Response //////////////////////////////
					old_tsim = new_tsim;
					new_tsim = 0;
					#pragma omp parallel for
					for(int MC_sim_it = 0; MC_sim_it < monte_carlo_iterations_tsim; MC_sim_it++)
						new_tsim += g.calculate_tsim(P);
					new_tsim = new_tsim/(float)monte_carlo_iterations_tsim;

					if (exp3){
						final_sim.push_back(new_tsim);
						DSS.push_back(total_DSS);
						sim_Sj.push_back(avg_similarity_Sj);
						sim_S_S_Sj.push_back(avg_tot_sim - avg_similarity_S_Sj);
					}
					
					cout << "\r" << "Round " << round << ". (alg: " << algorithm << " |C| = " << P << " k = " << k << ") Total Similarity: " << new_tsim;
					if (new_tsim != 0 && old_tsim != 0)
						cout << " (" << setprecision(2) << fixed << ((new_tsim - old_tsim)/old_tsim)*100 << "%)" << endl;
					else
						cout << endl;

					if (cond == "seeds")
						cout <<  " Strategy changes : " << changedSeedsets << endl;
					
					if (cond == "seeds"){
						condition = !sameSeedsets;
					}
					else if (cond == "perc"){
						condition = abs((new_tsim - old_tsim)/old_tsim) > percent;
					}
					else if (cond == "rounds"){
						condition = (round < rounds && !(round == 1 && algorithm != "EG"));
					}

					if (exp2){
						round_time.push_back(tot_time);
						round_score.push_back(new_tsim);
					}
					if (rounds > 1){
						// map<string, vector<double>>::iterator it1 = data_tsim_vec.find(algorithm + "." + to_string(round));
						// map<string, vector<double>>::iterator it2 = data_time_vec.find(algorithm + "." + to_string(round));
						// map<string, vector<double>>::iterator it3 = data_memo_vec.find(algorithm + "." + to_string(round));
						data_tsim_vec[algorithm + "." + to_string(round)].push_back(new_tsim);
						data_time_vec[algorithm + "." + to_string(round)].push_back(tot_time);
						data_memo_vec[algorithm + "." + to_string(round)].push_back(getCurrentMemoryUsage());
					}
					// cout << "*********** 5" << endl;
				}
			}
			///////////////////////// Basic Stats ////////////////////////////
			if (exp1 or exp5){
				data_tsim[a/2][event] = new_tsim;
				data_time[a/2][event] = tot_time;
				data_memo[a/2][event] = getCurrentMemoryUsage();
			}


			//////////////////////// Round Stats /////////////////////////////
			if (exp2){
				vector<double> copy_time(round_time);
				vector<double> copy_score(round_score);
				round_data_time[a/2] = copy_time;
				round_data_score[a/2] = copy_score;
			}
			// for (int i = 0; i < P; i++){
			// 	for (int seed : g.seedSets[i]) vec_seeds[seed]++;
			// }
		}
		delete [] eventsRandomPermutation;
	}

	// for (int i = 0; i < P; i++){
	// 	cout << i << ":";
	// 	for (int seed : g.seedSets[i]){
	// 		cout << " " << seed;
	// 	}
	// 	cout << endl;
	// }
	///////////////////////////////////////////////////////////////////////////////////////////

	gettimeofday(&total_end, NULL);

	cout << "Time: " << print_time_main(total_start, total_end)/1000 << "s" << endl;
	cout << "Memory: " << getCurrentMemoryUsage() << " MB" << endl;

	double ** data[] = {data_tsim, data_time, data_memo};
	if (exp1){
		cout << "###############1" << endl;
		output_stats_k_and_p(output, data, variable, algorithms, reps, eventsArray, budgetsArray);
	}
	if (exp2){
		cout << "###############2" << endl;
		output_stats_q_and_t(output, algorithms, round_data_time, round_data_score, rounds);
	}
	if (exp3){
		cout << "###############3" << endl;
		output_basic_utility(output, sim_Sj, sim_S_S_Sj, DSS, final_sim, variable, reps, eventsArray, budgetsArray);
	}
	if (exp4){
		cout << "###############4" << endl;
		output_mult_rounds(output, data_tsim_vec, data_time_vec, data_memo_vec, variable, reps, eventsArray, budgetsArray, rounds);
	}
	if (exp5){
		cout << "###############5" << endl;
		output_clust(output, data, variable, algorithms);
	}

	cout << "Done" << endl;

	// g.example(usersLocations, eventsLocations, P);
	// delete usersLocations;
	// delete eventsLocations;
}
