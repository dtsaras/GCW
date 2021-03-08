#include "option.h"

OptionParser::OptionParser(int argc, char ** argv){
	for (int i = 1; i < argc; i += 2){
		if (argv[i][0] == '-'){
			if (i + 1 <= argc){
				std::string tmp(argv[i]);
				paras.insert(std::make_pair(tmp, argv[i + 1]));
			}else{
				isValid = false;
				return;
			}
		} else {
			isValid = false;
			return;
		}
	}
	isValid = true;
}

char * OptionParser::getPara(const char * str){
	std::string tmp(str);
	std::map<std::string, char *>::iterator it = paras.find(tmp);
	if (it == paras.end())
		return NULL;
	else
		return it->second;
}

bool OptionParser::validCheck(){
	return isValid;
}

// void OptionParser::setPara(std::string &algorithms, std::string &dataset, std::string &usersLocationsFile, std::string &eventsLocationsFile, double &epsilon, double &delta, int &budget, int &P, int &num_of_monte_carlo_iterations, int &monte_carlo_iterations_tsim, int &t, double &percent, int &reps, std::string &variable, std::string &cond, int &rounds, float &pr){
void OptionParser::setPara(std::string &algorithms, std::string &dataset, std::string &usersLocationsFile, std::string &eventsLocationsFile, std::string &similarities_file, std::string &dt_type, float &epsilon, float &delta, int &budget, int &P, int &num_of_monte_carlo_iterations, int &monte_carlo_iterations_tsim, int &t, float &percent, int &reps, std::string &variable, std::string &cond, int &rounds, int &m, bool &exp1, bool &exp2, bool &exp3, bool &exp4, bool &exp5, std::string &fileName, std::string &ris_alg, bool &sim, bool &rand_edges, float &mu, float &sigma){
	// Algorith (EG, CG, RG, NN, RD)
	char *  algorithm = getPara("-alg");
	if (algorithm == NULL){
		algorithms = "EG";
	}
	else if (strcmp(algorithm, "all") == 0){
		algorithms = "EGNANNRD";
	}
	else {
		algorithms = algorithm;
	}

	
	// Dataset // 0: Foursquare, 1: Gowalla
	char *  dt = getPara("-dt");
	char * prob = getPara("-pr");
	char * com = getPara("-clust");
	char * box = getPara("-box");
	std::string path = "./data/";
	if (dt == NULL || strcmp(dt, "Gowalla") == 0 || strcmp(dt, "G") == 0 || strcmp(dt, "g") == 0){
		path += "Gowalla/";
		if (com == NULL || strcmp(com, "original") == 0){
			path += "original/";
		}
		else if(strcmp(com, "kmeans") == 0 || strcmp(com, "Kmeans") == 0 || strcmp(com, "k-means") == 0){
			path += "kmeans/";
		}
		else if(strcmp(com, "dgcd") == 0 || strcmp(com, "dgcd") == 0){
			path += "dgcd/";
		}
		dataset = path + "gowalla.bin";
		usersLocationsFile = path + "users_locations_gowalla.txt";
		eventsLocationsFile = "./data/Gowalla/events_locations_gowalla.txt";
		if (box != NULL && (strcmp(box, "T") == 0 || strcmp(box, "True") == 0 || strcmp(box, "true") == 0)){
			eventsLocationsFile = "./data/Gowalla/box_events_locations_gowalla.txt";
		}
		dt_type = "gsn";
	}
	else if (strcmp(dt, "Foursquare") == 0 || strcmp(dt, "F") == 0 || strcmp(dt, "f") == 0){
		path += "Foursquare/";
		if (com == NULL || strcmp(com, "original") == 0){
			path += "original/";
		}
		else if(strcmp(com, "kmeans") == 0 || strcmp(com, "Kmeans") == 0 || strcmp(com, "k-means") == 0){
			path += "kmeans/";
		}
		else if(strcmp(com, "dgcd") == 0 || strcmp(com, "DGCD") == 0){
			path += "dgcd/";
		}
		else if(strcmp(com, "forrest_fire") == 0 || strcmp(com, "ff") == 0){
			path += "ff/";
		}
		dataset = path + "foursquare.bin";
		usersLocationsFile = path + "users_locations_foursquare.txt";
		eventsLocationsFile = "./data/Foursquare/events_locations_foursquare.txt";

		if (box != NULL && (strcmp(box, "T") == 0 || strcmp(box, "True") == 0 || strcmp(box, "true") == 0)){
            eventsLocationsFile = "./data/Foursquare/box_events_locations_foursquare.txt";
        }
        dt_type = "gsn";
	}
	else if (dt == NULL || strcmp(dt, "Flixster") == 0 || strcmp(dt, "X") == 0 || strcmp(dt, "x") == 0){
		path += "Flixster/";
		if (com == NULL || strcmp(com, "original") == 0){
			path += "original/";
		}
		dataset = path + "flixster.bin";
		similarities_file = path + "flixster_similarities.txt";
		dt_type = "social";
		usersLocationsFile = "";
		eventsLocationsFile = "";
	}
	else if (dt == NULL || strcmp(dt, "Lastfm") == 0 || strcmp(dt, "L") == 0 || strcmp(dt, "l") == 0){
		path += "Lastfm/";
		if (com == NULL || strcmp(com, "original") == 0){
			path += "original/";
		}
		dataset = path + "lastfm.bin";
		similarities_file = path + "lastfm_similarities.txt";
		dt_type = "social";
		usersLocationsFile = "";
		eventsLocationsFile = "";
	}
	else if (strcmp(dt, "Dallas") == 0){
		path += "Dallas/";
		dataset = "./data/Dallas/dallas.bin";
		usersLocationsFile = "./data/Dallas/dallas_users_locations.txt";
		eventsLocationsFile = "./data/Dallas/dallas_events_locations.txt";
		dt_type = "gsn";
	}
	else if (strcmp(dt, "test") == 0){
		path += "Test/";
		dataset = "./data/Test/test.bin";
		usersLocationsFile = "./data/Test/test_users_locations.txt";
		eventsLocationsFile = "./data/Test/test_events_locations.txt";
		dt_type = "gsn";
	}

	// Epsilon
	char *  eps = getPara("-e");
	if (eps == NULL){
		epsilon = 0.1;
	}
	else{
		epsilon = std::stof(eps);
	}

	// Delta
	char *  del = getPara("-d");
	if (del == NULL){
		delta = 0.01;
	}
	else{
		delta = atof(del);
	}

	// Budget
	char *  bud = getPara("-k");
	if (bud == NULL){
		budget = 10;
	}
	else{
		budget = atoi(bud);
	}

	// Number of Events
	char *  ev = getPara("-c");
	if (ev == NULL){
		P = 16;
	}
	else{
		P = atoi(ev);
	}

	// Monte Carlo iterations for updating AP tables
	char *  MC_its = getPara("-mc1");
	if (MC_its == NULL){
		num_of_monte_carlo_iterations = 100;
	}
	else{
		num_of_monte_carlo_iterations = atoi(MC_its);
	}

	// Monte Carlo iterations for calculating total similarity
	char * MC_its_tsim = getPara("-mc2");
	if (MC_its_tsim == NULL){
		monte_carlo_iterations_tsim = 100;
	}
	else{
		monte_carlo_iterations_tsim = atoi(MC_its_tsim);
	}

	// Number of threads
	char * num_threads = getPara("-t");
	if (num_threads == NULL){
		t = 1;
	}
	else{
		t = atoi(num_threads);
	}

	// Quality Percent
	char * perc = getPara("-p");
	if (perc == NULL){
		percent = 0.03;
	}
	else{
		percent = atof(perc);
	}

	// Variable in the experiments
	char * var = getPara("-v");
	if (var == NULL){
		variable = "";
	}
	else{
		variable = var;
	}

	// Repetitions (number of different experments)
	char * repetitions = getPara("-r");
	if (repetitions == NULL){
		if (variable == "")
			reps = 1;
		else if (variable == "c")
			reps = 5;
		else if (variable == "k")
			reps = 5;
	}
	else{
		reps = atoi(repetitions);
	}

	// Which condition to stop the reps
	char * c = getPara("-cnd");
	if (c == NULL){
		cond = "rounds";
	}
	else{
		cond = c;
	}

	// If the condition is rounds, how many till it stops
	char * rnd = getPara("-rnd");
	if (rnd == NULL){
		rounds = 1;
	}
	else{
		rounds = atoi(rnd);
	}

	// Select model (IC or LT)
	char * mod = getPara("-m");
	if (mod == NULL || atoi(mod) == 1){
		m = 1;
	}
	else{
		m = 0;
	}

	// Select what experiment to run (exp1: quality (similarity, time, memory), exp2: similarity in each round, exp3: basic utility)
	char * exper = getPara("-exp");
	if (exper == NULL){
		exp1 = true;
		exp2 = false;
		exp3 = false;
		exp4 = false;
		exp5 = false;
	}
	else if (strcmp(exper, "1") == 0){
		exp1 = true;
		exp2 = false;
		exp3 = false;
		exp4 = false;
		exp5 = false;
	}
	else if (strcmp(exper, "2") == 0){
		exp1 = false;
		exp2 = true;
		exp3 = false;
		exp4 = false;
		exp5 = false;
	}
	else if (strcmp(exper, "3") == 0){
		exp1 = false;
		exp2 = false;
		exp3 = true;
		exp4 = false;
		exp5 = false;
	}
	else if (strcmp(exper, "4") == 0){
		exp1 = false;
		exp2 = false;
		exp3 = false;
		exp4 = true;
		exp5 = false;
	}
	else if (strcmp(exper, "5") == 0){
		exp1 = false;
		exp2 = false;
		exp3 = false;
		exp4 = false;
		exp5 = true;
	}


	
	// else{
	// 	pr = atof(prob);
	// }


	char *  ris_algorithm = getPara("-ris_alg");
	if (ris_algorithm == NULL){
		ris_alg = "DSSA";
	}
	else {
		ris_alg = ris_algorithm;
	}

	char * similarity = getPara("-sim");
	if (similarity == NULL){
		sim = true;
	}
	else if (strcmp(similarity, "no") == 0 || strcmp(similarity, "n") == 0 || strcmp(similarity, "no_sim") == 0){
		sim = false;
	}
	else if (strcmp(similarity, "yes") == 0 || strcmp(similarity, "y") == 0 || strcmp(similarity, "yes_sim") == 0){
		sim = true;
	}

	char * re = getPara("-rand_edges");
	if (re == NULL){
		rand_edges = false;
	}
	else if (strcmp(re, "t") == 0 || strcmp(re, "T") == 0 || strcmp(re, "true") == 0){
		rand_edges = true;
	}
	else if (strcmp(re, "f") == 0 || strcmp(re, "F") == 0 || strcmp(re, "false") == 0){
		rand_edges = false;
	}

	char * mean = getPara("-mu");
	if (mean == NULL){
		mu = 0.1;
	}
	else{
		mu = atof(mean);
	}

	char * sig = getPara("-sigma");
	if (sig == NULL){
		sigma = 0.1;
	}
	else{
		sigma = atof(sig);
	}

	// Defaults the name of the file
	char * outFile = getPara("-o");
	if (outFile != NULL && strcmp(outFile, "def") == 0){
		fileName += "results/";
		if (com == NULL || strcmp(com, "original") == 0){
			fileName += "original/";
		}
		else if(strcmp(com, "kmeans") == 0 || strcmp(com, "Kmeans") == 0 || strcmp(com, "k-means") == 0){
			fileName += "kmeans/";
		}
		else if(strcmp(com, "dgcd") == 0 || strcmp(com, "dgcd") == 0){
			fileName += "dgcd/";
		}
		else if(strcmp(com, "forrest_fire") == 0 || strcmp(com, "ff") == 0){
			fileName += "ff/";
		}
		if (dt == NULL || strcmp(dt, "Gowalla") == 0 || strcmp(dt, "G") == 0 || strcmp(dt, "g") == 0)
			fileName += "g";
		else if (strcmp(dt, "Foursquare") == 0 || strcmp(dt, "F") == 0 || strcmp(dt, "f") == 0)
			fileName += "f";
		else
			fileName += dt;
		if (m == 1)
			fileName += "_ic";
		else
			fileName += "_lt";
		if (sim)
			fileName += "_ysim";
		else
			fileName += "_nsim";
		if (rand_edges)
			fileName += "_rand_" + (std::string) mean + "_" + (std::string) sig;
		if (var != NULL)
			fileName += "_" + variable;
		if (exp1)
			fileName += ".dat";
		else if (exp2)
			fileName += "_rounds.dat";
		else if (exp3)
			fileName += "_utility.dat";
		else if (exp4)
			fileName += "_mult_rounds.dat";
		else if (exp5)
			fileName = "tmp_clust_exper.dat";
		std::cout << fileName << std::endl;
	}
	else if (outFile != NULL){
		fileName = std::string(outFile);
	}

	// // Output file
	// std::ofstream file;
	// char * outFile = getPara("-o");
	// if (outFile == NULL){
	// 	output = &std::cout;
	// }
	// else{
	// 	file.open(outFile);
	// 	output = &file;
	// }
}
