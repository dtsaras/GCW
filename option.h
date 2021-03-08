#ifndef _OPTION_H_
#define _OPTION_H_
#include <sstream>
#include <map>
#include <string.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>

class OptionParser{

private:
	std::map<std::string, char *> paras;
	bool isValid;
public:
	OptionParser(int argc, char ** argv);
	// void setPara(std::string &algorithms, std::string &dataset, std::string &usersLocationsFile, std::string &eventsLocationsFile, double &epsilon, double &delta, int &budget, int &P, int &num_of_monte_carlo_iterations, int &monte_carlo_iterations_tsim, int &t, double &percent, int &reps, std::string &variable, std::string &cond, int &rounds, float &pr);
	void setPara(std::string &algorithms, std::string &dataset, std::string &usersLocationsFile, std::string &eventsLocationsFile, std::string &similarities_file, std::string &dt_type, float &epsilon, float &delta, int &budget, int &P, int &num_of_monte_carlo_iterations, int &monte_carlo_iterations_tsim, int &t, float &percent, int &reps, std::string &variable, std::string &cond, int &rounds, int &m, bool &exp1, bool &exp2, bool &exp3, bool &exp4, bool &exp5, std::string &fileName, std::string &ris_alg, bool &sim, bool &rand_edges, float &mu, float &sigma);
	char * getPara(const char * str);
	bool validCheck();
};


#endif
