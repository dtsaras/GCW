#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <list>
#include <map>
#include <set>
#include <cstdlib>
#include <sys/time.h>
#include <sstream>
#include <cmath>
#include <fstream>
#include <math.h>
#include <string>
#include <ctime>
#include <float.h>
#include <random>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <queue>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <ctime>
#include <set>


#include "Utilities.h"


using namespace std;

#define pi 3.14159265358979323846


//  REAL DATASET

#define LAT_SIZE 200 
#define LNG_SIZE 200

#define MIN_LAT 30.045322
#define MIN_LNG -98.259187
#define MAX_LAT 33.192731
#define MAX_LNG -96.515322


/// @brief The usual PI/180 constant
#define DEG_TO_RAD 0.017453292519943295769236907684886
/// @brief Earth's quatratic mean radius for WGS-84
#define EARTH_RADIUS_IN_KILOMETERS 6371


#define DELTA_LAT ((MAX_LAT - MIN_LAT)/ (LAT_SIZE-1))
#define DELTA_LNG ((MAX_LNG - MIN_LNG)/ (LNG_SIZE-1))

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_set_num_threads(int t) { return 1;}
inline omp_int_t omp_get_thread_num() { return 0;}
#endif

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

Utilities::Utilities(){}
/*
void Utilities::createFilesForIBFS(const char* costs, const char* socialgraph, int numOfUsers, int numOfSocialEdges, int numOfTerminals){
     int totalNodes = numOfUsers+numOfTerminals;
     
     for(int i=0; i < numOfTerminals; i++){
     
     } 

}
*/







int* Utilities::readBudgetsFromFile(int V, int P, int k, int run_program)
{   
    int* budgets;
    budgets = new int[P];
    char filename[256]; 
    
    
    
    if (run_program==0)
    snprintf(filename, sizeof(char) * 256, "/home/lefteris/stcpu2/entaflos/CIM_Create_Budgets/Test/P=%i/file%i.txt", P, k);


    if (run_program==1)
        snprintf(filename, sizeof(char) * 256, "/home/lefteris/stcpu2/entaflos/CIM_Create_Budgets/Gowalla/P=%i/file%i.txt", P, k);


    if (run_program==2)
        snprintf(filename, sizeof(char) * 256, "/home/lefteris/stcpu2/entaflos/CIM_Create_Budgets/Foursquare/P=%i/file%i.txt", P, k);

    
    ifstream fin(filename);
    string line;   

    if (! fin)
        cout << "Cannot open budgets file " << filename << endl;
        
    
    int p=0;
    while(getline(fin, line))
    { 
        if(p==P) break;
        fin.clear();
        std::istringstream iss(line);
        iss.clear();
        iss >> budgets[p];
        if (! fin.good())
        {
            cout << "fin is not good: event = " << p << endl;
            continue;
        }
        p++;
    }
   
   
   
   return budgets;
}


double Utilities::computeMinimumDistance(double lat1, double lon1, double lat2, double lon2, char unit){
    double R = 6371; // Radius of the earth in km
    double dLat = deg2rad(lat2-lat1);  // deg2rad below
    double dLon = deg2rad(lon2-lon1); 
    double a = sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2) * sin(dLon/2); 
    double c = 2 * atan2(sqrt(a), sqrt(1-a)); 
    double d = R * c; // Distance in km
    return d;
}


double Utilities::computeSimilarity(double v_lat, double v_lon, double p_lat, double p_lon, double max_dist)
{
    double sim = 0;
    double dist = 0; 
    double alpha = 1;
    dist = computeMinimumDistance(v_lat, v_lon, p_lat, p_lon, 'K');

    // sim = max_dist * exp(- alpha * dist);
    sim = 1 - dist/(double)max_dist; 
    
    return sim;
 }


double Utilities::deg2rad(double deg) {
	return (deg * pi / 180);
}


double Utilities::rad2deg(double rad) {
	return (rad * 180 / pi);
}




double Utilities::print_time(struct timeval &start, struct timeval &end)
{
    double usec;

    usec = (end.tv_sec*1000 + (end.tv_usec/1000)) - (start.tv_sec*1000 + start.tv_usec/1000);
    return usec;
}
  
