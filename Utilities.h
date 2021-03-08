#include <set>

using namespace std;
class Utilities{

    public:	
        Utilities();
        int* readBudgetsFromFile(int V, int P, int k, int run_program);
        double print_time(struct timeval &start, struct timeval &end);
        double computeMinimumDistance(double lat1, double lon1, double lat2, double lon2, char unit);
        double computeSimilarity(double v_lat, double v_lon, double p_lat, double p_lon, double max_dist);
        double deg2rad(double deg);
        double rad2deg(double rad);
        

};









      
