#include <utility>
#include <map>
#include <vector>
#include <string>

using namespace std;

class Locations{
	public:
		double** locations;
		int size;
		map<pair<double, double>, vector<int>> positions;

		Locations(const string file);
		~Locations();
		void get_coordinates(int id, pair <double, double> &coord);
		vector<int>* get_users(pair<double, double> &coord);
		bool same_location(int id1, int id2);
};