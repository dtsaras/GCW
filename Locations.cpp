#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <list>
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
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <ctime>
#include <set>
#include "Locations.h"
using namespace std;

Locations::Locations(const string file){
    if (file != ""){
        string line;
        ifstream fin(file);
        vector<int> *users_list;
        pair <double, double> coord;
        int id=0;
        double lat, lng;
        fin >> size;
        locations = new double *[size + 1];
        map< pair <double, double>, vector<int> > positions;
        for(int i = 0; i <= size; i++){
            locations[i] = new double[2];
        }
        
        while(fin >> id >> setprecision(12) >> lat >> setprecision(12) >> lng)
        { 
            locations[id][0] = lat;
            locations[id][1] = lng;
            coord = make_pair(lat, lng);
            users_list = &positions[coord];
            users_list->push_back(id);
        }
        fin.close();
    }
}

Locations::~Locations(){
    for (int i = 0; i <= size; i++){
        delete [] locations[i];
    }
    delete [] locations;
}

void Locations::get_coordinates(int id, pair <double, double> &coord){
    coord = make_pair(locations[id][0], locations[id][1]);
}

bool Locations::same_location(int id1, int id2){
    pair <double, double> coord1;
    vector<int> users_list;

    get_coordinates(id1, coord1);
    users_list = *get_users(coord1);
    for (int i = 0; i < users_list.size(); ++i){
        if (users_list[i] == id2)
            return true;
    }
    return false;

}

vector<int>* Locations::get_users(pair<double, double> &coord){
    return &positions[coord];
}


