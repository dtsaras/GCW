/*
* Functionality: convert from a graph file in weighted edge list to a binary file
* Syntax:
	./el2bin <graph file input> <binary graph output>

* The graph file input must follow the following format:
	<number of nodes> <number of edges>
	<first node of edge 1> <second node of edge 1> <weight of edge 1>
	...
	<first node of the last edge> <second node of the last edge> <weight of the last edge>

* The binary graph output will be used in our SSA/DSSA algorithms for fast reading

* Author: Hung T. Nguyen (hungnt@vcu.edu)
*/

#include <cstdio>
#include <fstream>
#include <random>
#include "GLib.hpp"
#include <cmath>
#include <cstring>

using namespace std;

int main(int argc, char ** argv)
{
	ifstream in(argv[1]);
	srand(time(NULL));
	int n,u,v;
	long long m;
	float w;
	in >> n >> m;
	printf("%d %lld\n", n, m);
	vector<int> in_degree(n+1,0);
	vector<int> out_degree(n+1,0);
	vector<vector<int> > eList(n+1);
	vector<vector<int> > outList(n+1);
	vector<vector<float> > weight(n+1);
	vector<vector<float> > outWeight(n+1);
	vector<float> weightR(n+1,0);

	printf("Reading the graph!\n");

	for (long long i = 0; i < m; ++i){
		in >> u >> v >> w;
		out_degree[u]++;
		in_degree[v]++;
		eList[v].push_back(u);
		outList[u].push_back(v);
		outWeight[u].push_back(w);
		weight[v].push_back(w);
		weightR[u] += 1;
	}
	
	in.close();

    vector<size_t> idx(n);

	FILE * pFile;
	pFile = fopen(argv[2],"wb");
	fwrite(&n, sizeof(int), 1, pFile);
	fwrite(&m, sizeof(long long), 1, pFile);

    for (int i = 0; i < n; ++i){
		idx[i] = i;
	}
	vector<int> inv_idx(n);
	for (int i = 0; i < n; ++i){
		inv_idx[idx[i]]	= i;
	}
	
	vector<int> iTmp(n);
	
	for (int i = 0; i < n; ++i){
		iTmp[i] = in_degree[idx[i]+1];
	}
	
	// Write node in-degrees for DSSA
	fwrite(&iTmp[0], sizeof(int), n, pFile);

	// Write node out-degrees
	fwrite(&out_degree[1], sizeof(int), n, pFile);

	
	for (int i = 1; i <= n; ++i){
		// Write neighbors
		for (unsigned int j = 0; j < eList[idx[i-1]+1].size(); ++j){
			iTmp[j] = inv_idx[eList[idx[i-1]+1][j]-1]+1;
		}
		fwrite(&iTmp[0], sizeof(int), eList[idx[i-1]+1].size(), pFile);
		fwrite(&outList[i][0], sizeof(int), outList[i].size(), pFile);
	}

	for (int i = 1; i <= n; ++i){
		// Write weights
        fwrite(&weight[idx[i-1] + 1][0], sizeof(float), weight[idx[i-1]+1].size(), pFile);
        fwrite(&outWeight[i][0], sizeof(float), outWeight[i].size(), pFile);

    }

	fclose(pFile);
	printf("Done!\n");
	return 1;
}
