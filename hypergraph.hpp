#ifndef _HYPERGRAPH_H_
#define _HYPERGRAPH_H_

#include "rwgraph.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include "mappedHeap.hpp"
#include "HeapData.hpp"

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_set_num_threads(int t) { return 1;}
inline omp_int_t omp_get_thread_num() { return 0;}
#endif

using namespace std;

/*
* building the hypergraph procedure which generates hyperedges following LT model
*/
long long addHyperedge(Graph & g, HyperGraph & hg, int t, long long num, bool lt, sampler<int> &samp, int p)
{
	int numNodes = g.getSize();

	omp_set_num_threads(t);
	
	long long iter = 0;
	int c = 100;
	// cout << "@@@@@@@@@@ 1" << endl;
    #pragma omp parallel
	{
		vector<int> visit_mark(numNodes+1,0);
		vector<bool> visit(numNodes+1,false);
		vector<unsigned int> link;
		// cout << "@@@@@@@@@@ 2" << endl;     
		if (lt == 0){
			while (iter < num){
				// cout << "@@@@@@@@@@ 3" << endl;
            	for (int i = 0; i < c; ++i){
					vector<int> he;
           	        hg.pollingLT1(g,visit,visit_mark, samp, p);
               	}
               	// cout << "@@@@@@@@@@ 4" << endl;
				#pragma omp atomic
				iter += c;
			}
		} else {
			while (iter < num){
                for (int i = 0; i < c; ++i){
                    vector<int> he;
                    hg.pollingIC1(g,visit,visit_mark, samp, p);
                }
                #pragma omp atomic
                iter += c;
        	}
		}
	}
	// cout << "@@@@@@@@@@ 5" << endl;
	hg.updateDeg();
	// cout << "@@@@@@@@@@ 6" << endl;
	return hg.getNumEdge();
}

/*
* find seed nodes procedure using greedy algorithm
*/
void buildSeedSet(HyperGraph & hg, vector<int> & seeds, unsigned int n, int k, vector<double> &degree)
{	
	long long i;
	unsigned int j,l,maxInd;
	vector<int> e, nList;

	vector<int> nodeDegree(n,0);
	vector<int> indx(n,0);
	for (j = 0; j < n; ++j){
		indx[j] = j;
		nodeDegree[j] = hg.getNode(j+1).size();
	}

	InfCost<int> hd(&nodeDegree[0]);
	MappedHeap<InfCost<int> > heap(indx,hd);
	long long numEdge = hg.getNumEdge();

	// check if an edge is removed
	vector<bool> edgeMark(numEdge, false);
	vector<bool> nodeMark(n+1, true);
	
	double totalCost = 0;

	i=1;
	// building each seed at a time
	while(totalCost < k && !heap.empty()){
		maxInd = heap.pop()+1;
		nodeMark[maxInd] = false;

		totalCost++;

		e = hg.getNode(maxInd);
		
		degree[i] = degree[i-1]+nodeDegree[maxInd-1];
		
		seeds.push_back(maxInd);
		for (j = 0; j < e.size(); ++j){
			if (edgeMark[e[j]]){
				continue;
			}
	
			nList = hg.getEdge(e[j]);
			for (l = 0; l < nList.size(); ++l){
				nodeDegree[nList[l]-1]--;
				if (nodeMark[nList[l]]){
					heap.heapify(nList[l]-1);
				}
			}
			edgeMark[e[j]] = true;
		}
		i++;
	}

	vector<int>().swap(nodeDegree);
	vector<int>().swap(e);
	vector<int>().swap(nList);
	vector<int>().swap(indx);
	vector<bool>().swap(edgeMark);
}

#endif
