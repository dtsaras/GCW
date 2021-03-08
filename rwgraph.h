#ifndef _RWGRAPH_H
#define _RWGRAPH_H
#include <vector>
#include <random>
#include "sfmt/SFMT.h"
#include "sampler.h"
#include "Utilities.h"
#include "BinaryHeap.h"
#include "Locations.h"
#include <unordered_map>
#include <unordered_set>

typedef uint32_t UI;
typedef uint64_t ULL;


// Graph class defines fundamental operators on graph
class Graph
{
	friend class HyperGraph;
	private:
		sfmt_t sfmtSeed;

		// number of nodes
		unsigned int numNodes;
		// number of edges
		unsigned int numEdges;
		int numEvents;
		// adjacency list
		std::vector<std::vector<int> > adjList; // O(E)
		// std::vector<std::vector<int> > invAdjList; // O(E)
		std::vector<std::vector<float> > probabilities; // O(E)
		std::vector<std::vector<float> > invProbabilities; // O(E)
		std::vector<int> node_deg;  // O(V)
		// std::vector<int> inv_node_deg;  // O(V)
		// std::vector<std::vector<float> > weights;  // O(E)
		// std::vector<std::vector<UI> > weights;  // O(E)
		std::vector<bool> dart;
		int model;
	
	public:
		UI UI_MAX = 4294967295U;
		ULL ULL_MAX = 18446744073709551615ULL;
		// std::vector<UI> LT_node; // Thresholds for each node O(V)

		std::vector<int> inv_node_deg;  // O(V)
		std::vector<std::vector<int> > invAdjList; // O(E)
		std::vector<std::vector<float> > weights;  // O(E)

		std::vector<float> LT_node; // Thresholds for each node O(V)
		std::set<int>* seedSets; // Seedsets for each event O(P*k)
		std::set<int>* seedSets_before; // Seedsets for each event at previous round O(P*k)
		std::vector<std::vector<double> > similarities; // list of similarities between users O(V*P)
		std::multimap<double, int>* A; // O(V*P)
		// std::set<int>** A;  // A[v][p] set of events p' such that sim(v,p')>sim(v,p) // The events that have better weight than the event that plays best response for a user v
		std::vector<double> tsim_V_p; // sum of the similarities of all the users in V for the current best response event p O(P)
		// double* IP; // O(V)
		float** AP; // O(V*P)
		bool sim; 
		// Constructor
		Graph(std::string dataset, int m, bool s, bool rand_edges, float mu, float sigma);
		// Destructor
		~Graph();
		// Adjusts attributes for the proper number of competitors
		void change_events(int &P);

		void set_events(int &P);

		void create_AV_sets(Locations* usersLocations, Locations* eventsLocations, int P, Utilities util);

		void read_similarities(string filename, int P);

		// get a vector of neighbours of node u
		const std::vector<int> & operator [] (int u) const;
		// return weights of neighbours of node u
		// const std::vector<UI> & getWeight(int u) const;
		const std::vector<float> & getWeight(int u) const;

		// get a vector of neighbours of node u
		const std::vector<int> & operator [] (int u);
		// return weights of neighbours of node u
		// const std::vector<UI> & getWeight(int u);
		const std::vector<float> & getWeight(int u);

		void set_nodes(unsigned n);

		float get_AP_IC(int u, int v, int idx, int p);
		float get_inv_AP_IC(int u, int v, int idx, int p);
		float get_AP_LT(float active_prob, int u, int p);
		// get degree of node u
		int getDegree(int u) const;
		// get size of the graph
		int getSize() const;
		// get number of edges
		int getEdge() const;
		// read graph from a file under IC
		void readGraph_IC(const std::string filename, bool rand_edges, float m, float s);
		// read graoh from a file under LT
		void readGraph_LT(const std::string filename, bool rand_edges, float m, float s);
		// randomize the probabilities of the edges
		void rand_graph(float m, float s);
		// write the graph to file
		void writeToFile(const char * filename);
		// build seed set
		void initial_seedset(std::string algorithm, Utilities &util, Locations * usersLocations, Locations * eventsLocations, int P, int k);
		// The furthers distance between a user and a competitor
		void rand_LT();
		double find_max_dist(int P,  Locations* usersLocations, Locations* eventsLocations, Utilities util);
		void BFS_IC(int p, float** APc);
		void BFS_LT(int p, float** APc);
		double calculate_tsim(int P);
		void calculate_sim_stats(int P, int ev, double * stats);
		void example(Locations* usersLocations, Locations* eventsLocations, int P);
		
};

class HyperGraph
{
	private:
		// store the edges that a node is incident to
		std::vector<std::vector<int> > node_edge;
		// store hyperedges
		std::vector<std::vector<int> > edge_node;
		unsigned int curEdge;
		unsigned int maxDegree;
		unsigned int numNodes;
		sfmt_t sfmtSeed;
		// randomly selects a node in a list according to the weight distribution
		// inline int randIndex_bin(Graph &g, const vector<UI> &w, unsigned int si, int cur, int p);
		// inline int randIndex_lin(Graph &g, const vector<UI> &w, unsigned int si, int cur, int p);
		// inline int randIndex_dart(const std::vector<UI> &w, unsigned int si);
		inline int randIndex_bin(Graph &g, const vector<float> &w, unsigned int si, int cur, int p);
		inline int randIndex_lin(Graph &g, const vector<float> &w, unsigned int si, int cur, int p);
		inline int randIndex_dart(const std::vector<float> &w, unsigned int si);
	public:
		HyperGraph(unsigned int n);
		void updateDeg();
		void updateEdge();
		void addEdge(std::vector<int> & edge);
		void addEdgeD(std::vector<int> & edge);
		int getMaxDegree();
		const std::vector<int> & getEdge(int e) const;
		const std::vector<int> & getEdge(int e);
		const std::vector<int> & getNode(int n) const;
		const std::vector<int> & getNode(int n);
        int getNumEdge() const;
		void clearEdges();

		// different polling processes for different purposes: model and returned output
		void pollingLT1(Graph &g, std::vector<bool> &visit, std::vector<int> &mark_visit, sampler<int> &samp, int p);
		bool pollingLT2(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &visit_mark, sampler<int> &samp, int p);
		bool pollingLT(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit, sampler<int> &samp, int p);
		void pollingIC1(Graph &g, std::vector<bool> &visit, std::vector<int> &visit_mark, sampler<int> &samp, int p);
		bool pollingIC2(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &visit_mark, sampler<int> &samp, int p);
		bool pollingIC(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &visit_mark, sampler<int> &samp, int p);
};

float getCurrentMemoryUsage();

#endif
