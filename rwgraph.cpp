#include "rwgraph.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <random>
#include <cstring> 

using namespace std;

Graph::Graph(string dataset, int m, bool s, bool rand_edges, float mu, float sigma){
    model = m;
    sim = s;
    if (model == 0){
        readGraph_LT(dataset, rand_edges, mu, sigma);
        LT_node.resize(getSize() + 1);
    }
    else
        readGraph_IC(dataset, rand_edges, mu, sigma);
    sfmt_init_gen_rand(&sfmtSeed, rand());
    A = new multimap<double, int>[getSize() + 1];
    similarities.resize(getSize() + 1);
    numEvents = 0;
}

Graph::~Graph(){
    delete [] seedSets;
    delete [] seedSets_before;
    for (int p = 0; p < numEvents; ++p){
        delete [] AP[p];
    }
    delete [] AP;
}

void Graph::set_nodes(unsigned n){
    numNodes = n;
}

void Graph::set_events(int &P){
    numEvents = P;
    seedSets = new set<int>[P];
    seedSets_before = new set<int>[P];
    tsim_V_p.clear();
    tsim_V_p.resize(P, 0);
    for (int i = 1; i < getSize() + 1; i++){
        similarities[i].resize(P);
    }

    AP = new float*[P]; // Awareness probability of node from competitor
    for (int p = 0; p < P; p++)
        AP[p] = new float[getSize() + 1];
}

void Graph::change_events(int &P){
    delete [] seedSets;
    delete [] seedSets_before;
    for (int p = 0; p < numEvents; ++p){
        delete [] AP[p];
    }
    delete [] AP;
    set_events(P);
}

const vector<int> & Graph::operator [] (int u) const
{
	return invAdjList[u];
}


const vector<int> & Graph::operator [] (int u)
{
	return invAdjList[u];
}


// const vector<UI> & Graph::getWeight (int u) const
const vector<float> & Graph::getWeight (int u) const
{
    return weights[u];
}

// const vector<UI> & Graph::getWeight (int u)
const vector<float> & Graph::getWeight (int u)
{
    return weights[u];
}

/*
* get degree of node u
*/
int Graph::getDegree(int u) const
{
	return inv_node_deg[u];
}

/*
* get the number of nodes
*/
int Graph::getSize() const
{
	return numNodes;
}

/*
* get the number of edges
*/
int Graph::getEdge() const
{
	return numEdges;
}

float Graph::get_AP_IC(int u, int v, int idx, int p){
    float prob = probabilities[u][idx]; // (u) ---probabilities[u][idx]---> (v) (idx of (v) in adjacency list of (u))
    if (sim)
        prob *= similarities[v][p]; // similarity of the node (v) for the product p
    return prob;
}

float Graph::get_inv_AP_IC(int u, int v, int idx, int p){
    float prob = weights[u][idx]; // (u) ---probabilities[u][idx]---> (v) (idx of (v) in adjacency list of (u))
    if (sim)
        prob *= similarities[u][p]; // similarity of the node (u) for the product p
    return prob;
}

float Graph::get_AP_LT(float active_prob, int u, int p){
    if (sim)
        return active_prob * similarities[u][p];
    else
        return active_prob;
}

void Graph::create_AV_sets(Locations* usersLocations, Locations* eventsLocations, int P, Utilities util){
    double max_dist = find_max_dist(P, usersLocations, eventsLocations, util);
    for (int v = 1; v <= getSize(); v++){
        A[v].clear();
        for (int p = 0; p < P; p++){
            similarities[v][p] = util.computeSimilarity(usersLocations->locations[v][0], usersLocations->locations[v][1], eventsLocations->locations[p][0], eventsLocations->locations[p][1], max_dist);
            A[v].insert(pair<double, int>(-similarities[v][p], p));
            tsim_V_p[p] += similarities[v][p];
        }
    }
}

void Graph::read_similarities(string filename, int P){
    ifstream infile(filename);
    int num_nodes;
    int max_P;
    double tmp;
    stringstream tmp_s;
    int node_id;
    string a;
    infile >> num_nodes >> max_P;
    for (int v = 1; v <= getSize(); v++){
        A[v].clear();
        // infile.getline(tmp_s, 256, "\n");
        infile >> node_id;
        for (int p = 0; p < max_P; p++){
            if (p < P){
                infile >> setprecision(16) >> tmp;
                similarities[v][p] = tmp;
                A[v].insert(pair<double, int>(-similarities[v][p], p));
                tsim_V_p[p] += similarities[v][p];
            }
            else{
                infile >> setprecision(16) >> tmp;
            }
        }
        // infile.seekg(sizeof(double) * (max_P - P) * 2, infile.cur);
        // infile.ignore(max_P - P, "\n");
        // for (int i = 0; i < max_P - P; i++)
        //     infile.ignore(max_P - P, "\n");
    }
}

void Graph::rand_graph(float m, float s){
    default_random_engine generator;
    normal_distribution<float> distribution(m, s);
    float rand;
    int v;
    for (int i = 1; i <= numNodes; i++){
        probabilities[i].resize(node_deg[i]);
        for (int j = 0; j < node_deg[i]; j++){
            rand = distribution(generator);
            if (rand < 0)
                rand = 0;
            else if (rand > 1)
                rand = 1;
            probabilities[i][j] = rand;
            v = adjList[i][j];
            weights[v].push_back(rand);
        }
    }

    for (int i = 1; i <= numNodes; i++){
        weights[i].resize(inv_node_deg[i]);
    }
}

void Graph::readGraph_IC(const string filename, bool rand_edges, float m, float s){
    ifstream infile(filename, ios::in|ios::binary|ios::ate);
    int v;
    long long u;
    float prob;
    infile.seekg(0, ios::beg);

    // Read Num Nodes + Num Edges
    infile.read(reinterpret_cast<char*>(&v), sizeof(int));
    infile.read(reinterpret_cast<char*>(&u), sizeof(long long));

    numNodes = v;
    numEdges = u;

    // Read Node degrees
    inv_node_deg.resize(numNodes + 1);
    node_deg.resize(numNodes + 1);
    for (int i = 1; i <= numNodes; i++){
        infile.read(reinterpret_cast<char*>(&v), sizeof(int));
        inv_node_deg[i] = v;
    }
    for (int i = 1; i <= numNodes; i++){
        infile.read(reinterpret_cast<char*>(&v), sizeof(int));
        node_deg[i] = v;
    }


    // Read Neighbors
    adjList.resize(numNodes + 1);
    invAdjList.resize(numNodes + 1);
    for (int i = 1; i <= numNodes; i++){
        adjList[i].resize(node_deg[i]);
        invAdjList[i].resize(inv_node_deg[i]);
        // Incoming neighbors
        for (int j = 0; j < inv_node_deg[i]; ++j){
            infile.read(reinterpret_cast<char*>(&v), sizeof(int));
            invAdjList[i][j] = v;
        }
        // Outgoing neighbors
        for (int j = 0; j < node_deg[i]; ++j){
            infile.read(reinterpret_cast<char*>(&v), sizeof(int));
            adjList[i][j] = v;
        }
    }

    weights.resize(numNodes + 1);
    probabilities.resize(numNodes + 1);

    if (rand_edges)
        rand_graph(m, s);
    else{
        float tmp;
        // vector<UI> b;
        // invProbabilities.resize(numNodes + 1);
        for (unsigned int i = 1; i <= numNodes; ++i){
            weights[i].resize(inv_node_deg[i]);
            // invProbabilities[i].resize(inv_node_deg[i]);
            probabilities[i].resize(node_deg[i]);

            // vector<UI> tmp1(inv_node_deg[i] + 1, 0);
            for(int j = 0; j < inv_node_deg[i]; ++j){
                infile.read(reinterpret_cast<char*>(&tmp), sizeof(float));
                weights[i][j] = tmp;
                // weights[i][j] = tmp*UI_MAX;
                // invProbabilities[i][j] = tmp;
                // if (weights[i][j] <= 0)
                //     weights[i][j] = UI_MAX;
            }

            for(int j = 0; j < node_deg[i]; ++j){
                infile.read(reinterpret_cast<char*>(&tmp), sizeof(float));
                probabilities[i][j] = tmp;
            }
        }

    }
    

    infile.close();
}

void Graph::readGraph_LT(const string filename, bool rand_edges, float m, float s){
    ifstream infile(filename, ios::in|ios::binary|ios::ate);
    int v;
    long long u;
    float prob;
    infile.seekg(0, ios::beg);

    // Read Num Nodes + Num Edges
    infile.read(reinterpret_cast<char*>(&v), sizeof(int));
    infile.read(reinterpret_cast<char*>(&u), sizeof(long long));

    numNodes = v;
    numEdges = u;

    // Read Node degrees
    inv_node_deg.resize(numNodes + 1);
    node_deg.resize(numNodes + 1);
    for (int i = 1; i <= numNodes; i++){
        infile.read(reinterpret_cast<char*>(&v), sizeof(int));
        inv_node_deg[i] = v;
    }
    for (int i = 1; i <= numNodes; i++){
        infile.read(reinterpret_cast<char*>(&v), sizeof(int));
        node_deg[i] = v;
    }


    // Read Neighbors
    adjList.resize(numNodes + 1);
    invAdjList.resize(numNodes + 1);
    for (int i = 1; i <= numNodes; i++){
        adjList[i].resize(node_deg[i]);
        invAdjList[i].resize(inv_node_deg[i]);
        // Incoming neighbors
        for (int j = 0; j < inv_node_deg[i]; ++j){
            infile.read(reinterpret_cast<char*>(&v), sizeof(int));
            invAdjList[i][j] = v;
        }
        // Outgoing neighbors
        for (int j = 0; j < node_deg[i]; ++j){
            infile.read(reinterpret_cast<char*>(&v), sizeof(int));
            adjList[i][j] = v;
        }
    }
    weights.resize(numNodes + 1);
    probabilities.resize(numNodes + 1);

    if (rand_edges)
        rand_graph(m, s);
    else{
        float tmp;
        // vector<UI> b;
        // invProbabilities.resize(numNodes + 1);
        for (unsigned int i = 1; i <= numNodes; ++i){
            weights[i].resize(inv_node_deg[i]);
            // invProbabilities[i].resize(inv_node_deg[i]);
            probabilities[i].resize(node_deg[i]);

            // vector<UI> tmp1(inv_node_deg[i] + 1, 0);
            for(int j = 0; j < inv_node_deg[i]; ++j){
                infile.read(reinterpret_cast<char*>(&tmp), sizeof(float));
                weights[i][j] = tmp;
                // weights[i][j] = tmp*UI_MAX;
                // invProbabilities[i][j] = tmp;
                // if (weights[i][j] <= 0)
                //     weights[i][j] = UI_MAX;
            }

            for(int j = 0; j < node_deg[i]; ++j){
                infile.read(reinterpret_cast<char*>(&tmp), sizeof(float));
                probabilities[i][j] = tmp;
            }
        }

    }

    float total_weight;
    for (int i = 1; i <= numNodes; ++i){
        total_weight = 0;
        for (int j = 0; j < inv_node_deg[i]; ++j){
            total_weight += weights[i][j];
            if (total_weight >= 1){
                weights[i][j] = 1;
            }
            else {
                weights[i][j] = total_weight;
            }
        }
    }


    // size_t probs_pos = infile.tellg();
    // vector<UI> b;
    // weights.push_back(b);
    // for (unsigned int i = 1; i <= numNodes; ++i){
    //     vector<float> tmp(node_deg[i] + 1, 0);
    //     vector<UI> tmp1(node_deg[i] + 1, 0);
    //     for(int j = 1;j < node_deg[i] + 1; ++j){
    //         infile.read(reinterpret_cast<char*>(&tmp[j]), sizeof(float));
    //     }
    //     for(int j = 1;j < node_deg[i] + 1; ++j){
    //         tmp[j] += tmp[j-1];
    //         if (tmp[j] >= 1){
    //             tmp1[j] = UI_MAX;
    //         } else {
    //             tmp1[j] = tmp[j]*UI_MAX;
    //         }
    //     }
    //     weights.push_back(tmp1);
    //     node_deg[i]++;
    // }

    // // Read Probabilities
    // infile.seekg(probs_pos);
    // probabilities.resize(numNodes + 1);
    // for (int i = 1; i <= numNodes; i++){
    //     // probabilities[i].resize(node_deg[i]);
    //     for (int j = 0; j < node_deg[i] - 1; ++j){
    //         infile.read(reinterpret_cast<char*>(&prob), sizeof(float));
    //         // probabilities[i][j] = prob;
    //         probabilities[adjList[i][j]].push_back(prob);
    //     }
    // }
    // for (int i = 1; i <= numNodes; i++)
    //     probabilities[i].resize(node_deg[i] - 1);

    infile.close();

    // for (int i = 1; i <= numNodes; i++){
    //     cout << "In deg: " << inv_node_deg[i] << " w_size: " << weights[i].size() << endl;
    // }

    
    
    // cout << "  " << adjList.size() << "  " << probabilities.size() << "  " << weights.size() << endl;
    // for (int v = 1; v <= numNodes - 1; v++){
    //     for (int u = 0; u < node_deg[v]; u++){
    //         cout << v << "  " << adjList[v][u] << "  " << weights[v][u + 1] << endl;
    //     }
    // }

    // cout << " Read the graph *********" << endl;
}

void Graph::writeToFile(const char * filename)
{/*
	ofstream output(filename);
	for (unsigned int i = 0; i < numNodes; ++i){
		for (unsigned int j = 0; j < adjList[i].size(); ++j){
			if (adjList[i][j] > i){
				output << adjList[i][j] << " " << i << " " << weights[i][j] << endl;
			}
		}
	}
	output.close();
*/	
}


// inline int HyperGraph::randIndex_dart(const vector<UI> &w, unsigned int si)
inline int HyperGraph::randIndex_dart(const vector<float> &w, unsigned int si)
{
    int prob = 0;
    while (prob == 0){
        UI ran = sfmt_genrand_uint32(&sfmtSeed)%si;
        UI ran2 = sfmt_genrand_uint32(&sfmtSeed);
        if (w[ran] >= ran2)
            prob = ran;
    }
    return prob;
}

HyperGraph::HyperGraph(unsigned int n)
{
	sfmt_init_gen_rand(&sfmtSeed, rand());
	node_edge = vector<vector<int> >(n+1);
	maxDegree = 0;
	numNodes = n;
	curEdge=0;
}

void HyperGraph::updateDeg(){
	unsigned int num=edge_node.size();
	for (unsigned int i = curEdge; i < num; ++i){
		unsigned int num2 = edge_node[i].size();
		for (unsigned int j=0;j<num2;++j){
			node_edge[edge_node[i][j]].push_back(i);
		}
	}
	curEdge = edge_node.size();
}

void HyperGraph::updateEdge(){
	curEdge = edge_node.size();
}

/*
* Add a hyperedge into the hypergraph
*/
void HyperGraph::addEdge(vector<int> & edge)
{
	edge_node.push_back(edge);
	unsigned int ind = edge_node.size() - 1;
	for (unsigned int i = 0; i < edge.size(); ++i)
		node_edge[edge[i]].push_back(ind);
}

/*
* Add a hyperedge into the hypergraph while keeping track of the node with max degree
*/
void HyperGraph::addEdgeD(vector<int> & edge)
{
    edge_node.push_back(edge);
    int ind = edge_node.size() - 1;
    for (unsigned int i = 0; i < edge.size(); ++i){
        node_edge[edge[i]].push_back(ind);
    	if (node_edge[edge[i]].size() > maxDegree)
    		maxDegree = node_edge[edge[i]].size();
	}
}

/*
* get an edge from the hypergraph
*/
const vector<int> & HyperGraph::getEdge(int e) const{
	return edge_node[e];
}

const vector<int> & HyperGraph::getEdge(int e){
	return edge_node[e];
}

/*
* get the list of hyperedges incident to node n
*/
const vector<int> & HyperGraph::getNode(int n) const{
	return node_edge[n];
}

const vector<int> & HyperGraph::getNode(int n){
	return node_edge[n];
}

/*
* get the number of hyperedges
*/
int HyperGraph::getNumEdge() const
{
    return edge_node.size();
}

/*
* get the maximum degree
*/
int HyperGraph::getMaxDegree()
{
	return maxDegree;
}

/*
* remove all the hyperedges
*/
void HyperGraph::clearEdges()
{
    edge_node.clear();
    node_edge.clear();
    cout << "clear edges!" << endl;
    maxDegree = 0;
}


// choose a random edge in LT model based on linear search
// inline int HyperGraph::randIndex_lin(Graph &g, const vector<UI> &w, unsigned int si, int cur, int p)
inline int HyperGraph::randIndex_lin(Graph &g, const vector<float> &w, unsigned int si, int cur, int p)
{
    // UI ranNum = g.LT_node[cur]; 
    // float ranNum = g.LT_node[cur];
    float ranNum = sfmt_genrand_real1(&sfmtSeed);
    // UI ranNum = sfmt_genrand_uint32(&sfmtSeed);
    // if (si <= 1 || ranNum > w[si - 1] * g.similarities[cur][p])
    if (si <= 0 || ranNum > w[si - 1] * g.similarities[cur][p]){
        return -1;
    }

    for (unsigned int i = 0; i < si; ++i){
        if (ranNum <= w[i] * g.similarities[cur][p])
            return i;
    }
    return -1;
}

// choose a random live edge in LT model based on binary search
// inline int HyperGraph::randIndex_bin(Graph &g, const vector<UI> &w, unsigned int si, int cur, int p)
inline int HyperGraph::randIndex_bin(Graph &g, const vector<float> &w, unsigned int si, int cur, int p)
{
    // UI ran = g.LT_node[cur];
    // float ran = g.LT_node[cur];
    float ran = sfmt_genrand_real1(&sfmtSeed);
    // UI ran = sfmt_genrand_uint32(&sfmtSeed);
    if (si <= 0 || ran > w[si - 1] * g.similarities[cur][p]){
        return -1;
    }
    int left = 0;
    int right = si - 1;
    int prob;
    for (unsigned int i = 0; i < si; ++i){
        prob = (left + right)/2;
        // if (w[prob - 1] * g.similarities[cur][p] > ran){
        if (w[prob] * g.similarities[cur][p] > ran){
            right = prob - 1;
            continue;
        }
        if (w[prob] * g.similarities[cur][p] <= ran){
            left = prob + 1;
            continue;
        }
        break;
    }
    return prob;
}


bool HyperGraph::pollingLT(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark, sampler<int> &samp, int p)
{
    unsigned int i;
    bool t = false;
    unsigned int gSize = g.getSize();
    // unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
    unsigned cur = samp();
    unsigned int num_marked = 0;
    // cout << "startLT: " << cur << endl;
    for (i = 0; i < gSize; ++i){
    	if (link[cur] < k){
            t=true;
    		break;
        }
        if (visit[cur] == true) break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
    	num_marked++;
    	int ind;

        ind = randIndex_lin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);
        if (g.weights[cur].size() >= 32)
        //     ind = randIndex_bin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);
        // else
        //     ind = randIndex_lin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);
        
        if (ind == -1)
            break;

        // cur = g.invAdjList[cur][ind-1];
        cur = g.invAdjList[cur][ind];
        // cout << cur << endl;
    }
    for (i = 0; i < num_marked; ++i){
        visit[visit_mark[i]]=false;
    }
    return t;
}


void HyperGraph::pollingLT1(Graph &g, vector<bool> &visit, vector<int> &visit_mark, sampler<int> &samp, int p)
{
    unsigned int i;
    unsigned int gSize = g.getSize();
    // unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
    unsigned cur = samp();
    unsigned int num_marked = 0;
    // cout << "startLT1: " << cur << endl;
    // cout << "^^^^^^^^^^^^^^^ 1" << endl;
    for (i = 0; i < gSize; ++i){
        // cout << "^^^^^^^^^^^^^^^ 1" << endl;
        if (visit[cur] == true) break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
    	num_marked++;
        const vector<int> &neigh = g[cur];
        int ind;
        // cout << "^^^^^^^^^^^^^^^ 2" << endl;
        ind = randIndex_lin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);
        // if (g.weights[cur].size() >= 32)
        //     ind = randIndex_bin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);
        // else
        //     ind = randIndex_lin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);

        // cout << "^^^^^^^^^^^^^^^ 3" << endl;
        if (ind == -1)
            break;
        // cur = neigh[ind-1];
        cur = neigh[ind];
        // cout << cur << endl;
    }
    // cout << "^^^^^^^^^^^^^^^ 4" << endl;
	edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
    for (i = 0; i < num_marked; ++i){
        visit[visit_mark[i]]=false;
    }
    // cout << "^^^^^^^^^^^^^^^ 5" << endl;
}


/*
* polling process under LT model
*/ 
bool HyperGraph::pollingLT2(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark, sampler<int> &samp, int p)
{   
    // cout << "&&&&&&&&&&&&&&&&&& 1" << endl;
    unsigned int i;
    bool t = false;
    unsigned int gSize = g.getSize();
    // cout << "&&&&&&&&&&&&&&&&&& 2" << endl;
    // unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1; // rand node
    unsigned cur = samp();
    unsigned int num_marked = 0;
    // cout << "startLT2: " << cur << endl;
    for (i = 0; i < gSize; ++i){ // Goes over all nodes
        // cout << cur << " &&&&&&&&&&&&&&&&&& 3" << endl;
        if (visit[cur] == true) break; // Only unvisited
        // cout << "&&&&&&&&&&&&&&&&&& 3" << endl;
        visit[cur] = true;  // marks it visited
        // cout << "&&&&&&&&&&&&&&&&&& 3" << endl;
        visit_mark[num_marked] = cur; // tracks which one was marked at each iteration
        // cout << "&&&&&&&&&&&&&&&&&& 3" << endl;
        num_marked++; // How many were visited
        // cout << "&&&&&&&&&&&&&&&&&& 4" << endl;
        if (link[cur] < k) // Has been influenced
            t=true;
        int ind;
        // cout << "&&&&&&&&&&&&&&&&&& 5" << endl;
        ind = randIndex_lin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);
        // if (g.weights[cur].size() >= 32)
        //     ind = randIndex_bin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);
        // else
        //     ind = randIndex_lin(g, g.weights[cur],g.inv_node_deg[cur], cur, p);

        if (ind == -1)
            break;
        // cout << "&&&&&&&&&&&&&&&&&& 6" << endl;
        // cur = g.invAdjList[cur][ind-1];
        cur = g.invAdjList[cur][ind];
        // cout << "&&&&&&&&&&&&&&&&&& 7" << endl;
        // cout << cur << endl;
    }
    // cout << "&&&&&&&&&&&&&&&&&& 8" << endl;
    edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
    for (i = 0; i < num_marked; ++i){
        visit[visit_mark[i]]=false;
    }
    // cout << "&&&&&&&&&&&&&&&&&& 9" << endl;
    return t;
}


// void HyperGraph::pollingIC1(Graph &g, vector<bool> &visit, vector<int> &visit_mark, sampler<int> &samp, int p)
// {
//     int i;
//     unsigned cur = samp();
//     int num_marked=1;
//     int curPos=0;
//     visit[cur] = true;
//     visit_mark[0] = cur;
//     float prob;
//     while(curPos < num_marked){
//         cur = visit_mark[curPos];
//         curPos++;
//         // const vector<UI> &w=g.getWeight(cur);
//         const vector<float> &w=g.getWeight(cur);
//         const vector<int> &neigh = g[cur];
//         for (i = 0; i < g.inv_node_deg[cur]; ++i){
//             // if (sfmt_genrand_real1(&sfmtSeed) <  (w[i]) * g.similarities[cur][p]){ //Check here is it supposed to be g.similarities[cur][p] or g.similarities[g[cur][i]][p]
//             if (sfmt_genrand_real1(&sfmtSeed) <  g.get_inv_AP_IC(cur, v, i, p)){ //Check here is it supposed to be g.similarities[cur][p] or g.similarities[g[cur][i]][p]
//                 if (!visit[neigh[i]]){
//                     visit[neigh[i]] = true;
//                     visit_mark[num_marked]=neigh[i];
//                     num_marked++;
//                 }
//             }
//         }
//     }
//     edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
//     for(i = 0; i < num_marked;++i){
//         visit[visit_mark[i]]=false;
//     }
// }

void HyperGraph::pollingIC1(Graph &g, vector<bool> &visit, vector<int> &visit_mark, sampler<int> &samp, int p)
{
    int i, v;
    unsigned cur = samp();
    int num_marked=1;
    int curPos=0;
    visit[cur] = true;
    visit_mark[0] = cur;
    float prob;
    while(curPos < num_marked){
        cur = visit_mark[curPos];
        curPos++;
        for (i = 0; i < g.inv_node_deg[cur]; ++i){
            v = g.invAdjList[cur][i];
            if (sfmt_genrand_real1(&sfmtSeed) <  g.get_inv_AP_IC(cur, v, i, p)){ //Check here is it supposed to be g.similarities[cur][p] or g.similarities[g[cur][i]][p]
                if (!visit[v]){
                    visit[v] = true;
                    visit_mark[num_marked]=v;
                    num_marked++;
                }
            }
        }
    }
    edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
    for(i = 0; i < num_marked;++i){
        visit[visit_mark[i]]=false;
    }
}

/*
* polling process under IC model
*/
bool HyperGraph::pollingIC2(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark, sampler<int> &samp, int p)
{
    int i, v;
    unsigned cur = samp();
	int curPos=0;
    int num_marked=1;
	visit[cur] = true;
	visit_mark[0] = cur;
	bool t = false;
    // cout << "Starting Cur Node: " << cur << endl; // del
    while(curPos < num_marked){
		cur = visit_mark[curPos];
		curPos++;
		if (link[cur] < k)
            t=true;
		for (i = 0; i < g.inv_node_deg[cur]; ++i){
            v = g.invAdjList[cur][i];
            if (sfmt_genrand_real1(&sfmtSeed) <  g.get_inv_AP_IC(cur, v, i, p)){
                if (!visit[v]){
					visit[v] = true;
					visit_mark[num_marked]=v;
					num_marked++;
				}
			}
		}
    }
    edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));

	for(i = 0; i < num_marked;++i){
		visit[visit_mark[i]]=false;
	}
	return t;
}

// bool HyperGraph::pollingIC2(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark, sampler<int> &samp, int p)
// {
//     int i;
//     unsigned cur = samp();
//     int curPos=0;
//     int num_marked=1;
//     visit[cur] = true;
//     visit_mark[0] = cur;
//     bool t = false;
//     // cout << "Starting Cur Node: " << cur << endl; // del
//     while(curPos < num_marked){
//         cur = visit_mark[curPos];
//         curPos++;
//         if (link[cur] < k)
//            t=true;
//         // const vector<UI> &w=g.getWeight(cur);
//         const vector<float> &w=g.getWeight(cur);
//         const vector<int> &neigh = g[cur];
//         // cout << "Explore " << cur  << " incoming neighbors"<< endl; //del
//         for (i = 0; i < g.inv_node_deg[cur]; ++i){
//             // cout << "Weight: " << w[i] << "Sim: " << g.similarities[cur][p] << endl; 
//             // cout << "Attemp to visit " << neigh[i] << " with probability " << (w[i]) * g.similarities[cur][p] << endl; //del
//             // if (sfmt_genrand_real1(&sfmtSeed) <  (w[i]) * g.similarities[cur][p]){
//             if (sfmt_genrand_real1(&sfmtSeed) <  (w[i])){
//             // if (sfmt_genrand_uint32(&sfmtSeed) <  (w[i+1]) * g.similarities[cur][p] * g.UI_MAX){ //Check here is it supposed to be g.similarities[cur][p] or g.similarities[g[cur][i]][p]
//                 if (!visit[neigh[i]]){
//                     visit[neigh[i]] = true;
//                     visit_mark[num_marked]=neigh[i];
//                     num_marked++;
//                 }
//             }
//         }
//     }
//     edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));

//     for(i = 0; i < num_marked;++i){
//         visit[visit_mark[i]]=false;
//     }
//     return t;
// }

bool HyperGraph::pollingIC(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark, sampler<int> &samp, int p)
{
    int i;
    // unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(g.getSize())+1;
    unsigned cur = samp();
    int curPos=0;
    int num_marked=1;
    visit[cur] = true;
    visit_mark[0] = cur;
    bool t = false;
    while(curPos < num_marked){
        cur = visit_mark[curPos];
        curPos++;
        if (link[cur] < k){
            t=true;
            break;
        }
        // const vector<UI> &w=g.getWeight(cur);
        const vector<float> &w=g.getWeight(cur);
        const vector<int> &neigh = g[cur];
        for (i = 0; i < g.inv_node_deg[cur]; ++i){
            if (sfmt_genrand_real1(&sfmtSeed) <  w[i]){
            // if (sfmt_genrand_real1(&sfmtSeed) <  w[i] * g.similarities[cur][p]){
                if (!visit[neigh[i]]){
                    visit[neigh[i]] = true;
                    visit_mark[num_marked]=neigh[i];
                    num_marked++;
                }
            }
        }
    }
    for(i = 0; i < num_marked;++i){
        visit[visit_mark[i]]=false;
    }
    return t;
}

/*
* measure the consumed memory
*/
float getCurrentMemoryUsage() {

    // string pid = to_string(unsigned(getpid()));
    // string outfile = "tmp_" + pid + ".txt";
    // string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
    // system(command.c_str());

    // string mem_str;
    // ifstream ifs(outfile.c_str());
    // std::getline(ifs, mem_str);
    // ifs.close();

    // mem_str = mem_str.substr(0, mem_str.size()-1);
    // float mem = stof(mem_str);

    // command = "rm " + outfile;
    // system(command.c_str());

    // return mem/1024;

    return 0;
}

void Graph::initial_seedset(string alg, Utilities &util, Locations * usersLocations, Locations * eventsLocations, int P, int k){
    for (int p = 0; p < P; ++p){
        seedSets[p].clear();
        seedSets_before[p].clear();
    }
    

    ///////////////////////// Closest seedset (CG) (NN)////////////////////////
    if (alg == "CG" || alg == "NN" ){
        double dist;
        BinaryHeap heap(k);
        for (int p = 0; p < P; ++p){
            for (int v = 1; v <= getSize(); ++v){
                dist = util.computeMinimumDistance(usersLocations->locations[v][0], usersLocations->locations[v][1], eventsLocations->locations[p][0], eventsLocations->locations[p][1], 'K');
                if (dist < heap.getMaxDistance()){
                    heap.insertNewMax(dist, v);
                }
            }
            for (int i = 0; i < k; ++i){
                seedSets[p].insert(heap.getVertices()[i]);
                seedSets_before[p].insert(heap.getVertices()[i]);
            }
            heap.init();
        }
    }

    ///////////////////////// Random seedset (RG) (RD) ////////////////////////
    else if (alg == "RG" || alg == "RD"){
        int randomSeed;
        double dist;
        for (int p = 0; p < P; ++p){
            while (seedSets[p].size() < k + 1){
                randomSeed = rand()%(getSize() + 1) + 1;
                if (randomSeed != 0){
                    seedSets[p].insert(randomSeed);
                    seedSets_before[p].insert(randomSeed);
                }
            }
        }
    }

    // ///////////////////////// Empty seedset (EG) ////////////////////////
    // The seedset is already empty 
}

double Graph::find_max_dist(int P,  Locations* usersLocations, Locations* eventsLocations, Utilities util){
    double max_dist = 0;
    double dist;
    for (int p = 0; p < P; ++p){
        for (int v = 1; v <= getSize(); ++v){
            dist = util.computeMinimumDistance(usersLocations->locations[v][0], usersLocations->locations[v][1], eventsLocations->locations[p][0], eventsLocations->locations[p][1], 'K');
            if (dist > max_dist){
                max_dist = dist;
            }
        }
    }
    return max_dist;
}

void Graph::rand_LT(){
    // Random thresholds for each node
    for (int i = 1; i <= numNodes; i++){
        LT_node[i] = sfmt_genrand_real1(&sfmtSeed);
    }
}

void Graph::BFS_IC(int p, float** APc){
    queue<int> myqueue;
    vector<bool> visited(getSize() + 1,false);
    for (int seed : seedSets[p]){
        myqueue.push(seed);
        APc[p][seed] += 1;
        visited[seed] = true;
    }
    // (u) ---> (v)
    int u; //current node
    int v; //current node's friend
    while (!myqueue.empty()){
        u = myqueue.front();
        myqueue.pop();
        for (int i = 0; i < node_deg[u]; ++i){
            v = adjList[u][i];
            if (!visited[v] and sfmt_genrand_real1(&sfmtSeed) <= get_AP_IC(u, v, i, p)){
            // if (!visited[v] and sfmt_genrand_real1(&sfmtSeed) < probabilities[u][i] * similarities[v][p]){
                myqueue.push(v);
                APc[p][v] += 1;
                visited[v] = true;
            }
        }
    }
}

void Graph::BFS_LT(int p, float** APc){
    // (u) ---> (v)
    queue<int> myqueue;
    vector<bool> active(getSize() + 1,false); //tracks the users that are active
    vector<double> active_prob(getSize() + 1,0); //tracks the probability of the active incoming edges for each node
    int user_front, user_front_friend;
    int u, v;

    for (int seed : seedSets[p]){ // Sets seeds as active and counts that they are Aware of p
        active[seed] = true;
        APc[p][seed]++;
    }

    for (int seed : seedSets[p]){
        for (int i = 0; i < node_deg[seed]; i++){ // Inserts all the users that are not active to the queue and could be activated (seed) ---> (v)
            v = adjList[seed][i];
            active_prob[v] += probabilities[seed][i];
            if (active[v] == false){
                myqueue.push(v); // Now we need to check if (v) will become active since the seeds were active
            }
        }
    }

    while (!myqueue.empty()){ // Continues till none can get aware
        user_front = myqueue.front(); // Gets the first guy from the queue
        u = myqueue.front(); // Gets the first guy from the queue
        myqueue.pop();
        // if (active[user_front] == false && similarities[user_front][p] >= sfmt_to_real1(LT_node[user_front]) && active_prob[user_front] * similarities[user_front][p] >= sfmt_to_real1(LT_node[user_front])){ // Making sure the first guy is not active

        // if (active[u] == false && get_AP_LT(active_prob[u], u, p)  >= sfmt_to_real1(LT_node[user_front])){ // Making sure the first guy is not active
        if (active[u] == false && get_AP_LT(active_prob[u], u, p)  >= LT_node[user_front]){ // Making sure the first guy is not active
            active[u] = true;
            APc[p][u]++;
            for (int i = 0; i < node_deg[u]; i++){ // Inserts all the outgoing neighbors that are not active since they can become aware
                v = adjList[u][i];
                if (active[v] == false)
                    active_prob[v] += probabilities[u][i];
                    myqueue.push(v);
            }
        }
    }  
}

double Graph::calculate_tsim(int P){
    float** aware_all_ev = new float*[P];
    for (int p = 0; p < P; ++p){
        aware_all_ev[p] = new float[getSize() + 1];
        for(int v = 0; v <= getSize(); v++){
            aware_all_ev[p][v] = 0;
        }
    }
    
    for (int ev = 0; ev < P; ev++){
        if (model == 0){
            rand_LT();
            BFS_LT(ev, aware_all_ev);
        }
        else{
            BFS_IC(ev, aware_all_ev);
        }
    }

    double total_sim = 0;
    double max_weight;
    for (int v = 1; v <= getSize(); ++v){
        max_weight = 0;
        for (int p = 0; p < P; ++p){
            if (aware_all_ev[p][v] != 0 && similarities[v][p] > max_weight){
                max_weight = similarities[v][p];
            }
        }
        total_sim += max_weight;
    }

    for (int p = 0; p < P; ++p){
        delete [] aware_all_ev[p];
    }
    delete [] aware_all_ev;

    return total_sim;
}

void Graph::calculate_sim_stats(int P, int ev, double * stats){
    float** aware_all_ev = new float*[P];
    for (int p = 0; p < P; ++p){
        aware_all_ev[p] = new float[getSize() + 1];
        for(int v = 0; v <= getSize(); v++){
            aware_all_ev[p][v] = 0;
        }
    }

    for (int ev = 0; ev < P; ev++){
        if (model == 0){
            rand_LT();
            BFS_LT(ev, aware_all_ev);
        }
        else
            BFS_IC(ev, aware_all_ev);
    }
    double similarity_Sj = 0;
    double similarity_S_Sj = 0;
    double total_sim  = 0;
    double max_weight;
    double max_weight_not_ev;
    int max_event;
    for (int v = 1; v <= getSize(); ++v){
        // if (aware_all_ev[ev][v] == 1){
        //     similarity_Sj += similarities[v][ev];
        // }
        max_weight = 0;
        max_weight_not_ev = 0;
        for (int p = 0; p < P; ++p){
            if (aware_all_ev[p][v] == 1 && similarities[v][p] > max_weight_not_ev && p != ev){
                max_weight_not_ev = similarities[v][p];
            }
            if (aware_all_ev[p][v] == 1 && similarities[v][p] > max_weight){
                max_weight = similarities[v][p];
                max_event = p;
            }
        }
        if (max_event == ev){
            similarity_Sj += max_weight;
        }
        total_sim += max_weight;
        similarity_S_Sj += max_weight_not_ev;
    }
    for (int p = 0; p < P; ++p){
        delete [] aware_all_ev[p];
    }
    delete [] aware_all_ev;

    stats[0] = total_sim;
    stats[1] = similarity_Sj;
    stats[2] = similarity_S_Sj;
}

void Graph::example(Locations* usersLocations, Locations* eventsLocations, int P)
{
    ofstream events;
    ofstream all_but_seeds;
    ofstream not_influenced;
    ofstream influenced;
    ofstream graph;
    ofstream * first = new ofstream[P];
    ofstream * second = new ofstream[P];
    vector<int> all_seeds;


    influenced.open("influenced.txt");
    graph.open("graph.txt");

    events.open("example_events.csv");
    events << "id,latitude,longitude" << endl;
    all_but_seeds.open("example_all_but_seeds.csv");
    all_but_seeds << "id,latitude,longitude" << endl;
    not_influenced.open("example_not_influenced.csv");
    not_influenced << "id,latitude,longitude" << endl;



    vector<vector<int>> influenced_ev(P);

    int eds = 0;
    graph << numNodes << " " << numEdges << endl;
    for (int i = 1; i <= numNodes; ++i){
        eds += node_deg[i];
        for (int j = 1; j < node_deg[i]; ++j){
            cout << i << " " << adjList[i][j] << " " << probabilities[i][j] << endl;
            graph << i << " " << adjList[i][j] << " " << probabilities[i][j] << endl;
        }     
    }
    cout << "Number of edges: " << eds << endl;

    for (int p = 0; p < P; ++p){
        first[p].open("example_seeds_event" + to_string(p + 1) + ".csv");
        first[p] << "id,latitude,longitude" << endl;
        second[p].open("example_inf_event" + to_string(p + 1) + ".csv");
        second[p] << "id,latitude,longitude" << endl;
        events << p + 1 << "," << setprecision(12) << eventsLocations->locations[p][0] << "," << setprecision(12) << eventsLocations->locations[p][1] << endl;
        for (int seed: seedSets[p]){
            first[p] << seed << "," << setprecision(12) << usersLocations->locations[seed][0] << "," << setprecision(12) << usersLocations->locations[seed][1] << endl;
            all_seeds.push_back(seed);
        }
    }

    float** aware_all_ev = new float*[P];
    for (int p = 0; p < P; ++p){
        aware_all_ev[p] = new float[getSize() + 1];
        for(int v = 0; v <= getSize(); v++){
            aware_all_ev[p][v] = 0;
        }
    }
                    
    for (int p = 0; p < P; p++){
        if (model == 0)
            BFS_LT(p, aware_all_ev);
        else
            BFS_IC(p, aware_all_ev);
    }


    double max_weight;
    int max_event;
    vector<int>::iterator it;
    for (int v = 1; v <= getSize(); ++v){
        it = find(all_seeds.begin(), all_seeds.end(), v);
        if (all_seeds.end() == it){
            all_but_seeds << v << "," << setprecision(12) << usersLocations->locations[v][0] << "," << setprecision(12) << usersLocations->locations[v][1] << endl;
        }
        max_weight = 0;
        max_event = -1;
        for (int p = 0; p < P; ++p){
            if (aware_all_ev[p][v] == 1 && similarities[v][p] > max_weight){
                max_weight = similarities[v][p];
                max_event = p;
            }
        }
        if (max_event == -1){
            not_influenced << v << "," << setprecision(12) << usersLocations->locations[v][0] << "," << setprecision(12) << usersLocations->locations[v][1] << endl;
        }
        else{
            influenced_ev[max_event].push_back(v);
            second[max_event] << v << "," << setprecision(12) << usersLocations->locations[v][0] << "," << setprecision(12) << usersLocations->locations[v][1] << endl;
        }
    }

    for (int p = 0; p < P; ++p){
        for (int i = 0; i < influenced_ev[p].size(); ++i)
            influenced << influenced_ev[p][i] << " ";
        influenced << endl;
    }

    for (int p = 0; p < P; ++p){
        delete [] aware_all_ev[p];
    }
    delete [] aware_all_ev;
    delete [] first;
    delete [] second;
}