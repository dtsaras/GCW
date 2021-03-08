#include "BinaryHeap.h"
#include <iostream>
#include <limits>

using namespace std;

BinaryHeap::BinaryHeap(int numItems){
	size = numItems;
	distances = new double[size];
	vertices = new int[size];
	for (int i = 0; i < size; ++i){
		distances[i] = numeric_limits<double>::max();
		vertices[i] = 0;
	}

	for (int i = parent(size - 1); i >= 0; --i){
		hipify(i);
	}
}

int BinaryHeap::left(int idx) const{
	return (2*idx + 1);
}

int BinaryHeap::right(int idx) const{
	return (2*idx + 2);
}

int BinaryHeap::parent(int idx) const{
	return ((idx - 1)/2);
}

double BinaryHeap::priority(int idx) const{
	return distances[idx];
}

void BinaryHeap::swap(int idx1, int idx2){
	double tempDist = distances[idx1];
	int tempVert = vertices[idx1];
	distances[idx1] = distances[idx2];
	vertices[idx1] = vertices[idx2];
	distances[idx2] = tempDist;
	vertices[idx2] = tempVert;
}

void BinaryHeap::hipify(int idx){
	int l = left(idx);
	int r = right(idx);

	int max = idx;
	if (l < size && priority(l) > priority(idx)){
		max = l;
	}

	if (r < size && priority(r) > priority(max)){
		max = r;
	}

	if (max != idx){
		swap(max, idx);
		hipify(max);
	}
}

void BinaryHeap::percolate(int idx){
	int curr = idx;
	int p = parent(idx);
	while (curr > 0 && priority(p) < priority(curr)){
		swap(p, curr);
		curr = p;
		p = parent(curr);
	}
}

BinaryHeap::~BinaryHeap(){
	delete[] distances;
	delete[] vertices;
}

int BinaryHeap::getMax() const{
	return vertices[0];
}

double BinaryHeap::getMaxDistance() const{
	return distances[0];
}

void BinaryHeap::popMax(){
	distances[0] = -1;
	vertices[0] = -1;
	swap(0, size - 1);
	--size;
	hipify(0);
}

void BinaryHeap::insertNewMax(double dist, int vert){
	distances[0] = dist;
	vertices[0] = vert;
	hipify(0);
}

int* BinaryHeap::getVertices(){
	return vertices;
}

double* BinaryHeap::getDistances(){
	return distances;
}

void BinaryHeap::print(){
	cout << "[";
	for (int i = 0; i < size; ++i){
		cout << "(" << distances[i] << ", " << vertices[i] << "), ";
	}
	cout << "]" << endl;
}

void BinaryHeap::init(){
	for (int i = 0; i < size; ++i){
		distances[i] = numeric_limits<double>::max();
		vertices[i] = 0;
	}
}







