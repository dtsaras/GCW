#ifndef BINARY_HEAP_H
#define BINARY_HEAP_H

class BinaryHeap
{
	private:
		double* distances;
		int* vertices;
		int size;

	protected:
		void hipify(int idx);

		int right(int idx) const;

		int left(int idx) const;

		int parent(int idx) const;

		double priority(int idx) const;

		void swap(int idx1, int idx2);

		void percolate(int idx);

	public:
		BinaryHeap(int numItems);

		~BinaryHeap();

		int getMax() const;

		double getMaxDistance() const;

		void popMax();

		void insertNewMax(double dist, int vert);

		int* getVertices();

		double* getDistances();

		void print();

		void init();
};



#endif