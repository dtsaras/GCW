/*
 * 	Mapped Priority Queue
 *
 *  Created on: May 17, 2010
 *      Author: Thang N. Dinh
 *      Version: 1.0 	May 17, 2010
 *      		 1.03	Nov 12, 2010    Fix heap creation problem, fix down heap problem
 */

#ifndef MAPPED_PRIORITY_QUEUE_H_
#define MAPPED_PRIORITY_QUEUE_H_

#include	<iostream>
#include	<vector>
#include	<algorithm>

using namespace std;
//using namespace std::tr1;


template <class T>
struct lessArray {
	T *value;
public:
	lessArray(T* vl=NULL):value(vl){};
	bool operator()(int &x, int &y) const {
		return value[x] < value[y];
	}
};

template <class T>
struct greaterArray {
	T *value;
public:
	greaterArray(T* vl=NULL):value(vl){};
	bool operator()(int &x, int &y) const {
		return value[x] > value[y];
	}
};



/**
 * @class MappedHeap
 * @brief A priority queue that allows to update/modify the value of inserted elements.
 * The queue store only keys of objects, keys are expected to be consecutive
 * and non-negative. The max value of keys should have same order with the
 * number of objects in the priority queue.
 * @author  Thang N. Dinh
 *
 */
template<class CompareLess>
class MappedHeap {
private:
	vector<int> heap;
	vector<int> position;
	int n; // Heap size
	int up_heap(int k);
	int down_heap(int k);
	CompareLess compLess;
public:
	/**
	 * Default constructor that create an empty heap
	 */
	MappedHeap();

	/**
	 * Default constructor that create an empty heap
	 */
	MappedHeap(CompareLess cl):compLess(cl){};


	/**
	 * Create the heap from an array of indices
	 * Time complexity: O(n) where n is the number of indices
	 */
	MappedHeap(const vector<int> &v);

	/**
	 * Create the heap from an array of indices
	 * Time complexity: O(n) where n is the number of indices
	 */
	MappedHeap(const vector<int> &v, CompareLess cl);

	/**
	 * Adjust the heap when the values associated with
	 * the index 'key' is modified/updated
	 */
	int heapify(int key);

	/**
	 * Alias for heapify
	 */
	int modify(int key);

	/**
	 * Insert an index into heap
	 */
	int push(int key);

	/**
	 * Remove the element with the highest priority
	 */
	int pop();

	/**
	 * Get the index of the element with the highest priority
	 */
	int top();

	/**
	 * Check if the heap is empty or not
	 */
	bool empty();

	/**
	 * Number of elements in the heap
	 */
	size_t size();


	void printOut(){
		printf("\nHeap: ");
		for (int i = 0; i < heap.size(); ++i)
			printf("%d ", heap[i]);
		printf("\n");
	}
};

template<class CompareLess>
int MappedHeap<CompareLess>::up_heap(int k) {
	while (k > 1) {
		int p = k / 2;
		if (compLess(heap[p], heap[k])) {
			swap(position[heap[p]], position[heap[k]]);
			swap(heap[p], heap[k]);
		}
		k = p;
	}
	return k;
}

template<class CompareLess>
int MappedHeap<CompareLess>::down_heap(int k) {
	while (k * 2 <= n) {
		int child = k * 2, right_child = child + 1;
		if ((right_child <= n)
				&& (compLess(heap[child], heap[right_child])))
                    child = right_child;
		if ((compLess(heap[k], heap[child]))) {
			swap(position[heap[k]], position[heap[child]]);
            swap(heap[k], heap[child]);			
			k = child;
		} else
			break;
	}
	return k;
}

template<class CompareLess>
MappedHeap<CompareLess>::MappedHeap(const vector<int> &v) {
	if (v.empty())
		n = 0;
	else {
		n = v.size();
		heap.resize(n + 1);
		position.resize(max<int> (*max_element(v.begin(), v.end()), n) + 3, -1);
		for (int i = 0; i < n; i++) {
			heap[i + 1] = v[i];
			position[v[i]] = i + 1;
		}
		for (size_t i = n/ 2; i >= 1; i--)
			down_heap(i);
	}
}

template<class CompareLess>
MappedHeap<CompareLess>::MappedHeap(const vector<int> &v,CompareLess cl):compLess(cl) {
	if (v.empty())
		n = 0;
	else {
		n = v.size();
		heap.resize(n + 1);
		position.resize(max<int> (*max_element(v.begin(), v.end()), n) + 3, -1);
		for (int i = 0; i < n; i++) {
			heap[i + 1] = v[i];
			position[v[i]] = i + 1;
		}
		for (size_t i = n/ 2; i >= 1; i--)
			down_heap(i);
	}
}


template<class CompareLess>
MappedHeap<CompareLess>::MappedHeap() {
	heap.push_back(0); // Ignore heap[0]
	position.clear();
	n = 0;
}

template<class CompareLess>
inline int MappedHeap<CompareLess>::heapify(int key) {
	//DCHECK(key < (int)position.size(),"Error:MappedHeap:heapify: index out of range " << key <<endl )
	//DCHECK((1 <= position[key]) && (position[key] <= n), "Error:MappedHeap:heapify: invalid position " << position[key] << " key: "<<key<< endl);
	up_heap(position[key]);    
	return down_heap(position[key]);
}

template<class CompareLess>
inline int MappedHeap<CompareLess>::modify(int key) {
	return heapify(key);
}

template<class CompareLess>
int MappedHeap<CompareLess>::push(int key) {
	heap.push_back(key);
	n++;
	if (key > (int)position.size() - 1)
		position.resize(key + 3, -1);
	position[key] = n;
	return up_heap(n);
}

template<class CompareLess>
int MappedHeap<CompareLess>::pop() {
	if (!empty()) {
		int top = heap[1];
		position[heap[1]] = 0;
		heap[1] = heap[n];
		position[heap[1]] = 1;
		n--;
		down_heap(1);
		return top;
	} else {
		cerr << "Error: Cannot pop from an empty heap!" << endl;
		return -1;
	}
}

template<class CompareLess>
inline int MappedHeap<CompareLess>::top() {
	return heap[1];
}

template<class CompareLess>
inline bool MappedHeap<CompareLess>::empty() {
	return n <= 0;
}

template<class CompareLess>
inline size_t MappedHeap<CompareLess>::size() {
	return n;
}

#ifdef	_DEBUG
static vector<int> value;

class _CheckHeap {
	struct compareCheck {
		bool operator()(int &x, int &y) const {
			return value[y] < value[x];
		}

	};

public:
	_CheckHeap() {
		value.clear();
		for(int i=0; i < 10; i++)
		value.push_back(i);

		vector<int> v;
		v.push_back(5);
		v.push_back(2);
		v.push_back(4);
		v.push_back(6);
		MappedHeap<compareCheck> mh(v);

		cout <<"Testing heap...";
		CHECK( mh.top() == 2,"Error: MappedHeap constructor test failed!" )
		mh.push(1);
		CHECK( mh.top() == 1,"Error: MappedHeap push test failed!" )
		value[4] = 0;
		mh.heapify(4);
		CHECK( mh.top() == 4,"Error: MappedHeap heapify test failed!" )
		mh.pop();
		CHECK( mh.top() == 1,"Error: MappedHeap pop test failed!" )
		cout <<"..done"<<endl;
	}
};
//_CheckHeap _checkHeap;
#endif

#endif /* MAPPED_PRIORITY_QUEUE_H_ */
