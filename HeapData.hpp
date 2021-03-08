#ifndef _HEAP_DATA_H_
#define _HEAP_DATA_H_

template<class T1>
struct InfCost{
	T1 *v1;
public:
	InfCost(T1 *u1):v1(u1){};
	
	bool operator() (int &i, int &j) const{
		return v1[i] < v1[j];
	}
};

#endif
