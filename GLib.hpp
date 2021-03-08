#ifndef _GLIB_H_
#define _GLIB_H_

#include <vector>
#include <algorithm>

template <typename T>
void sort_indexes(const std::vector<T> &v, std::vector<size_t> &idx)
{
        // initialize origianl index locations
        for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
        
        // sort indexes based on comparing values in v
         sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2){return v[i1] > v[i2];});
}

#endif
