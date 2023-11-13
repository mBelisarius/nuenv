#ifndef NUENV_ALGORITHM_SEARCH_H
#define NUENV_ALGORITHM_SEARCH_H

#include "nuenv/src/Core/Container.h"
#include "nuenv/src/Core/Exception.h"

namespace nuenv {

/**
 * @brief Search through a sorted array and returns the index of the element
 *  that is immediately lower than 'val'.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param arr Array to search in.
 * @param val Value to search for.
 *
 * @return Position in the array where the value is immediately lower than 'val'.
 */
template<typename Scalar>
size_t binarySearch(const VectorX<Scalar>& arr, Scalar val)
{
    if (arr.size() == 0) { throw invalid_argument("Array must not be empty"); }

    size_t lower_index = 0;
    size_t upper_index = arr.size() - 1;

    while (upper_index - lower_index > 1)
    {
        size_t mid_index = lower_index + (upper_index - lower_index) / 2;
        Scalar mid = arr[mid_index];

        if (val == mid) { return mid_index; }
        if (val > mid) { lower_index = mid_index; }
        else { upper_index = mid_index; }
    }

    return lower_index;
}

}

#endif
