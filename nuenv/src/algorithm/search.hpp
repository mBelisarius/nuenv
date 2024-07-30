#ifndef NUENV_ALGORITHM_SEARCH_H_
#define NUENV_ALGORITHM_SEARCH_H_

#include "nuenv/core"

namespace nuenv {

namespace internal {

template<typename T>
T bchoice(bool condition, T v_true, T v_false) {
  return (condition * v_true) | (!condition * v_false);
}

template<typename Scalar>
Index linearSearch(const VectorX<Scalar>& arr,
				   Scalar val,
				   Index begin = 0,
				   Index end = 0) {
  assert((arr.size() > 0) && "Array must not be empty");

  Index size = arr.size();
  if (end == 0) { end = size - 1; }

  for (; begin <= end; ++begin) {
	if (arr[begin] == val) { break; }
  }

  // Search successful
  if ((begin < end) || (arr[end] == val)) { return begin; }

  return size;
}

template<typename Scalar>
Index linearSortedSearch(const VectorX<Scalar>& arr,
						 Scalar val,
						 Index begin = 0,
						 Index end = 0) {
  assert((arr.size() > 0) && "Array must not be empty");

  if (end == 0) { end = arr.size() - 1; }

  Index cnt = 0;
  for (Index i = end + 1; i > begin; --i) {
	cnt += (arr[i - 1] > val);
  }

  return end - bchoice(cnt > end, end, cnt);
}

template<typename Scalar>
Index binarySearch(const VectorX<Scalar>& arr,
				   Scalar val,
				   Index begin = 0,
				   Index end = 0) {
  assert((arr.size() > 0) && "Array must not be empty");

  Index size = arr.size();
  if (end == 0) { end = size; }

  while (end - begin > 1) {
	Index mid = (begin + end) / 2;
	bool condition = val < arr[mid];

	begin = bchoice(!condition, mid, begin);
	end = bchoice(condition, mid, end);
  }

  // Bounds checking
  begin = bchoice<Index>(val < arr[0], 0, begin);
  begin = bchoice<Index>(val > arr[size - 1], size - 1, begin);

  return begin;
}

template<typename Scalar, Index skip = 8>
Index skiplistSearch(const VectorX<Scalar>& arr,
					 Scalar val,
					 Index begin = 0,
					 Index end = 0) {
  assert((arr.size() >= skip) && "Array must have at least 'skip' elements");

  Index size = arr.size();
  VectorX<Scalar> skplst = VectorX<Scalar>(size / skip);
  for (Index i = 0; i < size / skip; ++i) {
	skplst[i] = arr[skip * i];
  }

  Index j = binarySearch(skplst, val, begin, end);
  if (j < size / skip - 1) { end = (j + 1) * skip; }
  begin = linearSortedSearch(arr, val, skip * j, end);

  return begin;
}

} // namespace internal

template<typename Scalar>
Index SearchSorted(const VectorX<Scalar>& arr, Scalar val) {
  assert((arr.size() > 0) && "Array must not be empty");

  if (arr.size() > 64) { return internal::skiplistSearch(arr, val); }

  return internal::linearSortedSearch(arr, val);
}

template<typename Scalar>
Index SearchUnsorted(const VectorX<Scalar>& arr, Scalar val) {
  assert((arr.size() > 0) && "Array must not be empty");

  return internal::linearSearch<Scalar>(arr, val);
}

} // namespace nuenv

#endif
