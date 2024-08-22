#ifndef NUENV_ALGORITHM_SPACE_H_
#define NUENV_ALGORITHM_SPACE_H_

#include "nuenv/src/core/container.hpp"
#include "nuenv/src/core/ctypes.hpp"
#include "nuenv/src/core/math.hpp"

namespace nuenv {

/**
 * @brief Generate an array with evenly spaced scalars over a specified interval.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param start Starting value of the sequence.
 * @param stop End value of the sequence.
 * @param num Number of samples to generate. Must be non-negative.
 *
 * @return Array of evenly spaced scalars.
 */
template<typename Scalar>
VectorX<Scalar> LinearSpace(Scalar start, Scalar stop, Index num) {
  return VectorX<Scalar>::LinSpaced(num, start, stop);
}

/**
 * @brief Generate an array with evenly spaced scalars on a log scale over a
 *  specified interval.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param start Starting value of the sequence. Must be positive.
 * @param stop End value of the sequence.
 * @param num Number of samples to generate. Must be non-negative.
 *
 * @return Array of evenly spaced scalars.
 */
template<typename Scalar>
VectorX<Scalar> LogarithmicSpace(Scalar start, Scalar stop, Index num) {
  assert((start > 0.0 && stop > 0.0) && "Arguments must not be 0");

  constexpr Scalar base = 2.0;

  Scalar a_log = log(start) / log(base);
  Scalar b_log = log(stop) / log(base);
  VectorX<Scalar> arr = LinearSpace(a_log, b_log, num);

  for (Index i = 0; i < num; i++) {
	arr[i] = pow(base, arr[i]);
  }

  // Returns with copy elision
  return arr;
}

/**
 * @brief Generate an array with a geometric progression of scalars over a
 *  specified interval.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param start Starting value of the sequence. Must be positive.
 * @param stop End value of the sequence.
 * @param num Number of samples to generate. Must be non-negative.
 *
 * @return Array of evenly spaced scalars.
 */
template<typename Scalar>
VectorX<Scalar> geometricSpace(Scalar start, Scalar stop, Index num) {
  assert((start > 0.0 && stop > 0.0) && "Arguments must not be 0");

  Scalar r = pow(stop / start, 1.0 / (static_cast<Scalar>(num) - 1.0));

  VectorX<Scalar> arr(num);

  for (Index i = 0; i < num; i++) {
	arr[i] = start * pow(r, i);
  }

  // Returns with copy elision
  return arr;
}

}

#endif
