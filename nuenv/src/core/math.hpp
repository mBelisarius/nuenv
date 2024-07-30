#ifndef NUENV_CORE_MATHUTILS_H
#define NUENV_CORE_MATHUTILS_H

#include <cmath>
#include <numbers>

namespace nuenv {

using namespace std::numbers;

using std::ceil;

using std::min;

using std::max;

using std::sqrt;

using std::pow;

using std::exp;

using std::log;

using std::log2;

using std::log10;

using std::sin;

/**
 * @brief Calculates the square of a given scalar value.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param x Value to be squared.
 *
 * @return Squared value.
 *
 * @note This function does not handle complex numbers. If the input scalar
 *       value is a complex number, the returned value will be the square of
 *       the real part of that number.
 */
template<typename Scalar>
Scalar pow2(Scalar x) {
  return x * x;
}

} // namespace nuenv

#endif
