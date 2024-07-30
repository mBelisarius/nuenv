#ifndef NUENV_INTERPOLATE_INTERP1D_H
#define NUENV_INTERPOLATE_INTERP1D_H

#include "nuenv/src/algorithm/search.hpp"
#include "nuenv/src/core/container.hpp"
#include "nuenv/src/core/math.hpp"

namespace nuenv {

/**
 * @class Interp1d
 *
 * @brief Interpolate a 1-dimensional function.
 *
 * This class implements methods whose call uses interpolation to find the value
 * of new points of some function f: 'y = f(x)'.
 *
 * @tparam Scalar Scalar type of the numbers.
 */
template<typename Scalar>
class Interp1d {
 public:
  Interp1d(const VectorX<Scalar>& x,
		   const VectorX<Scalar>& y,
		   bool check_bounds = true);

  Scalar linear(Scalar x);

  Scalar exponential(Scalar x);

 private:
  VectorX<Scalar> x_;
  VectorX<Scalar> y_;
  size_t size_;
  bool check_bounds_;
};

/**
 * Constructs the interpolator.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param x Array of x-values representing the independent variable,
 *  must be increasing.
 * @param y Array of y-values representing the dependent variable.
 * @param check_bounds Indicates whether to check if new points are within the
 *  domain bounds of 'x'.
 */
template<typename Scalar>
Interp1d<Scalar>::Interp1d(const VectorX<Scalar>& x,
						   const VectorX<Scalar>& y,
						   const bool check_bounds)
	: x_(x), y_(y), size_(x.size()), check_bounds_(check_bounds) {
  assert((x.size() > 0 && y.size() > 0) && "Arrays must not be empty");
  assert((x.size() == y.size()) && "Arrays 'x' and 'y' must have same size");
}

/**
 * @brief Linear interpolation.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param x Point to be interpolated.
 *
 * @return Interpolated value at 'x'.
 */
template<typename Scalar>
Scalar Interp1d<Scalar>::linear(Scalar x) {
  assert(!(check_bounds_ && x < x_[0] && x > x_[size_ - 1]) && "'x' is out of bounds");

  if (x < x_[0]) { return y_[0]; }
  else if (x > x_[size_ - 1]) { return y_[size_ - 1]; }

  size_t index = SearchSorted<Scalar>(x_, x);

  return y_[index]
	  + ((y_[index + 1] - y_[index]) / (x_[index + 1] - x_[index]))
		  * (x - x_[index]);
}

/**
 * @brief Exponential interpolation.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param x Point to be interpolated.
 *
 * @return Interpolated value at 'x'.
 */
template<typename Scalar>
Scalar Interp1d<Scalar>::exponential(Scalar x) {
  assert(!(check_bounds_ && x < x_[0] && x > x_[size_ - 1]) && "'x' is out of bounds");

  if (x < x_[0]) { return y_[0]; }
  else if (x > x_[size_ - 1]) { return y_[size_ - 1]; }

  size_t index = SearchSorted<Scalar>(x_, x);

  Scalar zeta = log(y_[index + 1] / y_[index])
	  / (x_[index + 1] - x_[index]);

  return y_[index] * exp(zeta * (x - x_[index]));
}

}

#endif
