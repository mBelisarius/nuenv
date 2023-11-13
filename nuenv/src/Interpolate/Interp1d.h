#ifndef NUENV_INTERPOLATE_INTERP1D_H
#define NUENV_INTERPOLATE_INTERP1D_H

#include "nuenv/src/Algorithm/Search.h"
#include "nuenv/src/Core/Container.h"
#include "nuenv/src/Core/Math.h"

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
class Interp1d
{
public:
    Interp1d(const VectorX<Scalar>& x,
             const VectorX<Scalar>& y,
             bool check_bounds = true,
             bool assume_sorted = false);

    Scalar linear(Scalar x);

    Scalar exponential(Scalar x);

private:
    VectorX<Scalar> m_x;
    VectorX<Scalar> m_y;
    size_t m_size;
    bool m_check_bounds;
};

/**
 * Constructs the interpolator.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param x Array of x-values representing the independent variable.
 * @param y Array of y-values representing the dependent variable.
 * @param check_bounds Indicates whether to check if new points are within the
 *  domain bounds of 'x'.
 * @param assume_sorted Indicates whether 'x' is assumed to be pre-sorted.
 */
template<typename Scalar>
Interp1d<Scalar>::Interp1d(const VectorX<Scalar>& x,
                           const VectorX<Scalar>& y,
                           const bool check_bounds,
                           const bool assume_sorted)
    : m_x(x), m_y(y), m_size(x.size()), m_check_bounds(check_bounds)
{
    if (x.size() != y.size())
    {
        throw invalid_argument("Arrays 'x' and 'y' must have same size");
    }

    // TODO: Implement sorting.
    if (!assume_sorted)
    {
        throw invalid_argument("Array 'x' must be pre-sorted");
    }
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
Scalar Interp1d<Scalar>::linear(Scalar x)
{
    if (m_check_bounds && ((x < m_x(0)) || (x > m_x(m_size - 1))))
    {
        throw invalid_argument("'x' is out of bounds");
    }

    size_t index = binarySearch<Scalar>(m_x, x);

    return m_y[index]
           + ((m_y[index + 1] - m_y[index]) / (m_x[index + 1] - m_x[index]))
           * (x - m_x[index]);
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
Scalar Interp1d<Scalar>::exponential(Scalar x)
{
    if (m_check_bounds && ((x < m_x(0)) || (x > m_x(m_size - 1))))
    {
        throw invalid_argument("'x' is out of bounds");
    }

    size_t index = binarySearch<Scalar>(m_x, x);

    Scalar zeta = log(m_y[index + 1] / m_y[index])
                  / (m_x[index + 1] - m_x[index]);

    return m_y[index] * exp(zeta * (x - m_x[index]));
}

}

#endif
