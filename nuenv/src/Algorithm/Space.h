#ifndef NUENV_ALGORITHM_SPACE_H
#define NUENV_ALGORITHM_SPACE_H

#include "nuenv/src/Core/Container.h"
#include "nuenv/src/Core/Exception.h"
#include "nuenv/src/Core/Math.h"

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
VectorX<Scalar> linearSpace(Scalar start, Scalar stop, size_t num)
{
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
VectorX<Scalar> logarithmicSpace(Scalar start, Scalar stop, size_t num)
{
    if (start == 0.0 || stop == 0.0)
    {
        throw invalid_argument("Argument must not be 0");
    }

    const Scalar base = 10.0;

    Scalar a_log = log(start) / log(base);
    Scalar b_log = log(stop) / log(base);

    VectorX<Scalar> vector = linearSpace(a_log, b_log, num);

    for (size_t i = 0; i < num; i++)
    {
        vector[i] = pow(base, vector[i]);
    }

    return vector;
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
VectorX<Scalar> geometricSpace(Scalar start, Scalar stop, size_t num)
{
    if (start == 0.0 || stop == 0.0)
    {
        throw invalid_argument("Argument must not be 0");
    }

    Scalar r = pow(stop / start, 1.0 / (static_cast<Scalar>(num) - 1.0));

    VectorX<Scalar> vector(num);

    for (size_t i = 0; i < num; i++)
    {
        vector[i] = start * pow(r, i);
    }

    return vector;
}

}

#endif
