#ifndef NUENV_INTEGRATE_SOLUTION_H
#define NUENV_INTEGRATE_SOLUTION_H

#include "nuenv/src/core/container.hpp"

namespace nuenv {

/**
 * @brief Collection of solutions of a differential equation.
 *
 * Collection of solutions of a differential equation containing the
 * time at which the ODE where evaluated 't', its result 'x' and the
 * size of the solution 'size'.
 *
 * @tparam Scalar Scalar type of the numbers.
 * @tparam ScalarField Scalar field type.
 */
template<typename Scalar, typename ScalarField>
struct OdeSolution
{
    OdeSolution(const VectorX<Scalar>& t_,
                const VectorX<ScalarField>& x_,
                const size_t size_)
        : t(t_), x(x_), size(size_) {}

    const VectorX<Scalar> t;
    const VectorX<ScalarField> x;
    const size_t size;
};

}

#endif
