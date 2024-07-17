#ifndef NUENV_INTEGRATE_ODESOLVER_H
#define NUENV_INTEGRATE_ODESOLVER_H

#include "nuenv/src/integrate/ode_solution.hpp"

namespace nuenv {

/**
 * @class OdeSolver
 * @brief Abstract base class for solving ordinary differential equations.
 *
 * Provides methods to iterate and solve a differential equation defined by a
 * given function and initial conditions.
 *
 * @tparam Scalar Scalar type of the numbers.
 * @tparam ScalarField Scalar field type.
 */
template<typename Scalar, typename ScalarField>
class OdeSolver
{
public:
    /**
        * @brief Iterate one step over the differential equation.
        *
        * @param t0 Initial time.
        * @param x0 Initial state.
        * @param step Step size.
        *
        * @return Solution of the differential equation at time 't0 + step'.
        */
    virtual ScalarField iter(Scalar t0,
                             ScalarField x0,
                             Scalar step) = 0;

    /**
     * @brief Solves a differential equation defined by the given function and
     *  initial state.
     *
     * Solves a differential equation defined by the given equations
     * 'OdeSystem' and initial state 'x0'.
     * It provides an interface to different solvers that implement the
     * 'Solver' interface.
     *
     * @param t_eval Time values at which to evaluate the solution.
     * @param x0 Initial state.
     * @param stopEvent Lambda function that returns 'true' if an event
     *  to stop the solver has occurred, 'false' otherwise.
     *
     * @return Solution to the differential equation at the specified times.
     */
    virtual OdeSolution<Scalar, ScalarField> solve(
        const VectorX<Scalar>& t_eval,
        ScalarField x0,
        Lambda<bool(ScalarField)> stopEvent) = 0;

    virtual ~OdeSolver() = default;
};

}

#endif
