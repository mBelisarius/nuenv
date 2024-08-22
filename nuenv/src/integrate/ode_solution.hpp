#ifndef NUENV_INTEGRATE_SOLUTION_H_
#define NUENV_INTEGRATE_SOLUTION_H_

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
struct OdeSolution {
  OdeSolution(const VectorX<Scalar>& _t,
			  const VectorX<ScalarField>& _x,
			  const size_t _size)
	  : t(_t), x(_x), size(_size) {}

  const VectorX<Scalar> t;
  const VectorX<ScalarField> x;
  const size_t size;
};

}

#endif
