#ifndef NUENV_INTEGRATE_RK4_H_
#define NUENV_INTEGRATE_RK4_H_

#include "nuenv/src/core/container.hpp"
#include "nuenv/src/core/lambda.hpp"
#include "nuenv/src/integrate/ode_solver.hpp"
#include "nuenv/src/integrate/ode_solution.hpp"

namespace nuenv {

#define RK4_TEMPLATE template<typename Scalar, typename ScalarField, class ODESystem>
#define RK4_EXTENSION Rk4<Scalar, ScalarField, ODESystem>
#define RK4_STATIC_CONST_SCALAR RK4_TEMPLATE const Scalar RK4_EXTENSION

/**
 * @class Rk4
 *
 * @brief Fourth-order Runge-Kutta method to solve ordinary differential
 *  equations.
 */
RK4_TEMPLATE
class Rk4 final : public OdeSolver<Scalar, ScalarField> {
 public:
  explicit Rk4(const ODESystem& ode);

  ScalarField iter(Scalar t0, ScalarField x0, Scalar step);

  OdeSolution<Scalar, ScalarField> solve(
	  const VectorX<Scalar>& t_eval,
	  ScalarField x0,
	  Lambda<bool(ScalarField)> stopEvent = [](ScalarField /*x*/) {
		return false;
	  });

 private:
  static const Scalar
	  c2, c3,
	  a21,
	  a31, a32,
	  a41, a42, a43,
	  b1, b2, b3, b4;

  mutable ScalarField k1, k2, k3, k4;

  ODESystem m_ode;
};

RK4_STATIC_CONST_SCALAR::c2 = 1.0 / 2.0;
RK4_STATIC_CONST_SCALAR::c3 = 1.0 / 2.0;

RK4_STATIC_CONST_SCALAR::a21 = 1.0 / 2.0;

RK4_STATIC_CONST_SCALAR::a31 = 0.0;
RK4_STATIC_CONST_SCALAR::a32 = 1.0 / 2.0;

RK4_STATIC_CONST_SCALAR::a41 = 0.0;
RK4_STATIC_CONST_SCALAR::a42 = 0.0;
RK4_STATIC_CONST_SCALAR::a43 = 1.0;

RK4_STATIC_CONST_SCALAR::b1 = 1.0 / 6.0;
RK4_STATIC_CONST_SCALAR::b2 = 1.0 / 3.0;
RK4_STATIC_CONST_SCALAR::b3 = 1.0 / 3.0;
RK4_STATIC_CONST_SCALAR::b4 = 1.0 / 6.0;

RK4_TEMPLATE
RK4_EXTENSION::Rk4(const ODESystem& ode)
	: m_ode(ode) {}

RK4_TEMPLATE
ScalarField RK4_EXTENSION::iter(Scalar t0,
								ScalarField x0,
								Scalar step) {
  k1 = m_ode(t0, x0);
  k2 = m_ode(t0 + c2 * step, x0 + step * (a21 * k1));
  k3 = m_ode(t0 + c3 * step, x0 + step * (a31 * k1 + a32 * k2));
  k4 = m_ode(t0 + step, x0 + step * (a41 * k1 + a42 * k2 + a43 * k3));

  ScalarField x1 = x0 + step * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4);

  return x1;
}

RK4_TEMPLATE
OdeSolution<Scalar, ScalarField>
RK4_EXTENSION::solve(const VectorX<Scalar>& t_eval,
					 ScalarField x0,
					 Lambda<bool(ScalarField)> stopEvent) {
  Scalar step;
  size_t size = t_eval.size();
  VectorX<ScalarField> x(size);
  x[0] = x0;

  size_t i = 1;
  for (; i < size; i++) {
	step = t_eval[i] - t_eval[i - 1];
	x[i] = iter(t_eval[i - 1], x[i - 1], step);

	if (stopEvent(x[i])) {
	  i++;
	  break;
	}
  }

  OdeSolution<Scalar, ScalarField> sol(t_eval.segment(0, i),
									   x.segment(0, i),
									   i);
  return sol;
}

}

#endif
