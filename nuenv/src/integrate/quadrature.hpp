#ifndef NUENV_INTEGRATE_QUADRATURE_H_
#define NUENV_INTEGRATE_QUADRATURE_H_

#include "nuenv/core"

namespace nuenv {

namespace internal {

template<typename Scalar>
struct ConstsG10K21 {
  static constexpr Index kNg = 10;
  static constexpr Index kNk = 21;

  static const VectorX_s<Scalar, kNg> kWg;
  static const VectorX_s<Scalar, kNg> kXg;
  static const VectorX_s<Scalar, kNk> kWk;
  static const VectorX_s<Scalar, kNk> kXk;
};

template<typename Scalar>
const VectorX_s<Scalar, ConstsG10K21<Scalar>::kNg> ConstsG10K21<Scalar>::kWg =
	{
		6.667134430868813759356880989333179e-02,
		1.494513491505805931457763396576973e-01,
		2.190863625159820439955349342281632e-01,
		2.692667193099963550912269215694694e-01,
		2.955242247147528701738929946513383e-01,
		2.955242247147528701738929946513383e-01,
		2.692667193099963550912269215694694e-01,
		2.190863625159820439955349342281632e-01,
		1.494513491505805931457763396576973e-01,
		6.667134430868813759356880989333179e-02
	};

template<typename Scalar>
const VectorX_s<Scalar, ConstsG10K21<Scalar>::kNg> ConstsG10K21<Scalar>::kXg =
	{
		-9.739065285171717200779640120844521e-01,
		-8.650633666889845107320966884234930e-01,
		-6.794095682990244062343273651148736e-01,
		-4.333953941292471907992659431657842e-01,
		-1.488743389816312108848260011297200e-01,
		1.488743389816312108848260011297200e-01,
		4.333953941292471907992659431657842e-01,
		6.794095682990244062343273651148736e-01,
		8.650633666889845107320966884234930e-01,
		9.739065285171717200779640120844521e-01
	};

template<typename Scalar>
const VectorX_s<Scalar, ConstsG10K21<Scalar>::kNk> ConstsG10K21<Scalar>::kWk =
	{
		1.169463886737187427806439606219205e-02,
		3.255816230796472747881897245938976e-02,
		5.475589657435199603138130024458018e-02,
		7.503967481091995276704314091619001e-02,
		9.312545458369760553506546508336634e-02,
		1.093871588022976418992105903258050e-01,
		1.234919762620658510779581098310742e-01,
		1.347092173114733259280540017717068e-01,
		1.427759385770600807970942731387171e-01,
		1.477391049013384913748415159720680e-01,
		1.494455540029169056649364683898212e-01,
		1.477391049013384913748415159720680e-01,
		1.427759385770600807970942731387171e-01,
		1.347092173114733259280540017717068e-01,
		1.234919762620658510779581098310742e-01,
		1.093871588022976418992105903258050e-01,
		9.312545458369760553506546508336634e-02,
		7.503967481091995276704314091619001e-02,
		5.475589657435199603138130024458018e-02,
		3.255816230796472747881897245938976e-02,
		1.169463886737187427806439606219205e-02
	};

template<typename Scalar>
const VectorX_s<Scalar, ConstsG10K21<Scalar>::kNk> ConstsG10K21<Scalar>::kXk =
	{
		-9.956571630258080807355272806890028e-01,
		-9.739065285171717200779640120844521e-01,
		-9.301574913557082260012071800595083e-01,
		-8.650633666889845107320966884234930e-01,
		-7.808177265864168970637175783450424e-01,
		-6.794095682990244062343273651148736e-01,
		-5.627571346686046833390000992726941e-01,
		-4.333953941292471907992659431657842e-01,
		-2.943928627014601981311266031038656e-01,
		-1.488743389816312108848260011297200e-01,
		0.000000000000000000000000000000000e+00,
		1.488743389816312108848260011297200e-01,
		2.943928627014601981311266031038656e-01,
		4.333953941292471907992659431657842e-01,
		5.627571346686046833390000992726941e-01,
		6.794095682990244062343273651148736e-01,
		7.808177265864168970637175783450424e-01,
		8.650633666889845107320966884234930e-01,
		9.301574913557082260012071800595083e-01,
		9.739065285171717200779640120844521e-01,
		9.956571630258080807355272806890028e-01
	};

}

/**
 * @brief Compute a definite integral.
 *
 * Integrate func from 'a' to 'b' (possibly infinite interval) using the
 * adaptive Gauss-Kronrod 10-21 technique.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param func Function or method to integrate. It should take a single scalar
 *  argument and return a scalar value.
 * @param a Lower limit of integration.
 * @param b Upper limit of integration.
 * @param tol Absolute error tolerance. Default is 6e-6.
 *
 * @return Integral of 'func' from 'a' to 'b'.
 */
template<typename Scalar>
Scalar quadratureA(Lambda<Scalar(Scalar)> func,
				   Scalar a,
				   Scalar b,
				   Scalar tol = 6e-6) {
  using internal::ConstsG10K21;

  Scalar aux, integral = 0.0;

  Scalar integral_g = 0.0;
  // TODO: Vectorize
  for (Index i = 0; i < ConstsG10K21<Scalar>::kNg; i++) {
	aux = ((b - a) * ConstsG10K21<Scalar>::kXg[i] + (b + a)) / 2.0;
	integral_g += ConstsG10K21<Scalar>::kWg[i] * func(aux);
  }

  integral_g *= (b - a) / 2.0;

  Scalar integral_k = 0.0;
  // TODO: Vectorize
  for (Index i = 0; i < ConstsG10K21<Scalar>::kNk; i++) {
	aux = ((b - a) * ConstsG10K21<Scalar>::kXk[i] + (b + a)) / 2.0;
	integral_k += ConstsG10K21<Scalar>::kWk[i] * func(aux);
  }

  integral_k *= (b - a) / 2.0;

  Scalar error = abs(integral_k - integral_g);

  if (error < tol) {
	integral += integral_k;
  } else {
	Index n = ceil(1.0 + log2(error / tol));
	Scalar h = (b - a) / static_cast<Scalar>(n);

	// TODO: Parallelize
	for (Index i = 0; i < n; i++) {
	  integral += quadratureA(func, a + i * h, a + (i + 1) * h, tol);
	}
  }

  return integral;
}

/**
 * @brief Performs the Gauss-Legendre quadrature with 1 point for a given
 *  function.
 *
 * This function computes the integral of a given function 'func' using
 * the Gauss-Legendre quadrature method with 1 point.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param func Function or method to integrate. It should take a single
 *  scalar argument and return a scalar value.
 * @param a Lower limit of integration.
 * @param b Upper limit of integration.
 *
 * @return Integral of 'func' over the interval [a, b].
 */
template<typename Scalar>
Scalar quadratureG1(Lambda<Scalar(Scalar)> func, Scalar a, Scalar b) {
  Scalar h2 = b - a;
  Scalar integral = h2 * func(a + h2 / 2.0);

  return integral;
}

/**
 * @brief Performs the Gauss-Legendre quadrature with 2 points for a given
 *  function.
 *
 * This function computes the integral of a given function 'func' using
 * the Gauss-Legendre quadrature method with 2 points.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param func Function or method to integrate. It should take a single
 *  scalar argument and return a scalar value.
 * @param a Lower limit of integration.
 * @param b Upper limit of integration.
 *
 * @return Integral of 'func' over the interval [a, b].
 */
template<typename Scalar>
Scalar quadratureG2(Lambda<Scalar(Scalar)> func, Scalar a, Scalar b) {
  static constexpr Index n = 2;
  static const VectorX_s<Scalar, n> w = {1.0, 1.0};
  static const VectorX_s<Scalar, n> x = {-0.5773502692, 0.5773502692};

  Scalar aux, integral = 0.0;
  for (Index i = 0; i < n; i++) {
	aux = ((b - a) * x[i] + (b + a)) / 2.0;
	integral += w[i] * func(aux);
  }

  integral *= (b - a) / 2.0;

  return integral;
}

/**
 * @brief Performs the Gauss-Legendre quadrature with 3 points for a given
 *  function.
 *
 * This function computes the integral of a given function 'func' using
 * the Gauss-Legendre quadrature method with 3 points.
 *
 * @tparam Scalar Scalar type of the numbers.
 *
 * @param func Function or method to integrate. It should take a single
 *  scalar argument and return a scalar value.
 * @param a Lower limit of integration.
 * @param b Upper limit of integration.
 *
 * @return Integral of 'func' over the interval [a, b].
 */
template<typename Scalar>
Scalar quadratureG3(Lambda<Scalar(Scalar)> func, Scalar a, Scalar b) {
  static constexpr Index n = 3;
  static const VectorX_s<Scalar, n> w = {
	  0.5555555556,
	  0.8888888889,
	  0.5555555556
  };
  static const VectorX_s<Scalar, n> x = {
	  -0.7745966692,
	  0.0,
	  0.7745966692
  };

  Scalar aux, integral = 0.0;
  for (Index i = 0; i < n; i++) {
	aux = ((b - a) * x[i] + (b + a)) / 2.0;
	integral += w[i] * func(aux);
  }

  integral *= (b - a) / 2.0;

  return integral;
}

}

#endif
