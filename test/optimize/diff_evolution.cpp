#include "nuenv/src/optimize/diff_evolution.hpp"

#include "nuenv/src/core/math.hpp"

#include <gtest/gtest.h>

namespace nuenv::test {

constexpr size_t REPEAT = 1000;

double ackley(const Vector2X<double>& x) {
  const double arg1 = -0.2 * sqrt(0.5 * (Pow2(x[0]) + Pow2(x[1])));
  const double arg2 = 0.5 * (cos(2.0 * pi * x[0]) + cos(2.0 * pi * x[1]));
  return -20.0 * exp(arg1) - exp(arg2) + 20.0 + e;
}

double schaffer2(const Vector2X<double>& x) {
  const double arg1 = Pow2(sin(Pow2(x[0]) - Pow2(x[1]))) - 0.5;
  const double arg2 = Pow2(1.0 + 0.001 * (Pow2(x[0]) + Pow2(x[1])));
  return 0.5 + arg1 / arg2;
}

template<size_t N = 2>
double rosenbrock(const VectorX_s<double, N>& x) {
  double sum = 0.0;
  for (size_t i = 0; i < N - 1; i++) {
	sum += 100.0 * Pow2(x[i + 1] - Pow2(x[i])) + Pow2(1 - x[i]);
  }

  return sum;
}

TEST(DiffEvolutionTest, FitnessAckley) {
  const VectorT<Vector2X<double>> bounds = {{-5.0, 5.0},
											{-5.0, 5.0}};

  auto fitness = [](const Vector2X<double>& x) {
	return ackley(x);
  };

  DiffEvolution<double, Vector2X<double>> opt(fitness, bounds);

  for (size_t i = 0; i < REPEAT; i++) {
	auto sol = opt.optimize();

	EXPECT_NEAR(sol[0], 0.0, 1e-4);
	EXPECT_NEAR(sol[1], 0.0, 1e-4);
  }
}

TEST(DiffEvolutionTest, FitnessSchaffer2) {
  const VectorT<Vector2X<double>> bounds = {{-100.0, 100.0},
											{-100.0, 100.0}};

  auto fitness = [](const Vector2X<double>& x) {
	return schaffer2(x);
  };

  DiffEvolution<double, Vector2X<double>> opt(fitness, bounds);

  for (size_t i = 0; i < REPEAT; i++) {
	auto sol = opt.optimize();

	EXPECT_NEAR(sol[0], 0.0, 1e-4);
	EXPECT_NEAR(sol[1], 0.0, 1e-4);
  }
}

TEST(DiffEvolutionTest, FitnessRosenbrock10) {
  const VectorT<Vector2X<double>> bounds = {{-100.0, 100.0},
											{-100.0, 100.0},
											{-100.0, 100.0},
											{-100.0, 100.0},
											{-100.0, 100.0},
											{-100.0, 100.0},
											{-100.0, 100.0},
											{-100.0, 100.0},
											{-100.0, 100.0},
											{-100.0, 100.0}};

  auto fitness = [](const VectorX_s<double, 10>& x) {
	return rosenbrock<10>(x);
  };

  DiffEvolution<double, VectorX_s<double, 10>> opt(fitness, bounds);

  for (size_t i = 0; i < REPEAT; i++) {
	auto sol = opt.optimize();

	for (const auto s : sol) {
	  EXPECT_NEAR(s, 1.0, 1e-4);
	}
  }
}

TEST(DiffEvolutionTest, FitnessRosenbrock10Constraint) {
  VectorT<Vector2X<double>> bounds = {{-100.0, 100.0},
									  {-100.0, 100.0},
									  {-100.0, 100.0},
									  {-100.0, 100.0},
									  {-100.0, 100.0},
									  {-100.0, 100.0},
									  {-100.0, 100.0},
									  {-100.0, 100.0},
									  {-100.0, 100.0},
									  {-100.0, 100.0}};

  auto fitness = [](const VectorX_s<double, 10>& x) {
	return rosenbrock<10>(x);
  };

  auto constraint = [](const VectorX_s<double, 10>& x) {
	return x.sum() - 10.0;
  };

  DiffEvolution<double, VectorX_s<double, 10>> opt(fitness, bounds,
												   constraint);

  for (size_t i = 0; i < REPEAT; i++) {
	auto sol = opt.optimize();

	for (auto s : sol) {
	  EXPECT_NEAR(s, 1.0, 1e-4);
	}
  }
}

} // namespace nuenv::test
