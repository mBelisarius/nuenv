#include "nuenv/src/integrate/quadrature.hpp"

#include "nuenv/src/core/math.hpp"

#include <gtest/gtest.h>

namespace nuenv::test {

TEST(QuadratureTest, QuadratureATest)
{
    // Define a lambda function for the test
    const Lambda<double(double)> func = [](const double x)
    {
        return -exp(-sqrt2 * x) - 0.05 * exp(0.5 * cos(20 * pi * x)) + 1 + e / 20;
    };

    constexpr double a = 0.0;
    constexpr double b = 1.0;
    constexpr double tol = 1e-8;

    const double result = quadratureA(func, a, b, tol);
    constexpr double expected_result = 0.54754263323770047;

    EXPECT_NEAR(result, expected_result, tol);
}

TEST(QuadratureTest, QuadratureG1Test)
{
    const Lambda<double(double)> func = [](const double x) { return x; };

    constexpr double a = 0.0;
    constexpr double b = 1.0;

    const double result = quadratureG1(func, a, b);
    constexpr double expected_result = 1.0 / 2.0;

    EXPECT_NEAR(result, expected_result, 1e-8);
}

TEST(QuadratureTest, QuadratureG2Test)
{
    const Lambda<double(double)> func = [](const double x) { return x * x; };

    constexpr double a = 0.0;
    constexpr double b = 1.0;

    const double result = quadratureG2(func, a, b);
    constexpr double expected_result = 1.0 / 3.0;

    EXPECT_NEAR(result, expected_result, 1e-8);
}

TEST(QuadratureTest, QuadratureG3Test)
{
    const Lambda<double(double)> func = [](const double x) { return x * x * x; };

    constexpr double a = 0.0;
    constexpr double b = 1.0;

    const double result = quadratureG3(func, a, b);
    constexpr double expected_result = 1.0 / 4.0;

    EXPECT_NEAR(result, expected_result, 1e-8);
}

} // namespace nuenv::test
