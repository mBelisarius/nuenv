#include "nuenv/src/integrate/rk4.hpp"

#include "nuenv/src/core/container.hpp"
#include "nuenv/src/core/math.hpp"
#include "nuenv/src/integrate/ode_solution.hpp"

#include <gtest/gtest.h>

namespace nuenv::test {

TEST(RK4Test, FirstOrderIter)
{
    class ODESystem
    {
    public:
        double operator()(const double t, const double x) const
        {
            return x - pow2(t) + 1;
        }
    };

    constexpr ODESystem ode_system;
    Rk4<double, double, ODESystem> rk4(ode_system);

    constexpr double t0 = 0.0;
    constexpr double x0 = 0.5;
    constexpr double step = 0.1;

    const auto result = rk4.iter(t0, x0, step);
    constexpr double expected_result = 0.6574144;

    EXPECT_NEAR(result, expected_result, 1e-4);
}

TEST(RK4Test, FirstOrderSolve)
{
    class ODESystem
    {
    public:
        double operator()(const double t, const double x) const
        {
            return x - pow2(t) + 1;
        }
    };

    constexpr ODESystem ode_system;
    Rk4<double, double, ODESystem> rk4(ode_system);

    VectorX_s<double, 6> t_eval = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
    constexpr double x0 = 0.5;

    const auto result = rk4.solve(t_eval, x0);
    VectorX_s<double, 6> expected_result = { 0.5000000, 0.6574144, 0.8292983,
                                             1.0150701, 1.2140869, 1.4256384 };

    EXPECT_EQ(result.t.size(), t_eval.size());

    for (int i = 0; i < result.t.size(); i++)
    {
        EXPECT_NEAR(result.t[i], t_eval[i], 1e-4);
        EXPECT_NEAR(result.x[i], expected_result[i], 1e-4);
    }
}

TEST(RK4Test, HigherOrderIter)
{
    /** y'' - 2 y' + 2 y = exp(2 t) sin (t) */
    class ODESystem
    {
    public:
        Vector2X<double> operator()(const double t, const Vector2X<double>& x) const
        {
            double u1 = x[1];
            double u2 = exp(2.0 * t) * sin(t) - 2.0 * x[0] + 2.0 * x[1];
            return { u1, u2 };
        }
    };

    constexpr ODESystem ode_system;
    Rk4<double, Vector2X<double>, ODESystem> rk4(ode_system);

    constexpr double t0 = 0.0;
    const Vector2X<double> x0 = { -0.4, -0.6 };
    constexpr double step = 0.1;

    auto result = rk4.iter(t0, x0, step);
    Vector2X<double> expected_result = { -0.46173297, -0.6316304 };

    for (size_t i = 0, max = x0.size(); i < max; i++)
    {
        EXPECT_NEAR(result[i], expected_result[i], 1e-4);
    }
}

TEST(RK4Test, HigherOrderSolve)
{
    class ODESystem
    {
    public:
        Vector2X<double> operator()(const double t, const Vector2X<double>& x) const
        {
            double u1 = x[1];
            double u2 = exp(2.0 * t) * sin(t) - 2.0 * x[0] + 2.0 * x[1];
            return { u1, u2 };
        }
    };

    constexpr ODESystem ode_system;
    Rk4<double, Vector2X<double>, ODESystem> rk4(ode_system);

    VectorX_s<double, 6> t_eval = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
    const Vector2X<double> x0 = { -0.4, -0.6 };

    const auto result = rk4.solve(t_eval, x0);
    VectorT<Vector2X<double>> expected_result =
            {
                    { -0.40000000, -0.6000000 },
                    { -0.46173297, -0.6316304 },
                    { -0.52555905, -0.6401478 },
                    { -0.58860005, -0.6136630 },
                    { -0.64661028, -0.5365821 },
                    { -0.69356395, -0.3887395 },
            };

    EXPECT_EQ(result.t.size(), t_eval.size());

    for (int i = 0; i < result.t.size(); i++)
    {
        EXPECT_NEAR(result.t[i], t_eval[i], 1e-4);
        for (size_t j = 0, max = x0.size(); j < max; j++)
        {
            EXPECT_NEAR(result.x[i][j], expected_result[i][j], 1e-4);
        }
    }
}

TEST(RK4Test, SystemIter)
{
    /**
     * y1' = -4 y1 + 3 y2 + 6
     * y2' = -2.4 y1 + 1.6 y2 + 3.6
     */
    class ODESystem
    {
    public:
        Vector2X<double> operator()(double t, const Vector2X<double>& x) const
        {
            double u1 = -4.0 * x[0] + 3 * x[1] + 6.0;
            double u2 = -2.4 * x[0] + 1.6 * x[1] + 3.6;
            return { u1, u2 };
        }
    };

    constexpr ODESystem ode_system;
    Rk4<double, Vector2X<double>, ODESystem> rk4(ode_system);

    constexpr double t0 = 0.0;
    const Vector2X<double> x0 = { 0.0, 0.0 };
    constexpr double step = 0.1;

    auto result = rk4.iter(t0, x0, step);
    Vector2X<double> expected_result = { 0.53826391, 0.31963204 };

    for (size_t i = 0, max = x0.size(); i < max; i++)
    {
        EXPECT_NEAR(result[i], expected_result[i], 1e-4);
    }
}

TEST(RK4Test, SystemSolve)
{
    class ODESystem
    {
    public:
        Vector2X<double> operator()(double t, const Vector2X<double>& x) const
        {
            double u1 = -4.0 * x[0] + 3 * x[1] + 6.0;
            double u2 = -2.4 * x[0] + 1.6 * x[1] + 3.6;
            return { u1, u2 };
        }
    };

    constexpr ODESystem ode_system;
    Rk4<double, Vector2X<double>, ODESystem> rk4(ode_system);

    VectorX_s<double, 6> t_eval = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
    const Vector2X<double> x0 = { 0.0, 0.0 };

    const auto result = rk4.solve(t_eval, x0);
    auto expected_result1 = [](const double t)
    {
        return -3.375 * exp(-2.0 * t) + 1.875 * exp(-0.4 * t) + 1.5;
    };
    auto expected_result2 = [](const double t)
    {
        return -2.25 * exp(-2.0 * t) + 2.25 * exp(-0.4 * t);
    };


    EXPECT_EQ(result.t.size(), t_eval.size());

    for (int i = 0; i < result.t.size(); i++)
    {
        EXPECT_NEAR(result.t[i], t_eval[i], 1e-4);
        EXPECT_NEAR(result.x[i][0], expected_result1(t_eval[i]), 1e-4);
        EXPECT_NEAR(result.x[i][1], expected_result2(t_eval[i]), 1e-4);
    }
}

} // namespace nuenv::test
