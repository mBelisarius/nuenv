#include "nuenv/src/Interpolate/Interp1d.h"

#include <gtest/gtest.h>

namespace nuenv::test {

TEST(Interp1dTest, LinearSorted)
{
    const VectorX_s<double, 3> x = { 1.0, 2.0, 3.0 };
    const VectorX_s<double, 3> y = { 10.0, 20.0, 30.0 };
    Interp1d<double> interp(x, y, true, true);

    EXPECT_NEAR(interp.linear(1.5), 15.0, 1e-8);
    EXPECT_NEAR(interp.linear(2.5), 25.0, 1e-8);
}

TEST(Interp1dTest, LinearSortedOOB)
{
    const VectorX_s<double, 3> x = { 1.0, 2.0, 3.0 };
    const VectorX_s<double, 3> y = { 10.0, 20.0, 30.0 };
    Interp1d<double> interp(x, y, false, true);

    EXPECT_NEAR(interp.linear(4.0), 40.0, 1e-8);
}

TEST(Interp1dTest, LinearNonSorted)
{
    const VectorX_s<double, 3> x = { 3.0, 2.0, 1.0 };
    const VectorX_s<double, 3> y = { 30.0, 20.0, 10.0 };
    Interp1d<double> interp(x, y, true, false);

    EXPECT_NEAR(interp.linear(1.5), 15.0, 1e-8);
    EXPECT_NEAR(interp.linear(2.5), 25.0, 1e-8);
}

TEST(Interp1dTest, ExponentialSorted)
{
    const VectorX_s<double, 2> x = { 1.0, 2.0 };
    const VectorX_s<double, 2> y = { 10.0, 20.0 };
    Interp1d<double> interp(x, y, true, true);

    EXPECT_NEAR(interp.exponential(1.5), sqrt(200.0), 1e-8);
}

TEST(Interp1dTest, ExponentialSortedOOB)
{
    const VectorX_s<double, 2> x = { 1.0, 2.0 };
    const VectorX_s<double, 2> y = { 10.0, 20.0 };
    Interp1d<double> interp(x, y, false, true);

    EXPECT_NEAR(interp.exponential(3.0), 40.0, 1e-8);
}

TEST(Interp1dTest, ExponentialNonSorted)
{
    const VectorX_s<double, 2> x = { 2.0, 1.0 };
    const VectorX_s<double, 2> y = { 20.0, 10.0 };
    Interp1d<double> interp(x, y, true, false);

    EXPECT_NEAR(interp.exponential(1.5), sqrt(200.0), 1e-8);
}

TEST(Interp1dTest, ArraysDifferentSize)
{
    const VectorX_s<double, 3> x = { 1.0, 2.0, 3.0 };
    const VectorX_s<double, 2> y = { 10.0, 20.0 };

    EXPECT_THROW(Interp1d<double> interp(x, y, true, true),
                 invalid_argument);
}

} // namespace nuenv::test
