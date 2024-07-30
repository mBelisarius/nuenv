#include "nuenv/src/interpolate/interp1d.hpp"

#include <gtest/gtest.h>

namespace nuenv::test {

TEST(Interp1dTest, LinearSorted) {
  const VectorX_s<double, 3> x = {1.0, 2.0, 3.0};
  const VectorX_s<double, 3> y = {10.0, 20.0, 30.0};
  Interp1d<double> interp(x, y, true);

  EXPECT_NEAR(interp.linear(1.5), 15.0, 1e-8);
  EXPECT_NEAR(interp.linear(2.5), 25.0, 1e-8);
}

TEST(Interp1dTest, LinearSortedOOB) {
  const VectorX_s<double, 3> x = {1.0, 2.0, 3.0};
  const VectorX_s<double, 3> y = {10.0, 20.0, 30.0};
  Interp1d<double> interp(x, y, false);

  EXPECT_NEAR(interp.linear(4.0), y[2], 1e-8);
}

TEST(Interp1dTest, ExponentialSorted) {
  const VectorX_s<double, 3> x = {1.0, 2.0, 3.0};
  const VectorX_s<double, 3> y = {10.0, 20.0, 30.0};
  Interp1d<double> interp(x, y, true);

  EXPECT_NEAR(interp.exponential(1.5), sqrt(200.0), 1e-8);
}

TEST(Interp1dTest, ExponentialSortedOOB) {
  const VectorX_s<double, 3> x = {1.0, 2.0, 3.0};
  const VectorX_s<double, 3> y = {10.0, 20.0, 30.0};
  Interp1d<double> interp(x, y, false);

  EXPECT_NEAR(interp.exponential(4.0), y[2], 1e-8);
}

} // namespace nuenv::test
