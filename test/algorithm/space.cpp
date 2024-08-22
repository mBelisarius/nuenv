#include "nuenv/src/algorithm/space.hpp"

#include "nuenv/src/core/container.hpp"

#include "gtest/gtest.h"

namespace nuenv::test {

TEST(SpaceTest, LinearSpaceTest) {
  const VectorX<double> expVector = VectorX<double>::LinSpaced(5, 1.0, 5.0);
  const VectorX<double> vector = LinearSpace(1.0, 5.0, 5);

  ASSERT_TRUE(vector.isApprox(expVector));

  // Check for zero number of elements
  EXPECT_EQ(LinearSpace(1.0, 5.0, 0).size(), 0);
}

TEST(SpaceTest, GeometricSpaceTest) {
  VectorX<double> vector = geometricSpace(1.0, 100.0, 5);

  EXPECT_NEAR(vector[0], sqrt(1.0), 1e-8);
  EXPECT_NEAR(vector[1], sqrt(10.0), 1e-8);
  EXPECT_NEAR(vector[2], sqrt(100.0), 1e-8);
  EXPECT_NEAR(vector[3], sqrt(1000.0), 1e-8);
  EXPECT_NEAR(vector[4], sqrt(10000.0), 1e-8);

  // Check for zero number of elements
  EXPECT_EQ(geometricSpace(1.0, 5.0, 0).size(), 0);
}

TEST(SpaceTest, LogarithmicSpaceTest) {
  VectorX<double> vector = LogarithmicSpace(1.0, 100.0, 5);

  EXPECT_NEAR(vector[0], 1.0, 1e-8);
  EXPECT_NEAR(vector[4], 100.0, 1e-8);

  //Check for zero number of elements
  EXPECT_EQ(LogarithmicSpace(1.0, 5.0, 0).size(), 0);
}

} // namespace nuenv::test
