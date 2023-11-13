#include "nuenv/src/Algorithm/Space.h"

#include "nuenv/src/Core/Container.h"
#include "nuenv/src/Core/Exception.h"

#include "gtest/gtest.h"

namespace nuenv::test {

TEST(SpaceTest, LinearSpaceTest)
{
    const VectorX<double> expVector = VectorX<double>::LinSpaced(5, 1.0, 5.0);
    const VectorX<double> vector = linearSpace(1.0, 5.0, 5);

    ASSERT_TRUE(vector.isApprox(expVector));

    // Check for zero number of elements
    EXPECT_EQ(linearSpace(1.0, 5.0, 0).size(), 0);
}

TEST(SpaceTest, GeometricSpaceTest)
{
    VectorX<double> vector = geometricSpace(1.0, 100.0, 5);

    EXPECT_NEAR(vector[0], sqrt(1.0), 1e-8);
    EXPECT_NEAR(vector[1], sqrt(10.0), 1e-8);
    EXPECT_NEAR(vector[2], sqrt(100.0), 1e-8);
    EXPECT_NEAR(vector[3], sqrt(1000.0), 1e-8);
    EXPECT_NEAR(vector[4], sqrt(10000.0), 1e-8);

    // Check for invalid arguments
    EXPECT_THROW(geometricSpace(0.0, 5.0, 5), invalid_argument);
    EXPECT_THROW(geometricSpace(5.0, 0.0, 5), invalid_argument);
    EXPECT_THROW(geometricSpace(0.0, 0.0, 5), invalid_argument);

    // Check for zero number of elements
    EXPECT_EQ(geometricSpace(1.0, 5.0, 0).size(), 0);
}

TEST(SpaceTest, LogarithmicSpaceTest)
{
    VectorX<double> vector = logarithmicSpace(1.0, 100.0, 5);

    EXPECT_NEAR(vector[0], 1.0, 1e-8);
    EXPECT_NEAR(vector[4], 100.0, 1e-8);

    // Check for invalid arguments
    EXPECT_THROW(logarithmicSpace(0.0, 5.0, 5), invalid_argument);
    EXPECT_THROW(logarithmicSpace(5.0, 0.0, 5), invalid_argument);
    EXPECT_THROW(logarithmicSpace(0.0, 0.0, 5), invalid_argument);

    //Check for zero number of elements
    EXPECT_EQ(logarithmicSpace(1.0, 5.0, 0).size(), 0);
}

} // namespace nuenv::test
