#include "nuenv/src/Algorithm/Search.h"

#include "nuenv/src/Core/Container.h"
#include "nuenv/src/Core/Exception.h"

#include <gtest/gtest.h>

// TODO: Test the performance with extremely large arrays
// TODO: Test with duplicates

namespace nuenv::test {

VectorX_s<int, 10> sortedVectorInt = { -5, -4, -3, -2, -1, 0, 2, 3, 4, 5 };

TEST(SearchTest, BinarySearchFoundElementInt)
{
    const auto result = binarySearch<int>(sortedVectorInt, 0);
    constexpr size_t expected = 5;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchLowerBoundInt)
{
    const auto result = binarySearch<int>(sortedVectorInt, -5);
    constexpr size_t expected = 0;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchUpperBoundInt)
{
    const auto result = binarySearch<int>(sortedVectorInt, 5);
    constexpr size_t expected = 9 - 1;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchElementNotInArrayLowerInt)
{
    const auto result = binarySearch<int>(sortedVectorInt, -6);
    constexpr size_t expected = 0;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchElementNotInArrayUpperInt)
{
    const auto result = binarySearch<int>(sortedVectorInt, 6);
    constexpr size_t expected = 9 - 1;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchElementNotInArrayMiddleInt)
{
    const auto result = binarySearch<int>(sortedVectorInt, 1);
    constexpr size_t expected = 5;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchEmptyArrayInt)
{
    const VectorX<int> emptyVector;

    EXPECT_THROW(binarySearch<int>(emptyVector, 0), invalid_argument);
}

VectorX_s<float, 10> sortedVectorFloat = { -5.5f, -4.4f, -3.3f, -2.2f, -1.1f,
                                           0.0f, 2.2f, 3.3f, 4.4f, 5.5f };

TEST(SearchTest, BinarySearchFoundElementFloat)
{
    const auto result = binarySearch<float>(sortedVectorFloat, 0.0f);
    constexpr size_t expected = 5;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchLowerBoundFloat)
{
    const auto result = binarySearch<float>(sortedVectorFloat, -5.5f);
    constexpr size_t expected = 0;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchUpperBoundFloat)
{
    const auto result = binarySearch<float>(sortedVectorFloat, 5.5f);
    constexpr size_t expected = 9 - 1;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchElementNotInArrayLowerFloat)
{
    const auto result = binarySearch<float>(sortedVectorFloat, -6.6f);
    constexpr size_t expected = 0;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchElementNotInArrayUpperFloat)
{
    const auto result = binarySearch<float>(sortedVectorFloat, 6.6f);
    constexpr size_t expected = 9 - 1;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchElementNotInArrayMiddleFloat)
{
    const auto result = binarySearch<float>(sortedVectorFloat, 1.1f);
    constexpr size_t expected = 5;

    EXPECT_EQ(result, expected);
}

TEST(SearchTest, BinarySearchEmptyArrayFloat)
{
    const VectorX<float> emptyVector;

    EXPECT_THROW(binarySearch<float>(emptyVector, 0.0f), invalid_argument);
}

} // namespace nuenv::test
