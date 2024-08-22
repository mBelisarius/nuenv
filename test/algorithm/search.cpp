#include "nuenv/src/algorithm/search.hpp"

#include "nuenv/src/core/container.hpp"

#include <gtest/gtest.h>

// TODO: Test the performance with extremely large arrays
// TODO: Test with duplicates

namespace nuenv::test {

const VectorX_s<int, 10> kVectorIntSorted = {-5, -4, -3, -2, -1, 0, 2, 3, 4, 5};

TEST(SearchTest, SearchSortedFoundElementInt) {
  const auto result = SearchSorted<int>(kVectorIntSorted, 0);
  constexpr Index expected = 5;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedLowerBoundInt) {
  const auto result = SearchSorted<int>(kVectorIntSorted, -5);
  constexpr Index expected = 0;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedUpperBoundInt) {
  const auto result = SearchSorted<int>(kVectorIntSorted, 5);
  constexpr Index expected = 9;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedElementNotInArrayLowerInt) {
  const auto result = SearchSorted<int>(kVectorIntSorted, -6);
  constexpr Index expected = 0;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedElementNotInArrayUpperInt) {
  const auto result = SearchSorted<int>(kVectorIntSorted, 6);
  constexpr Index expected = 9;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedElementNotInArrayMiddleInt) {
  const auto result = SearchSorted<int>(kVectorIntSorted, 1);
  constexpr Index expected = 5;

  EXPECT_EQ(result, expected);
}

const VectorX_s<float, 10> kVectorFloatSorted = {-5.5f, -4.4f, -3.3f, -2.2f, -1.1f,
												 0.0f, 2.2f, 3.3f, 4.4f, 5.5f};

TEST(SearchTest, SearchSortedFoundElementFloat) {
  const auto result = SearchSorted<float>(kVectorFloatSorted, 0.0f);
  constexpr Index expected = 5;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedLowerBoundFloat) {
  const auto result = SearchSorted<float>(kVectorFloatSorted, -5.5f);
  constexpr Index expected = 0;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedUpperBoundFloat) {
  const auto result = SearchSorted<float>(kVectorFloatSorted, 5.5f);
  constexpr Index expected = 9;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedElementNotInArrayLowerFloat) {
  const auto result = SearchSorted<float>(kVectorFloatSorted, -6.6f);
  constexpr Index expected = 0;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedElementNotInArrayUpperFloat) {
  const auto result = SearchSorted<float>(kVectorFloatSorted, 6.6f);
  constexpr Index expected = 9;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchSortedElementNotInArrayMiddleFloat) {
  const auto result = SearchSorted<float>(kVectorFloatSorted, 1.1f);
  constexpr Index expected = 5;

  EXPECT_EQ(result, expected);
}

VectorX_s<int, 10> kVectorIntUnsorted = {3, 0, -3, 5, -4, -1, 2, -5, 4, -2};

TEST(SearchTest, SearchUnsortedFoundElementInt) {
  const auto result = SearchUnsorted<int>(kVectorIntUnsorted, 0);
  constexpr Index expected = 1;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchUnsortedFirstElemInt) {
  const auto result = SearchUnsorted<int>(kVectorIntUnsorted, 3);
  constexpr Index expected = 0;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchUnsortedLastElemInt) {
  const auto result = SearchUnsorted<int>(kVectorIntUnsorted, -2);
  constexpr Index expected = 9;

  EXPECT_EQ(result, expected);
}

TEST(SearchTest, SearchUnsortedElementNotInArrayInt) {
  const auto result = SearchUnsorted<int>(kVectorIntUnsorted, 1);
  constexpr Index expected = 10;

  EXPECT_EQ(result, expected);
}


} // namespace nuenv::test
