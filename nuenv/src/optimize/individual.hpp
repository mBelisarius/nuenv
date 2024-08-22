#ifndef NUENV_OPTIMIZE_INDIVIDUAL_H_
#define NUENV_OPTIMIZE_INDIVIDUAL_H_

#include "nuenv/src/core/ctypes.hpp"

namespace nuenv {

/**
 * @class Individual
 *
 * @brief Individual contained in population-based optimization algorithms.
 *
 * Stores information about an individual's value, fitness, and constraint.
 * It provides methods to manipulate and compare individuals.
 *
 * @tparam Scalar Scalar type of the numbers.
 * @tparam ScalarField Scalar field type.
 */
template<typename Scalar, typename ScalarField>
class Individual {
 public:
  Individual();

  explicit Individual(ScalarField value);

  Individual(ScalarField value, Scalar fitness, Scalar constraint);

  Individual& operator=(const Individual& other);

  Scalar& operator[](size_t pos);

  bool isBetter(const Individual& other);

  ScalarField value;
  Scalar fitness;
  Scalar constraint;
};

/**
 * @brief Default constructor.
 */
template<typename Scalar, typename ScalarField>
Individual<Scalar, ScalarField>::Individual()
	: value(),
	  fitness(numeric_limits<Scalar>::max()),
	  constraint(numeric_limits<Scalar>::max()) {}

/**
 * @brief Constructor.
 *
 * @param value Value of the Individual.
 */
template<typename Scalar, typename ScalarField>
Individual<Scalar, ScalarField>::Individual(ScalarField value)
	: value(value),
	  fitness(numeric_limits<Scalar>::max()),
	  constraint(numeric_limits<Scalar>::max()) {}

/**
 * @brief Constructor.
 *
 * @param value Value of the Individual.
 * @param fitness Fitness value.
 * @param constraint Constraint value.
 */
template<typename Scalar, typename ScalarField>
Individual<Scalar, ScalarField>::Individual(ScalarField value,
											Scalar fitness,
											Scalar constraint)
	: value(value),
	  fitness(fitness),
	  constraint(constraint) {}

/**
 * @brief Assignment operator.
 */
template<typename Scalar, typename ScalarField>
Individual<Scalar, ScalarField>&
Individual<Scalar, ScalarField>::operator=(const Individual& other) {
  value = other.value;
  fitness = other.fitness;
  constraint = other.constraint;

  return *this;
}

/**
 * @brief Bracket operator to access the Individual's value.
 *
 * @param pos Position index.
 *
 * @return Value at 'pos'.
 */
template<typename Scalar, typename ScalarField>
Scalar& Individual<Scalar, ScalarField>::operator[](size_t pos) {
  return value[pos];
}

/**
 * @brief Checks if the current Individual is better than the given
 *  Individual.
 *
 * @param other Individual to compare with.
 *
 * @return True if the current Individual is better, false otherwise.
 */
template<typename Scalar, typename ScalarField>
bool Individual<Scalar, ScalarField>::isBetter(const Individual& other) {
  if (((other.constraint <= 0.0) && (other.fitness <= fitness)) ||
	  ((other.constraint <= 0.0) && (constraint > 0.0)) ||
	  ((other.constraint > 0.0) && (other.constraint < constraint))) {
	return false;
  }

  return true;
}

}

#endif
