#ifndef NUENV_OPTIMIZE_POPULATION_H_
#define NUENV_OPTIMIZE_POPULATION_H_

#include "nuenv/src/core/container.hpp"
#include "nuenv/src/optimize/individual.hpp"

namespace nuenv {

/**
 * @class Population
 *
 * @brief Represents a population of individuals in population-based
 *  optimization algorithms.
 *
 * This class provides methods to manage a population of individuals,
 * including creation, access and selection of the fittest.
 *
 * @tparam Scalar Scalar type of the numbers.
 * @tparam ScalarField Scalar field type.
 */
template<typename Scalar, typename ScalarField>
class Population {
 public:
  explicit Population(size_t popsize);

  explicit Population(
	  VectorT<Individual<Scalar, ScalarField>> individuals);

  Individual<Scalar, ScalarField>& operator[](size_t pos);

  Individual<Scalar, ScalarField>& best();

  size_t popsize;
  VectorT<Individual<Scalar, ScalarField>> individuals;
  Individual<Scalar, ScalarField> best_individual;
};

/**
 * @brief Constructor.
 *
 * @param popsize Population size.
 */
template<typename Scalar, typename ScalarField>
Population<Scalar, ScalarField>::Population(size_t popsize)
	: popsize(popsize),
	  individuals(popsize),
	  best_individual(individuals[0]) {}

/**
 * @brief Constructor.
 *
 * @param individuals Array of individuals of the population.
 */
template<typename Scalar, typename ScalarField>
Population<Scalar, ScalarField>::Population(
	VectorT<Individual<Scalar, ScalarField>> individuals)
	: popsize(individuals.size()),
	  individuals(individuals),
	  best_individual(best()) {}

/**
 * @brief Access an individual.
 *
 * @param pos Position of the individual to access.
 *
 * @return Reference to the individual at the given position.
 *
 * @throw std::out_of_range if the position is out of range.
 */
template<typename Scalar, typename ScalarField>
Individual<Scalar, ScalarField>& Population<Scalar, ScalarField>::operator[](size_t pos) {
  return individuals[pos];
}

/**
 * @brief Get a reference to the best individual in the population.
 *
 * @return Reference to the best individual.
 */
template<typename Scalar, typename ScalarField>
Individual<Scalar, ScalarField>& Population<Scalar, ScalarField>::best() {
  for (size_t i = 0; i < popsize; i++) {
	if (individuals[i].isBetter(best_individual)) {
	  best_individual = individuals[i];
	}
  }

  return best_individual;
}

}

#endif
