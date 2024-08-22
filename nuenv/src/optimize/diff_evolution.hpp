#ifndef NUENV_OPTIMIZE_DIFFEVOLUTION_H_
#define NUENV_OPTIMIZE_DIFFEVOLUTION_H_

#include "nuenv/src/core/container.hpp"
#include "nuenv/src/core/lambda.hpp"
#include "nuenv/src/core/math.hpp"
#include "nuenv/src/core/random.hpp"
#include "nuenv/src/optimize/population.hpp"

namespace nuenv {

/**
 * @class DiffEvolution
 *
 * @brief Finds the global minimum of a multivariate function using the
 *  Differential Evolution (DE) method.
 *
 * The differential evolution method is a stochastic population based method
 * that is useful for global optimization problems. It does not use gradient
 * methods to find the minimum, and can SearchSorted large areas of candidate space,
 * but often requires larger numbers of function evaluations than conventional
 * gradient-based techniques.
 *
 * The constraints are handled using the approach by [Lampinen, J., 2002].
 *
 * @tparam Scalar Scalar type of the numbers.
 * @tparam ScalarField Scalar field type.
 *
 * @see Lampinen, J., A constraint handling approach for the differential evolution algorithm. Proceedings of the 2002 Congress on Evolutionary Computation. CEC’02 (Cat. No. 02TH8600). Vol. 2. IEEE, 2002.
 */
template<typename Scalar, typename ScalarField>
class DiffEvolution {
 public:
  DiffEvolution(Lambda<Scalar(ScalarField)> fitness,
				VectorT<Vector2X<Scalar>> bounds,
				Lambda<Scalar(ScalarField)> constraints =
				[](ScalarField /*x*/) { return 0.0; },
				size_t maxiter = 1000,
				size_t breakafter = 100,
				size_t popsize = 16,
				Vector2X<Scalar> mutation = Vector2X<Scalar> {0.5, 0.8},
				Scalar recombination = 0.8);

  ScalarField optimize();

 protected:
  void ensureBounds(ScalarField& candidate);

  VectorX<size_t> generateSamples(Index n, size_t candidate_index);

  void generate();

  void mutateRand1(ScalarField& candidate, const VectorX<size_t>& samples);

  void mutateBest1(ScalarField& candidate, const VectorX<size_t>& samples);

  ScalarField crossover(size_t pos);

  bool iterate();

 private:
  Lambda<Scalar(ScalarField)> m_fitness;
  VectorT<Vector2X<Scalar>> m_bounds;
  size_t m_dim;
  Lambda<Scalar(ScalarField)> m_constraints;
  size_t m_maxiter;
  size_t m_breakafter;
  size_t m_popsize;
  Scalar m_recombination;

  Scalar m_mutation;
  Population<Scalar, ScalarField> m_population;
  random_device m_rand_device;
  minstd_rand m_gen;
  uniform_int_distribution<size_t> m_rand_dim;
  uniform_int_distribution<size_t> m_rand_popsize;
  VectorT<uniform_real_distribution<Scalar>> m_rand_bounds;
  uniform_real_distribution<Scalar> m_rand_mutation;
  uniform_real_distribution<Scalar> m_rand_recomb;
};

/**
 * Constructor.
 *
 * @param fitness Lambda function that represents the objective function to be
 *  minimized.
 * @param bounds Array of pairs representing the bounds for each dimension of
 *  the optimization space.
 * @param constraints Lambda function for constraints, should return a
 *  non-positive number is the constraints are satisfied.
 *  Default: no constraints.
 * @param maxiter Maximum number of iterations. Default is 1000.
 * @param breakafter Number of iterations without improvement to break early.
 *  Default is 100.
 * @param popsize Population size given in the form as 'popsize * dim'.
 *  Default is 16.
 * @param mutation Range for the mutation factor. Default is [0.5, 0.8].
 * @param recombination Recombination probability. Default is 0.8.
 */
template<typename Scalar, typename ScalarField>
DiffEvolution<Scalar, ScalarField>::DiffEvolution(
	Lambda<Scalar(ScalarField)> fitness,
	VectorT<Vector2X<Scalar>> bounds,
	Lambda<Scalar(ScalarField)> constraints,
	const size_t maxiter,
	const size_t breakafter,
	const size_t popsize,
	Vector2X<Scalar> mutation,
	Scalar recombination)
	: m_fitness(fitness),
	  m_bounds(bounds),
	  m_dim(bounds.size()),
	  m_constraints(constraints),
	  m_maxiter(maxiter),
	  m_breakafter(breakafter),
	  m_popsize(popsize * m_dim),
	  m_recombination(recombination),
	  m_population(m_popsize),
	  m_gen(m_rand_device()),
	  m_rand_dim(0, m_dim),
	  m_rand_popsize(0, m_popsize),
	  m_rand_bounds(m_dim),
	  m_rand_mutation(mutation[0], mutation[1]),
	  m_rand_recomb(0.0, 1.0) {
  // TODO: Throw exception
  assert(m_dim > 0);
  assert(m_popsize > 3);

  for (size_t i = 0; i < m_dim; i++) {
	m_rand_bounds[i] = uniform_real_distribution<Scalar>(m_bounds[i][0],
														 m_bounds[i][1]);
  }
}

/**
 * @brief Optimize the objective function.
 *
 * @return Scalar field that minimizes the objective function.
 */
template<typename Scalar, typename ScalarField>
ScalarField DiffEvolution<Scalar, ScalarField>::optimize() {
  generate();

  size_t last_best = 0;

  for (size_t iter = 0; iter < m_maxiter; iter++) {
	// Iterate to next generation
	if (iterate()) { last_best = iter; }

	if (iter - last_best >= m_breakafter) { break; }
  }

  return m_population.best_individual.value;
}

/**
 * @brief Ensure that the candidate solution remains within the specified bounds.
 *
 * Checks and enforces that each component of the candidate solution remains
 * within the specified bounds for each dimension, modifying the value in-place.
 * It provides three strategies for handling out-of-bounds values:
 *  1. Truncate: The value is truncated to the bounds if it exceeds them.
 *  2. Reflection: The value is reflected back into the bounds if it exceeds them.
 *  3. Random: A random value within the bounds is assigned if the value is out of bounds.
 *
 * @param candidate Candidate solution to be checked and adjusted.
 */
template<typename Scalar, typename ScalarField>
void DiffEvolution<Scalar, ScalarField>::ensureBounds(ScalarField& candidate) {
  // TODO: Select method
  // TODO: Vectorize
  for (size_t i = 0; i < m_dim; i++) {
	// Truncate
	// candidate[i] = max(candidate[i], m_bounds[i][0]);
	// candidate[i] = min(candidate[i], m_bounds[i][1]);

	// Reflection
	// candidate[i] = candidate[i] < m_bounds[i][0] ?
	//         2.0 * m_bounds[i][0] - candidate[i] : candidate[i];
	// candidate[i] = candidate[i] > m_bounds[i][1] ?
	//         2.0 * m_bounds[i][1] - candidate[i] : candidate[i];

	// Random
	if (candidate[i] < m_bounds[i][0] || candidate[i] > m_bounds[i][1]) {
	  candidate[i] = m_rand_bounds[i](m_gen);
	}
  }
}

/**
 * @brief Generate a set of unique random indexes.
 *
 * Generates 'n' unique random indexes used for selecting individuals from the
 * population. The candidate at 'candidate_index' should not be included in the
 * selection.
 *
 * @param n Number of unique random indexes to generate.
 * @param candidate_index Index of the candidate to exclude from the selection.
 *
 * @return Array of 'n' unique random indexes.
 */
template<typename Scalar, typename ScalarField>
VectorX<size_t> DiffEvolution<Scalar, ScalarField>::generateSamples(
	const Index n,
	const size_t candidate_index) {
  VectorX<size_t> indexes(n);

  bool success;

  for (Index i = 0; i < n; i++) {
	do {
	  // Use min(m_popsize - 1, rand) to fix C++ random where it may
	  // occasionally return b
	  indexes[i] = min(m_popsize - 1, m_rand_popsize(m_gen));

	  success = (indexes[i] != candidate_index);
	  for (Index j = 0; j < i; j++) {
		success &= (indexes[i] != indexes[j]);
	  }
	} while (!success);
  }

  return indexes;
}

/**
 * @brief Initialize the population and evaluate constraints and fitness.
 *
 * Initializes the population of candidate solutions by randomly generating
 * values within the specified bounds for each dimension. It then evaluates the
 * constraints and fitness values for each candidate solution. The method takes
 * into account only the candidate solutions that are feasible (satisfy
 * constraints) when computing fitness values.
 */
template<typename Scalar, typename ScalarField>
void DiffEvolution<Scalar, ScalarField>::generate() {
  // TODO: Change generation method to Halton
  // TODO: Parallelize
  for (size_t i = 0; i < m_popsize; i++) {
	for (size_t j = 0; j < m_dim; j++) {
	  m_population[i][j] = m_rand_bounds[j](m_gen);
	}

	m_population[i].constraint = m_constraints(m_population[i].value);

	// Only need to work out population energies for those that are feasible
	if (m_population[i].constraint <= 0.0) {
	  m_population[i].fitness = m_fitness(m_population[i].value);
	}
  }

  m_population.best();
}

/**
 * @brief Mutate a candidate solution using the "rand/1" strategy.
 *
 * Mutates the candidate solution using three random candidate solutions by
 * adding a scaled difference between two of them to the third. The indexes of
 * the three random candidate solutions to be used in the mutation are specified
 * by 'samples'.
 *
 * @param candidate Candidate solution to be mutated.
 * @param samples Array of indexes specifying the three candidate solutions to
 *  be used in the mutation.
 */
template<typename Scalar, typename ScalarField>
void DiffEvolution<Scalar, ScalarField>::mutateRand1(ScalarField& candidate,
													 const VectorX<size_t>& samples) {
  candidate = m_population[samples[2]].value
	  + m_mutation * (m_population[samples[0]].value
		  - m_population[samples[1]].value);
}

/**
 * @brief Mutate a candidate solution using the "best/1" strategy.
 *
 * Mutates the candidate solution by adding a scaled difference between two
 * random candidate solutions to the current best candidate solution. The
 * indexes of the two random candidate solutions to be used in the mutation are
 * specified by 'samples.'
 *
 * @param candidate The candidate solution to be mutated.
 * @param samples A vector of indexes specifying the three candidate solutions
 *  to be used in the mutation.
 */
template<typename Scalar, typename ScalarField>
void DiffEvolution<Scalar, ScalarField>::mutateBest1(ScalarField& candidate,
													 const VectorX<size_t>& samples) {
  candidate = m_population.best_individual.value
	  + m_mutation * (m_population[samples[0]].value
		  - m_population[samples[1]].value);
}

/**
 * @brief Generate a trial candidate solution through crossover and mutation.
 *
 * Generates a trial candidate solution for the 'i-th' candidate in the
 * population by applying crossover and mutation. The crossover strategy depends
 * on the feasibility of the candidate solution. For feasible solutions
 * (constraint <= 0), the "best/1" mutation strategy is used, and for infeasible
 * solutions, the "rand/1" strategy is applied. The method also checks a
 * recombination probability to decide whether a dimension should be replaced by
 * the mutated value.
 *
 * @param pos Index of the candidate solution in the population.
 *
 * @return Trial candidate solution generated through crossover and mutation.
 */
template<typename Scalar, typename ScalarField>
ScalarField DiffEvolution<Scalar, ScalarField>::crossover(size_t pos) {
  ScalarField trial = m_population[pos].value;

  if (m_population[pos].constraint <= 0.0) {
	const VectorX_s<size_t, 2> samples = generateSamples(2, pos);
	mutateBest1(trial, samples);
  } else {
	const VectorX_s<size_t, 3> samples = generateSamples(3, pos);
	mutateRand1(trial, samples);
  }

  // TODO: Optimize - dont calculate the mutation if its not gonna be replaced
  for (size_t j = 0; j < m_dim; j++) {
	if (!(m_rand_recomb(m_gen) < m_recombination || j == m_rand_dim(m_gen))) {
	  trial[j] = m_population[pos][j];
	}
  }

  ensureBounds(trial);

  return trial;
}

/**
 * @brief Perform one iteration of the algorithm.
 *
 * Performs one iteration of the algorithm. It generates trial candidate
 * solutions for each candidate in the population using the crossover and
 * mutation process. It checks the feasibility of trial solutions, computes
 * their fitness values, and compares them with the current population. If a
 * better solution is found, it updates the population.
 *
 * @return True if a better solution is found during the iteration,
 *  false otherwise.
 */
template<typename Scalar, typename ScalarField>
bool DiffEvolution<Scalar, ScalarField>::iterate() {
  bool found_better = false;

  m_mutation = m_rand_mutation(m_gen);

  Individual<Scalar, ScalarField> trial;

  // TODO: Parallelize
  for (size_t i = 0; i < m_popsize; i++) {
	trial.value = crossover(i);
	trial.constraint = m_constraints(trial.value);

	// Only need to work out population energies for those that are feasible
	if (trial.constraint <= 0.0) {
	  trial.fitness = m_fitness(trial.value);
	}

	if (trial.isBetter(m_population[i])) {
	  m_population[i] = trial;
	  found_better = true;
	}
  }

  m_population.best();

  return found_better;
}

}

#endif
