#ifndef NUENV_CORE_CONTAINER_H_
#define NUENV_CORE_CONTAINER_H_

#include "Eigen/Dense"
#include <vector>

namespace nuenv {

template<class T>
using VectorT = std::vector<T>;

template<typename Scalar>
using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

template<typename Scalar, size_t size>
using VectorX_s = Eigen::Matrix<Scalar, size, 1>;

template<typename Scalar>
using Vector2X = Eigen::Matrix<Scalar, 2, 1>;

template<typename Scalar>
using Vector3X = Eigen::Matrix<Scalar, 3, 1>;

template<typename Scalar>
using MatrixSQX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

}

#endif
