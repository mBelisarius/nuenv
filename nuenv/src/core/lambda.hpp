#ifndef NUENV_CORE_LAMBDA_H_
#define NUENV_CORE_LAMBDA_H_

#include <functional>

namespace nuenv {

template<typename Signature>
using Lambda = std::function<Signature>;

}

#endif
