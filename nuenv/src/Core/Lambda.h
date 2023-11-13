#ifndef NUENV_CORE_LAMBDA_H
#define NUENV_CORE_LAMBDA_H

#include <functional>

namespace nuenv {

template<typename Signature>
using Lambda = std::function<Signature>;

}

#endif
