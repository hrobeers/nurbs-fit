#pragma once

#ifndef NURBSFIT_CURVEFIT_HPP
#define NURBSFIT_CURVEFIT_HPP

#include <cassert>

#include "hrlib/io/vertexio.hpp"

namespace nurbsfit
{
  inline
  std::vector<hrlib::vertex<2>> fit_qb(const std::vector<hrlib::vertex<2>> &target) {
    assert(target.size()==4); // TODO warn or switch to least squares?

    auto A = target[0];
    auto B = target[1];
    auto C = target[2];
    auto D = target[3];

    decltype(A) P = {.0,.0};

    return {A, P, D};
  }
}

#endif //NURBSFIT_CURVEFIT_HPP
