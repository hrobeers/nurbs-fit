#pragma once

#ifndef NURBSFIT_PROPS_HPP
#define NURBSFIT_PROPS_HPP

namespace nurbsfit
{
  struct props {
    // Solver props
    double relax       = 0.3;
    double tol         = 0.01;
    size_t max_it      = 1024;
  };
}

#endif // NURBSFIT_PROPS_HPP
