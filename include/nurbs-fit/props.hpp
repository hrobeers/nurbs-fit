#pragma once

#ifndef NURBSFIT_PROPS_HPP
#define NURBSFIT_PROPS_HPP

#include <vector>
#include <string>

namespace nurbsfit
{
  enum spline_t { QB, CB };

  struct props {
    // Fit props
    spline_t spline    = CB;
    std::vector<std::vector<std::string>> constraints;

    // Solver props
    double relax       = 0.05;
    double tol         = 0.001;
    size_t max_it      = 1024;

    // Export props
    bool origin        = false;
  };
}

#endif // NURBSFIT_PROPS_HPP
