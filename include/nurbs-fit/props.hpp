#pragma once

#ifndef NURBSFIT_PROPS_HPP
#define NURBSFIT_PROPS_HPP

namespace nurbsfit
{
  enum spline_t { QB, CB };

  struct props {
    // Fit props
    spline_t spline    = CB;

    // Solver props
    //double relax       = 0.3;
    double tol         = 0.001;
    size_t max_it      = 1024;

    // Export props
    bool origin        = false;
  };
}

#endif // NURBSFIT_PROPS_HPP
