/*
 * This file is part of the nurbs-fit project (https://github.com/hrobeers/nurbs-fit).
 * Copyright (c) 2021 hrobeers.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#ifndef NURBSFIT_CURVEPROC_HPP
#define NURBSFIT_CURVEPROC_HPP

#include <functional>
#include <boost/math/special_functions/pow.hpp> 

#include "hrlib/io/vertexio.hpp"

namespace bm = boost::math;

namespace nurbsfit
{
  inline
  std::function<hrlib::vertex<2>(double)> to_func(const std::vector<hrlib::vertex<2>> &curve) {
    switch (curve.size()) {
    case 3:
      // Quadratic Bezier: (1-t)^2 P0 + 2t(1-t) Pc + u^2 P1 = P(t)
      return [curve](double t) -> hrlib::vertex<2> { return {
          bm::pow<2>(1-t)*curve[0][0] + 2*t*(1-t)*curve[1][0] + bm::pow<2>(t)*curve[2][0],
          bm::pow<2>(1-t)*curve[0][1] + 2*t*(1-t)*curve[1][1] + bm::pow<2>(t)*curve[2][1]
        };};
    case 4:
      // TODO
      assert(false);
      break;
    default:
      assert(false);
      break;
    }
  }
  inline
  std::function<hrlib::vertex<2>(double)> to_deriv(const std::vector<hrlib::vertex<2>> &curve) {
    switch (curve.size()) {
    case 3:
      // Quadratic Bezier: 2(1-t)(P1-P0) + 2t(P2-P1) = P'(t)
      return [curve](double t) -> hrlib::vertex<2> { return {
          2*(1-t)*(curve[1][0]-curve[0][0]) + 2*t*(curve[2][0]-curve[1][0]),
          2*(1-t)*(curve[1][1]-curve[0][1]) + 2*t*(curve[2][1]-curve[1][1])
        };};
    case 4:
      // TODO
      assert(false);
      break;
    default:
      assert(false);
      break;
    }
  }

  inline
  double t_prec_arclen(const std::function<hrlib::vertex<2>(double)> &f_curve, double perc_acrlen) {
    const double RES = 1024;

    std::vector<hrlib::vertex<2>> pts(RES+1);
    for (size_t i=0; i<=RES; i++)
      pts[i] = f_curve(i/RES);

    std::vector<double> arclen(RES+1);
    arclen[0] = 0;
    for (size_t i=1; i<=RES; i++)
      arclen[i] = std::sqrt(bm::pow<2>(pts[i][0]-pts[i-1][0]) + bm::pow<2>(pts[i][1]-pts[i-1][1])) + arclen[i-1];

    double length = arclen.back();
    double target = perc_acrlen * length;
    size_t idx=0;
    while (arclen[idx]<target) idx++;

    return idx/RES;
  }

  inline
  std::vector<hrlib::vertex<2>> center_origin(const std::vector<hrlib::vertex<2>> &curve) {
    auto f_curve = to_func(curve);
    auto f_deriv = to_deriv(curve);
    double t_center = t_prec_arclen(f_curve, 0.5);
    auto p_center = f_curve(t_center);

    // translate center to origin
    std::vector<hrlib::vertex<2>> to(curve.size());
    std::transform(curve.cbegin(), curve.cend(), to.begin(),
                   [&p_center](auto p) -> hrlib::vertex<2> { return { p[0]-p_center[0], p[1]-p_center[1] }; });

    // rotate
    auto d = f_deriv(t_center);
    auto a = std::atan(d[1]/d[0]);
    auto s = std::sin(-a);
    auto c = std::cos(-a);
    std::vector<hrlib::vertex<2>> ro(curve.size());
    std::transform(to.cbegin(), to.cend(), ro.begin(),
                   [s,c](auto p) -> hrlib::vertex<2> { return { p[0]*c-p[1]*s, p[0]*s+p[1]*c }; });


    // translate x coordinate back
    std::transform(ro.cbegin(), ro.cend(), ro.begin(),
                   [&p_center](auto p) -> hrlib::vertex<2> { return { p[0]+p_center[0], p[1] }; });

    return ro;
  }

  inline
  std::vector<hrlib::vertex<2>> to_cubic(const std::vector<hrlib::vertex<2>> &curve) {
    switch (curve.size()) {
    case 3:
      {
        std::vector<hrlib::vertex<2>> cubic(4);
        cubic[0] = curve[0];
        cubic[1] = { 1./3*(curve[0][0]+2.*curve[1][0]), 1./3*(curve[0][1]+2.*curve[1][1]) };
        cubic[2] = { 1./3*(curve[2][0]+2.*curve[1][0]), 1./3*(curve[2][1]+2.*curve[1][1]) };
        cubic[3] = curve[2];
        return cubic;
      }
    case 4:
      return curve;
    default:
      assert(false);
      break;
    }
  }
}

#endif //NURBSFIT_CURVEPROC_HPP
