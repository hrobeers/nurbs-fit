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

#ifndef NURBSFIT_CURVEFIT_HPP
#define NURBSFIT_CURVEFIT_HPP

#include <cassert>
#include <optional>
#include <math.h>

#include "hrlib/io/vertexio.hpp"

namespace nurbsfit
{
  std::vector<hrlib::vertex<2>> fit_qb(const std::vector<hrlib::vertex<2>> &P,
                                       double relax=0.05, double tol=0.001, size_t max_it=1024);
  std::vector<hrlib::vertex<2>> fit_cb(const std::vector<hrlib::vertex<2>> &P,
                                       std::optional<std::array<double,2>> tangents = std::nullopt,
                                       double relax=0.05, double tol=0.001, size_t max_it=1024);

  template<typename Tprops, int size>
  std::optional<std::array<double,size>> get_tangents(const Tprops &props) {
    if (props.constraints.size()==0)
      return std::nullopt;

    std::array<double,size> tangents;
    for (size_t i=0; i<size; i++)
      tangents[i] = NAN;

    size_t i=0;
    for (auto constraint : props.constraints) {
      switch(constraint.front()[0]) {
      case 't':
        tangents[i] = NAN;
        break;
      case 'a':
        tangents[i] = atof(constraint[1].c_str()) * M_PI / 180.0;
        break;
      }
      if (++i == size)
        break;
    }
    return tangents;
  }

  template<typename Tprops>
  std::vector<hrlib::vertex<2>> fit_qb(const std::vector<hrlib::vertex<2>> &P, const Tprops &props) {
    return fit_qb(P, props.relax, props.tol, props.max_it);
  }

  template<typename Tprops>
  std::vector<hrlib::vertex<2>> fit_cb(const std::vector<hrlib::vertex<2>> &P, const Tprops &props) {
    return fit_cb(P, get_tangents<Tprops,2>(props), props.relax, props.tol, props.max_it);
  }
}

#endif //NURBSFIT_CURVEFIT_HPP
