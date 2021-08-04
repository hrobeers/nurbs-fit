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

#include "nurbs-fit/curveproc.hpp"
#include "hrlib/io/vertexio.hpp"
#include "invert-matrix.hpp"

#include <boost/math/special_functions/pow.hpp>
namespace bm = boost::math;

namespace nurbsfit
{
  template<typename Tprops>
  std::vector<hrlib::vertex<2>> fit_qb(const std::vector<hrlib::vertex<2>> &P, const Tprops &props) {
    typedef hrlib::vertex<2> vertex;

    const size_t Dim = 2;
    size_t pcnt = P.size();
    size_t ucnt = pcnt-2;
    const size_t ccnt = 1;

    assert(pcnt>=4);

    // Order of rows:
    // P-equations, linearizations, extra equation
    size_t rows = Dim*ucnt + 2*ucnt + 1;
    // Order of cols:
    // u, u^3, (1-u)^3, Pc
    size_t cols = 3*ucnt + ccnt*Dim;

    const double init_relax = 16./props.max_it; // props.relax;
    double relax = init_relax;

    // Bezier points
    auto P0 = P.front();
    auto P1 = P.back();
    std::array<vertex, ccnt> Pc;

    std::vector<double> u(ucnt);
    for (size_t i=0; i<ucnt; i++)
      u[i] = (i+1.)/pcnt;
    std::vector<double> u_prev = u;

    size_t it = 0;
    while (it++<props.max_it) {
      using namespace boost::numeric::ublas;

      // Ax=b
      matrix<double> A = zero_matrix<double>(rows, cols);
      vector<double> b = zero_vector<double>(rows);

      size_t r = 0;

      // linearization of t^2 around u
      // Taylor: t^2 ~ u^2 + 2u (t-u)
      // Matrix: 2u*t -1*t^2 = u^2

      // t^3
      for (size_t i=0; i<ucnt; i++) {
        A(r,i) = 2*u[i];
        A(r,i+ucnt) = -1;
        b(r) = bm::pow<2>(u[i]);
        r++;
      }

      // linearization of (1-t)^2 around u
      // Taylor: (1-t)^2 ~ (1-u)^2 -2(1-u)(t-u)
      // Matrix: (2u-2)*t -1*(1-t)^2 = (u^2-1)

      // (1-t)^2
      for (size_t i=0; i<ucnt; i++) {
        A(r,i) = 2*u[i]-2;
        A(r,i+2*ucnt) = -1;
        b(r) = bm::pow<2>(u[i]) - 1;
        r++;
      }

      // Pu equations
      for (size_t i=0; i<ucnt; i++)
        for (size_t d=0; d<Dim; d++) {
          A(r,ucnt+i) = P1[d];
          A(r,2*ucnt+i) = P0[d];
          A(r,3*ucnt+d) = 2*(1-u[i])*u[i]; // Pc
          b(r) = P[i+1][d];
          r++;
        }

      decltype(A) At = identity_matrix<double>(rows);
      if (rows!=cols) At = trans(A);

      // Solve
      decltype(A) Ainv(cols,cols);
      A = prod(At,A);
      if (invert(A, Ainv)) {
        u_prev = u;
      } else {
        // We need to relax more
        relax *= 0.75;
        u = u_prev;
        continue;
      }

      Ainv = prod(Ainv,At);
      auto x = prod(Ainv,b);

      for (size_t i=0; i<ucnt; i++) {
        u[i] = (x(i)-u[i])*relax + u[i];
        u[i] = std::max(0., std::min(u[i], 1.));
      }
      relax = std::min(init_relax, relax*1.01); // carefully regenerate the relaxation

      Pc[0] = { x(3*ucnt), x(3*ucnt+1) };

      // Check the residuals
      std::vector<double> residuals;
      for (size_t i=0; i<ucnt; i++)
        residuals.push_back(std::abs(u[i]-x(i)));
      if (std::all_of(residuals.cbegin(), residuals.cend(), [&props](double r){ return r<props.tol; }))
          break;
    }

    if (it>props.max_it)
      std::cerr << "[NOT CONVERGED] ";
    else
      std::cerr << "[SUCCESS] ";
    std::cerr << "iterations: " << it << std::endl;
    return {P0, Pc[0], P1};
  }

  template<typename Tprops>
  std::vector<hrlib::vertex<2>> fit_cb(const std::vector<hrlib::vertex<2>> &P, const Tprops &props) {
    typedef hrlib::vertex<2> vertex;

    const size_t Dim = 2;
    size_t pcnt = P.size();
    size_t ucnt = pcnt-2;
    const size_t ccnt = 2;

    assert(pcnt>=5);

    // Order of rows:
    // P-equations, linearizations, extra equation
    size_t rows = Dim*ucnt + 2*ucnt + 1;
    // Order of cols:
    // u, u^3, (1-u)^3, Pc1, Pc2
    size_t cols = 3*ucnt + ccnt*Dim;

    // Hack to stop in time before exploding
    // TODO re-evaluate the initial relaxation
    const double init_relax = 16./props.max_it; // props.relax;
    double relax = init_relax;

    // Bezier points
    auto P0 = P.front();
    auto P1 = P.back();
    std::array<vertex, ccnt> Pc;

    std::vector<double> u(ucnt);
    for (size_t i=0; i<ucnt; i++)
      u[i] = (i+1.)/pcnt;
    std::vector<double> u_prev = u;

    size_t it = 0;
    while (it++<props.max_it) {
      using namespace boost::numeric::ublas;

      // Ax=b
      matrix<double> A = zero_matrix<double>(rows, cols);
      vector<double> b = zero_vector<double>(rows);

      size_t r = 0;

      // linearization of t^3 around u
      // Taylor: t^3 ~ u^3 + 3u^2 (t-u)
      // Matrix: 3u^2*t -1*t^3 = 2u^3

      // t^3
      for (size_t i=0; i<ucnt; i++) {
        A(r,i) = 3*bm::pow<2>(u[i]);
        A(r,i+ucnt) = -1;
        b(r) = 2*bm::pow<3>(u[i]);
        r++;
      }

      // linearization of (1-t)^3 around u
      // Taylor: (1-t)^3 ~ (1-u)^3 -3(1-u)^2(t-u)
      // Matrix: (3u^2-6u+3)*t +1*(1-t)^3 = 2u^3-3u^2+1

      // (1-t)^3
      for (size_t i=0; i<ucnt; i++) {
        A(r,i) = 3*bm::pow<2>(u[i])-6*u[i]+3;
        A(r,i+2*ucnt) = 1;
        b(r) = 2*bm::pow<3>(u[i]) - 3*bm::pow<2>(u[i]) + 1;
        r++;
      }

      // Pu equations
      for (size_t i=0; i<ucnt; i++)
        for (size_t d=0; d<Dim; d++) {
          A(r,ucnt+i) = P1[d];
          A(r,2*ucnt+i) = P0[d];
          A(r,3*ucnt+d) = 3*bm::pow<2>(1-u[i])*u[i]; // Pc1
          A(r,3*ucnt+Dim+d) = 3*(1-u[i])*bm::pow<2>(u[i]); // Pc2
          b(r) = P[i+1][d];
          r++;
        }

      // Enforce Pc1+Pc2 = P0+P1 (extra equation)
      // Avoids a control point explosion
      A(r,3*ucnt) = 1;
      A(r,3*ucnt+Dim) = 1;
      A(r,3*ucnt+1) = 1;
      A(r,3*ucnt+Dim+1) = 1;
      b(r) = P0[0]+P0[1]+P1[0]+P1[1];

      decltype(A) At = identity_matrix<double>(rows);
      if (rows!=cols) At = trans(A);

      // Solve
      decltype(A) Ainv(cols,cols);
      A = prod(At,A);
      if (invert(A, Ainv)) {
        u_prev = u;
      } else {
        // We need to relax more
        relax *= 0.75;
        u = u_prev;
        continue;
      }

      Ainv = prod(Ainv,At);
      auto x = prod(Ainv,b);

      for (size_t i=0; i<ucnt; i++) {
        u[i] = (x(i)-u[i])*relax + u[i];
        u[i] = std::max(0., std::min(u[i], 1.));
      }
      relax = std::min(init_relax, relax*1.01); // carefully regenerate the relaxation

      Pc[0] = { x(3*ucnt), x(3*ucnt+1) };
      Pc[1] = { x(3*ucnt+Dim), x(3*ucnt+Dim+1) };

      // Check the residuals
      std::vector<double> residuals;
      for (size_t i=0; i<ucnt; i++)
        residuals.push_back(std::abs(u[i]-x(i)));
      if (std::all_of(residuals.cbegin(), residuals.cend(), [&props](double r){ return r<props.tol; }))
        break;
    }

    if (it>props.max_it)
      std::cerr << "[NOT CONVERGED] ";
    else
      std::cerr << "[SUCCESS] ";
    std::cerr << "iterations: " << it << std::endl;

    return {P0, Pc[0], Pc[1], P1};
  }
}

#endif //NURBSFIT_CURVEFIT_HPP
