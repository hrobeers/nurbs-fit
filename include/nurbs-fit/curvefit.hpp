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
    const size_t pcnt = 4;
    const size_t ucnt = pcnt-2;
    const size_t ccnt = ucnt/2;

    assert(P.size()==pcnt); // TODO warn or switch to least squares?

    double relax = props.relax;

    // Bezier points
    vertex P0 = P.front();
    vertex P1 = P.back();
    std::array<vertex, ccnt> Pc;

    std::array<double,ucnt> u = {.25,.75};
    std::array<double,ucnt> u_prev = u;

    size_t it = 0;
    while (it++<props.max_it) {
      using namespace boost::numeric::ublas;
      // Quadratic Bezier: (1-t)^2 P0 + 2t(1-t) Pc + u^2 P1 = P(t)

      // Ax=b
      // Order of x variables:
      // u, v, u^2, v^2, (1-u)^2, (1-v)^2, Pcx, Pcy
      matrix<double> A = zero_matrix<double>(8, 8);
      vector<double> b = zero_vector<double>(8);

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

      // Solve
      decltype(A) Ainv(8,8);
      if (invert(A, Ainv)) {
        u_prev = u;
      } else {
        // We need to relax more
        relax *= 0.75;
        u = u_prev;
        continue;
      }

      auto x = prod(Ainv,b);

      for (size_t i=0; i<ucnt; i++) {
        u[i] = (x(i)-u[i])*relax + u[i];
        u[i] = std::max(0., std::min(u[i], 1.));
      }
      relax = std::min(props.relax, relax*1.01); // carefully regenerate the relaxation

      Pc[0] = { x(3*ucnt), x(3*ucnt+1) };

      // Check the residuals
      // u, u^2, (1-u)^2, Pc
      if (std::abs(u[0]-x(0))<props.tol &&
          std::abs(u[1]-x(1))<props.tol &&
          std::abs(bm::pow<2>(u[0])-x(2))<(props.tol*props.tol) &&
          std::abs(bm::pow<2>(u[1])-x(3))<(props.tol*props.tol)) { // &&
          /*
          (1-u)*(1-u)-x(4)<props.tol &&
          (1-v)*(1-v)-x(5)<props.tol) {
          */
        break;
      }
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
    const size_t pcnt = 6;
    const size_t ucnt = pcnt-2;
    const size_t ccnt = ucnt/2;

    assert(P.size()==pcnt); // TODO warn or switch to least squares?

    // Hack to stop in time before exploding
    // TODO re-evaluate the stop condition
    const double init_relax = 1./props.max_it; // props.relax;
    double relax = init_relax;

    // Bezier points
    auto P0 = P.front();
    auto P1 = P.back();
    std::array<vertex, ccnt> Pc;

    std::array<double,ucnt> u = {.2,.4,.6,.8};
    std::array<double,ucnt> u_prev = u;

    size_t it = 0;
    while (it++<props.max_it) {
      using namespace boost::numeric::ublas;
      // Quadratic Bezier: (1-t)^2 P0 + 2t(1-t) Pc + u^2 P1 = P(t)

      // Ax=b
      // Order of x variables:
      // u, u^3, (1-u)^3, Pc1, Pc2 -> 4+4+4+2+2 = 16
      matrix<double> A = zero_matrix<double>(16, 16);
      vector<double> b = zero_vector<double>(16);

      size_t r = 0;
      size_t c = 0;

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

      // Solve
      decltype(A) Ainv(16,16);
      if (invert(A, Ainv)) {
        u_prev = u;
      } else {
        // We need to relax more
        relax *= 0.75;
        u = u_prev;
        continue;
      }

      auto x = prod(Ainv,b);

      for (size_t i=0; i<ucnt; i++) {
        u[i] = (x(i)-u[i])*relax + u[i];
        u[i] = std::max(0., std::min(u[i], 1.));
      }
      relax = std::min(init_relax, relax*1.01); // carefully regenerate the relaxation

      Pc[0] = { x(3*ucnt), x(3*ucnt+1) };
      Pc[1] = { x(3*ucnt+Dim), x(3*ucnt+Dim+1) };

      // Check the residuals
      // u, u^3, (1-u)^3, Pc1, Pc2
      size_t in_tol = 0;
      for (size_t i; i<ucnt; i++)
        if (std::abs(u[i]-x(i))<props.tol)
          in_tol++;
      if (in_tol>=ucnt-1)
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
