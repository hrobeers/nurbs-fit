#pragma once

#ifndef NURBSFIT_CURVEFIT_HPP
#define NURBSFIT_CURVEFIT_HPP

#include <cassert>

#include "hrlib/io/vertexio.hpp"
#include "invert-matrix.hpp"

namespace nurbsfit
{
  template<typename Tprops>
  std::vector<hrlib::vertex<2>> fit_qb(const std::vector<hrlib::vertex<2>> &target, const Tprops &props) {
    assert(target.size()==4); // TODO warn or switch to least squares?

    double relax = props.relax;

    auto P0 = target[0];
    auto Pu = target[1];
    auto Pv = target[2];
    auto P1 = target[3];

    decltype(P0) Pc = {.0,.0};

    double u = 0.25;
    double v = 0.75;
    double u_prev = u;
    double v_prev = v;

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

      // linearization of t^2 around a
      // Taylor: t^2 ~ a^2 + 2a (t-a)
      // Matrix: 2a*t -1*t^2 = a^2

      // u^2
      A(r,0) = 2*u;
      A(r,2) = -1;
      b(r) = u*u;
      r++;
      // v^2
      A(r,1) = 2*v;
      A(r,3) = -1;
      b(r) = v*v;
      r++;

      // linearization of (1-t)^2 around a
      // Taylor: (1-t)^2 ~ (1-a)^2 -2(1-a)(t-a)
      // Matrix: (2a-2)*t -1*(1-t)^2 = (a^2-1)

      // (1-u)^2
      A(r,0) = (2*u-2);
      A(r,4) = -1;
      b(r) = u*u - 1;
      r++;
      // (1-v)^2
      A(r,1) = (2*v-2);
      A(r,5) = -1;
      b(r) = v*v - 1;
      r++;

      // Pu equations
      A(r,2) = P1[0];
      A(r,4) = P0[0];
      A(r,6) = 2*u*(1-u);
      b(r) = Pu[0];
      r++;
      A(r,2) = P1[1];
      A(r,4) = P0[1];
      A(r,7) = 2*u*(1-u);
      b(r) = Pu[1];
      r++;

      // Pv equations
      A(r,3) = P1[0];
      A(r,5) = P0[0];
      A(r,6) = 2*v*(1-v);
      b(r) = Pv[0];
      r++;
      A(r,3) = P1[1];
      A(r,5) = P0[1];
      A(r,7) = 2*v*(1-v);
      b(r) = Pv[1];
      r++;

      // Solve
      decltype(A) Ainv(8,8);
      if (invert(A, Ainv)) {
        u_prev = u;
        v_prev = v;
      } else {
        // We need to relax more
        relax *= 0.75;
        u = u_prev;
        v = v_prev;
        continue;
      }

      auto x = prod(Ainv,b);

      u = (x(0)-u)*relax + u;
      v = (x(1)-v)*relax + v;
      u = std::max(0., std::min(u, 1.));
      v = std::max(0., std::min(v, 1.));
      relax = std::min(props.relax, relax*1.01); // carefully regenerate the relaxation

      Pc[0] = x(6);
      Pc[1] = x(7);

      // Check the residuals
      /*
      std::cerr << "res u:  " << std::abs(u-x(0)) << std::endl;
      std::cerr << "res v:  " << std::abs(v-x(1)) << std::endl;
      std::cerr << "res u2: " << std::abs(u*u-x(2)) << std::endl;
      std::cerr << "res v2: " << std::abs(v*v-x(3)) << std::endl << std::endl;
      */
      // u, v, u^2, v^2, (1-u)^2, (1-v)^2, Pcx, Pcy
      if (std::abs(u-x(0))<props.tol &&
          std::abs(v-x(1))<props.tol &&
          std::abs(u*u-x(2))<(props.tol*props.tol) &&
          std::abs(v*v-x(3))<(props.tol*props.tol)) { // &&
          /*
          (1-u)*(1-u)-x(4)<props.tol &&
          (1-v)*(1-v)-x(5)<props.tol) {
          */
        std::cerr << "iterations: " << it << std::endl;
        break;
      }
    }

    return {P0, Pc, P1};
  }
}

#endif //NURBSFIT_CURVEFIT_HPP
