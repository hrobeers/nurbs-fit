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

#include <nurbs-fit/curvefit.hpp>

#include <list>

#include "nurbs-fit/curveproc.hpp"
#include "nurbs-fit/invert-matrix.hpp"

#include <boost/math/special_functions/pow.hpp>
namespace bm = boost::math;

// TODO move ublas includes here?

using namespace nurbsfit;

namespace nurbsfit {
  typedef std::function<void(std::vector<double>,
                             std::list<std::vector<double>>&,
                             std::vector<double>&)> f_Ab;

  // TODO templated Dim
  template<int ccnt>
  std::array<hrlib::vertex<2>, ccnt> solve2d(const f_Ab &fAb, const size_t ucnt,
                                             double relax, double tol, size_t max_it) {
    typedef hrlib::vertex<2> vertex;

    const size_t Dim = 2;
    const size_t pcnt = ucnt+2;

    const double init_relax = relax;

    // Result array
    std::array<vertex, ccnt> Pr;

    // First guess unknowns array
    std::vector<double> u(ucnt);
    for (size_t i=0; i<ucnt; i++)
      u[i] = (i+1.)/pcnt;
    std::vector<double> u_prev = u;

    // Iterative loop
    size_t it = 0;
    while (it++<max_it) {
      using namespace boost::numeric::ublas;

      std::list<std::vector<double>> vA;
      std::vector<double> vb;
      fAb(u, vA, vb);

      // Ax=b
      size_t rows = vA.size();
      size_t cols = vA.front().size();
      matrix<double> A = matrix<double>(rows, cols);
      vector<double> b = vector<double>(rows);

      // Fill the ublas matrices
      size_t ri=0;
      for (auto it=vA.cbegin(); it!=vA.cend(); ++it, ++ri) {
        b(ri) = vb[ri];
        for (size_t ci=0; ci<cols; ci++)
          A(ri,ci) = (*it)[ci];
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

      for (size_t i=0; i<ccnt; i++)
        Pr[i] = { x(3*ucnt+i*Dim), x(3*ucnt+i*Dim+1) };

      // Check the residuals
      std::vector<double> residuals;
      for (size_t i=0; i<ucnt; i++)
        residuals.push_back(std::abs(u[i]-x(i)));
      if (std::all_of(residuals.cbegin(), residuals.cend(), [&tol](double r){ return r<tol; }))
        break;
    }

    if (it>max_it)
      std::cerr << "[NOT CONVERGED] ";
    else
      std::cerr << "[SUCCESS] ";
    std::cerr << "iterations: " << it << std::endl;

    // TODO return std::array?
    return Pr;
  }
}

std::vector<hrlib::vertex<2>> nurbsfit::fit_qb(const std::vector<hrlib::vertex<2>> &P,
                                               double relax, double tol, size_t max_it) {
    typedef hrlib::vertex<2> vertex;

    const size_t Dim = 2;
    size_t pcnt = P.size();
    size_t ucnt = pcnt-2;
    const size_t ccnt = 1;

    assert(pcnt>=4);

    // Order of cols:
    // u, u^3, (1-u)^3, Pc1, Pc2
    size_t cols = 3*ucnt + ccnt*Dim;


    f_Ab fAb = [ucnt,cols,Dim,&P](auto u, auto &vA, auto &vb) {
      auto P0 = P.front();
      auto P1 = P.back();

      // linearization of t^2 around u
      // Taylor: t^2 ~ u^2 + 2u (t-u)
      // Matrix: 2u*t -1*t^2 = u^2

      // t^2
      for (size_t i=0; i<ucnt; i++) {
        std::vector<double> a(cols, 0);
        a[i] = 2*u[i];
        a[i+ucnt] = -1;
        vb.push_back(bm::pow<2>(u[i]));
        vA.push_back(std::move(a));
      }

      // linearization of (1-t)^2 around u
      // Taylor: (1-t)^2 ~ (1-u)^2 -2(1-u)(t-u)
      // Matrix: (2u-2)*t -1*(1-t)^2 = (u^2-1)

      // (1-t)^2
      for (size_t i=0; i<ucnt; i++) {
        std::vector<double> a(cols, 0);
        a[i] = 2*u[i]-2;
        a[i+2*ucnt] = -1;
        vb.push_back(bm::pow<2>(u[i]) - 1);
        vA.push_back(std::move(a));
      }

      // Pu equations
      for (size_t i=0; i<ucnt; i++)
        for (size_t d=0; d<Dim; d++) {
          std::vector<double> a(cols, 0);
          a[ucnt+i] = P1[d];
          a[2*ucnt+i] = P0[d];
          a[3*ucnt+d] = 2*(1-u[i])*u[i]; // Pc
          vb.push_back(P[i+1][d]);
          vA.push_back(std::move(a));
        }
    };

    std::vector<vertex> Pr;
    auto Pc = solve2d<1>(fAb, ucnt, relax, tol, max_it);
    Pr.push_back(P.front());
    Pr.insert(Pr.end(), Pc.cbegin(), Pc.cend());
    Pr.push_back(P.back());

    // TODO return std::array?
    return Pr;
}

std::vector<hrlib::vertex<2>> nurbsfit::fit_cb(const std::vector<hrlib::vertex<2>> &P,
                                               std::optional<std::array<double,2>> tangents,
                                               double relax, double tol, size_t max_it) {
    typedef hrlib::vertex<2> vertex;

    const size_t Dim = 2;
    size_t pcnt = P.size();
    size_t ucnt = pcnt-2;
    const size_t ccnt = 2;

    assert(pcnt>=5);

    // Order of cols:
    // u, u^3, (1-u)^3, Pc1, Pc2
    size_t cols = 3*ucnt + ccnt*Dim;


    f_Ab fAb = [ucnt,cols,Dim,tangents,&P](auto u, auto &vA, auto &vb) {
      auto P0 = P.front();
      auto P1 = P.back();

      // linearization of t^3 around u
      // Taylor: t^3 ~ u^3 + 3u^2 (t-u)
      // Matrix: 3u^2*t -1*t^3 = 2u^3

      // t^3
      for (size_t i=0; i<ucnt; i++) {
        std::vector<double> a(cols, 0);
        a[i] = 3*bm::pow<2>(u[i]);
        a[i+ucnt] = -1;
        vb.push_back(2*bm::pow<3>(u[i]));
        vA.push_back(std::move(a));
      }

      // linearization of (1-t)^3 around u
      // Taylor: (1-t)^3 ~ (1-u)^3 -3(1-u)^2(t-u)
      // Matrix: (3u^2-6u+3)*t +1*(1-t)^3 = 2u^3-3u^2+1

      // (1-t)^3
      for (size_t i=0; i<ucnt; i++) {
        std::vector<double> a(cols, 0);
        a[i] = 3*bm::pow<2>(u[i])-6*u[i]+3;
        a[i+2*ucnt] = 1;
        vb.push_back(2*bm::pow<3>(u[i]) - 3*bm::pow<2>(u[i]) + 1);
        vA.push_back(std::move(a));
      }

      // Pu equations
      for (size_t i=0; i<ucnt; i++)
        for (size_t d=0; d<Dim; d++) {
          std::vector<double> a(cols, 0);
          a[ucnt+i] = P1[d];
          a[2*ucnt+i] = P0[d];
          a[3*ucnt+d] = 3*bm::pow<2>(1-u[i])*u[i]; // Pc1
          a[3*ucnt+Dim+d] = 3*(1-u[i])*bm::pow<2>(u[i]); // Pc2
          vb.push_back(P[i+1][d]);
          vA.push_back(std::move(a));
        }

      if (ucnt < 8 && !tangents) {
        // Enforce Pc1x+Pc2x = P0x+P1x (extra equation)
        // Improves stability when X as main axis
        {
          std::vector<double> a(cols, 0);
          a[3*ucnt] = 1;
          a[3*ucnt+Dim] = 1;
          vb.push_back(P0[0]+P1[0]);
          vA.push_back(std::move(a));
        }
      }
      if (ucnt >= 8 || tangents) {
        // Fix tangents
        {
          std::vector<double> a(cols, 0);
          // f*dX = (1-f)*dY
          double f0 = (!tangents || std::isnan((*tangents)[0]))?
            (P[1][1]-P[0][1])/(P[1][0]-P[0][0]+P[1][1]-P[0][1]) :
            std::sin((*tangents)[0])/(std::sin((*tangents)[0])+std::cos((*tangents)[0]));
          a[3*ucnt] = f0;
          a[3*ucnt+1] = -(1-f0);
          vb.push_back(f0*P0[0]-(1-f0)*P0[1]);
          vA.push_back(std::move(a));
        }
        {
          std::vector<double> a(cols, 0);
          // f*dX = (1-f)*dY
          double f1 = (!tangents || std::isnan((*tangents)[1]))?
            (P[ucnt][1]-P[ucnt+1][1])/(P[ucnt][0]-P[ucnt+1][0]+P[ucnt][1]-P[ucnt+1][1]) :
            std::sin((*tangents)[1])/(std::sin((*tangents)[1])+std::cos((*tangents)[1]));
          a[3*ucnt+Dim] = f1;
          a[3*ucnt+Dim+1] = -(1-f1);
          vb.push_back(f1*P1[0]-(1-f1)*P1[1]);
          vA.push_back(std::move(a));
        }
      }
      /*
      // Sets the last point's tangent vertical
      A(r,3*ucnt+Dim) = 1;
      b(r) = 0;
      r++;
      // Sets 5*Pc1y = -Pc0x
      A(r,3*ucnt) = 1;
      A(r,3*ucnt+Dim+1) = 5;
      b(r) = P1[1]+P0[0];
      r++;
      */
    };

    std::vector<vertex> Pr;
    auto Pc = solve2d<2>(fAb, ucnt, relax, tol, max_it);
    Pr.push_back(P.front());
    Pr.insert(Pr.end(), Pc.cbegin(), Pc.cend());
    Pr.push_back(P.back());

    // TODO return std::array?
    return Pr;
}

/* Old formulation below
  // Different formulation of equations to solve.
  // Seems slightly less stable on the tested cases.
  template<typename Tprops>
  std::vector<hrlib::vertex<2>> fit_cb2(const std::vector<hrlib::vertex<2>> &P, const Tprops &props) {
    typedef hrlib::vertex<2> vertex;

    const size_t Dim = 2;
    size_t pcnt = P.size();
    size_t ucnt = pcnt-2;
    const size_t ccnt = 2;

    assert(pcnt>=5);

    // Order of rows:
    // P-equations, linearizations, extra equation
    size_t rows = Dim*ucnt + 3*ucnt + 1;
    // Order of cols:
    // u, u^3, u^2, Pc1, Pc2
    size_t cols = 3*ucnt + ccnt*Dim;

    // Hack to stop in time before exploding
    // TODO re-evaluate the initial relaxation
    const double init_relax = props.relax; // 64./props.max_it; // props.relax;
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

      // linearization of t^2 around u
      // Taylor: t^2 ~ u^2 + 2u (t-u)
      // Matrix: 2u*t -1*t^2 = u^2

      // t^2
      for (size_t i=0; i<ucnt; i++) {
        A(r,i) = 2*u[i];
        A(r,i+2*ucnt) = -1;
        b(r) = bm::pow<2>(u[i]);
        r++;
      }

      // u, u^3, u^2, Pc1, Pc2
      for (size_t i=0; i<ucnt; i++)
        for (size_t d=0; d<Dim; d++) {
          A(r,i) = -3*P0[d];         // u
          A(r,ucnt+i) = P1[d]-P0[d]; // u^3
          A(r,2*ucnt+i) = 3*P0[d];   // u^2
          A(r,3*ucnt+d) = 3*bm::pow<3>(u[i])-6*bm::pow<2>(u[i])+3*(u[i]); // Pc1
          A(r,3*ucnt+Dim+d) = -3*bm::pow<3>(u[i])+3*bm::pow<2>(u[i]); // Pc2
          b(r) = P[i+1][d]-P0[d];
          r++;
        }

      // Enforce Pc1x+Pc2x = P0x+P1x (extra equation)
      // Avoids a control point explosion
      A(r,3*ucnt) = 1;
      A(r,3*ucnt+Dim) = 1;
      b(r) = P0[0]+P1[0];

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
*/
