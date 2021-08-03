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

#include <cstdlib>
#include <vector>
#include <list>

#include "../version_autogen.hpp"

#include "nurbs-fit/main.hpp"
#include "nurbs-fit/curvefit.hpp"
#include "nurbs-fit/curveproc.hpp"
#include "nurbs-fit/props.hpp"
#include "nurbs-fit/colormap.hpp"

#include "hrlib/io/vertexio.hpp"

using namespace nurbsfit;
using namespace hrlib;

inline void showHelp(const po::options_description &cmdline_options, std::ostream &stream = std::cout)
{
    stream << cmdline_options << "\n";
}

int main(int argc, char *argv[])
{
  try
  {
    //
    // Positional options
    //
    po::positional_options_description positional_options;
    positional_options.add("input-file", -1);

    po::options_description positional_params;
    positional_params.add_options()
      ("input-file", po::value<std::string>(), "Optional input file (can also be provided through stdin)");


    //
    // Generic options
    //
    po::options_description generic_params("Generic options");
    generic_params.add_options()
      ("version,v", "Print version string")
      ("about,a", "About this application")
      ("help,h", "Print help");

    //
    // Curve-fit options
    //
    po::options_description curve_params("Curve options");
    curve_params.add_options()
      //("type,t", "Type of curve to fit. [qb, cb]")
      //("dimensions,d", "Nurber of dimensions. Currently only 2 is supported (option is ignored)")
      ("origin,o", "Shifts and rotates the curve to have the middle at the origin with an horizontal tangent")
      ;

    //
    // Process the actual command line arguments given by the user
    //
    po::options_description cmdline_options;
    cmdline_options.add(positional_params).add(generic_params).add(curve_params);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(positional_options).run(), vm);
    po::notify(vm);

    // Process the generic options
    bool do_exit = false;
    if (vm.count("version")) {
      std::cout << MAJOR_VERSION << '.'
                << MINOR_VERSION << '.'
                << REVISION << '.'
                << BUILD_NUMBER << '\n'
                << COMMIT_HASH << '\n';
      do_exit = true;
    }
    if (vm.count("help")) {
      showHelp(cmdline_options);
      do_exit = true;
    }
    if (do_exit) exit(EXIT_SUCCESS);


    //
    // Process cmd params
    //
    props p;
    p.origin = vm.count("origin");


    //
    // Run the application
    //

    if (true)
    {
      auto input = getInputStream(vm);
      std::istream* is = input.get();

      std::cout << "<svg>" << std::endl;

      std::list<std::vector<vertex<2>>> inputs;
      vertex<2> vtx;
      while (utf8::read_next_vertex<2>(*is, vtx)) {
        if (vtx[0] == 0)
          inputs.push_back({});
        inputs.back().push_back(vtx);
      }

      size_t i=0;
      for (auto &vts : inputs) {
        if (vts.size() != 4)
          continue;

        auto fit = fit_qb(vts,p);

        if (p.origin)
          // ceter on origin
          fit = center_origin(fit);

        // Convert to cubic for more compatibility
        fit = to_cubic(fit);

        std::cout << "<g stroke-width=\".3\">" << std::endl;

        if (!p.origin) {
          std::cout << "<g fill=\"blue\" stroke=\"none\">" << std::endl;
          for (auto v : vts)
            std::cout << "<circle cx=\"" << v[0] << "\" cy=\"" << -v[1] << "\" r=\".5\" />" << std::endl;
          std::cout << "</g>" << std::endl;

          std::cout << "<g fill=\"red\" stroke=\"none\">" << std::endl;
          std::cout << "<circle cx=\"" << fit[1][0] << "\" cy=\"" << -fit[1][1] << "\" r=\".5\" />" << std::endl;
          std::cout << "</g>" << std::endl;
        }

        std::cout << "<g  fill=\"none\" stroke=\"" << to_html(get_color(double(i++)/(inputs.size()-1))) << "\">" << std::endl;
        std::cout << "<path d=\"";
        std::cout << "M" << fit[0][0] << " " << -fit[0][1];
        std::cout << "C" << fit[1][0] << " " << -fit[1][1];
        std::cout << " " << fit[2][0] << " " << -fit[2][1];
        std::cout << " " << fit[3][0] << " " << -fit[3][1];
        std::cout << "\"/>" << std::endl;
        std::cout << "</g>" << std::endl;

        std::cout << "</g>" << std::endl;
      }

      std::cout << "</svg>" << std::endl;
    }
    else
    {
      showHelp(cmdline_options, std::cerr);
    }

    exit(EXIT_SUCCESS);
  }
  catch (std::exception &ex)
  {
    std::cerr << ex.what() << std::endl;
    exit(EXIT_FAILURE);
  }
}
