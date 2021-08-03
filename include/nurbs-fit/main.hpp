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

#ifndef NURBSFIT_MAIN_HPP
#define NURBSFIT_MAIN_HPP

#include <memory>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace nurbsfit
{
  class input {
    std::unique_ptr<std::istream> file;
    std::istream* str = nullptr;
  public:
    input(std::istream* s) : str(s) {}
    input(std::string file_name) { file.reset(new std::ifstream(file_name)); }

    std::istream* get() {
      if(file)
        return file.get();
      return str;
    }
  };

  inline
  input getInputStream(const po::variables_map &vm)
  {
    if (vm.count("input-file")) {
      auto filePath = vm["input-file"].as<std::string>();
      return std::move(input(filePath));
    }
    else {
      return std::move(input(&std::cin));
    }
  }
}

#endif //NURBSFIT_MAIN_HPP
