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
