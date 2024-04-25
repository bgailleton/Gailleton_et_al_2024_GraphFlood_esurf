#ifndef wrap_helper_HPP
#define wrap_helper_HPP

// STL imports
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <stack>
#include <stdlib.h>
#include <string>
#include <thread>
#include <vector>

#include "utils.hpp"

// defines all the format_input depnding on the eventual wrapper
#ifdef DAGGER_FT_PYTHON
#include "wrap_helper_python.hpp"
#elif DAGGER_FT_JULIA
#include "wrap_helper_julia.hpp"
#elif DAGGER_FT_MATLAB
#include "wrap_helper_MATLAB.hpp" // DOES NOT WORK AT THE MOMENT
// #include "wrap_helper_cpp.hpp"
#else
#include "wrap_helper_cpp.hpp"
#endif

namespace DAGGER {

template <typename in_t> auto format_input(in_t &tin) {
  auto ret = _format_input(tin);
  return ret;
}

template <class in_t, class out_t> out_t format_output(in_t &tin) {
  out_t ret = _format_output(tin);
  return ret;
}

} // namespace DAGGER

#endif
