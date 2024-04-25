//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef wrap_helper_python_HPP
#define wrap_helper_python_HPP

#ifndef DAGGER_FT_PYTHON
#define DAGGER_FT_PYTHON
#endif

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

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace DAGGER {

template <class T> class numvec {
public:
  T *ptr = nullptr;
  int isize = 0;
  size_t usize = 0;

  numvec(){};
  numvec(py::array_t<T, 1> &arr) {
    auto boeuf_heure = arr.request();
    this->ptr = (T *)boeuf_heure.ptr;
    this->isize = arr.size();
    this->usize = arr.size();
  };

  T &operator[](int i) { return this->ptr[i]; }

  void set(int i, T val) { (*this)[i] = val; }
  T get(int i) { return (*this)[i]; }

  std::vector<T> to_vec() {
    std::vector<T> out(this->isize);
    for (int i = 0; i < this->isize; ++i)
      out[i] = (*this)[i];
  }

  // size_t size(){return this->usize;}
  const size_t size() { return this->usize; }
};

template <class T> py::array _format_output(std::vector<T> &yolo) {
  return py::array(yolo.size(), yolo.data());
}
template <class T> py::array _format_output(pvector<T> &yolo) {
  return py::array(yolo.data->size(), yolo.data->data());
}
template <class T> std::vector<T> _format_output(pvector<T> &yolo) {
  return yolo.to_vec();
}

// template<class T>
// std::vector<T> _format_output(std::vector<T>& yolo){return yolo;}

template <class T> py::array _format_output(numvec<T> &yolo) {
  auto vec = yolo.to_vec();
  return py::array(vec.size(), vec.data());
}
// template<class T>
py::array _format_output(py::array &yolo) { return yolo; }

template <class T> std::vector<T> to_vec(numvec<T> &in) {
  std::vector<T> out(in.size());
  for (size_t i = 0; i < in.size(); ++i)
    out[i] = in[i];
  return out;
}

// template<class T>
// std::vector<T> to_vec(numvec<T> in)
// {
// 	std::vector<T> out(in.size());
// 	for(size_t i=0;i<in.size(); ++i)
// 		out[i] = in[i];
// 	return out;
// }

template <class T> std::vector<T> to_vec(py::array_t<T, 1> &tin) {
  numvec<double> in(tin);
  std::vector<T> out(in.size());
  for (size_t i = 0; i < in.size(); ++i)
    out[i] = in[i];
  return out;
}

template <class T> numvec<T> _format_input(py::array_t<T, 1> &yolo) {
  numvec<T> gabul(yolo);
  return gabul;
}

template <class T> numvec<T> _format_input(numvec<T> &yolo) { return yolo; }

// end of namesapce DAGGER
}; // namespace DAGGER

#endif
