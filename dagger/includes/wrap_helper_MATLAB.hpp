//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// WIP
// MATLAB™®©
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef wrap_helper_MATLAB_HPP
#define wrap_helper_MATLAB_HPP

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

// MATLAB™®© (proprietay, license protected, not free) related includes
#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

// MATLAB™®© (proprietay, license protected, not free) engine stuff
using namespace matlab::engine;

namespace DAGGER {

template <class T> class matlvec {
public:
  // Create MATLAB ™®© (proprietay, license protected, not free) data array
  // factory
  matlab::data::ArrayFactory factory;

  int isize = 0;
  size_t usize = 0;
  std::vector<T> arr;

  matlvec(){};
  matlvec(matlab::data::TypedArray<T> &arr) {
    this->arr = std::vector<T>(arr.begin(), arr.end());
    this->isize = arr.getNumberOfElements();
    this->usize = arr.getNumberOfElements();
  };

  matlvec(matlvec<T> &tin) {
    this->arr = tin.arr;
    this->isize = tin.isize;
    this->usize = tin.usize;
  }

  matlvec(const matlvec<T> &tin) {
    this->arr = tin.arr;
    this->isize = tin.isize;
    this->usize = tin.usize;
  }

  T &operator[](int i) { return this->arr[i]; }

  void set(int i, T val) { this->arr[i] = val; }

  T get(int i) { return this->arr[i]; }

  std::vector<T> to_vec() {
    std::vector<T> out(this->isize);
    for (int i = 0; i < this->isize; ++i)
      out[i] = this->arr[i];
  }

  // size_t size(){return this->usize;}
  const size_t size() { return this->usize; }
};

template <class A, class T> matlab::data::TypedArray<T> vec2Marray(A &vec) {
  matlab::data::ArrayFactory factory;
  auto dim = vec.size();
  matlab::data::TypedArray<T> retr =
      factory.createArray({1, dim}, vec.begin(), vec.end());
  return retr;
}

template <class A, class T> matlab::data::TypedArray<T> pvec2Marray(A &vec) {
  matlab::data::ArrayFactory factory;
  auto dim = vec.data->size();
  return factory.createArray({1, dim}, vec.data->begin(), vec.data->end());
}

template <class T>
matlab::data::TypedArray<T> _format_output(std::vector<T> &yolo) {
  return vec2Marray<std::vector<T>, T>(yolo);
}
template <class T>
matlab::data::TypedArray<T> _format_output(pvector<T> &yolo) {
  return pvec2Marray<pvector<T>, T>(yolo);
}
template <class T> std::vector<T> _format_output(pvector<T> &yolo) {
  return yolo.to_vec();
}

// template<class T>
// std::vector<T> _format_output(std::vector<T>& yolo){return yolo;}

template <class T>
matlab::data::TypedArray<T> _format_output(matlvec<T> &yolo) {
  auto vec = yolo.to_vec();
  return vec2Marray(vec);
}

template <class T> matlvec<T> _format_output(matlvec<T> &yolo) { return yolo; }

template <class T>
matlab::data::TypedArray<T> _format_output(matlab::data::TypedArray<T> &yolo) {
  return yolo;
}

// template<class T>
// std::vector<T> to_vec(matlvec<T>& in)
// {
// 	std::vector<T> out(in.arr);
// 	return out;
// }

template <class T> std::vector<T> to_vec(matlvec<T> in) {
  std::vector<T> out(in.arr);
  return out;
}

template <class T> std::vector<T> to_vec(matlab::data::TypedArray<T> &tin) {
  std::vector<T> out(tin.begin(), tin.end());
  return out;
}

template <class A> std::vector<double> to_vec(A &tin) {
  matlvec<double> in(tin);
  std::vector<double> out(in.size());
  for (size_t i = 0; i < in.size(); ++i)
    out[i] = in[i];
  return out;
}

template <class T> matlvec<T> _format_input(matlab::data::TypedArray<T> &yolo) {
  matlvec<T> gabul(yolo);
  return gabul;
}

template <class T> matlvec<T> _format_input(matlvec<T> &yolo) { return yolo; }

template <class A> A _format_input(A &yolo) { return to_vec(yolo); }

// end of namesapce DAGGER
}; // namespace DAGGER

#endif
