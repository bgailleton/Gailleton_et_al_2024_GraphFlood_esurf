//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef wrap_helper_JULIA_HPP
#define wrap_helper_JULIA_HPP

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

#include "jlcxx/functions.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

namespace DAGGER {

template <class T> class julvec {
public:
  int isize = 0;
  size_t usize = 0;
  jlcxx::ArrayRef<T> *arr;

  julvec(){};
  julvec(jlcxx::ArrayRef<T> &arr) {
    this->arr = &arr;
    this->isize = arr.size();
    this->usize = arr.size();
  };

  T &operator[](int i) { return (*this->arr)[i]; }

  void set(int i, T val) { (*this->arr)[i] = val; }

  T get(int i) { return (*this->arr)[i]; }

  std::vector<T> to_vec() {
    std::vector<T> out(this->isize);
    for (int i = 0; i < this->isize; ++i)
      out[i] = (*this->arr)[i];
  }

  // size_t size(){return this->usize;}
  const size_t size() { return this->usize; }
};

template <class T> jlcxx::ArrayRef<T> _format_output(std::vector<T> &yolo) {
  return jlcxx::ArrayRef<T>(yolo.data());
}
template <class T> jlcxx::ArrayRef<T> _format_output(pvector<T> &yolo) {
  return jlcxx::ArrayRef<T>(yolo.data->data());
}
template <class T> std::vector<T> _format_output(pvector<T> &yolo) {
  return yolo.to_vec();
}

// template<class T>
// std::vector<T> _format_output(std::vector<T>& yolo){return yolo;}

template <class T> jlcxx::ArrayRef<T> _format_output(julvec<T> &yolo) {
  auto vec = yolo.to_vec();
  return jlcxx::ArrayRef<T>(vec.data());
}
template <class T> jlcxx::ArrayRef<T> _format_output(jlcxx::ArrayRef<T> &yolo) {
  return yolo;
}

template <class T> std::vector<T> to_vec(julvec<T> &in) {
  std::vector<T> out(in.size());
  for (size_t i = 0; i < in.size(); ++i)
    out[i] = in[i];
  return out;
}

// template<class T>
// std::vector<T> to_vec(julvec<T> in)
// {
// 	std::vector<T> out(in.size());
// 	for(size_t i=0;i<in.size(); ++i)
// 		out[i] = in[i];
// 	return out;
// }

template <class T> std::vector<T> to_vec(jlcxx::ArrayRef<T> &tin) {
  julvec<double> in(tin);
  std::vector<T> out(in.size());
  for (size_t i = 0; i < in.size(); ++i)
    out[i] = in[i];
  return out;
}

template <class T> julvec<T> _format_input(jlcxx::ArrayRef<T> &yolo) {
  julvec<T> gabul(yolo);
  return gabul;
}

template <class T> julvec<T> _format_input(julvec<T> &yolo) { return yolo; }

// end of namesapce DAGGER
}; // namespace DAGGER

#endif
