#ifndef WRAP_HELPER_CPP_HPP
#define WRAP_HELPER_CPP_HPP

namespace DAGGER {

template <class T> std::vector<T> &_format_output(std::vector<T> &in) {
  return in;
}

// template<class T>
// std::vector<T>& _format_input(std::vector<T>& in){return in;}

// template<class T>
// std::vector<T>& to_vec(std::vector<T>& in){return in;}

// end of namespace DAGGER
}; // namespace DAGGER

#endif
