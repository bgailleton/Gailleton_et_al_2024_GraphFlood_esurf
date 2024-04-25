#ifndef FASTFLOOD_RECORDER_HPP
#define FASTFLOOD_RECORDER_HPP

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

#ifdef OPENMP_YOLO
#include <omp.h>
#endif

// local includes
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> The connector classes
#include "D8connector.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

namespace DAGGER {

template <class float_t> class fastflood_recorder {
public:
  fastflood_recorder() { ; }
  fastflood_recorder(int nnodes, int nlinks) {
    this->nnodes = nnodes;
    this->nlinks = nlinks;
  }

  void init() {
    if (this->edot2record)
      this->edot2record_init();
    if (this->ddot2record)
      this->ddot2record_init();
    if (this->lateral_edot2record)
      this->lateral_edot2record_init();
    if (this->lateral_ddot2record)
      this->lateral_ddot2record_init();
    if (this->qs2record)
      this->qs2record_init();
    if (this->dhw2record)
      this->dhw2record_init();
    if (this->vmot2record)
      this->vmot2record_init();
    if (this->tau2record)
      this->tau2record_init();
  }

  void init_water() {
    if (this->dhw2record)
      this->dhw2record_init();
  }

  void init_geo() {
    if (this->edot2record)
      this->edot2record_init();
    if (this->ddot2record)
      this->ddot2record_init();
    if (this->lateral_edot2record)
      this->lateral_edot2record_init();
    if (this->lateral_ddot2record)
      this->lateral_ddot2record_init();
    if (this->qs2record)
      this->qs2record_init();
    if (this->tau2record)
      this->tau2record_init();
    if (this->vmot2record)
      this->vmot2record_init();
  }

  // nnodes for the recording
  int nnodes = 0;
  int nlinks = 0;

  // stuff to record, eventually
  bool edot2record = false;
  std::vector<float_t> edot;
  void enable_edot_recording() {
    this->edot2record = true;
    this->edot = std::vector<float_t>(this->nnodes, 0.);
  };
  void disable_edot_recording() {
    this->edot2record = false;
    this->edot = std::vector<float_t>();
  }
  template <class out_t> out_t get_edot() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->edot);
  }
  void edot2record_init() {
    this->edot = std::vector<float_t>(this->nnodes, 0.);
  }

  bool ddot2record = false;
  std::vector<float_t> ddot;
  void enable_ddot_recording() {
    this->ddot2record = true;
    this->ddot = std::vector<float_t>(this->nnodes, 0.);
  };
  void disable_ddot_recording() {
    this->ddot2record = false;
    this->ddot = std::vector<float_t>();
  }
  template <class out_t> out_t get_ddot() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->ddot);
  }
  void ddot2record_init() {
    this->ddot = std::vector<float_t>(this->nnodes, 0.);
  }

  bool lateral_edot2record = false;
  std::vector<float_t> lateral_edot;
  void enable_lateral_edot_recording() {
    this->lateral_edot2record = true;
    this->lateral_edot = std::vector<float_t>(this->nnodes, 0.);
  };
  void disable_lateral_edot_recording() {
    this->lateral_edot2record = false;
    this->lateral_edot = std::vector<float_t>();
  }
  template <class out_t> out_t get_lateral_edot() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(
        this->lateral_edot);
  }
  void lateral_edot2record_init() {
    this->lateral_edot = std::vector<float_t>(this->nnodes, 0.);
  }

  bool lateral_ddot2record = false;
  std::vector<float_t> lateral_ddot;
  void enable_lateral_ddot_recording() {
    this->lateral_ddot2record = true;
    this->lateral_ddot = std::vector<float_t>(this->nnodes, 0.);
  };
  void disable_lateral_ddot_recording() {
    this->lateral_ddot2record = false;
    this->lateral_ddot = std::vector<float_t>();
  }
  template <class out_t> out_t get_lateral_ddot() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(
        this->lateral_ddot);
  }
  void lateral_ddot2record_init() {
    this->lateral_ddot = std::vector<float_t>(this->nnodes, 0.);
  }

  bool qs2record = false;
  std::vector<float_t> qs;
  void enable_qs_recording() {
    this->qs2record = true;
    this->qs = std::vector<float_t>(this->nnodes, 0.);
  };
  void disable_qs_recording() {
    this->qs2record = false;
    this->qs = std::vector<float_t>();
  }
  template <class out_t> out_t get_qs() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->qs);
  }
  void qs2record_init() { this->qs = std::vector<float_t>(this->nnodes, 0.); }

  bool dhw2record = false;
  std::vector<float_t> dhw;
  void enable_dhw_recording() {
    this->dhw2record = true;
    this->dhw = std::vector<float_t>(this->nnodes, 0.);
  };
  void disable_dhw_recording() {
    this->dhw2record = false;
    this->dhw = std::vector<float_t>();
  }
  template <class out_t> out_t get_dhw() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->dhw);
  }
  void dhw2record_init() { this->dhw = std::vector<float_t>(this->nnodes, 0.); }

  bool tau2record = false;
  std::vector<float_t> tau;
  void enable_tau_recording() {
    this->tau2record = true;
    this->tau = std::vector<float_t>(this->nnodes, 0.);
  };
  void disable_tau_recording() {
    this->tau2record = false;
    this->tau = std::vector<float_t>();
  }
  template <class out_t> out_t get_tau() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->tau);
  }
  void tau2record_init() { this->tau = std::vector<float_t>(this->nnodes, 0.); }

  bool vmot2record = false;
  std::vector<float_t> vmot;
  void enable_vmot_recording() {
    this->vmot2record = true;
    this->vmot = std::vector<float_t>(this->nnodes, 0.);
  };
  void disable_vmot_recording() {
    this->vmot2record = false;
    this->vmot = std::vector<float_t>();
  }
  template <class out_t> out_t get_vmot() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->vmot);
  }
  void vmot2record_init() {
    this->vmot = std::vector<float_t>(this->nnodes, 0.);
  }
};

} // namespace DAGGER

#endif
