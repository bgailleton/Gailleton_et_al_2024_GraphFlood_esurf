//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef D4connector_HPP
#define D4connector_HPP

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

// local includes
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

namespace DAGGER {

// Useful namespace for my priority queue
using PQ_i_d = std::priority_queue<PQ_helper<int, double>,
                                   std::vector<PQ_helper<int, double>>,
                                   std::greater<PQ_helper<int, double>>>;

template <class T> class D4connector {
public:
  // General informations about the graph
  // #-> number of nodes (integer and unsized int for loops or list innit).
  int nnodes = 0;
  size_t nnodes_t = 0;
  // #-> number of nodes in x direction.
  int nx = 0;
  // #-> number of nodes in x direction.
  int ny = 0;
  // #-> length in the x direction.
  T dx;
  // #-> length in the y direction.
  T dy;
  T dxy;
  T dxmax;
  T dxmin;
  // #-> cell area
  T cellarea;
  T Xmin;
  T Ymin;
  T Xmax;
  T Ymax;
  int not_a_node;

  // Bearer of values
  // #->boundary: index: node index, value: boundary value
  // #-->The possible values are:
  // #---> -1 = no in, no out, not process, internal or external
  // #--->  0 = in, no out, all fluxes leave the system, internal or external
  // #--->  1 = "normal" cell, in, out, internal
  // #--->  2 = external periodic boundary
  std::vector<int> boundary;

  // Helpers for neighbouring operations
  // -> neighbourer holdsthe indices to loop through for each boundary condition
  std::vector<std::vector<int>> neighbourer;
  // -> lengthener is the dx on each directions
  std::vector<T> lengthener;

  // Coordinate stuff
  // Xs and Ys are vectors of nx and ny size converting row to y and col to X
  // Extents holds the cxmin,xmax,ymin,ymax (extent is an option in matplotlib
  // imshow plots)
  std::vector<T> Xs, Ys, extents;

  std::shared_ptr<easyRand> randu = std::make_shared<easyRand>();

  // Default constructor, empty
  D4connector(){};

  // construct directly from dimensions
  D4connector(int nx, int ny, T dx, T dy, T xmin, T ymin) {
    // initialisation is offset to dedicated function because some languages
    // like Julia are complicating non default constructor initialisation.
    this->init_dimensions(nx, ny, dx, dy, xmin, ymin);
  }

  // initialise with dimension
  void init_dimensions(int nx, int ny, T dx, T dy, T xmin, T ymin) {
    int nnodes = nx * ny;
    this->set_dimensions(nx, ny, nnodes, dx, dy, xmin, ymin);
    this->set_default_boundaries("4edges");
  }

  /// @Description: Initialising the "neighbourer", a data structure managing
  /// the iterations though the neighbours of a given node Function of boundary
  /// conditions, one can feed the neighbourer with an indice which returns the
  /// index to ADD to the node to get its neighbour. The neighbour order is (in
  /// table referential, be careful if Y axis is inverted) top-left, top,
  /// top-right, left, right, bottom-left, bottom, bottom-right
  /// @Authors: B.G.
  /// @Date: 2021
  void initialise_neighbourer() {

    this->dxmin = std::min(this->dx, this->dy);
    this->dxmax = std::max(this->dx, this->dy);

    this->lengthener = std::initializer_list<T>{dy, dx, dx, dy};
    this->neighbourer.clear();

    // these vectors are additioned to the node indice to test the neighbors
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, -1, 1, this->nx}); // internal node 0
    this->neighbourer.emplace_back(std::initializer_list<int>{
        (this->ny - 1) * this->nx, -1, 1, this->nx}); // periodic_first_row 1
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, -1, 1, -(this->ny - 1) * this->nx}); // periodic_last_row 2
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, (this->nx - 1), 1, this->nx}); // periodic_first_col 3
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, -1, -this->nx + 1, this->nx}); // periodic last_col 4
    this->neighbourer.emplace_back(std::initializer_list<int>{
        this->not_a_node, -1, 1, this->nx}); // normal_first_row 5
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, -1, 1, this->not_a_node}); // normal_last_row 6
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, this->not_a_node, 1, this->nx}); // normal_first_col 7
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, -1, this->not_a_node, this->nx}); // normal_last_col 8
    this->neighbourer.emplace_back(std::initializer_list<int>{
        this->not_a_node, this->not_a_node, 1, this->nx}); // normal_top_left 9
    this->neighbourer.emplace_back(
        std::initializer_list<int>{this->not_a_node, -1, this->not_a_node,
                                   this->nx}); // normal_top_right 10
    this->neighbourer.emplace_back(
        std::initializer_list<int>{-this->nx, this->not_a_node, 1,
                                   this->not_a_node}); // normal_bottom_left 11
    this->neighbourer.emplace_back(
        std::initializer_list<int>{-this->nx, -1, this->not_a_node,
                                   this->not_a_node}); // normal_bottom_right 12
    this->neighbourer.emplace_back(
        std::initializer_list<int>{(this->ny - 1) * this->nx, this->nx - 1, 1,
                                   this->nx}); // top_left_periodic 13
    this->neighbourer.emplace_back(
        std::initializer_list<int>{(this->ny - 1) * this->nx, -1, -this->nx + 1,
                                   this->nx}); // top_right_periodic 14
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, this->nx - 1, 1,
        -(this->ny - 1) * this->nx}); // periodic_bottom_left 15
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx, -1, 1 - this->nx + 1,
        -(this->ny - 1) * this->nx}); // periodic_bottom_right 16
  }

  /// @description: Sets the dimension of the graph and initialise the
  /// neighbourer
  void set_dimensions(int nx, int ny, int nnodes, T dx, T dy, T xmin, T ymin) {
    this->nx = nx;
    this->ny = ny;
    this->nnodes = nnodes;
    this->nnodes_t = size_t(nnodes);
    this->dx = dx;
    this->dy = dy;
    this->dxy = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
    this->cellarea = this->dx * this->dy;
    this->Xmin = xmin;
    this->Ymin = ymin;

    // Not a node is utilised to detect when a neighbouring operation returns
    // not a node
    this->not_a_node = -nx * ny - 10;

    this->initialise_neighbourer();

    // Initialise coordinate stuff
    this->Xs = std::vector<T>(this->nx);
    this->Ys = std::vector<T>(this->ny);
    // this->extents = std::vector<T>(4,0);

    for (int i = 0; i < this->nx; ++i)
      this->Xs[i] = this->Xmin + this->dx / 2 + i * this->dx;
    for (int i = this->ny - 1; i >= 0; --i) {
      this->Ys[i] = this->Ymin + this->dy / 2 + i * this->dy;
    }

    this->Xmax = this->Xs.back() + this->dx / 2;
    this->Ymax = this->Ys.back() + this->dy / 2;

    std::reverse(this->Ys.begin(), this->Ys.end());

    this->extents = {this->Xmin, this->Xmin + (this->nx + 1) * this->dx,
                     this->Ymin, this->Ymin + (this->ny + 1) * this->dy};
  }

  void set_default_boundaries(std::string bountype) {

    this->boundary = std::vector<int>(this->nnodes_t, 1);

    if (bountype == "4edges") {
      for (size_t i = 0; i < this->nnodes_t; ++i) {
        if (this->is_on_dem_edge(i))
          this->boundary[i] = 0;
      }
    } else if (bountype == "periodic_EW") {
      for (int i = 0; i < this->nnodes; ++i) {
        if (this->is_on_top_row(i) || this->is_on_bottom_row(i))
          this->boundary[i] = 0;
        else if (this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
          this->boundary[i] = 2;
      }
    } else if (bountype == "periodic_NS") {
      for (int i = 0; i < this->nnodes; ++i) {
        if (this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
          this->boundary[i] = 0;
        else if (this->is_on_top_row(i) || this->is_on_bottom_row(i))
          this->boundary[i] = 2;
      }
    } else {
      throw std::runtime_error("invalid periodic boundaries");
    }
  }

  template <class bou_t> void set_custom_boundaries(bou_t &tbound) {
    auto bound = format_input(tbound);
    this->boundary = std::vector<int>(this->nnodes, 0);
    for (int i = 0; i < this->nnodes; ++i) {
      this->boundary[i] = bound[i];
    }
  }

  // Set all the out boundaries to 3, meaning they can now give to lower
  // elevation neighbours
  void set_out_boundaries_to_permissive() {
    for (auto &v : this->boundary) {
      if (v == 0)
        v = 3;
    }
  }

  int get_boundary_at_node(int i) { return this->boundary[i]; }

  bool is_in_bound(int i) {
    return (i >= 0 && i < this->nnodes) ? true : false;
  }

  void fill_linknodes(std::vector<int> &linknodes) {

    for (int i = 0; i < this->nnodes; ++i) {
      int o = this->get_right_idx(i);
      linknodes[i * 4] = (this->is_in_bound(o) ? i : -1);
      linknodes[i * 4 + 1] = (this->is_in_bound(o) ? o : -1);
      o = this->get_bottom_idx(i);
      linknodes[i * 4 + 2] = (this->is_in_bound(o) ? i : -1);
      linknodes[i * 4 + 3] = (this->is_in_bound(o) ? o : -1);
    }
  }

  // std::vector<int> get_

  int get_right_idx_links(int i) { return i * 2; }
  int get_bottom_idx_links(int i) { return i * 2 + 1; }
  int get_left_idx_links(int i) {
    int yolo = this->get_left_idx(i);
    if (yolo >= 0) {
      return this->get_right_idx_links(yolo);
    } else {
      return this->not_a_node;
    }
  }
  int get_top_idx_links(int i) {
    int yolo = this->get_top_idx(i);
    if (yolo >= 0) {
      return this->get_bottom_idx_links(yolo);
    } else {
      return this->not_a_node;
    }
  }

  int get_right_idx_linknodes(int i) { return 2 * i * 2; }
  int get_bottom_idx_linknodes(int i) { return 2 * (i * 2 + 1); }
  int get_left_idx_linknodes(int i) {
    int yolo = this->get_left_idx(i);
    if (yolo >= 0) {
      return 2 * this->get_right_idx_links(yolo);
    } else {
      return this->not_a_node;
    }
  }
  int get_top_idx_linknodes(int i) {
    int yolo = this->get_top_idx(i);
    if (yolo >= 0) {
      return 2 * this->get_bottom_idx_links(yolo);
    } else {
      return this->not_a_node;
    }
  }

  std::vector<int> get_empty_neighbour() { return std::vector<int>(4, 0); }

  int get_neighbour_idx(int i, std::vector<int> &out) {

    int n = 0;

    int next = this->get_right_idx(i);
    if (this->is_in_bound(next)) {
      if (this->can_flow_even_go_there(next)) {
        out[n] = next;
        ++n;
      }
    }

    next = this->get_bottom_idx(i);
    if (this->is_in_bound(next)) {
      if (this->can_flow_even_go_there(next)) {
        out[n] = next;
        ++n;
      }
    }

    next = this->get_left_idx(i);
    if (this->is_in_bound(next)) {
      if (this->can_flow_even_go_there(next)) {
        out[n] = next;
        ++n;
      }
    }

    next = this->get_top_idx(i);
    if (this->is_in_bound(next)) {
      if (this->can_flow_even_go_there(next)) {
        out[n] = next;
        ++n;
      }
    }

    return n;
  }

  int get_neighbour_idx_links(int i, std::vector<int> &out) {
    int n = 0;
    int next = this->get_right_idx_links(i);
    if (next >= 0 && next < this->nnodes * 2) {
      out[n] = next;
      ++n;
    }

    next = this->get_bottom_idx_links(i);
    if (next >= 0 && next < this->nnodes * 2) {
      out[n] = next;
      ++n;
    }

    next = this->get_left_idx_links(i);
    if (next >= 0 && next < this->nnodes * 2) {
      out[n] = next;
      ++n;
    }

    next = this->get_top_idx_links(i);
    if (next >= 0 && next < this->nnodes * 2) {
      out[n] = next;
      ++n;
    }

    return n;
  }

  int get_neighbour_idx_linknodes(int i, std::vector<int> &out) {
    int n = 0;
    int next = this->get_right_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }

    next = this->get_bottom_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }

    next = this->get_left_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }

    next = this->get_top_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }

    return n;
  }

  int get_neighbour_idx_distance(int i, std::vector<std::pair<int, T>> &out) {
    int n = 0;
    int next = this->get_right_idx(i);

    if (this->can_flow_even_go_there(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }

    next = this->get_bottom_idx(i);
    if (this->can_flow_even_go_there(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }

    next = this->get_left_idx(i);
    if (this->can_flow_even_go_there(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }

    next = this->get_top_idx(i);
    if (this->can_flow_even_go_there(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }

    return n;
  }

  int get_neighbour_idx_distance_links(int i,
                                       std::vector<std::pair<int, T>> &out) {
    int n = 0;
    int next = this->get_right_idx_links(i);
    if (next >= 0 && next < this->nnodes * 2) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }

    next = this->get_bottom_idx_links(i);
    if (next >= 0 && next < this->nnodes * 2) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }

    next = this->get_left_idx_links(i);
    if (next >= 0 && next < this->nnodes * 2) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }

    next = this->get_top_idx_links(i);
    if (next >= 0 && next < this->nnodes * 2) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }

    return n;
  }

  int get_neighbour_idx_distance_linknodes(
      int i, std::vector<std::pair<int, T>> &out) {
    int n = 0;
    int next = this->get_right_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }

    next = this->get_bottom_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }

    next = this->get_left_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }

    next = this->get_top_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }

    return n;
  }

  template <class ii_t>
  ii_t get_neighbour_idx_from_normalised_dxdy(ii_t node, float_t normdx,
                                              float_t normdy) {

    if (std::abs(normdx) > std::abs(normdy)) {
      if (normdx > 0)
        return this->get_right_idx(node);
      else
        return this->get_left_idx(node);
    } else {
      if (normdy > 0)
        return this->get_bottom_idx(node);
      else
        return this->get_top_idx(node);
    }
  }

  template <class i_t, class ii_t>
  std::pair<T, T> get_directed_dxdy_from_links_idx(i_t li, ii_t original_node,
                                                   ii_t n1, ii_t n2) {
    int C_C = (original_node == n1) ? 1 : -1;
    if (li % 4 == 0)
      return std::make_pair(C_C * this->dx, C_C * 0.);
    else if (li % 4 == 1)
      return std::make_pair(0., C_C * this->dy);
    else
      return std::make_pair(C_C * this->dx, C_C * 0.);
  }

  template <class i_t, class ii_t>
  std::pair<T, T> get_dxdy_from_links_idx(i_t li) {
    if (li % 4 == 0)
      return std::make_pair(this->dx, 0.);
    else if (li % 4 == 1)
      return std::make_pair(0., this->dy);

    else
      return std::make_pair(this->dx, 0.);
  }

  // ------------------------------------------------
  //	                             	              __
  //                                             / _)
  //                                    _.----._/ /
  //   Conversion METHODS        /         /
  //                                __/ (  | (  |
  //                               /__.-'|_|--|_|

  // All the methods related to geometrical conversions
  // ------------------------------------------------
  // WILL NEED SOME WORK HERE DEPENDING ON THE GRAPH ORIENTATION AND ALL

  inline void rowcol_from_node_id(int node_index, int &row, int &col) {
    col = node_index % this->nx;
    row = (int)std::floor(node_index / this->nx);
  }

  inline int nodeid_from_row_col(int row, int col) {
    return row * this->nx + col;
  }

  inline int nodeid_from_XY(T X, T Y) {
    int col = floor((X - this->Xmin) / this->dx);
    int row = this->ny - floor((Y - this->Ymin) / this->dy);
    return this->nodeid_from_row_col(row, col);
  }
  // void rowcol_from_XY()

  inline void XY_from_nodeid(int node_index, T &tX, T &tY) {
    int row, col;
    this->rowcol_from_node_id(node_index, row, col);
    this->XY_from_rowcol(row, col, tX, tY);
  }

  inline void XY_from_rowcol(int row, int col, T &tX, T &tY) {
    tX = this->Xs[col];
    tY = this->Ys[row];
  }

  inline T get_X_from_i(int i) {
    int col = i % this->nx;
    return this->Xs[col];
  }

  inline T get_Y_from_i(int i) {
    int row = (int)std::floor(i / this->nx);
    return this->Ys[row];
  }

  // ------------------------------------------------

  //	                             	            __
  //                                             / _)
  //                                    _.----._/ /
  //   Neighbouring METHODS            /         /
  //                                __/ (  | (  |
  //                               /__.-'|_|--|_|

  // All the methods related to accessing and calculating neighbours
  // ------------------------------------------------

  std::vector<Neighbour<int, T>> get_neighbours(int i,
                                                bool ignore_nodata = false) {

    // preformatting the output
    std::vector<Neighbour<int, T>> neighbours;

    // Reserving size depending on the flow topology
    neighbours.reserve(4);

    size_t id_neighbourer = this->_get_neighbourer_id(i);

    // Now I need to determine the index of the neighbourer vector
    // Which provides adder to determine the neighbours f(boundary_conditions)
    // internal node 0
    // periodic_first_row 1
    // periodic_last_row 2
    // periodic_first_col 3
    // periodic last_col 4
    // normal_first_row 5
    // normal_last_row 6
    // normal_first_col 7
    // normal_last_col 8
    // normal_top_left 9
    // normal_top_right 10
    // normal_bottom_left 11
    // normal_bottom_right 12
    // top_left_periodic 13
    // top_right_periodic 14
    // periodic_bottom_left 15
    // periodic_bottom_right 16

    this->set_neighbours_in_vector(neighbours, i, id_neighbourer,
                                   ignore_nodata);

    // and return it
    return neighbours;
  }

  // This function is the helper of the neighbouring function: it fills the
  // neighbours vector with the actul values.
  void set_neighbours_in_vector(std::vector<Neighbour<int, T>> &neighbours,
                                int &i, size_t &id_neighbourer,
                                bool ignore_nodata) {
    // node index of the current neighbour
    int tn;

    // If you wonder why I am not iterating with a vector and everything here,
    // it is for the small perf gain of doing it this way
    tn = this->neighbourer[id_neighbourer][1];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && this->can_flow_even_go_there(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[1]));
    tn = this->neighbourer[id_neighbourer][3];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && this->can_flow_even_go_there(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[3]));
    tn = this->neighbourer[id_neighbourer][0];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && this->can_flow_even_go_there(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[0]));
    tn = this->neighbourer[id_neighbourer][2];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && this->can_flow_even_go_there(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[2]));
  }

  inline size_t _get_neighbourer_id(int i) {

    size_t id_neighbourer = -1;
    // internal node, so neighbourer is 0
    if (this->boundary[i] == 1)
      id_neighbourer = 0;
    else {
      // Case top left corner
      if (i == 0) {
        // Case Periodic
        if (this->boundary[i] == 2)
          id_neighbourer = 13;
        // Case Noraml
        else if (this->boundary[i] == 0 || this->boundary[i] == -1 ||
                 this->boundary[i] == 3 || this->boundary[i] == 4)
          id_neighbourer = 9;
      }
      // Case top right corner
      else if (i == this->nx - 1) {
        // Case Periodic
        if (this->boundary[i] == 2)
          id_neighbourer = 14;
        // Case Noraml
        else if (this->boundary[i] == 0 || this->boundary[i] == -1 ||
                 this->boundary[i] == 3 || this->boundary[i] == 4)
          id_neighbourer = 10;
      }
      // Case bottom left corner
      else if (i == (this->nnodes - this->nx)) {
        // Case Periodic
        if (this->boundary[i] == 2)
          id_neighbourer = 15;
        // Case Noraml
        else if (this->boundary[i] == 0 || this->boundary[i] == -1 ||
                 this->boundary[i] == 3 || this->boundary[i] == 4)
          id_neighbourer = 11;
      }
      // Case bottom right corner
      else if (i == (this->nnodes - 1)) {
        // Case Periodic
        if (this->boundary[i] == 2)
          id_neighbourer = 16;
        // Case Noraml
        else if (this->boundary[i] == 0 || this->boundary[i] == -1 ||
                 this->boundary[i] == 3 || this->boundary[i] == 4)
          id_neighbourer = 12;
      }
      // Cases first row (no corner)
      else if (i < this->nx - 1) {
        // Case periodic
        if (this->boundary[i] == 2)
          id_neighbourer = 1;
        // Case normal
        else if (this->boundary[i] == 0 || this->boundary[i] == -1 ||
                 this->boundary[i] == 3 || this->boundary[i] == 4)
          id_neighbourer = 5;
      }
      // Cases last row (no corner)
      else if (i > (this->ny - 1) * this->nx) {
        // Case periodic
        if (this->boundary[i] == 2)
          id_neighbourer = 2;
        // Case normal
        else if (this->boundary[i] == 0 || this->boundary[i] == -1 ||
                 this->boundary[i] == 3 || this->boundary[i] == 4)
          id_neighbourer = 6;
      }
      // Cases first col (no corner)
      else if (i % this->nx == 0) {
        // Case periodic
        if (this->boundary[i] == 2)
          id_neighbourer = 3;
        // Case normal
        else if (this->boundary[i] == 0 || this->boundary[i] == -1 ||
                 this->boundary[i] == 3 || this->boundary[i] == 4)
          id_neighbourer = 7;
      }
      // Cases last col (no corner)
      else if (i % this->nx == this->nx - 1) {
        // Case periodic
        if (this->boundary[i] == 2)
          id_neighbourer = 4;
        // Case normal
        else if (this->boundary[i] == 0 || this->boundary[i] == -1 ||
                 this->boundary[i] == 3 || this->boundary[i] == 4)
          id_neighbourer = 8;
      } else
        id_neighbourer = 0;
    }

    // if(id_neighbourer == -1)
    // {
    // 	int row,col;
    // 	this->rowcol_from_node_id(i,row,col);
    // 	std::cout << "boundary was " << this->boundary[i] << " i: " << i << "/"
    // << this->nnodes << " row: " << row << " col: " << col << std::endl;
    // throw std::runtime_error("neighbouring issue");
    // }

    return id_neighbourer;
  }

  // D4 single neighbour routines
  Neighbour<int, T> get_left_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][1] + i,
                             this->dx);
  }

  Neighbour<int, T> get_top_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][0] + i,
                             this->dy);
  }

  Neighbour<int, T> get_right_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][2] + i,
                             this->dx);
  }

  Neighbour<int, T> get_bottom_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][3] + i,
                             this->dy);
  }

  // D4 single neighbour routines
  int get_left_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][1] + i;
  }

  int get_top_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][0] + i;
  }

  int get_right_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][2] + i;
  }

  int get_bottom_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][3] + i;
  }

  std::vector<int> get_D4_neighbours_only_id(int i) {
    std::vector<int> neighs;
    neighs.reserve(4);
    int tn = this->get_left_idx(i);
    if (tn >= 0)
      neighs.emplace_back(tn);
    tn = this->get_top_idx(i);
    if (tn >= 0)
      neighs.emplace_back(tn);
    tn = this->get_right_idx(i);
    if (tn >= 0)
      neighs.emplace_back(tn);
    tn = this->get_bottom_idx(i);
    if (tn >= 0)
      neighs.emplace_back(tn);
    return neighs;
  }

  // Some of my algorithm are adapted from richdem and require iterating through
  // neighbours the way they do starting from left and going clockwise around
  // the node
  std::vector<int> get_neighbours_idx_richdemlike(int i) {
    std::vector<int> neighs;
    neighs.reserve(4);
    neighs.emplace_back(this->get_left_idx(i));
    neighs.emplace_back(this->get_top_idx(i));
    neighs.emplace_back(this->get_right_idx(i));
    neighs.emplace_back(this->get_bottom_idx(i));
    return neighs;
  }

  // Method to test whether a node can outlet flow OUT of the model
  inline bool can_flow_out_there(int i) {
    if (this->boundary[i] == 3 || this->boundary[i] == 0)
      return true;
    else {
      return false;
    }
  }

  // method to check if a node can even accept flow
  inline bool can_flow_even_go_there(int i) {
    if (this->boundary[i] < 0)
      return false;
    else {
      return true;
    }
  }

  inline bool is_active(int i) {
    if (this->can_flow_out_there(i) && this->boundary[i] != 3)
      return false;
    else if (this->boundary[i] == 3)
      return true;
    else if (this->boundary[i] == 4)
      return false;

    if (this->can_flow_even_go_there(i))
      return true;
    return false;
  }

  // Returns true if the node is on the dem edge
  // (i.e. one of the 4 lines)
  bool is_on_dem_edge(int i) {
    if (this->is_on_top_row(i) || this->is_on_bottom_row(i) ||
        this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
      return true;
    return false;
  }

  bool is_on_top_row(int i) {
    if (i < this->nx)
      return true;
    return false;
  }

  bool is_on_bottom_row(int i) {
    if (i >= this->nnodes - this->nx)
      return true;
    else
      return false;
  }

  bool is_on_leftest_col(int i) {
    if (i % this->nx == 0)
      return true;
    else
      return false;
  }

  bool is_on_rightest_col(int i) {
    if (i % this->nx == this->nx - 1)
      return true;
    else
      return false;
  }

  void print_dim() {
    std::cout << "nx:" << this->nx << std::endl;
    std::cout << "ny:" << this->ny << std::endl;
    std::cout << "nnodes:" << this->nnodes << std::endl;
    std::cout << "dx:" << this->dx << std::endl;
    std::cout << "dy:" << this->dy << std::endl;
  }

  template <class topo_t, class out_t> out_t get_HS(topo_t &ttopography) {
    auto topography = format_input(ttopography);
    double altitude = 45;
    double azimuth = 315;
    double z_factor = 1;
    double pi = 3.1415926;

    // std::vector<double> hillshade(ptr, ptr + this->nnodes);
    std::vector<double> hillshade(this->nnodes, 0.);

    // convert zenith and azimuth into radians for calculation
    double zenith_rad = (90 - altitude) * pi / 180.0;
    double azimuth_math = 360 - azimuth + 90;
    if (azimuth_math >= 360.0)
      azimuth_math = azimuth_math - 360;
    double azimuth_rad = azimuth_math * pi / 180.0;

    for (int i = 0; i < this->nnodes; ++i) {
      // Ignoring no data
      if (this->boundary[i] < 0)
        continue;

      double slope_rad = 0;
      double aspect_rad = 0;
      double dzdx = 0;
      double dzdy = 0;

      double ij = topography[i];
      double ijp1 = topography[this->get_right_neighbour(i).node];
      double ip1j = topography[this->get_bottom_neighbour(i).node];
      double im1j = topography[this->get_top_neighbour(i).node];
      double ijm1 = topography[this->get_left_neighbour(i).node];

      if (ij > 0) {
        dzdx = ((ijp1 + 2 * ip1j) - (2 * im1j)) / (4 * this->dx);
        dzdy = ((2 * ijp1) - (2 * ijm1)) / (4 * this->dy);
        slope_rad = atan(z_factor * sqrt((dzdx * dzdx) + (dzdy * dzdy)));
        if (dzdx != 0) {
          aspect_rad = std::atan2(dzdy, (dzdx * -1));
          if (aspect_rad < 0)
            aspect_rad = 2 * pi + aspect_rad;
        } else {
          if (dzdy > 0)
            aspect_rad = pi / 2;
          else if (dzdy < 0)
            aspect_rad = 2 * pi - pi / 2;
        }

        hillshade[i] = ((std::cos(zenith_rad) * std::cos(slope_rad)) +
                        (std::sin(zenith_rad) * std::sin(slope_rad) *
                         std::cos(azimuth_rad - aspect_rad)));
        // std::cout << hillshade[i] << "|";
        if (hillshade[i] < 0)
          hillshade[i] = 0;
      }
    }

    return format_output<decltype(hillshade), out_t>(hillshade);
  }

  // Largely adapted from RichDEM
  template <class topo_t>
  std::vector<double> PriorityFlood(topo_t &ttopography) {
    std::random_device rd;  // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_real_distribution<> distr(1e-7, 1e-6); // define the range

    auto tttopography = format_input(ttopography);
    std::vector<double> topography = to_vec(tttopography);

    PQ_i_d open;
    std::queue<PQ_helper<int, double>> pit;

    uint64_t processed_cells = 0;
    uint64_t pitc = 0;
    int false_pit_cells = 0;
    double PitTop = -9999;

    std::vector<bool> closed(this->nnodes, false);

    // RDLOG_PROGRESS<<"Adding cells to the priority queue...";
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->can_flow_out_there(i)) {
        closed[i] = true;
        open.emplace(i, topography[i]);
      }
    }

    // RDLOG_PROGRESS<<"Performing Priority-Flood+Epsilon...";
    auto neighbours = this->get_empty_neighbour();

    while (open.size() > 0 || pit.size() > 0) {
      PQ_helper<int, double> c;

      if (pit.size() > 0 && open.size() > 0 &&
          open.top().score == pit.front().score) {
        c = open.top();
        open.pop();
        PitTop = -9999;
      } else if (pit.size() > 0) {
        c = pit.front();
        pit.pop();
        if (PitTop == -9999)
          PitTop = topography[c.node];
      } else {
        c = open.top();
        open.pop();
        PitTop = -9999;
      }
      processed_cells++;

      int nn = this->get_neighbour_idx(c.node, neighbours);
      for (int tnn = 0; tnn < nn; ++nn) {
        int n = neighbours[tnn];
        if (this->is_active(n) == false)
          continue;

        if (closed[n])
          continue;

        closed[n] = true;
        if (topography[n] <=
            std::nextafter(c.score, std::numeric_limits<double>::infinity())) {
          if (PitTop != -9999 && PitTop < topography[n] &&
              std::nextafter(c.score,
                             std::numeric_limits<double>::infinity()) >=
                  topography[n])
            ++false_pit_cells;

          ++pitc;
          // topography[n]=std::nextafter(c.score,std::numeric_limits<double>::infinity());
          topography[n] = c.score + distr(gen);
          pit.emplace(n, topography[n]);
        } else
          open.emplace(n, topography[n]);
      }
    }

    return topography;
  }

  template <class topo_t>
  void InitPriorityQue(topo_t &topography, std::vector<bool> &flag,
                       PQ_i_d &priorityQueue) {

    std::queue<PQ_helper<int, double>> depressionQue;

    // push border cells into the PQ
    auto neighbours = this->get_empty_neighbour();
    for (int i = 0; i < this->nnodes; ++i) {
      if (flag[i])
        continue;

      if (this->can_flow_even_go_there(i) == false) {
        flag[i] = true;
        int nn = this->get_neighbour_idx(i, neighbours);

        for (int tnn = 0; tnn < nn; ++tnn) {

          int n = neighbours[tnn];

          if (flag[n])
            continue;
          if (this->can_flow_even_go_there(n)) {
            priorityQueue.emplace(n, topography[n]);
            flag[n] = true;
          }
        }
      } else if (this->can_flow_out_there(i)) {
        // on the DEM border
        priorityQueue.emplace(i, topography[i]);
        flag[i] = true;
      }
    }
  }

  template <class topo_t>
  void ProcessTraceQue(topo_t &topography, std::vector<bool> &flag,
                       std::queue<PQ_helper<int, double>> &traceQueue,
                       PQ_i_d &priorityQueue) {
    std::queue<PQ_helper<int, double>> potentialQueue;
    int indexThreshold = 2; // index threshold, default to 2
    auto neighbours = this->get_empty_neighbour();
    while (!traceQueue.empty()) {
      const auto node = traceQueue.front();
      traceQueue.pop();

      bool Mask[5][5] = {{false}, {false}, {false}, {false}, {false}};

      int nn = this->get_neighbour_idx(node.node, neighbours);

      for (int tnn = 0; tnn < nn; ++tnn) {
        int n = neighbours[tnn];

        if (flag[n])
          continue;

        if (topography[n] > node.score) {
          traceQueue.emplace(n, topography[n]);
          flag[n] = true;
        } else {
          // initialize all masks as false
          bool have_spill_path_or_lower_spill_outlet =
              false; // whether cell n has a spill path or a lower spill outlet
                     // than node if n is a depression cell
          auto nneighs = this->get_neighbours_idx_richdemlike(n);
          int row, col;
          this->rowcol_from_node_id(node.node, row, col);
          int nrow, ncol;
          this->rowcol_from_node_id(n, nrow, ncol);
          int incr = 0;
          for (auto nn : nneighs) {
            ++incr;
            // Checking node validity
            if (nn < 0 || nn >= this->nnodes)
              continue;

            int nnrow, nncol;
            this->rowcol_from_node_id(nn, nnrow, nncol);
            if ((Mask[nnrow - row + 2][nncol - col + 2]) ||
                (flag[nn] && topography[nn] < node.score)) {
              Mask[nrow - row + 2][ncol - col + 2] = true;
              have_spill_path_or_lower_spill_outlet = true;
              break;
            }
          }

          if (!have_spill_path_or_lower_spill_outlet) {
            if (incr < indexThreshold)
              potentialQueue.push(node);
            else
              priorityQueue.push(node);
            break; // make sure node is not pushed twice into PQ
          }
        }
      } // end of for loop
    }

    while (!potentialQueue.empty()) {
      const auto node = potentialQueue.front();
      potentialQueue.pop();

      // first case
      auto neigh = this->get_neighbours_idx_richdemlike(node.node);
      int incr = 0;
      for (auto n : neigh) {
        ++incr;

        if (flag[n])
          continue;

        priorityQueue.push(node);
        break;
      }
    }
  }

  template <class topo_t>
  void ProcessPit(topo_t &topography, std::vector<bool> &flag,
                  std::queue<PQ_helper<int, double>> &depressionQue,
                  std::queue<PQ_helper<int, double>> &traceQueue) {
    auto neighbours = this->get_empty_neighbour();
    while (!depressionQue.empty()) {
      auto node = depressionQue.front();
      depressionQue.pop();
      int nn = this->get_neighbour_idx(node.node, neighbours);
      for (int tnn = 0; tnn < nn; ++tnn) {
        int n = neighbours[tnn];
        if (flag[n])
          continue;

        const auto iSpill = topography[n];
        if (iSpill > node.score) { // Dlope cell
          flag[n] = true;
          traceQueue.emplace(n, iSpill);
          continue;
        }
        // Depression cell
        flag[n] = true;
        topography[n] = node.score + 1e-6 * this->randu->get();
        depressionQue.emplace(n, topography[n]);
      }
    }
  }

  template <class topo_t>
  std::vector<double> PriorityFlood_Wei2018(topo_t &ttopography) {
    std::queue<PQ_helper<int, double>> traceQueue;
    std::queue<PQ_helper<int, double>> depressionQue;

    std::vector<bool> flag(this->nnodes, false);
    std::vector<double> topography = to_vec(ttopography);

    PQ_i_d priorityQueue;

    this->InitPriorityQue(topography, flag, priorityQueue);

    auto neighbours = this->get_empty_neighbour();

    while (!priorityQueue.empty()) {
      const auto tmpNode = priorityQueue.top();
      priorityQueue.pop();

      int nn = this->get_neighbour_idx(tmpNode.node, neighbours);

      for (int tnn = 0; tnn < nn; ++tnn) {
        int n = neighbours[tnn];
        if (flag[n])
          continue;

        auto iSpill = topography[n];

        if (iSpill <= tmpNode.score) {
          // depression cell
          topography[n] = tmpNode.score + 1e-3 + 1e-6 * this->randu->get();
          flag[n] = true;
          depressionQue.emplace(n, topography[n]);
          this->ProcessPit(topography, flag, depressionQue, traceQueue);
        } else {
          // slope cell
          flag[n] = true;
          traceQueue.emplace(n, iSpill);
        }
        ProcessTraceQue(topography, flag, traceQueue, priorityQueue);
      }
    }

    return topography;
  }

  // this function makes the connector generic and allow other connector to have
  // varying cellarea
  template <class i_t> T get_area_at_node(i_t i) { return this->cellarea; }

  // template<>
  // void fill_neighbour_matrices()

  template <class i_t> T get_dx_from_links_idx(i_t i) {
    if (i % 2 == 0)
      return this->dx;
    else if (i % 2 == 1)
      return this->dy;
    else
      return this->dx;
  }

  template <class i_t> T get_traverse_dx_from_links_idx(i_t i) {
    if (i % 2 == 0)
      return this->dy;
    else if (i % 2 == 1)
      return this->dx;
    else
      return this->dy;
  }

  // return the orthogonal node from a pair of node / link indices
  template <class i_t>
  std::pair<i_t, i_t> get_orthogonal_nodes(i_t node, i_t link) {

    if (node % 4 == 0)
      return {this->get_top_idx(node), this->get_bottom_idx(node)};
    else if (node % 4 == 1)
      return {this->get_left_idx(node), this->get_right_idx(node)};
    else
      throw std::runtime_error(
          "Fatal error in DAGGER::D4connector::get_orthogonal_nodes");
  }

  template <class i_t> i_t linkidx_from_nodes(i_t n1, i_t n2) {
    if (n1 < n2)
      std::swap(n1, n2);
    int delta = n2 - n1;

    if (delta == 1)
      return n1 * 4;
    else if (delta == this->nnodes + 1)
      return n1 * 4 + 1;
    else
      throw std::runtime_error(
          "Fatal error in DAGGER::D8connector::linkidx_from_nodes");
  }

  template <class i_t> T get_dx_from_linknodes_idx(i_t i) {
    int j = std::floor(i / 2);
    return this->get_dx_from_links_idx(j);
  }

  T get_travers_dy_from_dx(T tdx) {
    if (tdx == this->dx)
      return this->dy;
    else if (tdx == this->dy)
      return this->dx;
    else
      return tdx;
  }

  std::vector<int> get_ilinknodes_from_node_v1(int i) {
    std::vector<int> out;
    out.reserve(4);
    int tn = this->get_right_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 2)
      out.emplace_back(tn);
    tn = this->get_bottom_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 2)
      out.emplace_back(tn);
    tn = this->get_left_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 2)
      out.emplace_back(tn);
    tn = this->get_top_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 2)
      out.emplace_back(tn);
    return out;
  }

  std::vector<std::pair<int, bool>> get_ilinknodes_from_nodev2(int i) {
    std::vector<std::pair<int, bool>> out;
    out.reserve(4);

    if (this->boundary[i] != 1) {
      int tn = this->get_right_idx_links(i);
      bool val = false;
      if (tn > 0 && tn < this->nnodes * 2)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_bottom_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 2)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_left_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 2)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_top_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 2)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
    } else
      out = {
          std::make_pair(this->get_right_idx_links(i), true),
          std::make_pair(this->get_bottom_idx_links(i), true),
          std::make_pair(this->get_left_idx_links(i), true),
          std::make_pair(this->get_top_idx_links(i), true),
      };
    return out;
  }

  void get_ilinknodes_from_nodev3(int i,
                                  std::vector<std::pair<int, bool>> &out) {

    int tn = this->get_right_idx_links(i);
    bool val = false;
    if (tn > 0 && tn < this->nnodes * 2)
      val = true;
    out[0].first = tn;
    out[0].second = val;

    tn = this->get_bottom_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 2)
      val = true;
    out[1].first = tn;
    out[1].second = val;

    tn = this->get_left_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 2)
      val = true;
    out[2].first = tn;
    out[2].second = val;

    tn = this->get_top_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 2)
      val = true;
    out[3].first = tn;
    out[3].second = val;
  }

  int get_ilinknodes_from_node(int i, std::vector<int> &out) {
    int nn = 0;
    int tn = this->get_right_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 2) {
      out[nn] = tn;
      ++nn;
    }

    tn = this->get_bottom_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 2) {
      out[nn] = tn;
      ++nn;
    }

    tn = this->get_left_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 2) {
      out[nn] = tn;
      ++nn;
    }

    tn = this->get_top_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 2) {
      out[nn] = tn;
      ++nn;
    }

    return nn;
  }

  std::vector<bool> get_mask_array() {
    std::vector<bool> mask(this->nnodes, true);
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->is_active(i) == false)
        mask[i] = false;
    }
    return mask;
  }

  template <class topo_t>
  void set_values_at_boundaries(topo_t &tarray, float_t val) {
    auto array = format_input(tarray);
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->can_flow_out_there(i))
        array[i] = val;
    }
  }

  int get_nneighbours() { return 4; }
};

// end of namespace
}; // namespace DAGGER

#endif
