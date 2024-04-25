//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef D8connector_HPP
#define D8connector_HPP

// STL imports
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <initializer_list>
#include <iomanip>
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

#include "boundary_conditions.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

namespace DAGGER {

// Useful namespace for my priority queue
using PQ_i_d = std::priority_queue<PQ_helper<int, double>,
                                   std::vector<PQ_helper<int, double>>,
                                   std::greater<PQ_helper<int, double>>>;

template <class T> class D8connector {
public:
  // General informations about the graph
  // #-> number of nodes (integer and unsized int for loops or list innit).
  int nnodes = 0;
  size_t nnodes_t = 0;

  static const int nneighbours = 8;
  static const size_t nneighbours_t = 8;

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

  bool stochastic_slope_on = false;
  T stochastic_slope_coeff = 1.; // (0 -> deactivates)
  void set_stochaticiy_for_SFD(T val) {
    if (val <= 0) {
      this->stochastic_slope_on = false;
      this->stochastic_slope_coeff = 1.;
    } else {
      this->stochastic_slope_on = true;
      this->stochastic_slope_coeff = val;
    }
  }

  // DEPRECATED BOUNDARY SYSTEM
  // // Bearer of values
  // // #->boundary: index: node index, value: boundary value
  // // #-->The possible values are:
  // // #---> -1 = no in, no out, not process, internal or external
  // // #--->  0 = in, no out, all fluxes leave the system, internal or external
  // // #--->  1 = "normal" cell, in, out, internal
  // // #--->  2 = external periodic boundary
  // // #--->	3 = in, can out but can also give to neighbours
  // std::vector<int> boundary;

  // NEW BOUNDARY SYSTEM
  CodeBC boundaries;

  // Helpers for neighbouring operations
  // -> neighbourer holdsthe indices to loop through for each boundary condition
  std::vector<std::vector<int>> neighbourer;
  // -> lengthener is the dx on each directions
  std::vector<T> lengthener;

  // uint8_t vector for each link indicating its directionality:
  // - 0 is inverse (linknode 2 --> linknode 1)
  // - 1 is normal (linknode 1 --> linknode 2)
  // - 3 is invalid (link index may exist, but is either inexsitant (index is
  // conserved for speed reason when the link index is calculated from node
  // index) ) or temporarily invalid due to dynamic boundary conditions The
  // index itself is calculated from a connector
  std::vector<std::uint8_t> links;

  // integer vector of 2*links size with the node indices of each link
  // for example, the nodes of link #42 would be indices 84 and 85
  std::vector<int> linknodes;

  // Single graph receivers
  // -> Sreceivers: steepest recervers (nnodes size),
  // -> number of donors (nnodes size),
  // -> Steepest donors (nnodes * nneighbours size)
  // --> Sdonors of node i are located from index i*nneighbours to index
  // i*nneighbours + nSdonors[i] not included
  std::vector<int> Sreceivers, nSdonors, Sdonors;

  // Single graph distance to receivers
  std::vector<T> Sdistance2receivers;

  // Steepest slope
  std::vector<T> SS;

  // Coordinate stuff
  // Xs and Ys are vectors of nx and ny size converting row to y and col to X
  // Extents holds the cxmin,xmax,ymin,ymax (extent is an option in matplotlib
  // imshow plots)
  std::vector<T> Xs, Ys, extents;

  std::shared_ptr<easyRand> randu = std::make_shared<easyRand>();

  // Default constructor, empty
  D8connector(){};

  // construct directly from dimensions
  D8connector(int nx, int ny, T dx, T dy, T xmin, T ymin) {
    // initialisation is offset to dedicated function because some languages
    // like Julia are complicating non default constructor initialisation.
    this->init_dimensions(nx, ny, dx, dy, xmin, ymin);
    this->_allocate_vectors();
    this->fill_linknodes();
  }

  template <class topo_t> void update_links_from_topo(topo_t &ttopo) {
    auto topo = format_input(ttopo);
    this->update_links(topo);
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
    T diag = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
    this->dxy = diag;
    this->dxmin = std::min(this->dx, this->dy);
    this->dxmax = std::max(this->dx, this->dy);

    this->lengthener =
        std::initializer_list<T>{diag, dy, diag, dx, dx, diag, dy, diag};
    this->neighbourer.clear();

    // these vectors are additioned to the node indice to test the neighbors
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx - 1, -this->nx, -this->nx + 1, -1, 1, this->nx - 1, this->nx,
        this->nx + 1}); // internal node 0
    this->neighbourer.emplace_back(std::initializer_list<int>{
        (this->ny - 1) * this->nx - 1, (this->ny - 1) * this->nx,
        (this->ny - 1) * this->nx + 1, -1, 1, this->nx - 1, this->nx,
        this->nx + 1}); // periodic_first_row 1
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx - 1, -this->nx, -this->nx + 1, -1, 1,
        -(this->ny - 1) * this->nx - 1, -(this->ny - 1) * this->nx,
        -(this->ny - 1) * this->nx + 1}); // periodic_last_row 2
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -1, -this->nx, -this->nx + 1, (this->nx - 1), 1, 2 * this->nx - 1,
        this->nx, this->nx + 1}); // periodic_first_col 3
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx - 1, -this->nx, -2 * this->nx + 1, -1, -this->nx + 1,
        this->nx - 1, this->nx, 1}); // periodic last_col 4
    this->neighbourer.emplace_back(std::initializer_list<int>{
        this->not_a_node, this->not_a_node, this->not_a_node, -1, 1,
        this->nx - 1, this->nx, this->nx + 1}); // normal_first_row 5
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx - 1, -this->nx, -this->nx + 1, -1, 1, this->not_a_node,
        this->not_a_node, this->not_a_node}); // normal_last_row 6
    this->neighbourer.emplace_back(std::initializer_list<int>{
        this->not_a_node, -this->nx, -this->nx + 1, this->not_a_node, 1,
        this->not_a_node, this->nx, this->nx + 1}); // normal_first_col 7
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx - 1, -this->nx, this->not_a_node, -1, this->not_a_node,
        this->nx - 1, this->nx, this->not_a_node}); // normal_last_col 8
    this->neighbourer.emplace_back(std::initializer_list<int>{
        this->not_a_node, this->not_a_node, this->not_a_node, this->not_a_node,
        1, this->not_a_node, this->nx, this->nx + 1}); // normal_top_left 9
    this->neighbourer.emplace_back(std::initializer_list<int>{
        this->not_a_node, this->not_a_node, this->not_a_node, -1,
        this->not_a_node, this->nx - 1, this->nx,
        this->not_a_node}); // normal_top_right 10
    this->neighbourer.emplace_back(std::initializer_list<int>{
        this->not_a_node, -this->nx, -this->nx + 1, this->not_a_node, 1,
        this->not_a_node, this->not_a_node,
        this->not_a_node}); // normal_bottom_left 11
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx - 1, -this->nx, this->not_a_node, -1, this->not_a_node,
        this->not_a_node, this->not_a_node,
        this->not_a_node}); // normal_bottom_right 12
    this->neighbourer.emplace_back(std::initializer_list<int>{
        this->ny * this->nx - 1, (this->ny - 1) * this->nx,
        (this->ny - 1) * this->nx + 1, this->nx - 1, 1, 2 * this->nx - 1,
        this->nx, this->nx + 1}); // top_left_periodic 13
    this->neighbourer.emplace_back(std::initializer_list<int>{
        (this->ny - 1) * this->nx - 1, (this->ny - 1) * this->nx,
        (this->ny - 2) * this->nx + 1, -1, -this->nx + 1, this->nx - 1,
        this->nx, 1}); // top_right_periodic 14
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -1, -this->nx, -this->nx + 1, this->nx - 1, 1,
        -(this->ny - 2) * this->nx - 1, -(this->ny - 1) * this->nx,
        -(this->ny - 1) * this->nx + 1}); // periodic_bottom_left 15
    this->neighbourer.emplace_back(std::initializer_list<int>{
        -this->nx - 1, -this->nx, -this->nx + 1, -1, 1 - this->nx + 1,
        -(this->ny - 1) * this->nx - 1, -(this->ny - 1) * this->nx,
        -(this->ny) * this->nx + 1}); // periodic_bottom_right 16
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
    this->cellarea = this->dx * this->dy;
    this->Xmin = xmin;
    this->Ymin = ymin;

    // Not a node is utilised to detect when a neighbouring operation returns
    // not a node
    this->not_a_node = -nx * ny * 2;

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

    this->boundaries.codes = std::vector<BC>(this->nnodes_t, BC::FLOW);

    if (bountype == "4edges") {
      for (size_t i = 0; i < this->nnodes_t; ++i) {
        if (this->is_on_dem_edge(i))
          this->boundaries.codes[i] = BC::FORCE_OUT;
      }
    } else if (bountype == "periodic_EW") {
      for (int i = 0; i < this->nnodes; ++i) {
        if (this->is_on_top_row(i) || this->is_on_bottom_row(i))
          this->boundaries.codes[i] = BC::FORCE_OUT;
        else if (this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
          this->boundaries.codes[i] = BC::PERIODIC_BORDER;
      }
    } else if (bountype == "periodic_NS") {
      for (int i = 0; i < this->nnodes; ++i) {
        if (this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
          this->boundaries.codes[i] = BC::FORCE_OUT;
        else if (this->is_on_top_row(i) || this->is_on_bottom_row(i))
          this->boundaries.codes[i] = BC::PERIODIC_BORDER;
      }
    } else {
      throw std::runtime_error("invalid periodic boundaries");
    }

    // Recomputing fillnodes
    // this->fill_linknodes();
    this->_allocate_vectors();
    this->fill_linknodes();
  }

  // Set all the out boundaries to 3, meaning they can now give to lower
  // elevation neighbours
  void set_out_boundaries_to_permissive() {
    for (auto &v : this->boundaries.codes) {
      if (v == BC::FORCE_OUT)
        v = BC::OUT;
    }
  }

  template <class bou_t> void set_custom_boundaries(bou_t &tbound) {
    auto bound = format_input(tbound);

    std::vector<BC> ubound(bound.size(), BC::FLOW);
    for (int i = 0; i < tbound.size(); ++i)
      ubound[i] = static_cast<BC>(bound[i]);

    this->boundaries.set_codes(ubound);
    this->_allocate_vectors();
    this->fill_linknodes();
  }

  BC get_boundary_at_node(int i) { return this->boundaries.codes[i]; }

  bool is_in_bound(int i) {
    return (i >= 0 && i < this->nnodes) ? true : false;
  }

  // void fill_linknodes()
  // {
  // 	bool peridic = false;
  // 	for(int i=0; i<this->nnodes; ++i)
  // 	{

  // 		if(this->boundaries.is_periodic(i))
  // 			peridic = true;

  // 		bool NDT = this->boundaries.can_create_link(i);
  // 		bool i_forcing = this->boundaries.forcing_io(i);

  // 		int o = this->get_right_idx(i);
  // 		bool linkvalid = (NDT &&
  // this->is_in_bound_and_can_create_link(o)); 		bool
  // both_same_forcing
  // =
  // i_forcing && this->boundaries.forcing_io(o);
  // if(both_same_forcing
  // &&
  // ((this->boundaries.force_giving(o) && this->boundaries.force_giving(i)) ||
  // 					 (this->boundaries.force_receiving(o) &&
  // this->boundaries.force_receiving(i)))) 			linkvalid =
  // false; 		this->linknodes[i*8] = (linkvalid) ? i:-1;
  // this->linknodes[i*8 + 1] = (linkvalid) ? o:-1;

  // 		o = this->get_bottomright_idx(i);
  // 		linkvalid = (NDT &&  this->is_in_bound_and_can_create_link(o));
  // 		both_same_forcing = i_forcing && this->boundaries.forcing_io(o);
  // 		if(both_same_forcing && ((this->boundaries.force_giving(o) &&
  // this->boundaries.force_giving(i)) ||
  // (this->boundaries.force_receiving(o)
  // && this->boundaries.force_receiving(i)))) 			linkvalid =
  // false; 		this->linknodes[i*8 + 2] = (linkvalid) ? i:-1;
  // this->linknodes[i*8 + 3] = (linkvalid) ? o:-1;

  // 		o = this->get_bottom_idx(i);
  // 		linkvalid = (NDT &&  this->is_in_bound_and_can_create_link(o));
  // 		both_same_forcing = i_forcing && this->boundaries.forcing_io(o);
  // 		if(both_same_forcing && ((this->boundaries.force_giving(o) &&
  // this->boundaries.force_giving(i)) ||
  // (this->boundaries.force_receiving(o)
  // && this->boundaries.force_receiving(i)))) 			linkvalid =
  // false; 		this->linknodes[i*8 + 4] = (linkvalid) ? i:-1;
  // this->linknodes[i*8 + 5] = (linkvalid) ? o:-1;

  // 		o = this->get_bottomleft_idx(i);
  // 		linkvalid = (NDT &&  this->is_in_bound_and_can_create_link(o));
  // 		both_same_forcing = i_forcing && this->boundaries.forcing_io(o);
  // 		if(both_same_forcing && ((this->boundaries.force_giving(o) &&
  // this->boundaries.force_giving(i)) ||
  // (this->boundaries.force_receiving(o)
  // && this->boundaries.force_receiving(i)))) 			linkvalid =
  // false; 		this->linknodes[i*8 + 6] = (linkvalid) ? i:-1;
  // this->linknodes[i*8 + 7] = (linkvalid) ? o:-1;

  // 	}

  // 	if(peridic == false)
  // 		return;

  // 	// for(int i=0; i<this->nnodes; ++i)
  // 	// {

  // 	// }

  // }

  void fill_linknodes() {
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->is_in_bound_and_can_create_link(i) == false)
        continue;

      int o = this->get_right_idx(i);
      this->_fill_linknode(i, o, 0);
      o = this->get_bottomright_idx(i);
      this->_fill_linknode(i, o, 2);
      o = this->get_bottom_idx(i);
      this->_fill_linknode(i, o, 4);
      o = this->get_bottomleft_idx(i);
      this->_fill_linknode(i, o, 6);

      if (this->boundaries.is_periodic(i))
      // if(true)
      {
        o = this->get_left_idx(i);
        this->_fill_linknode(o, i, 0);
        o = this->get_topleft_idx(i);
        this->_fill_linknode(o, i, 2);
        o = this->get_top_idx(i);
        this->_fill_linknode(o, i, 4);
        o = this->get_topright_idx(i);
        this->_fill_linknode(o, i, 6);
      }
    }
  }

  void _fill_linknode(int f, int o, int add) {

    // if(this->is_in_bound_and_can_create_link(f) == false ||
    // 	this->is_in_bound_and_can_create_link(o) == false ||
    // 	(this->boundaries.force_giving(f) && this->boundaries.force_giving(o))
    // || 	(this->boundaries.force_receiving(f) &&
    // this->boundaries.force_receiving(o))
    // 	)

    if (this->is_in_bound_and_can_create_link(f) == false ||
        this->is_in_bound_and_can_create_link(o) == false)
      return;

    if (this->linknodes[f * 8 + add] != -1)
      return;

    if ((this->boundaries.can_give(f) && this->boundaries.can_receive(o)) ||
        (this->boundaries.can_give(o) && this->boundaries.can_receive(f))) {
      this->linknodes[f * 8 + add] = f;
      this->linknodes[f * 8 + add + 1] = o;
    }

    else {
      this->linknodes[f * 8 + add] = -1;
      this->linknodes[f * 8 + add + 1] = -1;
    }

    // if( this->boundaries.force_giving(f) && this->boundaries.can_receive(o) )
  }

  template <class i_t>
  void get_nodes_from_linkidx_implicit(int li, int &n1, int &n2) {
    n1 = floor(li / 4);
    n2 = n1 + li % 4;

    bool nodes_can_create_links =
        (this->is_in_bound_and_can_create_link(n1) == false ||
         this->is_in_bound_and_can_create_link(n1) == false);
    bool both_forced = false;
    if (nodes_can_create_links) {
      if (this->boundaries.forcing_io(n1) && this->boundaries.forcing_io(n2))
        both_forced = true;
    }

    // Then not a valid link
    if (nodes_can_create_links || both_forced) {
      n1 = -1;
      n2 = -1;
    }
  }

  // Function updating ONLY the MFD receivers
  // This is useful in the cases where SFD recs are conditionned by an other
  // mean and cannot be touched (e.g. Cordonnier)
  template <class topo_t> void update_links_MFD_only(topo_t &topography) {
    // iterating though every links
    for (size_t i = 0; i < this->links.size(); ++i) {
      // Getting hte 2 nodes of the current link
      int from = this->linknodes[i * 2];
      int to = this->linknodes[i * 2 + 1];

      if (from == -1) {
        this->links[i] = 2;
        continue;
      }

      if (this->boundaries.forcing_io(from) ||
          this->boundaries.forcing_io(to)) {

        if (this->boundaries.force_giving(from) ||
            this->boundaries.force_receiving(to))
          this->links[i] = 1;
        else if (this->boundaries.force_giving(to) ||
                 this->boundaries.force_receiving(from))
          this->links[i] = 0;
        else {
          throw std::runtime_error("Should not happen 777");
          this->links[i] = 2;
        }

        continue;
      }

      // by convention true -> topo1 > topo2
      if (topography[from] > topography[to] &&
          this->boundaries.can_give(from) && this->boundaries.can_receive(to))
        this->links[i] = 1;
      else if (this->boundaries.can_give(to) &&
               this->boundaries.can_receive(from))
        this->links[i] = 0;
      // If the configuration cannot allow the link, it is temporarily disabled
      else
        this->links[i] = 2;
    }
    // done
  }

  // Updates all the link and the SFD info
  template <class topo_t> void update_links(topo_t &topography) {
    // std::cout << "BUNT" << std::endl;

    // am I using a stochastic adjustment for deciding on the steepest slope
    // iterating through all the nodes
    for (size_t i = 0; i < this->links.size(); ++i) {

      // Checking the validity of the link
      if (this->linknodes[i * 2] < 0) {
        this->links[i] = 2;
        continue;
      }

      // Getting ht etwo nodes of the links
      int from = this->linknodes[i * 2];
      int to = this->linknodes[i * 2 + 1];

      // std::cout << std::setprecision(18);

      // if( 808 * this->nx + 733 == from)
      // 	std::cout << i <<" TESTED NODE from:: " << topography[from] << "
      // v " << topography[to]  << "||" << (topography[from] > topography[to])
      // << std::endl; if( 808 * this->nx + 733 == to) 	std::cout << i <<"
      // TESTED NODE to:: " << topography[to] << " v " << topography[from] <<
      // "||" << (topography[to] > topography[from]) << std::endl;

      // std::cout << from << "|" << to << "||";

      // getting the link infos
      // -> dx
      T dx = this->get_dx_from_links_idx(i);
      // -> slope
      T slope = (topography[from] - topography[to]) / dx;

      if (this->boundaries.forcing_io(from) ||
          this->boundaries.forcing_io(to)) {

        if (this->boundaries.force_giving(from) ||
            this->boundaries.force_receiving(to))
          slope = 1 / std::abs(slope);

        else if (this->boundaries.force_giving(to) ||
                 this->boundaries.force_receiving(from))
          slope = -1 / std::abs(slope);
      }

      if (this->stochastic_slope_on)
        slope *= (this->randu->get() * this->stochastic_slope_coeff) + 1e-6;

      // if slope is positive, to is the receiver by convention
      if (slope > 0 && this->boundaries.can_give(from) &&
          this->boundaries.can_receive(to)) {
        // Conventional direction
        this->links[i] = 1;
        // if Steepest Slope is higher than the current recorded one
        if (this->SS[from] < slope) {
          // saving the Sreceivers info as temporary best choice
          this->Sreceivers[from] = to;
          this->Sdistance2receivers[from] = dx;
          this->SS[from] = slope;
        }
      } else if (this->boundaries.can_give(to) &&
                 this->boundaries.can_receive(from)) {
        // Otherwise the convention is inverted:
        // isrec is falese and to is giving to from
        this->links[i] = 0;
        // NOte that slope is absolute values
        slope = std::abs(slope);
        if (this->SS[to] < slope) {
          this->Sreceivers[to] = from;
          this->Sdistance2receivers[to] = dx;
          this->SS[to] = slope;
        }
      } else
        this->links[i] = 2;
    }

    // Finally inverting the Sreceivers into the Sdonors info
    // Required for several routines
    this->compute_SF_donors_from_receivers();
  }

  // Fucntion inverting the SFD receivers into donors
  void compute_SF_donors_from_receivers() {
    // Initialising the graph dimesions for the donors
    // All of thenm have the graph dimension
    this->Sdonors = std::vector<int>(this->nnodes * this->nneighbours, -1);
    this->nSdonors = std::vector<int>(this->nnodes, 0);

    // iterating through all the nodes
    for (int i = 0; i < this->nnodes; ++i) {
      // SF so rid == i cause there is only 1 rec
      int trec = this->Sreceivers[i];
      if (trec == i)
        continue;

      // feeding the Sdonors array at rec position with current node and...
      this->Sdonors[trec * this->nneighbours + this->nSdonors[trec]] = i;
      // ... incrementing the number of Sdonors
      this->nSdonors[trec] += 1;
    }
    // done
  }

  // Same function than above but without reallocating the memory (can save a
  // bit of time depending on the context)
  void recompute_SF_donors_from_receivers() {

    for (int i = 0; i < this->nnodes; ++i) {
      for (int j = 0; j < this->nneighbours; ++j)
        this->Sdonors[i * this->nneighbours + j] = -1;
      this->nSdonors[i] = 0;
    }

    for (int i = 0; i < this->nnodes; ++i) {
      // SF so rid == i cause there is only 1 rec
      int trec = this->Sreceivers[i];
      if (trec == i)
        continue;
      this->Sdonors[trec * this->nneighbours + this->nSdonors[trec]] = i;
      this->nSdonors[trec] += 1;
    }
  }

  // Helper functions to allocate and reallocate vectors when
  // computing/recomputing the graph
  void _allocate_vectors() {
    this->links =
        std::vector<std::uint8_t>(int(this->nnodes * this->nneighbours / 2), 2);
    this->linknodes =
        std::vector<int>(int(this->nnodes * this->nneighbours), -1);
    this->Sreceivers = std::vector<int>(this->nnodes, -1);
    for (int i = 0; i < this->nnodes; ++i)
      this->Sreceivers[i] = i;
    this->Sdistance2receivers = std::vector<T>(this->nnodes, -1);
    this->SS = std::vector<T>(this->nnodes, 0.);
  }

  void _reallocate_vectors() {
    for (int i = 0; i < this->nnodes; ++i) {
      this->Sreceivers[i] = i;
      this->Sdistance2receivers[i] = 0;
      this->SS[i] = 0;
    }
  }

  bool is_in_bound_and_can_create_link(int o) {
    bool linkvalid = (this->is_in_bound(o));
    if (linkvalid)
      linkvalid = this->boundaries.can_create_link(o);
    return linkvalid;
  }

  int get_right_idx_links(int i) { return i * 4; }
  int get_bottomright_idx_links(int i) { return i * 4 + 1; }
  int get_bottom_idx_links(int i) { return i * 4 + 2; }
  int get_bottomleft_idx_links(int i) { return i * 4 + 3; }
  int get_left_idx_links(int i) {
    int yolo = this->get_left_idx(i);
    if (yolo >= 0) {
      return this->get_right_idx_links(yolo);
    } else {
      return this->not_a_node;
    }
  }
  int get_topleft_idx_links(int i) {
    int yolo = this->get_topleft_idx(i);
    if (yolo >= 0) {
      return this->get_bottomright_idx_links(yolo);
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
  int get_topright_idx_links(int i) {
    int yolo = this->get_topright_idx(i);
    if (yolo >= 0) {
      return this->get_bottomleft_idx_links(yolo);
    } else {
      return this->not_a_node;
    }
  }
  int get_right_idx_linknodes(int i) { return 2 * i * 4; }
  int get_bottomright_idx_linknodes(int i) { return 2 * (i * 4 + 1); }
  int get_bottom_idx_linknodes(int i) { return 2 * (i * 4 + 2); }
  int get_bottomleft_idx_linknodes(int i) { return 2 * (i * 4 + 3); }
  int get_left_idx_linknodes(int i) {
    int yolo = this->get_left_idx(i);
    if (yolo >= 0) {
      return 2 * this->get_right_idx_links(yolo);
    } else {
      return this->not_a_node;
    }
  }
  int get_topleft_idx_linknodes(int i) {
    int yolo = this->get_topleft_idx(i);
    if (yolo >= 0) {
      return 2 * this->get_bottomright_idx_links(yolo);
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
  int get_topright_idx_linknodes(int i) {
    int yolo = this->get_topright_idx(i);
    if (yolo >= 0) {
      return 2 * this->get_bottomleft_idx_links(yolo);
    } else {
      return this->not_a_node;
    }
  }

  std::vector<int> get_empty_neighbour() { return std::vector<int>(8, 0); }

  int get_neighbour_idx(int i, std::vector<int> &out) {

    int n = 0;

    int next = this->get_right_idx(i);
    if (this->is_in_bound(next)) {
      if (!this->boundaries.no_data(next)) {
        out[n] = next;
        ++n;
      }
    }
    next = this->get_bottomright_idx(i);
    if (this->is_in_bound(next)) {
      if (!this->boundaries.no_data(next)) {
        out[n] = next;
        ++n;
      }
    }
    next = this->get_bottom_idx(i);
    if (this->is_in_bound(next)) {
      if (!this->boundaries.no_data(next)) {
        out[n] = next;
        ++n;
      }
    }
    next = this->get_bottomleft_idx(i);
    if (this->is_in_bound(next)) {
      if (!this->boundaries.no_data(next)) {
        out[n] = next;
        ++n;
      }
    }
    next = this->get_left_idx(i);
    if (this->is_in_bound(next)) {
      if (!this->boundaries.no_data(next)) {
        out[n] = next;
        ++n;
      }
    }
    next = this->get_topleft_idx(i);
    if (this->is_in_bound(next)) {
      if (!this->boundaries.no_data(next)) {
        out[n] = next;
        ++n;
      }
    }
    next = this->get_top_idx(i);
    if (this->is_in_bound(next)) {
      if (!this->boundaries.no_data(next)) {
        out[n] = next;
        ++n;
      }
    }
    next = this->get_topright_idx(i);
    if (this->is_in_bound(next)) {
      if (!this->boundaries.no_data(next)) {
        out[n] = next;
        ++n;
      }
    }

    return n;
  }

  int get_neighbour_idx_links(int i, std::vector<int> &out) {
    int n = 0;
    int next = this->get_right_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }
    next = this->get_bottomright_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }
    next = this->get_bottom_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }
    next = this->get_bottomleft_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }
    next = this->get_left_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }
    next = this->get_topleft_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }
    next = this->get_top_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }
    next = this->get_topright_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = next;
      ++n;
    }

    if (this->is_on_dem_edge(i)) {
    }

    return n;
  }

  int get_neighbour_idx_linknodes(int i, std::vector<int> &out) {
    int n = 0;
    int next = this->get_right_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = next;
      ++n;
    }
    next = this->get_bottomright_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = next;
      ++n;
    }
    next = this->get_bottom_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = next;
      ++n;
    }
    next = this->get_bottomleft_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = next;
      ++n;
    }
    next = this->get_left_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = next;
      ++n;
    }
    next = this->get_topleft_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = next;
      ++n;
    }
    next = this->get_top_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = next;
      ++n;
    }
    next = this->get_topright_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = next;
      ++n;
    }
    return n;
  }

  int get_neighbour_idx_distance(int i, std::vector<std::pair<int, T>> &out) {
    int n = 0;
    int next = this->get_right_idx(i);

    if (!this->boundaries.no_data(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }
    next = this->get_bottomright_idx(i);
    if (!this->boundaries.no_data(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_bottom_idx(i);
    if (!this->boundaries.no_data(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }
    next = this->get_bottomleft_idx(i);
    if (!this->boundaries.no_data(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_left_idx(i);
    if (!this->boundaries.no_data(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }
    next = this->get_topleft_idx(i);
    if (!this->boundaries.no_data(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_top_idx(i);
    if (!this->boundaries.no_data(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }
    next = this->get_topright_idx(i);
    if (!this->boundaries.no_data(next) && this->is_in_bound(next)) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }

    return n;
  }

  int get_neighbour_idx_distance_links(int i,
                                       std::vector<std::pair<int, T>> &out) {
    int n = 0;
    int next = this->get_right_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }
    next = this->get_bottomright_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_bottom_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }
    next = this->get_bottomleft_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_left_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }
    next = this->get_topleft_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_top_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }
    next = this->get_topright_idx_links(i);
    if (next >= 0 && next < this->nnodes * 4) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    return n;
  }

  int get_neighbour_idx_distance_linknodes(
      int i, std::vector<std::pair<int, T>> &out) {
    int n = 0;
    int next = this->get_right_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }
    next = this->get_bottomright_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_bottom_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }
    next = this->get_bottomleft_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_left_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = std::make_pair(next, this->dx);
      ++n;
    }
    next = this->get_topleft_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    next = this->get_top_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = std::make_pair(next, this->dy);
      ++n;
    }
    next = this->get_topright_idx_linknodes(i);
    if (next >= 0 && next < this->nnodes * 8) {
      out[n] = std::make_pair(next, this->dxy);
      ++n;
    }
    return n;
  }

  // Takes a node index as input and a normalised displacement in the X and Y
  // direction (between -1 and 1) Returns the neighbouring node in the given
  // direction
  template <class ii_t>
  ii_t get_neighbour_idx_from_normalised_dxdy(ii_t node, T normdx, T normdy) {
    if (normdx >= 0.5) {
      if (normdy <= -0.5)
        return this->get_topright_idx(node);
      else if (normdy <= 0.5)
        return this->get_right_idx(node);
      else
        return this->get_bottomright_idx(node);
    } else if (normdx >= -0.5) {
      if (normdy <= 0)
        return this->get_top_idx(node);
      else
        return this->get_bottom_idx(node);
    } else {
      if (normdy <= -0.5)
        return this->get_topleft_idx(node);
      else if (normdy < 0.5)
        return this->get_left_idx(node);
      else
        return this->get_bottomleft_idx(node);
    }
  }

  // get receivers of node i and put them in the recs vector fed in
  // It returns the number of receivers in the recs vector
  // THis whole process optimises repeated receiver fetching, by never
  // reallocating/initialising the vector recs
  int get_receivers_idx(int i, std::vector<int> &recs) {
    // getting the related links stroing them temporarily in the rec vec
    int nli = this->get_neighbour_idx_links(i, recs);

    // going through the linksß
    // The idx is the idx of insertion in the recs vectors
    // the idx of insertion is always <= of the index of reading (both receivers
    // and links-to-assess are stored in the recs vector)
    int idx = 0;
    // counter used to keep track of the number of receivers: starts at the
    // number of links related to the given node and get decremented at each
    // donor link
    int newli = nli;

    // Iterating through the links
    for (int ti = 0; ti < nli; ++ti) {
      // current link index
      int li = recs[ti];

      // this link is a rec if the node is the first of the linknode and links
      // is true
      if (i == this->linknodes[li * 2] && this->links[li] == 1) {
        // in which case the receiver of the current node is the +  1
        recs[idx] = this->linknodes[li * 2 + 1];
        ++idx;
      }
      // OR if the current node is the +1 and the links false
      else if (this->links[li] == 0 && i == this->linknodes[li * 2 + 1]) {
        // in which case the receivers is the 0 node
        recs[idx] = this->linknodes[li * 2];
        ++idx;
      } else {
        // otherwise, it's not a rec and we decrease the newli
        --newli;
      }
    }
    // recs is changed in place, and we return the number of recs newli
    return newli;
  }

  // Getting the id of the receivers in the links array
  // see get_receivers_idx for full comments about the section
  int get_receivers_idx_links(int i, std::vector<int> &recs) {
    // getting the related links
    int nli = this->get_neighbour_idx_links(i, recs);

    // going through the linksß
    int idx = 0;
    int newli = nli;

    for (int ti = 0; ti < nli; ++ti) {
      // checking the orientation
      int li = recs[ti];
      if (i == this->linknodes[li * 2] && this->links[li]) {
        recs[idx] = li;
        ++idx;
      } else if (this->links[li] == 0 && i == this->linknodes[li * 2 + 1]) {
        recs[idx] = li;
        ++idx;
      } else {
        --newli;
      }
    }
    return newli;
  }

  int get_nth_steepest_donors_of_node(int node, int j) {
    return this->Sdonors[node * this->nneighbours + j];
  }

  // Getting donor indicies
  // see get_receivers_idx for full comments about the section
  int get_donors_idx(int i, std::vector<int> &dons) {
    // getting the related links
    int nli = this->get_neighbour_idx_links(i, dons);

    // going through the links
    int idx = 0;
    int newli = nli;

    for (int ti = 0; ti < nli; ++ti) {
      // checking the orientation
      int li = dons[ti];
      if (i == this->linknodes[li * 2] && this->links[li] == 0) {
        dons[idx] = this->linknodes[li * 2 + 1];
        ++idx;
      } else if (this->links[li] == 1 && i == this->linknodes[li * 2 + 1]) {
        dons[idx] = this->linknodes[li * 2];
        ++idx;
      } else {
        --newli;
      }
    }
    return newli;
  }

  // getting links indices of hte donors (in the links array)
  // see get_receivers_idx for full comments about the section
  int get_donors_idx_links(int i, std::vector<int> &dons) {
    // getting the related links
    int nli = this->get_neighbour_idx_links(i, dons);

    // going through the linksß
    int idx = 0;
    int newli = nli;

    for (int ti = 0; ti < nli; ++ti) {
      // checking the orientation
      int li = dons[ti];
      if (i == this->linknodes[li * 2] && this->links[li] == 0) {
        dons[idx] = li;
        ++idx;
      } else if (this->links[li] == 1 && i == this->linknodes[li * 2 + 1]) {
        dons[idx] = li;
        ++idx;
      } else {
        --newli;
      }
    }
    return newli;
  }

  // Debug function printing to the prompt the single receiver of a node
  // Will probably get deprecated
  std::vector<int> get_rowcol_Sreceivers(int row, int col) {
    int node = this->nodeid_from_row_col(row, col);
    std::vector<int> out_receivers;
    int trow, tcol;
    this->rowcol_from_node_id(this->Sreceivers[node], trow, tcol);
    out_receivers = std::vector<int>{trow, tcol};

    std::cout << "Srec is " << this->Sreceivers[node] << " node was " << node
              << std::endl;
    return out_receivers;
  }

  template <class topo_t> void print_receivers(int i, topo_t &ttopography) {
    std::cout << std::setprecision(12);
    auto topography = format_input<topo_t>(ttopography);

    auto receivers = this->get_empty_neighbour();
    int nn = this->get_receivers_idx(i, receivers);

    std::cout << "Topography is " << topography[i] << "# receivers: " << nn
              << std::endl;
    for (int tr = 0; tr < nn; ++tr) {
      int r = receivers[tr];
      int row, col;
      this->rowcol_from_node_id(r, row, col);
      std::cout << "Rec " << r << " row " << row << " col " << col << " topo "
                << topography[r] << std::endl;
    }

    auto neighbours = this->get_empty_neighbour();
    nn = this->get_neighbour_idx(i, neighbours);
    std::cout << "Neighbours are :" << std::endl;

    for (int tr = 0; tr < nn; ++tr) {
      int r = neighbours[tr];
      int row, col;
      this->rowcol_from_node_id(r, row, col);
      std::cout << "Neighbour " << r << " row " << row << " col " << col
                << " topo " << topography[r] << std::endl;
    }
  }

  // Returns the number of links stored in the graph
  // Note that it comprises some unvalid linked!
  int get_rec_array_size() { return int(this->links.size()); }

  /// Takes an array of nnodes size and sum the values at the outlets
  /// This can be useful for checking mass balances for example
  /// if true, include_internal_pits allow the code to add internal unprocessed
  /// pits, wether they are on purpose or not
  template <class array_t, class U>
  U sum_at_outlets(array_t &tarray, bool include_internal_pits = true) {
    std::cout << "DEPRECATION WARNING::sum_at_outlets::should be moved as a "
                 "standalone algorithm"
              << std::endl;
    auto array = format_input<array_t>(tarray);
    U out = 0;
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->Sreceivers[i] == i) {
        if (include_internal_pits) {
          out += array[i];
        } else if (this->flow_out_model(i)) {
          out += array[i];
        }
      }
    }
    return out;
  }

  /// Takes an array of nnodes size and sum the values at the outlets
  /// This can be useful for checking mass balances for example
  /// if true, include_internal_pits allow the code to add internal unprocessed
  /// pits, wether they are on purpose or not
  template <class array_t, class out_t>
  out_t keep_only_at_outlets(array_t &tarray,
                             bool include_internal_pits = true) {
    auto array = format_input<array_t>(tarray);
    std::vector<T> out = std::vector<T>(this->nnodes, 0);
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->Sreceivers[i] == i) {
        if (include_internal_pits)
          out[i] = array[i];
        else if (this->flow_out_model(i))
          out[i] = array[i];
      }
    }
    return format_output<decltype(out), out_t>(out);
  }

  // Checks the validity of a link
  // Invalid links have a linknode value of -1
  template <class ti_t> bool is_link_valid(ti_t i) {
    return (this->links[i] == 2) ? false : true;
  }

  // return true is the node is active: i.e. flow transfer through it
  // reflects the connector version of the function, but adds extra graph
  // specific checks. For example a node with permissive outletting border will
  // be active in the connector sense, but can be inactive in the graph if it
  // has no downstream neighbours.
  template <class ti_t> bool flow_out_model(ti_t node) {

    bool con_val = this->boundaries.can_out(node);
    if (con_val && this->Sreceivers[node] == int(node))
      return true;
    return false;
  }

  template <class ti_t> bool is_pit(ti_t node) {

    bool con_val = this->boundaries.can_out(node);
    if (!con_val && this->Sreceivers[node] == int(node))
      return true;
    return false;
  }

  template <class ti_t> bool flow_out_or_pit(ti_t node) {
    if (this->Sreceivers[node] == int(node))
      return true;
    return false;
  }

  // Returns the pair of node making a link, starting from the donor to the
  // receiver
  template <class ti_t> std::pair<ti_t, ti_t> get_from_to_links(ti_t i) {
    if (this->links[i] == 1)
      return std::make_pair(this->linknodes[i * 2], this->linknodes[i * 2 + 1]);
    else if (this->links[i] == 0)
      return std::make_pair(this->linknodes[i * 2 + 1], this->linknodes[i * 2]);
    else
      return {-1, -1};
  }

  template <class ti_t>
  void get_from_to_links(ti_t i, std::pair<ti_t, ti_t> &fromto) {
    if (this->links[i] == 1) {
      fromto.first = this->linknodes[i * 2];
      fromto.second = this->linknodes[i * 2 + 1];
    } else if (this->links[i] == 2) {
      fromto.first = this->linknodes[i * 2 + 1];
      fromto.second = this->linknodes[i * 2];
    } else
      fromto = {-1, -1};
  }

  template <class ti_t> ti_t get_from_links(ti_t i) {
    if (this->links[i] == 1)
      return this->linknodes[i * 2];
    else if (this->links[i] == 0)
      return this->linknodes[i * 2 + 1];
    else
      return -1;
  }

  template <class ti_t> ti_t get_to_links(ti_t i) {
    if (this->links[i] == 0)
      return this->linknodes[i * 2];
    else if (this->links[i] == 1)
      return this->linknodes[i * 2 + 1];
    else
      return -1;
  }

  template <class ti_t> ti_t get_other_node_from_links(ti_t li, ti_t ni) {
    if (this->linknodes[li * 2 + 1] == ni)
      return this->linknodes[li * 2];
    else if (this->linknodes[li * 2] == ni)
      return this->linknodes[li * 2 + 1];
    else
      return -1;
  }

  std::vector<int> get_n_receivers() {
    std::vector<int> nrecs(this->nnodes, 0);
    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i)) {
        auto frto = this->get_from_to_links(i);
        ++nrecs[frto.first];
      }
    }
    return nrecs;
  }

  template <class out_t> out_t get_SFD_receivers() {
    return format_output<std::vector<int>, out_t>(this->Sreceivers);
  }

  template <class out_t> out_t get_SFD_dx() {
    return format_output<std::vector<T>, out_t>(this->Sdistance2receivers);
  }

  template <class out_t> out_t get_SFD_ndonors() {
    return format_output<std::vector<int>, out_t>(this->nSdonors);
  }

  template <class out_t> out_t get_SFD_donors_flat() {
    return format_output<std::vector<int>, out_t>(this->Sdonors);
  }

  template <class out_t> out_t get_SFD_donors_list() {
    std::vector<std::vector<int>> out(this->nnodes);
    for (int i = 0; i < this->nnodes; ++i) {
      std::vector<int> tvec;
      for (int j = 0; j < this->nSdonors[i]; ++j)
        tvec.emplace_back(this->Sdonors[i * this->nneighbours + j]);
      out[i] = tvec;
    }

    return out;
  }

  template <class out_t> out_t get_links() { return this->links; }

  template <class out_t> out_t get_linknodes_flat() {
    return format_output<std::vector<int>, out_t>(this->linknodes);
  }

  template <class out_t> out_t get_linknodes_flat_D4() {
    std::vector<int> temp(int(this->linknodes.size() / 2), -1);
    int j = 0;
    int counter = -1;
    for (int i = 0; i < int(this->linknodes.size()); i += 2) {
      ++counter;
      if (counter == 0 || counter == 2) {
        temp[j] = this->linknodes[i];
        ++j;
        temp[j] = this->linknodes[i + 1];
        ++j;
      }

      if (counter == 3)
        counter = -1;
    }

    return format_output<std::vector<int>, out_t>(temp);
  }

  template <class out_t> out_t get_linkdx_flat_D4() {
    std::vector<T> temp(int(this->links.size() / 2), -1);
    int j = 0;
    int counter = -1;
    for (int i = 0; i < int(this->links.size()); ++i) {
      ++counter;

      if (counter == 0 || counter == 2) {
        T dx = this->get_dx_from_links_idx(i);
        temp[j] = dx;
        ++j;
      }

      if (counter == 3)
        counter = -1;
    }

    return format_output<std::vector<T>, out_t>(temp);
  }

  template <class out_t> out_t get_linknodes_list() {
    std::vector<std::vector<int>> out(this->links.size());
    for (size_t i = 0; i < this->links.size(); ++i) {
      out[i] =
          std::vector<int>{this->linknodes[i * 2], this->linknodes[i * 2 + 1]};
    }
    return out;
  }

  template <class out_t> out_t get_linknodes_list_oriented() {
    std::vector<std::vector<int>> out(this->links.size());
    for (size_t i = 0; i < this->links.size(); ++i) {
      out[i] = (this->links[i] >= 1)
                   ? std::vector<int>{this->linknodes[i * 2],
                                      this->linknodes[i * 2 + 1]}
                   : std::vector<int>{this->linknodes[i * 2 + 1],
                                      this->linknodes[i * 2]};
    }
    return out;
  }

  int get_SFD_receivers_at_node(int i) { return this->Sreceivers[i]; }

  int get_SFD_dx_at_node(int i) { return this->Sdistance2receivers[i]; }

  int get_SFD_ndonors_at_node(int i) { return this->nSdonors[i]; }

  template <class out_t> out_t get_SFD_donors_at_node(int i) {
    std::vector<int> out(this->nneighbours);
    for (int j = 0; j < this->nSdonors[i]; ++j)
      out.emplace_back(this->Sdonors[i * this->nneighbours + j]);
    return out;
  }

  /*
          =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
          =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
          =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
                                . - ~ ~ ~ - .
        ..     _      .-~               ~-.
       //|     \ `..~                      `.
      || |      }  }              /       \  \
  (\   \\ \~^..'                 |         }  \
   \`.-~  o      /       }       |        /    \
   (__          |       /        |       /      `.
    `- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
                |     /          |     /     ~-.     ~- _
                |_____|          |_____|         ~ - . _ _~_-_

          =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
          =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
          Functions computing gradients
          =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
          =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
  */

  template <class out_t, class topo_t>
  out_t get_SFD_gradient(topo_t &ttopography) {
    auto topography = format_input<topo_t>(ttopography);
    auto gradient = this->_get_SFD_gradient(topography);
    return format_output<decltype(gradient), out_t>(gradient);
  }

  template <class topo_t> std::vector<T> _get_SFD_gradient(topo_t &topography) {
    std::vector<T> gradient(this->nnodes, 0.);
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->Sreceivers[i] != i)
        gradient[i] = (topography[i] - topography[this->Sreceivers[i]]) /
                      this->Sdistance2receivers[i];
    }
    return gradient;
  }

  template <class out_t, class topo_t>
  out_t get_links_gradient(topo_t &ttopography, T min_slope) {
    auto topography = format_input<topo_t>(ttopography);
    std::vector<T> gradient = this->_get_links_gradient(topography, min_slope);
    return format_output<decltype(gradient), out_t>(gradient);
  }

  template <class topo_t>
  std::vector<T> _get_links_gradient(topo_t &topography, T min_slope) {

    std::vector<T> gradient = std::vector<T>(this->links.size(), 0);

    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i)) {
        gradient[i] = std::max((topography[this->linknodes[i * 2]] -
                                topography[this->linknodes[i * 2 + 1]]) /
                                   this->get_dx_from_links_idx(i),
                               min_slope);
      }
    }

    return gradient;
  }

  template <class out_t, class topo_t>
  out_t get_MFD_mean_gradient(topo_t &ttopography) {
    auto topography = format_input<topo_t>(ttopography);
    auto gradient = this->_get_MFD_mean_gradient(topography);
    return format_output<decltype(gradient), out_t>(gradient);
  }

  template <class topo_t>
  std::vector<T> _get_MFD_mean_gradient(topo_t &topography) {

    std::vector<T> gradient = std::vector<T>(this->nnodes, 0);
    std::vector<int> ngradient = std::vector<int>(this->nnodes, 0);

    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i)) {
        T this_gradient = std::abs(topography[this->linknodes[i * 2] -
                                              this->linknodes[i * 2 + 1]]) /
                          this->get_dx_from_links_idx(i);
        auto frto = this->get_from_to_links(i);
        gradient[frto.first] += this_gradient;
        ++ngradient[frto.first];
      }
    }

    for (int i = 0; i < this->nnodes; ++i) {
      if (ngradient[i] > 0)
        gradient[i] = gradient[i] / ngradient[i];
    }

    return gradient;
  }

  template <class out_t, class topo_t>
  out_t get_MFD_weighted_gradient(topo_t &ttopography, topo_t &tweights) {
    auto topography = format_input<topo_t>(ttopography);
    auto weights = format_input<topo_t>(tweights);
    auto gradient = this->_get_MFD_weighted_gradient(topography, weights);
    return format_output<decltype(gradient), out_t>(gradient);
  }

  template <class topo_t>
  std::vector<T> _get_MFD_weighted_gradient(topo_t &topography,
                                            topo_t &weights) {

    std::vector<T> gradient = std::vector<T>(this->nnodes, 0);
    std::vector<T> wgradient = std::vector<T>(this->nnodes, 0.);

    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i)) {
        T this_gradient = std::abs(topography[this->linknodes[i * 2] -
                                              this->linknodes[i * 2 + 1]]) /
                          this->get_dx_from_links_idx(i);
        auto frto = this->get_from_to_links(i);
        gradient[frto.first] += this_gradient * weights[i];
        wgradient[frto.first] += weights[i];
      }
    }

    for (int i = 0; i < this->nnodes; ++i) {
      if (wgradient[i] > 0)
        gradient[i] = gradient[i] / wgradient[i];
    }

    return gradient;
  }

  /*
  =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
  =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
  =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
                        . - ~ ~ ~ - .
..     _      .-~               ~-.
//|     \ `..~                      `.
|| |      }  }              /       \  \
(\   \\ \~^..'                 |         }  \
\`.-~  o      /       }       |        /    \
(__          |       /        |       /      `.
`- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
        |     /          |     /     ~-.     ~- _
        |_____|          |_____|         ~ - . _ _~_-_

  =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
  =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
  Functions to calculate weights
  =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
  =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
  */

  template <class out_t, class topo_t>
  out_t get_link_weights(topo_t &tgradient, T exp) {
    std::vector<T> weights(this->links.size(), 0.);
    auto gradient = format_input<topo_t>(tgradient);

    if (exp <= 0) {
      this->_get_link_weights_f_nrecs(weights);
    } else if (exp == 1) {
      this->_get_link_weights_proposlope(weights, gradient);
    } else {
      this->_get_link_weights_exp(weights, gradient, exp);
    }

    return format_output<decltype(weights), out_t>(weights);
  }

  void _get_link_weights_f_nrecs(std::vector<T> &weights) {
    auto nrecs = this->get_n_receivers();
    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i)) {
        int nr = nrecs[this->get_from_links(i)];
        if (nr > 0) {
          weights[i] = 1. / nr;
        } else
          weights[i] = 1.;
      }
    }
  }

  template <class topo_t>
  void _get_link_weights_proposlope(std::vector<T> &weights, topo_t &gradient) {
    std::vector<T> sumgrad(this->nnodes, 0.);
    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i)) {
        sumgrad[this->get_from_links(i)] += gradient[i];
      }
    }

    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i))
        weights[i] = gradient[i] / sumgrad[this->get_from_links(i)];
    }
  }

  template <class topo_t>
  void _get_link_weights_exp(std::vector<T> &weights, topo_t &gradient, T exp) {
    std::vector<T> sumgrad(this->nnodes, 0.);
    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i)) {
        sumgrad[this->get_from_links(i)] += std::pow(gradient[i], exp);
      }
    }

    for (size_t i = 0; i < this->links.size(); ++i) {
      if (this->is_link_valid(i))
        weights[i] =
            std::pow(gradient[i], exp) / sumgrad[this->get_from_links(i)];
    }
  }

  // ------------------------------------------------
  //	                             	              __
  //                                             / _)
  //                                    _.----._/ /
  //   Conversion METHODS              /         /
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
    neighbours.reserve(8);

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
         (ignore_nodata && !this->boundaries.no_data(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[1]));
    tn = this->neighbourer[id_neighbourer][3];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && !this->boundaries.no_data(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[3]));
    tn = this->neighbourer[id_neighbourer][4];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && !this->boundaries.no_data(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[4]));
    tn = this->neighbourer[id_neighbourer][6];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && !this->boundaries.no_data(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[6]));
    tn = this->neighbourer[id_neighbourer][0];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && !this->boundaries.no_data(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[0]));
    tn = this->neighbourer[id_neighbourer][2];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && !this->boundaries.no_data(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[2]));
    tn = this->neighbourer[id_neighbourer][5];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && !this->boundaries.no_data(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[5]));
    tn = this->neighbourer[id_neighbourer][7];
    if (tn != this->not_a_node &&
        (ignore_nodata == false ||
         (ignore_nodata && !this->boundaries.no_data(tn))))
      neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[7]));
  }

  inline size_t _get_neighbourer_id(int i) {

    size_t id_neighbourer = -1;
    // internal node, so neighbourer is 0
    if (this->boundaries.is_normal_node(i))
      id_neighbourer = 0;
    else {
      // Case top left corner
      if (i == 0) {
        // Case Periodic
        if (this->boundaries.is_periodic(i))
          id_neighbourer = 13;
        // Case Noraml
        else if (this->boundaries.normal_neighbouring_at_boundary(i))
          id_neighbourer = 9;
      }
      // Case top right corner
      else if (i == this->nx - 1) {
        // Case Periodic
        if (this->boundaries.is_periodic(i))
          id_neighbourer = 14;
        // Case Noraml
        else if (this->boundaries.normal_neighbouring_at_boundary(i))
          id_neighbourer = 10;
      }
      // Case bottom left corner
      else if (i == (this->nnodes - this->nx)) {
        // Case Periodic
        if (this->boundaries.is_periodic(i))
          id_neighbourer = 15;
        // Case Noraml
        else if (this->boundaries.normal_neighbouring_at_boundary(i))
          id_neighbourer = 11;
      }
      // Case bottom right corner
      else if (i == (this->nnodes - 1)) {
        // Case Periodic
        if (this->boundaries.is_periodic(i))
          id_neighbourer = 16;
        // Case Noraml
        else if (this->boundaries.normal_neighbouring_at_boundary(i))
          id_neighbourer = 12;
      }
      // Cases first row (no corner)
      else if (i < this->nx - 1) {
        // Case periodic
        if (this->boundaries.is_periodic(i))
          id_neighbourer = 1;
        // Case normal
        else if (this->boundaries.normal_neighbouring_at_boundary(i))
          id_neighbourer = 5;
      }
      // Cases last row (no corner)
      else if (i > (this->ny - 1) * this->nx) {
        // Case periodic
        if (this->boundaries.is_periodic(i))
          id_neighbourer = 2;
        // Case normal
        else if (this->boundaries.normal_neighbouring_at_boundary(i))
          id_neighbourer = 6;
      }
      // Cases first col (no corner)
      else if (i % this->nx == 0) {
        // Case periodic
        if (this->boundaries.is_periodic(i))
          id_neighbourer = 3;
        // Case normal
        else if (this->boundaries.normal_neighbouring_at_boundary(i))
          id_neighbourer = 7;
      }
      // Cases last col (no corner)
      else if (i % this->nx == this->nx - 1) {
        // Case periodic
        if (this->boundaries.is_periodic(i))
          id_neighbourer = 4;
        // Case normal
        else if (this->boundaries.normal_neighbouring_at_boundary(i))
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

  inline size_t _get_neighbourer_id_assume_periodic(int i) {

    size_t id_neighbourer = -1;
    // internal node, so neighbourer is 0
    if (this->boundaries.is_normal_node(i))
      id_neighbourer = 0;
    else {
      // Case top left corner
      if (i == 0) {
        // Case Periodic
        id_neighbourer = 13;
      }
      // Case top right corner
      else if (i == this->nx - 1) {
        // Case Periodic

        id_neighbourer = 14;

      }
      // Case bottom left corner
      else if (i == (this->nnodes - this->nx)) {

        id_neighbourer = 15;

      }
      // Case bottom right corner
      else if (i == (this->nnodes - 1)) {

        id_neighbourer = 16;

      }
      // Cases first row (no corner)
      else if (i < this->nx - 1) {

        id_neighbourer = 1;

      }
      // Cases last row (no corner)
      else if (i > (this->ny - 1) * this->nx) {

        id_neighbourer = 2;

      }
      // Cases first col (no corner)
      else if (i % this->nx == 0) {

        id_neighbourer = 3;

      }
      // Cases last col (no corner)
      else if (i % this->nx == this->nx - 1) {
        id_neighbourer = 4;

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
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][3] + i,
                             this->dx);
  }

  Neighbour<int, T> get_top_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][1] + i,
                             this->dy);
  }

  Neighbour<int, T> get_right_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][4] + i,
                             this->dx);
  }

  Neighbour<int, T> get_bottom_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][6] + i,
                             this->dy);
  }

  // D8 single neighbour extra routines
  Neighbour<int, T> get_topleft_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][0] + i,
                             this->dxy);
  }

  Neighbour<int, T> get_topright_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][2] + i,
                             this->dxy);
  }

  Neighbour<int, T> get_bottomleft_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][5] + i,
                             this->dxy);
  }

  Neighbour<int, T> get_bottomright_neighbour(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return Neighbour<int, T>(this->neighbourer[id_neighbourer][7] + i,
                             this->dxy);
  }

  // D4 single neighbour routines
  int get_left_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][3] + i;
  }

  int get_top_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][1] + i;
  }

  int get_right_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][4] + i;
  }

  int get_bottom_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][6] + i;
  }

  // D8 single neighbour extra routines
  int get_topleft_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][0] + i;
  }

  int get_topright_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][2] + i;
  }

  int get_bottomleft_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][5] + i;
  }

  int get_bottomright_idx(int i) {
    size_t id_neighbourer = this->_get_neighbourer_id(i);
    return this->neighbourer[id_neighbourer][7] + i;
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

  // std::vector<int> get_neighbour_idx(int i)
  // {
  // 	std::vector<int> neighs;neighs.reserve(8);
  // 	int tn = this->get_left_idx(i);
  // 	if(tn >= 0 )
  // 		neighs.emplace_back(tn);
  // 	tn = this->get_top_idx(i);
  // 	if(tn >= 0 )
  // 		neighs.emplace_back(tn);
  // 	tn = this->get_right_idx(i);
  // 	if(tn >= 0 )
  // 		neighs.emplace_back(tn);
  // 	tn = this->get_bottom_idx(i);
  // 	if(tn >= 0 )
  // 		neighs.emplace_back(tn);
  // 	tn = this->get_topleft_idx(i);
  // 	if(tn >= 0 )
  // 		neighs.emplace_back(tn);
  // 	tn = this->get_topright_idx(i);
  // 	if(tn >= 0 )
  // 		neighs.emplace_back(tn);
  // 	tn = this->get_bottomright_idx(i);
  // 	if(tn >= 0 )
  // 		neighs.emplace_back(tn);
  // 	tn = this->get_bottomleft_idx(i);
  // 	if(tn >= 0 )
  // 		neighs.emplace_back(tn);

  // 	return neighs;
  // }

  // Some of my algorithm are adapted from richdem and require iterating through
  // neighbours the way they do starting from left and going clockwise around
  // the node
  int get_neighbours_idx_richdemlike(int i, std::vector<int> &neighs) {
    int j = 0;
    neighs[j] = (this->get_left_idx(i));
    ++j;
    neighs[j] = (this->get_topleft_idx(i));
    ++j;
    neighs[j] = (this->get_top_idx(i));
    ++j;
    neighs[j] = (this->get_topright_idx(i));
    ++j;
    neighs[j] = (this->get_right_idx(i));
    ++j;
    neighs[j] = (this->get_bottomright_idx(i));
    ++j;
    neighs[j] = (this->get_bottom_idx(i));
    ++j;
    neighs[j] = (this->get_bottomleft_idx(i));
    ++j;
    return 8;
  }

  // Method to test whether a node can outlet flow OUT of the model
  // DEPRECATED
  // inline bool can_flow_out_there(int i)
  // {
  // 	if(this->boundary[i] == 3 || this->boundary[i] == 0 )
  // 		return true;
  // 	else
  // 	{
  // 		return false;
  // 	}
  // }

  // // method to check if a node can even accept flow
  // inline bool can_flow_even_go_there(int i)
  // {
  // 	if(this->boundary[i] < 0)
  // 		return false;
  // 	else
  // 	{
  // 		return true;
  // 	}
  // }

  // inline bool is_active(int i)
  // {
  // 	if(this->can_flow_out_there(i) && this->boundary[i] != 3)
  // 		return false;
  // 	else if(this->boundary[i] == 3)
  // 		return true;
  // 	else if(this->boundary[i] == 4)
  // 		return false;

  // 	if(this->can_flow_even_go_there(i))
  // 		return true;
  // 	return false;
  // }

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
      if (this->boundaries.no_data(i))
        continue;

      double slope_rad = 0;
      double aspect_rad = 0;
      double dzdx = 0;
      double dzdy = 0;

      double ij = topography[i];
      double ijp1 = topography[this->get_right_neighbour(i).node];
      double ip1j = topography[this->get_bottom_neighbour(i).node];
      double ip1jp1 = topography[this->get_bottomright_neighbour(i).node];
      double im1jm1 = topography[this->get_topleft_neighbour(i).node];
      double im1j = topography[this->get_top_neighbour(i).node];
      double im1jp1 = topography[this->get_topright_neighbour(i).node];
      double ijm1 = topography[this->get_left_neighbour(i).node];
      double ip1jm1 = topography[this->get_bottomleft_neighbour(i).node];

      if (ij > 0) {
        dzdx = ((ijp1 + 2 * ip1j + ip1jp1) - (im1jm1 + 2 * im1j + im1jp1)) /
               (8 * this->dx);
        dzdy = ((im1jp1 + 2 * ijp1 + ip1jp1) - (im1jm1 + 2 * ijm1 + ip1jm1)) /
               (8 * this->dy);
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

  // this function makes the connector generic and allow other connector to have
  // varying cellarea
  template <class i_t> T get_area_at_node(i_t i) { return this->cellarea; }

  // template<>
  // void fill_neighbour_matrices()

  template <class i_t> T get_dx_from_links_idx(i_t i) {
    if (i % 4 == 0)
      return this->dx;
    else if (i % 4 == 1)
      return this->dxy;
    else if (i % 4 == 2)
      return this->dy;
    else if (i % 4 == 3)
      return this->dxy;
    else
      return this->dx;
  }

  template <class i_t> T get_traverse_dx_from_links_idx(i_t i) {
    if (i % 4 == 0)
      return this->dy;
    else if (i % 4 == 1)
      return this->dxy;
    else if (i % 4 == 2)
      return this->dx;
    else if (i % 4 == 3)
      return this->dxy;
    else
      return this->dy;
  }

  // return the orthogonal node from a pair of node / link indices
  template <class i_t>
  std::pair<i_t, i_t> get_orthogonal_nodes(i_t node_from, i_t node_to) {

    if (node_to == this->get_top_idx(node_from))
      return {this->get_left_idx(node_from), this->get_right_idx(node_from)};
    if (node_to == this->get_topright_idx(node_from))
      return {this->get_top_idx(node_from), this->get_right_idx(node_from)};
    if (node_to == this->get_right_idx(node_from))
      return {this->get_top_idx(node_from), this->get_bottom_idx(node_from)};
    if (node_to == this->get_bottomright_idx(node_from))
      return {this->get_right_idx(node_from), this->get_bottom_idx(node_from)};
    if (node_to == this->get_bottom_idx(node_from))
      return {this->get_left_idx(node_from), this->get_right_idx(node_from)};
    if (node_to == this->get_bottomleft_idx(node_from))
      return {this->get_left_idx(node_from), this->get_bottom_idx(node_from)};
    if (node_to == this->get_left_idx(node_from))
      return {this->get_top_idx(node_from), this->get_bottom_idx(node_from)};
    if (node_to == this->get_topleft_idx(node_from))
      return {this->get_top_idx(node_from), this->get_left_idx(node_from)};
    else
      return {-1, -1};

    // else
    // {
    // 	std::cout << "Node was " << node_from << " to was " << node_to << "
    // possible neighbours to node_from:" << std::endl; 	std::cout <<
    // this->get_top_idx(node_from) << std::endl; 	std::cout <<
    // this->get_topright_idx(node_from) << std::endl; 	std::cout <<
    // this->get_right_idx(node_from) << std::endl; 	std::cout <<
    // this->get_bottomright_idx(node_from) << std::endl; 	std::cout <<
    // this->get_bottom_idx(node_from) << std::endl; 	std::cout <<
    // this->get_bottomleft_idx(node_from) << std::endl; 	std::cout <<
    // this->get_left_idx(node_from) << std::endl; 	std::cout <<
    // this->get_topleft_idx(node_from) << std::endl; 	std::cout << "Boundary
    // was " << this->boundary[node_from] << std::endl; 	throw
    // std::runtime_error("Fatal error in
    // DAGGER::D8connector::get_orthogonal_nodes -> from: " +
    // std::to_string(node_from) + " and node to " + std::to_string(node_to));
    // }
  }

  template <class i_t, class ii_t>
  std::pair<T, T> get_directed_dxdy_from_links_idx(i_t li, ii_t original_node,
                                                   ii_t n1, ii_t n2) {
    int C_C = (original_node == n1) ? 1 : -1;
    if (li % 4 == 0)
      return std::make_pair(C_C * this->dx, C_C * 0.);
    else if (li % 4 == 1)
      return std::make_pair(C_C * this->dx, C_C * this->dy);
    else if (li % 4 == 2)
      return std::make_pair(C_C * 0., C_C * this->dy);
    else if (li % 4 == 3)
      return std::make_pair(-1 * C_C * this->dx, C_C * this->dy);
    else
      return std::make_pair(C_C * this->dx, C_C * 0.);
  }

  template <class i_t, class ii_t>
  std::pair<T, T> get_dxdy_from_links_idx(i_t li) {
    if (li % 4 == 0)
      return std::make_pair(this->dx, 0.);
    else if (li % 4 == 1)
      return std::make_pair(this->dx, this->dy);
    else if (li % 4 == 2)
      return std::make_pair(0., this->dy);
    else if (li % 4 == 3)
      return std::make_pair(this->dx, this->dy);
    else
      return std::make_pair(this->dx, 0.);
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
      return this->dxy;
  }

  // EXPERIMENTAL!!!! POSSIBLE ISSUES WITH BC
  template <class i_t> i_t linkidx_from_nodes(i_t n1, i_t n2) {
    if (n1 > n2)
      std::swap(n1, n2);
    int delta = n2 - n1;
    // std::cout << n1 << "|" << n2 << "|" << delta << "||" << nnodes <<
    // std::endl;

    if (delta == 1)
      return n1 * 4;
    else if (delta == this->nx + 1)
      return n1 * 4 + 1;
    else if (delta == this->nx)
      return n1 * 4 + 2;
    else if (delta == this->nx - 1)
      return n1 * 4 + 3;
    else
      throw std::runtime_error(
          "Fatal error in DAGGER::D8connector::linkidx_from_nodes");
  }

  std::vector<int> get_ilinknodes_from_node_v1(int i) {
    std::vector<int> out;
    out.reserve(8);
    int tn = this->get_right_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4)
      out.emplace_back(tn);
    tn = this->get_bottomright_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4)
      out.emplace_back(tn);
    tn = this->get_bottom_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4)
      out.emplace_back(tn);
    tn = this->get_bottomleft_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4)
      out.emplace_back(tn);
    tn = this->get_left_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4)
      out.emplace_back(tn);
    tn = this->get_topleft_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4)
      out.emplace_back(tn);
    tn = this->get_top_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4)
      out.emplace_back(tn);
    tn = this->get_topright_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4)
      out.emplace_back(tn);
    return out;
  }

  std::vector<std::pair<int, bool>> get_ilinknodes_from_nodev2(int i) {
    std::vector<std::pair<int, bool>> out;
    out.reserve(8);

    if (!this->boundaries.normal_neighbouring_at_boundary(i) &&
        !this->boundaries.is_normal_node(i)) {
      int tn = this->get_right_idx_links(i);
      bool val = false;
      if (tn > 0 && tn < this->nnodes * 4)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_bottomright_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 4)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_bottom_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 4)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_bottomleft_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 4)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_left_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 4)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_topleft_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 4)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_top_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 4)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
      tn = this->get_topright_idx_links(i);
      if (tn > 0 && tn < this->nnodes * 4)
        val = true;
      out.emplace_back(std::make_pair(tn, val));
    } else
      out = {std::make_pair(this->get_right_idx_links(i), true),
             std::make_pair(this->get_bottomright_idx_links(i), true),
             std::make_pair(this->get_bottom_idx_links(i), true),
             std::make_pair(this->get_bottomleft_idx_links(i), true),
             std::make_pair(this->get_left_idx_links(i), true),
             std::make_pair(this->get_topleft_idx_links(i), true),
             std::make_pair(this->get_top_idx_links(i), true),
             std::make_pair(this->get_topright_idx_links(i), true)};

    return out;
  }

  void get_ilinknodes_from_nodev3(int i,
                                  std::vector<std::pair<int, bool>> &out) {

    int tn = this->get_right_idx_links(i);
    bool val = false;
    if (tn > 0 && tn < this->nnodes * 4)
      val = true;
    out[0].first = tn;
    out[0].second = val;
    tn = this->get_bottomright_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 4)
      val = true;
    out[1].first = tn;
    out[1].second = val;
    tn = this->get_bottom_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 4)
      val = true;
    out[2].first = tn;
    out[2].second = val;
    tn = this->get_bottomleft_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 4)
      val = true;
    out[3].first = tn;
    out[3].second = val;
    tn = this->get_left_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 4)
      val = true;
    out[4].first = tn;
    out[4].second = val;
    tn = this->get_topleft_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 4)
      val = true;
    out[5].first = tn;
    out[5].second = val;
    tn = this->get_top_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 4)
      val = true;
    out[6].first = tn;
    out[6].second = val;
    tn = this->get_topright_idx_links(i);
    val = false;
    if (tn > 0 && tn < this->nnodes * 4)
      val = true;
    out[7].first = tn;
    out[7].second = val;
  }

  int get_ilinknodes_from_node(int i, std::vector<int> &out) {
    int nn = 0;
    int tn = this->get_right_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4) {
      out[nn] = tn;
      ++nn;
    }
    tn = this->get_bottomright_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4) {
      out[nn] = tn;
      ++nn;
    }
    tn = this->get_bottom_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4) {
      out[nn] = tn;
      ++nn;
    }
    tn = this->get_bottomleft_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4) {
      out[nn] = tn;
      ++nn;
    }
    tn = this->get_left_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4) {
      out[nn] = tn;
      ++nn;
    }
    tn = this->get_topleft_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4) {
      out[nn] = tn;
      ++nn;
    }
    tn = this->get_top_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4) {
      out[nn] = tn;
      ++nn;
    }
    tn = this->get_topright_idx_links(i);
    if (tn > 0 && tn < this->nnodes * 4) {
      out[nn] = tn;
      ++nn;
    }
    return nn;
  }

  // returns a 1D bool array, false where no data
  std::vector<bool> get_mask_array() {
    std::vector<bool> mask(this->nnodes, true);
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->boundaries.no_data(i))
        mask[i] = false;
    }
    return mask;
  }

  template <class topo_t> void set_values_at_boundaries(topo_t &tarray, T val) {
    auto array = format_input(tarray);
    for (int i = 0; i < this->nnodes; ++i) {
      if (this->boundaries.can_out(i))
        array[i] = val;
    }
  }

  int get_nneighbours() { return 8; }
};

//  DEPRECATED AND MOVED OUTSIDE OF THE CONNECTOR
// 	template <class topo_t>
// 	std::vector<double> PriorityFlood(topo_t& ttopography)
// 	{
// 		std::random_device rd; // obtain a random number from hardware
// 	  std::mt19937 gen(rd()); // seed the generator
// 	  std::uniform_real_distribution<> distr(1e-7, 1e-6); // define the
// range

// 		auto tttopography = format_input(ttopography);
// 		std::vector<double> topography = to_vec(tttopography);

// 	  PQ_i_d open;
// 	  std::queue<PQ_helper<int,double> > pit;

// 	  uint64_t processed_cells = 0;
// 	  uint64_t pitc            = 0;
// 	  int      false_pit_cells = 0;
// 	  double PitTop = -9999;

// 	  std::vector<bool> closed(this->nnodes,false);

// 	  // RDLOG_PROGRESS<<"Adding cells to the priority queue...";
// 	  for(int i = 0; i<this->nnodes; ++i)
// 	  {
// 	  	if(this->can_flow_out_there(i))
// 	  	{
// 	  		closed[i] = true;
// 	  		open.emplace(i,topography[i]);
// 	  	}
// 	  }

// 	  // RDLOG_PROGRESS<<"Performing Priority-Flood+Epsilon...";
// 	  auto neighbours = this->get_empty_neighbour();

// 	  while(open.size()>0 || pit.size()>0)
// 	  {
// 	    PQ_helper<int,double> c;

// 	    if(pit.size()>0 && open.size()>0 &&
// open.top().score==pit.front().score)
// 	    {
// 	      c=open.top(); open.pop();
// 	      PitTop=-9999;
// 	    }
// 	    else if(pit.size()>0)
// 	    {
// 	      c=pit.front(); pit.pop();
// 	      if(PitTop==-9999)
// 	        PitTop=topography[c.node];
// 	    }
// 	    else
// 	    {
// 	      c=open.top(); open.pop();
// 	      PitTop=-9999;
// 	    }
// 	    processed_cells++;

// 	    int nn = this->get_neighbour_idx(c.node, neighbours);
// 	    for(int tnn = 0; tnn<nn; ++nn)
// 	    {
// 	    	int n = neighbours[tnn];
// 	      if(this->is_active(n) == false)
// 	      	continue;

// 	      if(closed[n])
// 	        continue;

// 	      closed[n]=true;
// 	      if(topography[n] <=
// std::nextafter(c.score,std::numeric_limits<double>::infinity()))
// 	      {
// 	        if(PitTop!=-9999 && PitTop<topography[n] &&
// std::nextafter(c.score,std::numeric_limits<double>::infinity())>=topography[n])
// 	          ++false_pit_cells;

// 	        ++pitc;
// 	        //
// topography[n]=std::nextafter(c.score,std::numeric_limits<double>::infinity());
// 	        topography[n] = c.score + distr(gen);
// 	        pit.emplace(n,topography[n]);
// 	      } else
// 	        open.emplace(n,topography[n]);
// 	    }
// 	  }

// 	  return topography;
// 	}

// 	template<class topo_t>
// 	void InitPriorityQue(
// 	  topo_t& topography,
// 	  std::vector<bool>& flag,
// 	  PQ_i_d& priorityQueue
// 	)
// 	{

// 	  std::queue<PQ_helper<int,double> > depressionQue;

// 	  // push border cells into the PQ
// 	  auto neighbours = this->get_empty_neighbour();
// 	  for(int i=0; i < this->nnodes; ++i)
// 	  {
// 	    if (flag[i]) continue;

// 	    if (this->can_flow_even_go_there(i) == false)
// 	    {
// 	      flag[i] = true;
// 	      int nn = this->get_neighbour_idx(i,neighbours);

// 	      for (int tnn =0; tnn < nn; ++tnn)
// 	      {

// 	      	int n = neighbours[tnn];

// 	        if (flag[n])
// 	          continue;
// 	        if (this->can_flow_even_go_there(n))
// 	        {
// 	          priorityQueue.emplace(n,topography[n]);
// 	          flag[n]=true;
// 	        }
// 	      }
// 	    }
// 	    else if(this->can_flow_out_there(i))
// 	    {
// 	      //on the DEM border
// 	      priorityQueue.emplace(i,topography[i]);
// 	      flag[i]=true;
// 	    }
// 	  }
// 	}

// 	template<class topo_t>
// 	void ProcessTraceQue(
// 	  topo_t& topography,
// 	  std::vector<bool>& flag,
// 	  std::queue<PQ_helper<int,double> >& traceQueue,
// 	  PQ_i_d& priorityQueue
// 	)
// 	{
// 	  std::queue<PQ_helper<int,double>  > potentialQueue;
// 	  int indexThreshold=2;  //index threshold, default to 2
// 	  auto neighbours = this->get_empty_neighbour();
// 	  while (!traceQueue.empty())
// 	  {
// 	    const auto node = traceQueue.front();
// 	    traceQueue.pop();

// 	    bool Mask[5][5]={{false},{false},{false},{false},{false}};

// 	    int nn = this->get_neighbour_idx(node.node, neighbours);

// 	    for (int tnn=0; tnn<nn; ++tnn)
// 	    {
// 	    	int n = neighbours[tnn];

// 	      if(flag[n])
// 	        continue;

// 	      if (topography[n] > node.score)
// 	      {
// 	        traceQueue.emplace(n, topography[n]);
// 	        flag[n]=true;
// 	      }
// 	      else
// 	      {
// 	        //initialize all masks as false
// 	        bool have_spill_path_or_lower_spill_outlet=false; //whether cell
// n has a spill path or a lower spill outlet than node if n is a depression
// cell 	        auto nneighs = this->get_neighbours_idx_richdemlike(n);
// int row,col; this->rowcol_from_node_id(node.node,row,col); int nrow,ncol;
// this->rowcol_from_node_id(n,nrow,ncol); 	        int incr=0; for(auto
// nn:nneighs)
// 	        {
// 	        	++incr;
// 	        	// Checking node validity
// 	        	if(nn<0 || nn>= this->nnodes)
// 	        		continue;

// 	        	int nnrow,nncol;
// 	        	this->rowcol_from_node_id(nn,nnrow,nncol);
// 	          if( (Mask[nnrow-row+2][nncol-col+2]) || (flag[nn] &&
// topography[nn] < node.score) )
// 	          {
// 	            Mask[nrow-row+2][ncol-col+2]=true;
// 	            have_spill_path_or_lower_spill_outlet=true;
// 	            break;
// 	          }
// 	        }

// 	        if(!have_spill_path_or_lower_spill_outlet)
// 	        {
// 	          if (incr<indexThreshold) potentialQueue.push(node);
// 	          else
// 	            priorityQueue.push(node);
// 	          break; // make sure node is not pushed twice into PQ
// 	        }
// 	      }
// 	    }//end of for loop
// 	  }

// 	  while (!potentialQueue.empty())
// 	  {
// 	    const auto node = potentialQueue.front();
// 	    potentialQueue.pop();

// 	    //first case
// 	    auto neigh = this->get_neighbours_idx_richdemlike(node.node);
// 	    int incr = 0;
// 	    for (auto n:neigh)
// 	    {
// 	    	++incr;

// 	      if(flag[n])
// 	        continue;

// 	      priorityQueue.push(node);
// 	      break;
// 	    }
// 	  }
// 	}

// 	template<class topo_t>
// 	void ProcessPit(
// 	  topo_t& topography,
// 	  std::vector<bool>& flag,
// 	  std::queue<PQ_helper<int,double> >& depressionQue,
// 	  std::queue<PQ_helper<int,double> >& traceQueue
// 	)
// 	{
// 		auto neighbours = this->get_empty_neighbour();
// 	  while (!depressionQue.empty())
// 	  {
// 	    auto node = depressionQue.front();
// 	    depressionQue.pop();
// 	    int nn = this->get_neighbour_idx(node.node, neighbours);
// 	    for (int tnn = 0; tnn<nn;++tnn)
// 	    {
// 	    	int n = neighbours[tnn];
// 	      if (flag[n])
// 	        continue;

// 	      const auto iSpill = topography[n];
// 	      if (iSpill > node.score)
// 	      { // Dlope cell
// 	        flag[n]=true;
// 	        traceQueue.emplace(n,iSpill);
// 	        continue;
// 	      }
// 	      // Depression cell
// 	      flag[n] = true;
// 	      topography[n] = node.score + 1e-6 * this->randu->get();
// 	      depressionQue.emplace(n,topography[n]);
// 	    }
// 	  }
// 	}

// 	template<class topo_t>
// 	std::vector<double> PriorityFlood_Wei2018(topo_t &ttopography)
// 	{
// 	  std::queue<PQ_helper<int,double> > traceQueue;
// 	  std::queue<PQ_helper<int,double> > depressionQue;

// 	  std::vector<bool> flag(this->nnodes,false);
// 	  std::vector<double> topography = to_vec(ttopography);

// 	  PQ_i_d priorityQueue;

// 	  this->InitPriorityQue(topography,flag,priorityQueue);

// 	  auto neighbours = this->get_empty_neighbour();

// 	  while (!priorityQueue.empty())
// 	  {
// 	    const auto tmpNode = priorityQueue.top();
// 	    priorityQueue.pop();

// 	    int nn = this->get_neighbour_idx(tmpNode.node, neighbours);

// 	    for(int tnn=0; tnn< nn; ++tnn)
// 	    {
// 	    	int n = neighbours[tnn];
// 	      if(flag[n])
// 	        continue;

// 	      auto iSpill = topography[n];

// 	      if (iSpill <= tmpNode.score)
// 	      {
// 	        //depression cell
// 	        topography[n] = tmpNode.score + 1e-3 + 1e-6 *
// this->randu->get(); 	        flag[n] = true;
// depressionQue.emplace(n,topography[n]);
// 	        this->ProcessPit(topography,flag,depressionQue,traceQueue);
// 	      }
// 	      else
// 	      {
// 	        //slope cell
// 	        flag[n] = true;
// 	        traceQueue.emplace(n,iSpill);
// 	      }
// 	      ProcessTraceQue(topography,flag,traceQueue,priorityQueue);
// 	    }
// 	  }

// 	  return topography;
// 	}

// end of namespace
}; // namespace DAGGER

#endif
