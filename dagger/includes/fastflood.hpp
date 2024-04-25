#ifndef FASTFLOOD_HPP
#define FASTFLOOD_HPP

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
#include "fastflood_recorder.hpp"
#include "wrap_helper.hpp"

namespace DAGGER {

float_t GRAVITY = 9.81, FIVETHIRD = 5. / 3., RHO = 1000;

template <class float_t, class Graph_t, class Connector_t, class _topo_t>
class fastflood {
public:
  Graph_t *graph;
  Connector_t *connector;
  std::vector<double> topography, Qbase;

  std::vector<float_t> hw, weights, Qwin, Qwout, Qsin, Qsout, qs, Qlink, lastdt,
      lastdl, lastS;
  std::vector<int> linknodesD4;

  std::vector<int> Qs_entry_point;
  std::vector<float_t> Qs_entry;

  fastflood_recorder<float_t> rec;

  float_t mannings = 0.033, topological_number = 4. / 8, froude_number = 0.8,
          alpha = 0.5;

  float_t parting_coeff = 1, maxdt = 0, sensibility_to_flowdepth = 2,
          depth_threshold = 2, stochaslope = -1, boundary_slope = 1e-3;

  bool hflow = false, depth_limiter = true;

  DAGGER::easyRand randu;

  bool REPORT_N_PIXEL_FILLED = false;

  bool spatially_varying_dt = false;
  std::vector<float_t> dts = {5e-3}, saved_dts;
  float_t min_dt = 1e-4, max_dt = 10, k_dt = 1;
  bool Afdt = false;

  fastflood(){};
  template <class topo_t>
  fastflood(Graph_t &graph, Connector_t &connector, topo_t &topo,
            topo_t &Qbase) {
    // graph.compute_graph("none",topo, connector,true, true);
    this->graph = &graph;
    this->connector = &connector;
    this->rec.nnodes = connector.nnodes;
    this->rec.nlinks = connector.nnodes * 4;
    auto temp = DAGGER::format_input<topo_t>(topo);
    this->topography = DAGGER::to_vec(temp);
    temp = DAGGER::format_input<topo_t>(Qbase);
    this->Qbase = DAGGER::to_vec(temp);
    this->hw = std::vector<float_t>(this->connector->nnodes, 0.);
    this->Qwin = std::vector<float_t>(this->connector->nnodes, 0.);
    this->Qwout = std::vector<float_t>(this->connector->nnodes, 0.);
    this->weights = std::vector<float_t>(this->connector->links.size(), 0.);
  }

  std::vector<float_t> get_surface() {
    std::vector<float_t> surf(this->graph->nnodes, 0.);
    for (int i = 0; i < this->graph->nnodes; ++i)
      surf[i] = this->topography[i] + this->hw[i];
    return surf;
  }

  template <class in_t> void spatial_dt(in_t &DTS) {
    this->spatially_varying_dt = true;
    auto tdts = DAGGER::format_input(DTS);
    this->dts = DAGGER::to_vec(tdts);
  }
  void set_dt(float_t dt) {
    this->spatially_varying_dt = false;
    this->dts = {dt};
  }
  void set_min_dt(float_t dt) { this->min_dt = dt; }
  void set_max_dt(float_t dt) { this->max_dt = dt; }
  void set_k_dt(float_t dt) { this->k_dt = dt; }
  void set_boundary_slope(float_t tS) { this->boundary_slope = tS; }
  void config_Afdt(float_t min_dt, float_t max_dt, float_t k_dt) {
    this->enable_Afdt();
    this->min_dt = min_dt;
    this->max_dt = max_dt;
    this->k_dt = k_dt;
    this->saved_dts = std::vector<float_t>(this->graph->nnodes, 0.);
  }
  void enable_Afdt() { this->Afdt = true; }
  void disable_Afdt() { this->Afdt = false; }

  void set_stochaslope(float_t stosto) {
    this->stochaslope = stosto;
    this->connector->set_stochaticiy_for_SFD(stosto);
  }

  template <class tt0> float_t dt(tt0 i) {
    // float tdt = (this->spatially_varying_dt)?this->dts[i]:this->dts[0];
    // if(this->Afdt)
    // {
    // 	if(this->Qwin[i]>0)
    // 	{
    // 		float_t ratio_Q = this->Qwout[i] / this->Qwin[i];
    // 		tdt = this->k_dt *  tdt / ratio_Q;
    // 		if(tdt < this->min_dt)
    // 			tdt = this->min_dt;
    // 		else if (this->max_dt < tdt)
    // 			tdt = this->max_dt;
    // 	}
    // 	this->saved_dts[i] = tdt;
    // 	return tdt;
    // }
    // else if(this->spatially_varying_dt)
    // 	return this->dts[i];
    // else
    return this->dts[0];
    // return tdt;
  }

  void set_parting_coeff(float_t tpart) { this->parting_coeff = tpart; }
  void set_out_boundaries_to_permissive() {
    this->connector->set_out_boundaries_to_permissive();
  }
  // void set_edges_to_0()
  // {
  // 	for(int i = 0; i<this->graph->nnodes; ++i)
  // 	{
  // 		if(this->connector->can_flow_out_there(i))
  // 			this->hw[i] = 0;
  // 	}
  // }
  void enable_hflow() { this->hflow = true; }
  void disable_hflow() { this->hflow = false; }
  void set_sensibility_to_flowdepth(float_t exp) {
    this->sensibility_to_flowdepth = exp;
  }
  float_t get_sensibility_to_flowdepth() {
    return this->sensibility_to_flowdepth;
  }

  // void run_SFD(float_t dt)

  // void run_dynamic(float_t dt)
  // {
  // 	std::vector<float_t> surf(this->connector.nnodes,0);
  // 	for(int i=0; i< this->connector.nnodes; ++i)
  // 		surf[i] = this->topography[i] + this->hw[i];

  // 	this->graph.update_recs(surf, *this->connector);
  // 	auto gradient = this->graph->get_links_gradient( *this->connector,surf);
  // 	this->graph->_get_link_weights_exp(this->weights, gradient, 0.5);
  // 	auto gmax = this->graph->_get_max_val_link_array(gradient);
  // 	for(int i=0; i< this->connector->nnodes; ++i)
  // 	{
  // 		this->hw[i] += this->Qwin[i] + this->Qbase[i]
  // 	}

  void fill_up(float_t tdt) {
    std::cout << "fill_up is DEPRECATED" << std::endl;
    // std::vector<float_t> surf(this->connector->nnodes,0);

    // // float_t sumgrad_sqrt = 0, sumHW = 0;

    // for(int i=0; i< this->connector->nnodes; ++i)
    // 	surf[i] = this->topography[i] + this->hw[i];

    // // std::vector<float_t> surfpp(surf);
    // std::vector<float_t> surfpp(surf);
    // // std::cout << "DEBUGDEBUGA11" << std::endl;
    // (*this->graph)._compute_graph(surfpp, *(this->connector), true, true);
    // // py::array temp  = this->graph->compute_graph<Connector_t,
    // std::vector<double> , py::array>(depression_solver, surfpp,
    // *(this->connector), true, true);
    // // std::cout << "DEBUGDEBUGA2" << std::endl;
    // // float_t totaddedbystuff = 0;
    // int n_pixel_filled = 0;
    // for(int i=this->connector->nnodes-1; i>=0; --i)
    // {

    // 	if(surfpp[i] > surf[i])
    // 	{
    // 		this->hw[i] += surfpp[i] - surf[i];
    // 		n_pixel_filled++;
    // 		// totaddedbystuff += surfpp[i] - surf[i];
    // 	}

    // }

    // if(this->REPORT_N_PIXEL_FILLED)
    // 	std:: cout << n_pixel_filled << "where filled" << std::endl;
  }

  // Function precomputing the graph and returning a filled topography
  // only_SS determines if the graph only needs Steepest Descent (true) or
  // multiple flow (false)
  void graph_automator(bool only_SS) {
    std::vector<float_t> surf(this->connector->nnodes, 0);
    for (int i = 0; i < this->connector->nnodes; ++i)
      surf[i] = this->topography[i] + this->hw[i];
    std::vector<float_t> surfpp(surf);
    (*this->graph)._compute_graph(surfpp, only_SS, true);

    for (int i = 0; i < this->connector->nnodes; ++i) {

      if (surfpp[i] > surf[i]) {
        this->hw[i] += surfpp[i] - surf[i];
      }
    }
  }

  void run_SFD() {
    // std::cout << "DEBUGDEBUGA" << std::endl;
    std::vector<float_t> surf(this->connector->nnodes, 0);

    // float_t sumgrad_sqrt = 0, sumHW = 0;

    for (int i = 0; i < this->connector->nnodes; ++i)
      surf[i] = this->topography[i] + this->hw[i];

    // std::cout << "DEBUGDEBUGA1" << std::endl;
    double totqin = 0, totqinout = 0;
    ;

    // std::vector<float_t> surfpp(surf);
    std::vector<float_t> surfpp(surf);
    // std::cout << "DEBUGDEBUGA11" << std::endl;
    // this->graph->depression_resolver = DAGGER::DEPRES::priority_flood_2014 ;
    (*this->graph)._compute_graph(surfpp, true, true);
    // py::array temp  = this->graph->compute_graph<Connector_t,
    // std::vector<double> , py::array>(depression_solver, surfpp,
    // *(this->connector), true, true); std::cout << "DEBUGDEBUGA2" <<
    // std::endl; float_t totaddedbystuff = 0;

    int n_pixel_filled = 0;
    for (int i = 0; i < this->connector->nnodes; ++i) {

      if (surfpp[i] > surf[i]) {
        this->hw[i] += surfpp[i] - surf[i];
        n_pixel_filled++;
        // totaddedbystuff += surfpp[i] - surf[i];
      }
    }
    if (this->REPORT_N_PIXEL_FILLED)
      std::cout << n_pixel_filled << "where filled" << std::endl;

    // debugging process of checking the SD receivers
    // this->check_SD_val();

    this->Qwin = this->graph->_accumulate_variable_downstream_SFD(this->Qbase);
    std::vector<float_t> gradients = this->connector->_get_SFD_gradient(surfpp);

    for (auto &v : this->Qwout)
      v = 0;

    for (auto v : this->Qbase)
      totqin += v;

    // std::cout << "DEBUGDEBUGA5" << std::endl;
    int nneg = 0;
    double exp1 = 1. / 2.;
    double exp2 = 5. / 3.;
    std::vector<bool> done(this->connector->nnodes, false);
    for (int ti = this->connector->nnodes - 1; ti >= 0; --ti) {
      int i = int(this->graph->Sstack[ti]);
      int rec = this->connector->_Sreceivers[i];

      if (gradients[i] < 0) {
        ++nneg;
        gradients[i] = 0;
      }
      // sumgrad_sqrt += std::sqrt(gradients[i]);

      if (this->connector->flow_out_or_pit(i) == false) {
        float_t thw =
            (this->hflow) ? this->get_hflow_from_nodes(i, rec) : this->hw[i];

        this->Qwout[i] = this->connector->Sdistance2receivers[i] /
                         this->mannings * std::pow(thw, exp2) *
                         std::pow(gradients[i], exp1);

        this->hw[i] += dt(i) * (this->Qwin[i] - this->Qwout[i]) /
                       this->connector->get_area_at_node(i);

        if (this->hw[i] < 0)
          this->hw[i] = 0;
      } else {
        totqinout += this->Qwin[i];
      }
    }
    // std::cout << "In: " << totqin << " out: " << totqinout << " Balance::" <<
    // totqin - totqinout << std::endl; std::cout << "man::" << 1/mannings << "
    // hwpow::" << sumHW << " grad::" << sumgrad_sqrt << " sumgrad2 " <<
    // sumgrad2 << " nmneg = " << nneg<< std::endl; std::cout << "DEBUGDEBUGA5"
    // << std::endl;
    this->manage_hydrology_at_BC();
  }

  // This function runs one iteration of the static multiple flow solver for the
  // SWE
  void run_MFD() {

    std::cout << "run_MFD is temporarily unavailable " << std::endl;
    // this->rec.init_water();

    // // timer, to ignore
    // DAGGER::ocarina link;
    // link.tik();

    // // The surface (Hw + topo)
    // std::vector<float_t> surf(this->connector->nnodes,0);
    // for(int i=0; i< this->connector->nnodes; ++i)
    // 	surf[i] = this->topography[i] + this->hw[i];
    // // link.tok("surfcalc");

    // // Filling the surface (filling local minimas)
    // link.tik();
    // std::vector<float_t> surfpp(surf);
    // (*this->graph)._compute_graph(surfpp, *(this->connector), false, true);
    // // link.tok("graph");

    // // applying the filling to the water height
    // link.tik();
    // int n_pixel_filled = 0;
    // for(int i=0; i<this->connector->nnodes; ++i)
    // {
    // 	if(surfpp[i] > surf[i])
    // 	{
    // 		this->hw[i] += surfpp[i] - surf[i];
    // 		// Tracking the number of pixel filled as one indicator for
    // stability 		n_pixel_filled++;
    // 	}
    // }
    // if(this->REPORT_N_PIXEL_FILLED)
    // 	std:: cout << n_pixel_filled << "where filled" << std::endl;

    // // reinitialising the Qw in and out
    // for(auto& v: this->Qwin)
    // 	v=0;
    // for(auto& v: this->Qwout)
    // 	v=0;

    // // Getting things ready for the main loop
    // //# allocating the vector fetching the neighbours
    // auto reclink = this->connector->get_empty_neighbour();
    // //# allocating the local vector storing sqrt of the gradient
    // std::vector<float_t> sqrtgrads(reclink.size(),0);
    // //# iterating through all the nodes
    // for(int i = this->graph->nnodes-1; i>=0; --i)
    // {
    // 	//## Current node
    // 	int node = this->graph->stack[i];

    // 	//## is the node active?
    // 	if(!this->graph->flow_out_model(node, *this->connector) == false)
    // 		continue;

    // 	//## getting the receiver links (all where heaviside is 1)
    // 	int nn = this->graph->get_receivers_idx_links(node, *this->connector,
    // reclink);

    // 	//## Accumulating local discharge
    // 	this->Qwin[node] += this->Qbase[node];

    // 	//## (Rare) no receivers?? then I skip
    // 	if(nn == 0)
    // 	{
    // 		continue;
    // 	}

    // 	//## bunch of temporary variables to accumulate factors, reuse
    // characteristics and determine max slope 	float_t maxgrad_sqrt = 0,
    // sumgrad_sqrt = 0, facgrad = 0, cellarea =
    // this->connector->get_area_at_node(i);

    // 	//## iterating through neighbours
    // 	for(int ttl = 0; ttl<nn; ++ttl)
    // 	{
    // 		//### neighbouring link id
    // 		int li = reclink[ttl];

    // 		//### Not valid 4 some reasons? skip
    // 		if(this->graph->is_link_valid(li) == false)
    // 			continue;

    // 		//### local dx
    // 		float_t dx = this->connector->get_dx_from_links_idx(li);

    // 		//### receiving node
    // 		int to = this->graph->get_to_links(li);

    // 		//### local gradient, set to a numrically unsignificant minimum
    // 		float_t tgrad = (surfpp[node] - surfpp[to]) / dx;
    // 		tgrad = std::max(tgrad,1e-8);
    // 		//### its squareroot
    // 		float_t tsqrt = std::sqrt(tgrad);
    // 		//### Storing local slope to the power needed
    // 		sqrtgrads[ttl] = (parting_coeff == 0.5)?tsqrt:std::pow(tgrad,
    // this->parting_coeff); 		sqrtgrads[ttl] *= 1 - this->randu.get()
    // * this->stochaslope;
    // 		//### SUmming the slope to the power needed
    // 		sumgrad_sqrt += sqrtgrads[ttl];
    // 		//### accumulating the equation factor
    // 		if(this->hflow)
    // 		{

    // 			facgrad +=
    // std::pow(this->get_hflow_at_link(reclink[ttl]), FIVETHIRD)
    // * tgrad/dx;
    // 		}
    // 		else
    // 			facgrad += tgrad/dx;

    // 		//### Saving the max sqrt slope
    // 		if(tsqrt > maxgrad_sqrt)
    // 			maxgrad_sqrt = tsqrt;
    // 	}

    // 	//## Parting the Qin to the receivers
    // 	//## Going through neighbours
    // 	for(int ttl = 0; ttl<nn; ++ttl)
    // 	{
    // 		///### Getting the link id
    // 		int li = reclink[ttl];
    // 		if(this->graph->is_link_valid(li) == false)
    // 			continue;
    // 		///### the other node
    // 		int to = this->graph->get_to_links(li);
    // 		///### partitioning proportional to the slope to the power
    // this->parting_coeff 		this->Qwin[to] +=
    // sqrtgrads[ttl]/sumgrad_sqrt
    // * this->Qwin[node];
    // 	}

    // 	//## Finally calculating Qout
    // 	if(this->hflow)
    // 		this->Qwout[node] = this->topological_number *
    // facgrad/maxgrad_sqrt
    // * cellarea/this->mannings;// * (1/(4 * this->connector->dxy )+ 1/
    // (2*this->connector->dx) + 1/(2 * this->connector->dy)); 	else
    // 		this->Qwout[node] = this->topological_number *
    // facgrad/maxgrad_sqrt
    // * std::pow(this->hw[node],FIVETHIRD) * cellarea/this->mannings;// * (1/(4
    // * this->connector->dxy )+ 1/ (2*this->connector->dx) + 1/(2 *
    // this->connector->dy));

    // 	//## Divergence of Q to apply changes to water height
    // 	float_t tdhw = (this->Qwin[node] - this->Qwout[node]) * dt(node) /
    // this->connector->get_area_at_node(node); 	if(this->rec.dhw2record)
    // this->rec.dhw[node] = tdhw; 	this->hw[node] += tdhw;

    // 	//## water height cannot be 0
    // 	if(this->hw[node] < 0) this->hw[node] = 0;

    // }// end of the loop

    // this->manage_hydrology_at_BC();

    // Done
  }

  // This function runs one iteration of the static multiple flow solver for the
  // SWE
  void run_MFD_dynamic() {
    // timer, to ignore
    DAGGER::ocarina link;
    link.tik();

    // The surface (Hw + topo)
    std::vector<float_t> surf(this->connector->nnodes, 0);
    for (int i = 0; i < this->connector->nnodes; ++i)
      surf[i] = this->topography[i] + this->hw[i];
    // link.tok("surfcalc");

    // Filling the surface (filling local minimas)
    link.tik();
    std::vector<float_t> surfpp(surf);
    (*this->graph)._compute_graph(surfpp, false, true);
    // link.tok("graph");

    // NO FILLING IN DYNAMIC
    // // applying the filling to the water height
    // link.tik();
    // int n_pixel_filled = 0;
    // for(int i=0; i<this->connector->nnodes; ++i)
    // {
    // 	if(surfpp[i] > surf[i])
    // 	{
    // 		this->hw[i] += surfpp[i] - surf[i];
    // 		// Tracking the number of pixel filled as one indicator for
    // stability 		n_pixel_filled++;
    // 	}
    // }
    // if(this->REPORT_N_PIXEL_FILLED)
    // 	std:: cout << n_pixel_filled << "where filled" << std::endl;

    // reinitialising the Qw in and out
    for (auto &v : this->Qwin)
      v = 0;
    for (auto &v : this->Qwout)
      v = 0;

    // Getting things ready for the main loop
    // # allocating the vector fetching the neighbours
    auto reclink = this->connector->get_empty_neighbour();
    // # allocating the local vector storing sqrt of the gradient
    std::vector<float_t> sqrtgrads(reclink.size(), 0);
    // # iterating through all the nodes
    for (int i = this->graph->nnodes - 1; i >= 0; --i) {
      // ## Current node
      int node = this->graph->stack[i];

      // ## is the node active?
      if (this->connector->flow_out_model(node))
        continue;

      // ## getting the receiver links (all where heaviside is 1)
      int nn = this->connector->get_receivers_idx_links(node, reclink);

      // ## Accumulating local discharge
      this->Qwin[node] += this->Qbase[node];

      // ## (Rare) no receivers?? then I skip
      if (nn == 0) {
        this->hw[i] += this->Qwin[i] / this->connector->get_area_at_node(i);
        continue;
      }

      // ## bunch of temporary variables to accumulate factors, reuse
      //  characteristics and determine max slope
      float_t maxgrad_sqrt = 0, sumgrad_sqrt = 0, facgrad = 0,
              cellarea = this->connector->get_area_at_node(i);

      // ## iterating through neighbours
      for (int ttl = 0; ttl < nn; ++ttl) {
        // ### neighbouring link id
        int li = reclink[ttl];

        // ### Not valid 4 some reasons? skip
        if (this->connector->is_link_valid(li) == false)
          continue;

        // ### local dx
        float_t dx = this->connector->get_dx_from_links_idx(li);

        // ### receiving node
        int to = this->connector->get_to_links(li);

        // ### local gradient, set to a numrically unsignificant minimum
        float_t tgrad = (surfpp[node] - surfpp[to]) / dx;
        tgrad = std::max(tgrad, 1e-8);
        // ### its squareroot
        float_t tsqrt = std::sqrt(tgrad);
        // ### Storing local slope to the power needed
        sqrtgrads[ttl] = (parting_coeff == 0.5)
                             ? tsqrt
                             : std::pow(tgrad, this->parting_coeff);
        // ### SUmming the slope to the power needed
        sumgrad_sqrt += sqrtgrads[ttl];
        // ### accumulating the equation factor
        facgrad += tgrad / dx;
        // ### Saving the max sqrt slope
        if (tsqrt > maxgrad_sqrt)
          maxgrad_sqrt = tsqrt;
      }

      // ## Parting the Qin to the receivers
      // ## Going through neighbours
      for (int ttl = 0; ttl < nn; ++ttl) {
        /// ### Getting the link id
        int li = reclink[ttl];
        if (this->connector->is_link_valid(li) == false)
          continue;
        /// ### the other node
        int to = this->connector->get_to_links(li);
        /// ### partitioning proportional to the slope to the power
        ///  this->parting_coeff
        this->Qwin[to] += sqrtgrads[ttl] / sumgrad_sqrt * this->Qwout[node];
      }

      // ## Calculating Qout
      this->Qwout[node] = this->topological_number * facgrad / maxgrad_sqrt *
                          std::pow(this->hw[node], FIVETHIRD) * cellarea /
                          this->mannings; // * (1/(4 * this->connector->dxy )+
                                          // 1/ (2*this->connector->dx) + 1/(2 *
                                          // this->connector->dy));

      // ## Divergence of Q to apply changes to water height
      this->hw[node] += (this->Qwin[node] - this->Qwout[node]) * dt(node) /
                        this->connector->get_area_at_node(node);

      // ## water height cannot be 0
      if (this->hw[node] < 0)
        this->hw[node] = 0;

    } // end of the loop

    this->manage_hydrology_at_BC();

    // Done
  }

  void manage_hydrology_at_BC() {
    return;
    auto neighbours = this->connector->get_empty_neighbour();
    for (int i = 0; i < this->graph->nnodes; ++i) {
      if (!this->connector->flow_out_model(i) == false) {
        int nn = this->connector->get_neighbour_idx(i, neighbours);
        float_t tsurf = std::numeric_limits<float_t>::max();
        int tnode = -1;
        for (int j = 0; j < nn; ++j) {
          if (this->get_surface_at_node(neighbours[j]) < tsurf) {
            tsurf = this->get_surface_at_node(neighbours[j]);
            tnode = neighbours[j];
          }
        }
        if (tnode >= 0) {
          // std::cout << "recasting " << i << std::endl;
          this->topography[i] = this->topography[tnode];
          this->hw[i] = this->hw[tnode];
        }
      }
    }
  }

  void run_MFD_exp(int n_run_by_graph_cycle = 1) {
    DAGGER::ocarina link;
    link.tik();

    std::vector<float_t> surf(this->connector->nnodes, 0);

    for (int i = 0; i < this->connector->nnodes; ++i)
      surf[i] = this->topography[i] + this->hw[i];

    // link.tok("surfcalc");

    link.tik();
    std::vector<float_t> surfpp(surf);
    // std::cout << "DEBUGDEBUGA11" << std::endl;
    (*this->graph)._compute_graph(surfpp, false, true);

    link.tik();
    int n_pixel_filled = 0;
    for (int i = 0; i < this->connector->nnodes; ++i) {

      if (surfpp[i] > surf[i]) {
        this->hw[i] += surfpp[i] - surf[i];
        n_pixel_filled++;
        // totaddedbystuff += surfpp[i] - surf[i];
      }
    }
    if (this->REPORT_N_PIXEL_FILLED)
      std::cout << n_pixel_filled << "where filled" << std::endl;

    for (auto &v : this->Qwin)
      v = 0;

    for (auto &v : this->Qwout)
      v = 0;

    for (int n_opti = 0; n_opti < n_run_by_graph_cycle; ++n_opti) {

      // std::cout << "DEBUGDEBUGA4" << std::endl;
      auto reclink = this->connector->get_empty_neighbour();
      std::vector<float_t> sqrtgrads(reclink.size(), 0);
      for (int i = this->graph->nnodes - 1; i >= 0; --i) {
        // std::cout << "DEBUGDEBUGA41" << std::endl;
        int node = this->graph->stack[i];
        // std::cout << "DEBUGDEBUGA42" << std::endl;

        if (this->connector->flow_out_model(node))
          continue;

        int nn = this->connector->get_receivers_idx_links(node, reclink);
        this->Qwin[node] += this->Qbase[node];

        if (nn == 0) {
          continue;
        }

        float_t maxgrad_sqrt = 0, sumgrad_sqrt = 0, facgrad = 0;
        for (int ttl = 0; ttl < nn; ++ttl) {
          int li = reclink[ttl];
          if (this->connector->is_link_valid(li) == false)
            continue;

          float_t dx = this->connector->get_dx_from_links_idx(li);

          int to = this->connector->get_to_links(li);
          float_t tgrad = (surfpp[node] - surfpp[to]) / dx;
          tgrad = std::max(tgrad, 1e-8);
          float_t tsqrt = std::sqrt(tgrad);
          sqrtgrads[ttl] = (parting_coeff == 0.5)
                               ? tsqrt
                               : std::pow(tgrad, this->parting_coeff);
          sumgrad_sqrt += sqrtgrads[ttl];
          facgrad += dx * tgrad;
          if (tsqrt > maxgrad_sqrt)
            maxgrad_sqrt = tsqrt;
        }

        if (sumgrad_sqrt > 0) {
          for (int ttl = 0; ttl < nn; ++ttl) {
            int li = reclink[ttl];
            if (this->connector->is_link_valid(li) == false)
              continue;

            int to = this->connector->get_to_links(li);

            this->Qwin[to] += sqrtgrads[ttl] / sumgrad_sqrt * this->Qwin[node];
          }

          this->Qwout[node] = facgrad / maxgrad_sqrt *
                              std::pow(this->hw[node], 5. / 3) * 0.5 /
                              this->mannings;
        } else {
          throw std::runtime_error("Happens?");
        }

        this->hw[node] += (this->Qwin[node] - this->Qwout[node]) * dt(node) /
                          this->connector->get_area_at_node(node);
        this->Qwin[node] = 0;
        this->Qwout[node] = 0;

        if (this->hw[node] < 0)
          this->hw[node] = 0;
      }
    }
  }

  // Function getting the full effective width based on Thomas Bernard's work
  template <class out_t>
  out_t get_w_eff(float_t prec, bool recompute_graph = true) {
    std::vector<float_t> w_eff(this->graph->nnodes, 0.);
    this->_get_w_eff(w_eff, prec, recompute_graph);
    return DAGGER::format_output<std::vector<float_t>, out_t>(w_eff);
  }

  // internal function calculatuing the effectinve width
  void _get_w_eff(std::vector<float_t> &w_eff, float_t prec,
                  bool recompute_graph) {

    // I need first the effective drainage area
    std::vector<float_t> a_eff(this->graph->nnodes, 0.);
    this->_get_a_eff(a_eff, prec, recompute_graph);
    auto A = this->graph->_accumulate_constant_downstream_SFD(
        this->connector->get_area_at_node(0));

    for (int i = 0; i < this->connector->nnodes; ++i) {
      if (!this->connector->flow_out_model(i)) {
        w_eff[i] = A[i] / a_eff[i];
      }
    }
  }

  template <class out_t> out_t get_a_eff(float_t prec, bool recompute_graph) {
    std::vector<float_t> a_eff(this->graph->nnodes, 0.);
    this->_get_a_eff(a_eff, prec, recompute_graph);
    return DAGGER::format_output<std::vector<float_t>, out_t>(a_eff);
  }

  void _get_a_eff(std::vector<float_t> &a_eff, float_t prec,
                  bool recompute_graph) {

    if (recompute_graph) {
      std::vector<float_t> surf(this->connector->nnodes, 0);

      for (int i = 0; i < this->connector->nnodes; ++i)
        surf[i] = this->topography[i] + this->hw[i];
      std::vector<float_t> surfpp(surf);
      (*this->graph)._compute_graph(surfpp, false, true);
    }

    std::vector<float_t> S_h(this->connector->nnodes, 0);
    this->_get_hydraulic_slope_D8(S_h, false);

    for (int i = 0; i < this->graph->nnodes; ++i) {
      if (!this->connector->flow_out_model(i)) {
        a_eff[i] = (std::pow(this->hw[i], FIVETHIRD) * std::sqrt(S_h[i])) /
                   (this->mannings * prec);
      }
    }
  }

  template <class out_t> out_t get_hydraulic_slope_D8(bool recompute_graph) {
    std::vector<float_t> S_h(this->graph->nnodes, 0.);
    this->_get_hydraulic_slope_D8(S_h, recompute_graph);
    return DAGGER::format_output<std::vector<float_t>, out_t>(S_h);
  }

  void _get_hydraulic_slope_D8(std::vector<float_t> &S_h,
                               bool recompute_graph) {

    if (recompute_graph) {
      std::vector<float_t> surf(this->connector->nnodes, 0);

      for (int i = 0; i < this->connector->nnodes; ++i)
        surf[i] = this->topography[i] + this->hw[i];
      std::vector<float_t> surfpp(surf);
      (*this->graph)._compute_graph(surfpp, false, true);
    }

    for (int i = 0; i < this->graph->nnodes; ++i) {
      if (!this->connector->flow_out_model(i)) {
        int rec = this->connector->_Sreceivers[i];
        float_t dx = this->connector->Sdistance2receivers[i];
        S_h[i] = std::max((this->topography[i] + this->hw[i] -
                           this->topography[rec] - this->hw[rec]) /
                              dx,
                          0.);
      }
    }
  }

  // Implementing Caesar-Lisflood flooding algorithm in DAGGER's logic
  // For the physics see paper from Bates et al., 2010 and the implementation is
  // adapted and modified from HAIL-CAESAR implementation by Declan Valters
  float_t caesar_lisflood() {

    // FIrst checking if the Qlinks has been initialised
    if (this->Qlink.size() == 0) {
      this->Qlink = std::vector<float_t>(this->connector->links.size(), 0.);
    }

    // Initializing the dhw
    std::vector<float_t> dhw(this->graph->nnodes, 0.);

    std::vector<float_t> surf(this->connector->nnodes, 0);
    for (int i = 0; i < this->connector->nnodes; ++i)
      surf[i] = this->topography[i] + this->hw[i];

    // calculationg the maximum timestep respecting the CFD
    float_t maxhw =
        std::max(*std::max_element(this->hw.begin(), this->hw.end()), 1e-3);

    float_t tdt =
        this->connector->dxmin * this->alpha / std::sqrt(GRAVITY * maxhw);

    if (tdt == 0)
      tdt = 1e-4;

    // main loop
    for (size_t i = 0; i < this->Qlink.size(); ++i) {
      // valid link??
      if (this->connector->linknodes[i * 2] == -1)
        continue;

      // by convention n1 and n2 are in the arbitrary order of the linknode
      // array (not really important)
      int n1 = this->connector->linknodes[i * 2],
          n2 = this->connector->linknodes[i * 2 + 1];

      // if both hw are below 0 I do not compute and qlink is null
      if (this->hw[n1] <= 0 && this->hw[n2] <= 0) {
        this->Qlink[i] = 0;
        continue;
      }

      // calculating gradient
      float_t Sw =
          (surf[n1] - surf[n2]) / this->connector->get_dx_from_links_idx(i);

      // Getting hflow, the flow height between the two cell:
      float_t hflow = std::max(surf[n1], surf[n2]) -
                      std::max(this->topography[n1], this->topography[n2]);

      // no flow, no algo
      if (hflow <= 0) {
        this->Qlink[i] = 0;
        continue;
      }

      // Get the trasverse dx
      float_t tdx = this->connector->get_traverse_dx_from_links_idx(i) *
                    this->topological_number;

      // Get the discharge
      this->Qlink[i] =
          tdx * (this->Qlink[i] - GRAVITY * hflow * tdt * Sw) /
          (1 + (GRAVITY * hflow * tdt * std::pow(this->mannings, 2) *
                std::abs(this->Qlink[i]) / std::pow(hflow, 10. / 3.)));

      // // Discahrge checkers
      // if(std::abs(this->Qlink[i])/hflow / std::sqrt(GRAVITY * hflow) >
      // this->froude_number ) 	this->Qlink[i] = DAGGER::sgn(this->Qlink[i]) *
      // hflow * std::sqrt(GRAVITY * hflow) * this->froude_number;

      // if(std::abs(this->Qlink[i]) * tdt / tdx > this->hw[n2]/4 &&
      // DAGGER::sgn(this->Qlink[i]) > 0 ) 	this->Qlink[i] = this->hw[n2] *
      // this->connector->get_area_at_node(n2)/5/tdt;

      // else if (std::abs(this->Qlink[i]) * tdt / tdx >  this->hw[n1]/4 &&
      // DAGGER::sgn(this->Qlink[i]) < 0 ) 	this->Qlink[i] = -1 *
      // this->hw[n1] * this->connector->get_area_at_node(n1)/5/tdt;

      // if(this->Qlink[i] > 0)
      dhw[n1] += tdt * this->Qlink[i];
      dhw[n2] -= tdt * this->Qlink[i];
    }

    for (int i = 0; i < this->graph->nnodes; ++i) {
      if (!this->connector->flow_out_model(i)) {
        this->hw[i] += (dhw[i]) / this->connector->get_area_at_node(i);
        if (this->hw[i] < 0)
          this->hw[i] = 0;
      }
    }

    for (int i = 0; i < this->graph->nnodes; ++i) {
      if (!this->connector->flow_out_model(i))
        this->hw[i] +=
            (this->Qbase[i] * tdt) / this->connector->get_area_at_node(i);
    }

    return tdt;
  }

  // Simply copy Qbase in Qwin in order to have something for the first
  // iteration of the nograph version
  void init_Qwin_with_Qbase() {
    for (int i = 0; i < this->graph->nnodes; ++i)
      this->Qwin[i] = this->Qbase[i];
  }
  void increment_hw_from_Qbase(float tdt) {
    for (int i = 0; i < this->graph->nnodes; ++i)
      this->hw[i] +=
          tdt * this->Qbase[i] / this->connector->get_area_at_node(i);
  }

  void set_topological_number(float_t val) { this->topological_number = val; }

  void run_MFD_old(std::string depression_solver) {}

  template <class out_t> out_t get_hw() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->hw);
  }
  template <class out_t> out_t get_Qwin() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->Qwin);
  }
  template <class out_t> out_t get_Qs() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->Qsin);
  }
  template <class out_t> out_t get_Qwout() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->Qwout);
  }
  template <class out_t> out_t get_topography() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->topography);
  }

  template <class out_t> out_t get_spatial_dts() {
    return DAGGER::format_output<std::vector<float_t>, out_t>(this->saved_dts);
  }

  template <class topo_t> void set_Qbase(topo_t &Qbase) {
    auto tin = DAGGER::format_input(Qbase);
    std::vector<double> temp = DAGGER::to_vec(tin);
    this->Qbase = std::move(temp);
  }

  template <class topo_t, class int_t>
  void set_Qs_entry_points(int_t &nonodes, topo_t &Qsb) {
    auto tins = DAGGER::format_input(Qsb);
    auto tini = DAGGER::format_input(nonodes);

    std::vector<double> temp = DAGGER::to_vec(tins);
    this->Qs_entry = std::move(temp);
    std::vector<int> temp2 = DAGGER::to_vec(tini);
    this->Qs_entry_point = std::move(temp2);
  }

  void add_to_hw(float_t val) {
    for (int i = 0; i < this->connector->nnodes; ++i) {
      if (this->connector->flow_out_model(i) == false)
        this->hw[i] += val;
    }
  }

  // SIMPLE TEST WITH EROSION

  void init_erosion() {
    this->Qsin = std::vector<float_t>(this->connector->nnodes, 0.);
    this->Qsout = std::vector<float_t>(this->connector->nnodes, 0.);

    for (size_t i = 0; i < this->Qs_entry_point.size(); ++i)
      this->Qsin[this->Qs_entry_point[i]] = this->Qs_entry[i];
  }

  void force_entry_point_above() {
    std::vector<int> recidx = this->connector->get_empty_neighbour();
    for (size_t i = 0; i < this->Qs_entry_point.size(); ++i) {
      int node = this->Qs_entry_point[i];
      int nn = this->connector->get_neighbour_idx(node, recidx);
      for (int j = 0; j < nn; ++j) {
        int onode = recidx[j];
        if (this->topography[node] + this->hw[node] <
            this->topography[onode] + this->hw[onode]) {
          this->topography[node] = this->topography[onode] + 1e-6;
          this->hw[node] = this->hw[onode] + 1e-6;
        }
      }
    }
  }

  void run_SFD_with_erosion(float_t ke, float_t aexp, float_t tau_c,
                            float_t transport_length, float_t ke_lat,
                            float_t kd_lat) {
    // Reinit the sediment flux
    // this->init_erosion();
    this->qs = std::vector<float_t>(this->connector->nnodes, 0.);
    for (size_t i = 0; i < this->Qs_entry_point.size(); ++i) {
      // std::cout << this->Qs_entry_point[i] << " gets " << this->Qs_entry[i]
      // << std::endl;
      this->qs[this->Qs_entry_point[i]] = this->Qs_entry[i];
    }

    this->force_entry_point_above();

    std::vector<float_t> vmot(this->connector->nnodes, 0.);

    // Getting the filled topo and processing the graph in a SD-D8
    this->graph_automator(true);

    // Some debugging recorder
    float_t mean_edot = 0;
    int N = 0;

    // Iterating top to bottom
    for (int i = this->graph->nnodes - 1; i >= 0; --i) {

      // Next node in line (topological order downstram direction)
      int node = this->graph->Sstack[i];

      // if node is inactive or there is no water -> skip
      if (!this->connector->flow_out_model(node) == false ||
          this->hw[node] <= 0)
        continue;

      // Receiver node
      int rec = this->connector->_Sreceivers[node];

      // distance to the next node and to othogonal nodes
      float_t tdx = this->connector->Sdistance2receivers[node];
      float_t tdw = this->connector->get_travers_dy_from_dx(tdx);

      // Getting the height of flow
      float_t thf = this->get_hflow_from_nodes(node, rec);
      // Setting the water height for the calculation
      float_t thw = this->hflow ? thf : this->hw[node];

      // Setting up a regulator for high water depth
      float_t regulator = 1.;
      // if(thw > 1) thw = 1;
      // if hflow is activated, I check how disconnected it is from the next
      // node
      if (thf < this->hw[node])
        regulator =
            std::pow(thf / this->hw[node], this->sensibility_to_flowdepth);
      // if depth limiter is activated, I reduce erosion after a certain depth,
      // balancing the effect of over carving when
      if (this->depth_limiter && thw > this->depth_threshold)
        regulator *= std::pow(this->depth_threshold / thw,
                              this->sensibility_to_flowdepth);

      // if(regulator )

      // Getting the hydraulic slope
      float_t Sw = std::max(this->topography[node] + this->hw[node] -
                                this->topography[rec] - this->hw[rec],
                            1e-6) /
                   tdx;

      // Calculating tau
      float_t tau = RHO * GRAVITY * Sw * thw * regulator;

      // Debugging statement
      // if(tau > 4000)
      // {
      // 	std::cout << Sw << "||" << thw << std::endl;
      // 	std::cout << "S calc :: " << this->topography[node]<< " + " <<
      // this->hw[node] << " - " << this->topography[rec] << " - " <<
      // this->hw[rec] << " Node - rec were " << node << " and " << rec <<
      // std::endl;
      // }

      // rinit all to 0
      float_t edot = 0, ddot = 0, latedotA = 0, latedotB = 0, latddotA = 0,
              latddotB = 0;

      // if tau exceed the critical motion threshold: calculating edot
      if (tau > tau_c)
        edot = ke * std::pow(tau - tau_c, aexp);

      // DEBUGGING RECORDER
      mean_edot += edot;
      N++;

      // Calculating ddot whatsoever
      ddot = this->qs[node] / transport_length;

      // Accessing orthogonal nodes
      std::pair<int, int> orthonodes =
          this->connector->get_orthogonal_nodes(node, rec);

      // if the orthonode A is valid, run lateral erosion/deposition
      if (orthonodes.first >= 0 && orthonodes.first < this->connector->nnodes) {
        // getting other node
        int &to = orthonodes.first;
        // if lateral node has higher ground: I erode
        if (this->topography[to] > this->topography[node])
          latedotA = (this->topography[to] - this->topography[node]) / tdw *
                     ke_lat * edot;
        // if lateral node has lower ground: I deposit
        if (this->topography[to] < this->topography[node])
          latddotA = (this->topography[node] - this->topography[to]) / tdw *
                     kd_lat * ddot;
      }

      // Same with laterl node B
      if (orthonodes.second >= 0 &&
          orthonodes.second < this->connector->nnodes) {
        int &to = orthonodes.second;
        if (this->topography[to] > this->topography[node])
          latedotB = (this->topography[to] - this->topography[node]) / tdw *
                     ke_lat * edot;
        if (this->topography[to] < this->topography[node])
          latddotB = (this->topography[node] - this->topography[to]) / tdw *
                     kd_lat * ddot;
      }

      // SEt of checks I'll delete after understanding the model
      bool quit = false;

      if (std::isfinite(regulator) == false) {
        std::cout << "regulator" << std::endl;
        quit = true;
      }

      if (std::isfinite(edot) == false) {
        std::cout << "edot" << std::endl;
        quit = true;
      }

      if (std::isfinite(this->qs[node]) == false) {
        std::cout << "this->qs[node]" << std::endl;
        quit = true;
      }

      if (std::isfinite(ddot) == false) {
        std::cout << "ddot" << std::endl;
        quit = true;
      }

      if (std::isfinite(latedotA) == false) {
        std::cout << "latedotA" << std::endl;
        quit = true;
      }

      if (std::isfinite(latedotB) == false) {
        std::cout << "latedotB" << std::endl;
        quit = true;
      }

      if (std::isfinite(tau) == false) {
        std::cout << "tau" << std::endl;
        quit = true;
      }
      if (quit) {
        std::cout << tau << "|" << edot << "|" << latedotA << "|" << latedotB
                  << "|" << ddot << "|" << latddotA << "|" << latddotB << "|"
                  << tdx << "||" << tdw << "|" << this->qs[node] << std::endl;
        throw std::runtime_error("non finite erosion");
      }

      // Second regulator round in case deposition > qs
      regulator = 1.;

      // calculating the full d to be taken from qs
      float_t totd = ddot + latddotA + latddotB;

      // calculating the sediment flux to be consumed by deposition
      if (totd * tdx > this->qs[node] && totd > 0) {
        // regulator is then the max proportion of sediment that can be eaten by
        // dep processes
        regulator = (this->qs[node]) / (totd * tdx);
        if (regulator > 1)
          throw std::runtime_error("should not happen regulator > 1.");
      }

      // Applying the reducer to the  deposition rates
      ddot *= regulator;
      latddotA *= regulator;
      latddotB *= regulator;

      // registering the versical motions to the current node (dep add z and e
      // remove z)
      vmot[node] += (ddot - edot) * this->dt(node);

      // same thing with adjacent nodes A and B
      if (orthonodes.first >= 0 && orthonodes.first < this->connector->nnodes) {
        vmot[orthonodes.first] += (latddotA - latedotA) * this->dt(node);
      }

      if (orthonodes.second >= 0 &&
          orthonodes.first < this->connector->nnodes) {
        vmot[orthonodes.second] += (latddotB - latedotB) * this->dt(node);
      }

      // qs receives what has been removed from land minus regulation and what
      // remains
      this->qs[rec] +=
          (edot + latedotA + latedotB - ddot - latddotA - latddotB) * tdx +
          this->qs[node];

      if (this->qs[rec] < 0) {
        this->qs[rec] = 0;
        // std::cout << "WARNING QS NEG" << std::endl;
      }

      if (std::isfinite(
              (edot + latedotA + latedotB - ddot - latddotA - latddotB) *
              tdx) == false) {
        std::cout << tau << "||" << tau_c << "||" << edot << "||" << latedotA
                  << "||" << latedotB << "||" << ddot << "||" << latddotA
                  << "||" << latddotB << ":::" << regulator << "tt" << tdx
                  << " ii " << tdw << std::endl;
        throw std::runtime_error("non finite flux");
      }
      // if(ddot * tdx > this->qs[node]) ddot = this->qs[node]/tdx;
    }
    // std::cout << "mean edot was " << mean_edot << std::endl;

    for (int i = 0; i < this->graph->nnodes; ++i)
      this->topography[i] += vmot[i];
  }

  void run_MFD_erosion(float_t ke, float_t aexp, float_t tau_c,
                       float_t transport_length, float_t ke_lat, float_t alpha,
                       float_t D, float_t rho_s) {
    return;
    // // Reinit the sediment flux
    // this->init_erosion();

    // std::vector<float_t>vmot(this->connector->nnodes,0.);

    // // Getting the filled topo and processing the graph in MFD (only_SS =
    // false) this->graph_automator(false);

    // // pre-caching the receivers links
    // std::vector<int> reclink = this->connector->get_empty_neighbour();
    // // pre-caching the receivers weights
    // std::vector<float> recweight(reclink.size(),0.);
    // // pre-caching the receivers slopes
    // std::vector<float> recslope(reclink.size(),0.);

    // for(int i = this->graph->nnodes-1; i>=0; --i)
    // {

    // 	int node = this->graph->stack[i];
    // 	if(!this->graph->flow_out_model(node,*this->connector) == false)
    // continue;

    // 	if(this->hw[i] == 0) continue;

    // 	int nr = this->graph->get_receivers_idx_links(node, *this->connector,
    // reclink);

    // 	// First calculating the weights and slope
    // 	// # sumslopes sums the real slopes and wsumslopes is the weighted
    // equivalent (the parting coeff is not always one) 	float_t
    // wsumslopes = 0., sumslopes = 0.; 	for(int j = 0; j<nr; ++j)
    // 	{
    // 		//#getting the link ID
    // 		int lix = reclink[j];
    // 		//#getting the slope and casting it to a minimum numeraicla
    // value 		float_t tsw = std::max(this->get_Sw_at_link(lix), 1e-6);
    // 		recslope[j] = tsw;

    // 		//# Registering and summing the slopes

    // 		wsumslopes += (this->parting_coeff == 1) ? tsw :
    // std::pow(tsw,parting_coeff); 		recweight[j] =
    // (this->parting_coeff
    // == 1) ? tsw : std::pow(tsw,parting_coeff);; 		sumslopes +=
    // tsw;
    // 	}

    // 	// normalising weights
    // 	for(int j = 0; j<nr; ++j)
    // 		recweight[j] = recweight[j]/wsumslopes;

    // 	// Dealing with the divergence of sediments
    // 	this->Qsout[node] = this->Qsin[node];

    // 	float_t cA = this->connector->get_area_at_node(node);
    // 	// first deposition:
    // 	float_t d_dot = this->Qsin[node]/((this->connector->dx +
    // this->connector->dy)/2 * transport_length); 	if(d_dot * cA >
    // this->Qsin[node])
    // 	{
    // 		d_dot = this->Qsin[node]/cA;
    // 	}

    // 	this->Qsout[node] -= d_dot * cA;

    // 	float_t latQ = 0.;
    // 	for(int j=0; j<nr;++j)
    // 	{
    // 		float_t tau = RHO * GRAVITY * this->hw[node] * recslope[j];
    // 		if(tau > tau_c)
    // 		{
    // 			int to = this->graph->get_to_links(reclink[j]);
    // 			float_t tdl =
    // this->connector->get_traverse_dx_from_links_idx(reclink[j]);
    // float_t edot = ke * std::pow(tau - tau_c,aexp);
    // this->Qsout[node]
    // += this->topological_number * edot * cA; 			auto
    // orthonodes = this->connector->get_orthogonal_nodes(node, to);

    // 			if(orthonodes.first >= 0)
    // 			{
    // 				if(this->topography[orthonodes.first] >
    // this->topography[node])
    // 				{
    // 					float_t latedot =
    // (this->topography[orthonodes.first]
    // - this->topography[node])/tdl * edot * ke_lat;
    // latQ += this->topological_number * latedot * cA;
    // vmot[orthonodes.first] -= this->topological_number * latedot *
    // this->dt(node);
    // 				}
    // 				else
    // 				{
    // 					float_t latddot =
    // ((this->topography[node]
    // - this->topography[orthonodes.first])/tdl)/(alpha * std::sqrt(tau /
    // ((rho_s
    // - RHO) * GRAVITY * D ) )) * this->Qsin[node] ;
    // vmot[orthonodes.first] += this->topological_number * latddot *
    // this->dt(node); 					latQ -=
    // this->topological_number * latddot * cA;
    // 				}
    // 			}

    // 			if(orthonodes.second >= 0)
    // 			{
    // 				if(this->topography[orthonodes.second] >
    // this->topography[node])
    // 				{
    // 					float_t latedot =
    // (this->topography[orthonodes.second]
    // - this->topography[node])/tdl * edot * ke_lat;
    // latQ += this->topological_number * latedot * cA;
    // vmot[orthonodes.second] -= this->topological_number * latedot *
    // this->dt(node);
    // 				}
    // 				else
    // 				{
    // 					float_t latddot =
    // ((this->topography[node]
    // - this->topography[orthonodes.second])/tdl)/(alpha * std::sqrt(tau /
    // ((rho_s - RHO) * GRAVITY * D ) )) * this->Qsin[node] ;
    // 					vmot[orthonodes.second] +=
    // this->topological_number
    // * latddot
    // * this->dt(node); 					latQ -=
    // this->topological_number * latddot * cA;
    // 				}
    // 			}

    // 		}

    // 	}

    // 	vmot[node] += (this->Qsin[node] - this->Qsout[node])/cA *
    // this->dt(node);

    // 	this->Qsout[node] += latQ;
    // 	if(this->Qsout[node] < 0) this->Qsout[node] = 0;

    // 	for(int j=0; j<nr;++j)
    // 	{
    // 		int to = this->graph->get_to_links(reclink[j]);
    // 		this->Qsin[to] += this->Qsout[node] * recweight[j];
    // 	}

    // }

    // for(int i = 0; i<this->graph->nnodes; ++i)
    // 	this->topography[i] += vmot[i];

    // std::vector<int> dodons = this->connector->get_empty_neighbour();
    // for(int i=0; i<this->graph->nnodes;++i)
    // {
    // 	if(!this->graph->flow_out_model(i,*this->connector) == false)
    // 	{
    // 		this->graph->get_donors_idx(i, *this->connector, dodons);
    // 		for(auto dono:dodons)
    // 		{
    // 			if(this->topography[dono]<this->topography[i])
    // 				this->topography[i] = this->topography[dono] -
    // 1e-6;
    // 		}
    // 	}
    // }
  }

  void run_MFD_erosion_B(float_t ke, float_t aexp, float_t tau_c,
                         float_t transport_length, float_t ke_lat,
                         float_t kd_lat) {

    // Reinit the sediment fluxes
    this->init_erosion();
    this->rec.init_geo();

    // This vector records and apply the vertical motions at the end of the time
    // step THis reduce the stochasticity but ensure consistency in the sediment
    // flux
    std::vector<float_t> vmot(this->connector->nnodes, 0.);

    // Getting the filled topo and processing the graph in MFD (only_SS = false)
    this->graph_automator(false);

    // pre-caching the receivers links
    std::vector<int> reclink = this->connector->get_empty_neighbour();
    // pre-caching the receivers weights
    std::vector<float> recweight(reclink.size(), 0.);
    // pre-caching the receivers slopes
    std::vector<float> recslope(reclink.size(), 0.);

    // debugging variables recording the factors regulating sediment fluxes
    float_t max_fac = 0;
    float_t min_fac = 1;

    // MAIN LOOP
    for (int i = this->graph->nnodes - 1; i >= 0; --i) {
      // Getting ht enext upstreamest non processed node
      int node = this->graph->stack[i];

      // If the current node is not active: ignore
      if (!this->connector->flow_out_model(node) == false)
        continue;

      // No water? no flow
      if (this->hw[i] == 0)
        continue;

      // getting the receiving links
      int nr = this->connector->get_receivers_idx_links(node, reclink);

      // First calculating the weights and slope
      // # sumslopes sums the real slopes and wsumslopes is the weighted
      // equivalent (the parting coeff is not always one)
      float_t wsumslopes = 0., sumslopes = 0.;
      for (int j = 0; j < nr; ++j) {
        // #getting the link ID
        int lix = reclink[j];
        // #getting the slope and casting it to a minimum numeraicla value
        float_t tsw = std::max(this->get_Sw_at_link(lix), 1e-6);
        if (!this->connector->flow_out_model(
                this->connector->get_to_links(lix)) == false)
          tsw = 1e-6;
        recslope[j] = tsw;

        // # Registering and summing the slopes

        wsumslopes +=
            (this->parting_coeff == 1) ? tsw : std::pow(tsw, parting_coeff);
        recweight[j] =
            (this->parting_coeff == 1) ? tsw : std::pow(tsw, parting_coeff);
        ;
        sumslopes += tsw;
      }

      // normalising weights
      float_t sumsum = 0;
      for (int j = 0; j < nr; ++j) {
        recweight[j] = recweight[j] / wsumslopes;
        sumsum += recweight[j];
      }

      // Dealing with the divergence of sediments
      this->Qsout[node] = this->Qsin[node];

      float_t cA = this->connector->get_area_at_node(node);

      for (int j = 0; j < nr; ++j) {

        int to = this->connector->get_to_links(reclink[j]);
        float_t tdl =
            this->connector->get_traverse_dx_from_links_idx(reclink[j]);
        float_t tdx = this->connector->get_dx_from_links_idx(reclink[j]);
        // float_t tdx = this->connector->get_dx_from_links_idx(reclink[j]);
        float_t tQin = this->Qsin[node] * recweight[j];

        float_t thflow = this->get_hflow_at_link(reclink[j]);

        float_t thw = (this->hflow) ? thflow : this->hw[node];

        // float_t tau = RHO * GRAVITY * thw * recslope[j] * recweight[j];
        float_t tau = RHO * GRAVITY * thw * recslope[j];

        if (this->rec.tau2record) {
          this->rec.tau[node] += tau;
        }

        // // EXPERIMENTAL!!!!!!!!!!
        float_t regulator = 1.;
        // if(thflow < this->hw[node])
        // {
        // 	regulator=std::pow(thflow/this->hw[node],
        // this->sensibility_to_flowdepth);
        // 	// throw std::runtime_error("exp RAN");
        // 	tau = tau * regulator;
        // 	// if(std::pow(thflow/this->hw[node],
        // this->sensibility_to_flowdepth) > 1)
        // 	// 	throw std::runtime_error("!");
        // }

        // first deposition:
        // float_t d_dot = (regulator>0) ? recweight[j]/regulator * tQin/(tdl *
        // transport_length):tQin/cA;
        float_t d_dot = tQin / (transport_length);

        // if(thflow < this->hw[i]/2) tau * thflow/this

        if (true) {

          float_t edot = 0.;

          float_t latedotA = 0;
          float_t latddotA = 0;
          float_t latedotB = 0;
          float_t latddotB = 0;

          auto orthonodes = this->connector->get_orthogonal_nodes(node, to);

          bool crit = tau > tau_c;

          // if(crit) edot = ke * std::pow(tau - tau_c,aexp);
          if (crit)
            edot = ke * std::pow(tau - tau_c, aexp) * recweight[j];

          if (orthonodes.first >= 0) {
            if (!this->connector->flow_out_model(orthonodes.first)) {
              if (this->topography[orthonodes.first] > this->topography[node] &&
                  crit) {
                latedotA = regulator *
                           (this->topography[orthonodes.first] -
                            this->topography[node]) /
                           tdl * edot * ke_lat;
              } else if (tQin > 0 && kd_lat > 0) {
                latddotA = ((this->topography[node] -
                             this->topography[orthonodes.first]) /
                            tdl) *
                           kd_lat * d_dot;
              }
            }
          }

          if (orthonodes.second >= 0) {
            if (!this->connector->flow_out_model(orthonodes.second)) {
              if (this->topography[orthonodes.second] >
                      this->topography[node] &&
                  crit) {
                latedotB = regulator *
                           (this->topography[orthonodes.second] -
                            this->topography[node]) /
                           tdl * edot * ke_lat;
              } else if (tQin > 0 && kd_lat > 0) {
                latddotB = ((this->topography[node] -
                             this->topography[orthonodes.second]) /
                            tdl) *
                           kd_lat * d_dot;
                ;
              }
            }
          }

          float_t fac = 1.;
          // float endQ = tQin/cA + edot + latedotA + latedotB - d_dot -
          // latddotA - latddotB; if(endQ < 0 && (abs(endQ) + tQin/cA) > 0) fac
          // = tQin/(cA * abs(endQ) + tQin);
          float_t endQ = (d_dot + latddotA + latddotB) * tdx;
          if (endQ > tQin)
            fac = tQin / endQ;

          d_dot *= fac;
          latddotA *= fac;
          latddotB *= fac;

          bool quit = false;

          if (std::isfinite(fac) == false) {
            std::cout << "fac" << std::endl;
            quit = true;
          }

          if (std::isfinite(edot) == false) {
            std::cout << "edot" << std::endl;
            quit = true;
          }

          if (std::isfinite(tQin) == false) {
            std::cout << "tQin" << std::endl;
            quit = true;
          }

          if (std::isfinite(d_dot) == false) {
            std::cout << "d_dot" << std::endl;
            quit = true;
          }

          if (std::isfinite(latedotA) == false) {
            std::cout << "latedotA" << std::endl;
            quit = true;
          }

          if (std::isfinite(latedotB) == false) {
            std::cout << "latedotB" << std::endl;
            quit = true;
          }
          if (quit) {
            std::cout << fac << "|" << edot << "|" << latedotA << "|"
                      << latedotB << "|" << d_dot << "|" << latddotA << "|"
                      << latddotB << "|" << cA << "|" << tQin << std::endl;
            throw std::runtime_error("non finite erosion");
          }

          if (max_fac < fac)
            max_fac = fac;
          if (min_fac > fac)
            min_fac = fac;

          vmot[node] -= (edot - d_dot) * this->dt(node);
          if (orthonodes.first >= 0)
            vmot[orthonodes.first] -= (latedotA - latddotA) * this->dt(node);
          if (orthonodes.second >= 0)
            vmot[orthonodes.second] -= (latedotB - latddotB) * this->dt(node);

          if (std::isfinite(
                  (edot + latedotB + latddotA - d_dot - latddotA - latddotB) *
                      cA +
                  tQin) == false) {
            std::cout << fac *
                                 (edot + latedotB + latddotA - d_dot -
                                  latddotA - latddotB) *
                                 cA +
                             tQin
                      << std::endl;
            std::cout << fac << "|" << edot << "|" << latedotA << "|"
                      << latedotB << "|" << d_dot << "|" << latddotA << "|"
                      << latddotB << "|" << cA << "|" << tQin << std::endl;
            throw std::runtime_error("nonfinitestuff");
          }

          this->Qsin[to] +=
              (edot + latedotB + latddotA - d_dot - latddotA - latddotB) * tdx +
              tQin;
          if (this->Qsin[to] < 0)
            this->Qsin[to] = 0;

          if (this->rec.edot2record)
            this->rec.edot[node] += edot;
          if (this->rec.ddot2record)
            this->rec.ddot[node] += d_dot;
          if (this->rec.lateral_edot2record && orthonodes.first >= 0) {
            this->rec.lateral_edot[orthonodes.first] += latedotA;
          }
          if (this->rec.lateral_edot2record && orthonodes.second >= 0) {
            this->rec.lateral_edot[orthonodes.second] += latedotB;
          }
          if (this->rec.lateral_ddot2record && orthonodes.first >= 0) {
            this->rec.lateral_ddot[orthonodes.first] += latddotA;
          }
          if (this->rec.lateral_ddot2record && orthonodes.second >= 0) {
            this->rec.lateral_ddot[orthonodes.second] += latddotB;
          }
        }
      }

      // for(int j=0; j<nr;++j)
      // {
      // 	int to = this->graph->get_to_links(reclink[j]);
      // 	this->Qsin[to] += this->Qsin[node] * recweight[j];
      // }
    }

    // std::cout << max_fac << "<--- max fac||min fac ---->" << min_fac <<
    // std::endl;

    for (int i = 0; i < this->graph->nnodes; ++i) {
      this->topography[i] += vmot[i];
    }

    // this->softbound(true,true);
  }

  void run_MFD_static(float_t ke, float_t aexp, float_t tau_c,
                      float_t transport_length, float_t ke_lat, float_t kd_lat,
                      float_t dt_hydro, int n_erosion, float_t dt_erosion,
                      bool stochastic_slope, bool no_erosion_at_all) {
    bool erosion = false;
    this->stochaslope = 0;

    float_t vmottot = 0;

    for (int nit = 1; nit <= n_erosion; ++nit) {

      if (nit == n_erosion && !no_erosion_at_all) {
        if (stochastic_slope)
          this->stochaslope = 1.;
        erosion = true;
      }

      // Reinit the sediment fluxes
      this->init_erosion();
      if (erosion)
        this->rec.init_geo();
      this->rec.init_water();

      // This vector records and apply the vertical motions at the end of the
      // time step THis reduce the stochasticity but ensure consistency in the
      // sediment flux
      std::vector<float_t> vmot(this->connector->nnodes, 0.);

      // Getting the filled topo and processing the graph in MFD (only_SS =
      // false)
      this->graph_automator(false);

      // reinitialising the Qw in and out
      for (auto &v : this->Qwin)
        v = 0;
      for (auto &v : this->Qwout)
        v = 0;

      // pre-caching the receivers links
      std::vector<int> reclink = this->connector->get_empty_neighbour();
      // pre-caching the receivers weights
      std::vector<float> recweight(reclink.size(), 0.);
      // pre-caching the receivers slopes
      std::vector<float> recslope(reclink.size(), 0.);

      // debugging variables recording the factors regulating sediment fluxes
      float_t max_fac = 0;
      float_t min_fac = 1;

      // MAIN LOOP
      for (int i = this->graph->nnodes - 1; i >= 0; --i) {
        // Getting ht enext upstreamest non processed node
        int node = this->graph->stack[i];

        // If the current node is not active: ignore
        if (this->connector->flow_out_or_pit(node))
          continue;

        // No water? no flow
        // if(this->hw[i] == 0) continue;

        // getting the receiving links
        int nr = this->connector->get_receivers_idx_links(node, reclink);

        // ## (Rare) no receivers?? then I skip
        if (nr == 0) {
          continue;
        }

        // ## Accumulating local discharge
        this->Qwin[node] += this->Qbase[node];

        // SECTION OF THE CODE MANAGING WHO GIVES TO THE BOUNDARY
        // bool is_to_boundary = false;
        // for(int j = 0; j< nr; ++j)
        // {
        // 	int to = this->graph->get_to_links(reclink[j]);
        // 	if(!this->graph->flow_out_or_pit(to,*this->connector) &&
        // this->connector->boundary[to] != 3) 		is_to_boundary = true;
        // }

        // if(is_to_boundary)
        // {
        // 	int nj = 0;
        // 	for(int j = 0; j< nr; ++j)
        // 	{
        // 		int to = this->graph->get_to_links(reclink[j]);

        // 		if(!this->graph->flow_out_model(this->graph->get_to_links(reclink[j]),*this->connector)
        // && this->connector->boundary[to] != 3)
        // 		{
        // 			reclink[nj] = reclink[j];
        // 			++nj;
        // 		}

        // 	}
        // 	nr = nj;
        // 	//## (Rare) no receivers?? then I skip
        // 	if(nr == 0)
        // 	{
        // 		continue;
        // 	}

        // }

        // First calculating the weights and slope
        // # sumslopes sums the real slopes and wsumslopes is the weighted
        // equivalent (the parting coeff is not always one)
        float_t wsumslopes = 0., sumslopes = 0., maxslope = 0.;
        for (int j = 0; j < nr; ++j) {

          // #getting the link ID
          int lix = reclink[j];

          if (this->connector->is_link_valid(lix) == false)
            continue;

          // #getting the slope and casting it to a minimum numeraicla value
          //  float_t tsw = (is_to_boundary)? this->boundary_slope :
          //  std::max(this->get_Sw_at_link(lix), 1e-6);
          float_t tsw = std::max(this->get_Sw_at_link(lix), 1e-6);

          // float_t tsw = std::max(this->get_Sw_at_link_custom_dx(lix,
          // this->connector->dx), 1e-6);
          // if(!this->graph->flow_out_model(this->graph->get_to_links(lix),*this->connector)
          // == false) tsw = 1e-6;

          recslope[j] = tsw;

          if (tsw > maxslope)
            maxslope = tsw;

          // # Registering and summing the slopes
          float_t weisw =
              (this->parting_coeff == 1) ? tsw : std::pow(tsw, parting_coeff);

          if (this->stochaslope > 0)
            weisw *= this->randu.get() + 1e-6;

          wsumslopes += weisw;
          recweight[j] = weisw;
          sumslopes += tsw;
        }

        // normalising weights
        float_t sumsum = 0;
        for (int j = 0; j < nr; ++j) {
          if (this->connector->is_link_valid(reclink[j]) == false)
            continue;
          recweight[j] = recweight[j] / wsumslopes;
          sumsum += recweight[j];
        }

        // Dealing with the divergence of sediments
        this->Qsout[node] = this->Qsin[node];

        float_t cA = this->connector->get_area_at_node(node);
        float_t facgrad = 0;

        for (int j = 0; j < nr; ++j) {
          // ### Not valid 4 some reasons? skip
          if (this->connector->is_link_valid(reclink[j]) == false)
            continue;

          int to = this->connector->get_to_links(reclink[j]);

          // bool is_to_boundary =
          // !!this->graph->flow_out_model(to,*this->connector);

          float_t tdl =
              this->connector->get_traverse_dx_from_links_idx(reclink[j]);

          float_t tdx = this->connector->get_dx_from_links_idx(reclink[j]);
          // float_t tdx = (is_to_boundary)? this->connector->dx :
          // this->connector->get_dx_from_links_idx(reclink[j]); float_t tdx =
          // this->connector->dx;

          // float_t tdx = this->connector->get_dx_from_links_idx(reclink[j]);
          float_t tQin = this->Qsin[node] * recweight[j];

          float_t thflow = this->get_hflow_at_link(reclink[j]);

          // ### accumulating the equation factor
          if (this->hflow) {
            facgrad +=
                std::pow(this->get_hflow_at_link(reclink[j]), FIVETHIRD) *
                recslope[j] / tdx;
          } else
            facgrad += recslope[j] / tdx;

          this->Qwin[to] += this->Qwin[node] * recweight[j];

          // // EXPERIMENTAL!!!!!!!!!!
          float_t regulator = 1.;
          // if(thflow < this->hw[node])
          // {
          // 	regulator=std::pow(thflow/this->hw[node],
          // this->sensibility_to_flowdepth);
          // 	// throw std::runtime_error("exp RAN");
          // 	tau = tau * regulator;
          // 	// if(std::pow(thflow/this->hw[node],
          // this->sensibility_to_flowdepth) > 1)
          // 	// 	throw std::runtime_error("!");
          // }

          // first deposition:
          // float_t d_dot = (regulator>0) ? recweight[j]/regulator * tQin/(tdl
          // * transport_length):tQin/cA;

          // if(thflow < this->hw[i]/2) tau * thflow/this

          if (this->connector->boundaries.forcing_io(node)) {
            this->Qsin[to] += this->Qsin[node] * recweight[j];
          } else if (erosion && this->hw[node] > 0) {

            float_t thw = (this->hflow) ? thflow : this->hw[node];

            // float_t tau = RHO * GRAVITY * thw * recslope[j] * recweight[j];
            float_t tau = RHO * GRAVITY * thw * recslope[j];

            // if(tau > 80)
            // std::cout << "RHO::" << RHO << " | " << "GRAVITY::" << GRAVITY <<
            // " | " << "thw::" << thw << " | " << "recslope[j]::" <<
            // recslope[j] << " | " << std::endl;

            float_t d_dot = tQin / (transport_length);

            if (this->rec.tau2record) {
              this->rec.tau[node] += tau * recweight[j];
            }

            float_t edot = 0.;

            float_t latedotA = 0;
            float_t latddotA = 0;
            float_t latedotB = 0;
            float_t latddotB = 0;

            auto orthonodes = this->connector->get_orthogonal_nodes(node, to);

            bool crit = tau > tau_c;

            if (crit)
              edot = ke * std::pow(tau - tau_c, aexp);
            // if(crit) edot = ke * std::pow(tau - tau_c,aexp) * recweight[j];

            if (orthonodes.first >= 0) {
              if (this->connector->flow_out_model(orthonodes.first) == false) {
                if (this->topography[orthonodes.first] >
                        this->topography[node] &&
                    crit) {
                  latedotA = regulator *
                             (this->topography[orthonodes.first] -
                              this->topography[node]) /
                             tdl * edot * ke_lat;
                } else if (tQin > 0 && kd_lat > 0) {
                  latddotA = ((this->topography[node] -
                               this->topography[orthonodes.first]) /
                              tdl) *
                             kd_lat * d_dot;
                }
              }
            }

            if (orthonodes.second >= 0) {
              if (this->connector->flow_out_model(orthonodes.second) == false) {
                if (this->topography[orthonodes.second] >
                        this->topography[node] &&
                    crit) {
                  latedotB = regulator *
                             (this->topography[orthonodes.second] -
                              this->topography[node]) /
                             tdl * edot * ke_lat;
                } else if (tQin > 0 && kd_lat > 0) {
                  latddotB = ((this->topography[node] -
                               this->topography[orthonodes.second]) /
                              tdl) *
                             kd_lat * d_dot;
                  ;
                }
              }
            }

            float_t fac = 1.;
            // float endQ = tQin/cA + edot + latedotA + latedotB - d_dot -
            // latddotA - latddotB; if(endQ < 0 && (abs(endQ) + tQin/cA) > 0)
            // fac = tQin/(cA * abs(endQ) + tQin);
            float_t endQ = (d_dot * tdx + latddotA * tdl + latddotB * tdl);
            if (endQ > tQin)
              fac = tQin / endQ;

            d_dot *= fac;
            latddotA *= fac;
            latddotB *= fac;

            bool quit = false;

            if (std::isfinite(fac) == false) {
              std::cout << "fac" << std::endl;
              quit = true;
            }

            if (std::isfinite(edot) == false) {
              std::cout << "edot" << std::endl;
              quit = true;
            }

            if (std::isfinite(tQin) == false) {
              std::cout << "tQin" << std::endl;
              quit = true;
            }

            if (std::isfinite(d_dot) == false) {
              std::cout << "d_dot" << std::endl;
              quit = true;
            }

            if (std::isfinite(latedotA) == false) {
              std::cout << "latedotA" << std::endl;
              quit = true;
            }

            if (std::isfinite(latedotB) == false) {
              std::cout << "latedotB" << std::endl;
              quit = true;
            }
            if (quit) {
              std::cout << fac << "|" << edot << "|" << latedotA << "|"
                        << latedotB << "|" << d_dot << "|" << latddotA << "|"
                        << latddotB << "|" << cA << "|" << tQin << std::endl;
              throw std::runtime_error("non finite erosion");
            }

            if (max_fac < fac)
              max_fac = fac;
            if (min_fac > fac)
              min_fac = fac;

            vmot[node] -= (edot - d_dot) * dt_erosion;
            if (orthonodes.first >= 0)
              vmot[orthonodes.first] -= (latedotA - latddotA) * dt_erosion;
            if (orthonodes.second >= 0)
              vmot[orthonodes.second] -= (latedotB - latddotB) * dt_erosion;

            if (std::isfinite(
                    (edot + latedotB + latddotA - d_dot - latddotA - latddotB) *
                        cA +
                    tQin) == false) {
              std::cout << fac *
                                   (edot + latedotB + latddotA - d_dot -
                                    latddotA - latddotB) *
                                   cA +
                               tQin
                        << std::endl;
              std::cout << fac << "|" << edot << "|" << latedotA << "|"
                        << latedotB << "|" << d_dot << "|" << latddotA << "|"
                        << latddotB << "|" << cA << "|" << tQin << std::endl;
              throw std::runtime_error("nonfinitestuff");
            }

            this->Qsin[to] +=
                (edot + latedotB + latddotA - d_dot - latddotA - latddotB) *
                    tdx +
                tQin;
            if (this->Qsin[to] < 0) {
              std::cout << "happens" << std::endl;
              this->Qsin[to] = 0;
            }

            vmottot -=
                (edot + latedotB + latddotA - d_dot - latddotA - latddotB) *
                dt_erosion;

            if (this->rec.edot2record)
              this->rec.edot[node] += edot;
            if (this->rec.ddot2record)
              this->rec.ddot[node] += d_dot;
            if (this->rec.lateral_edot2record && orthonodes.first >= 0) {
              this->rec.lateral_edot[orthonodes.first] += latedotA;
            }
            if (this->rec.lateral_edot2record && orthonodes.second >= 0) {
              this->rec.lateral_edot[orthonodes.second] += latedotB;
            }
            if (this->rec.lateral_ddot2record && orthonodes.first >= 0) {
              this->rec.lateral_ddot[orthonodes.first] += latddotA;
            }
            if (this->rec.lateral_ddot2record && orthonodes.second >= 0) {
              this->rec.lateral_ddot[orthonodes.second] += latddotB;
            }
          }
        }

        // ## Finally calculating Qout
        if (this->hflow)
          this->Qwout[node] = this->topological_number * facgrad /
                              std::sqrt(maxslope) * cA /
                              this->mannings; // * (1/(4 * this->connector->dxy
                                              // )+ 1/ (2*this->connector->dx) +
                                              // 1/(2 * this->connector->dy));
        else
          this->Qwout[node] = this->topological_number * facgrad /
                              std::sqrt(maxslope) *
                              std::pow(this->hw[node], FIVETHIRD) * cA /
                              this->mannings; // * (1/(4 * this->connector->dxy
                                              // )+ 1/ (2*this->connector->dx) +
                                              // 1/(2 * this->connector->dy));

        // ## Divergence of Q to apply changes to water height
        float_t tdhw = (this->Qwin[node] - this->Qwout[node]) * dt_hydro /
                       this->connector->get_area_at_node(node);
        if (this->rec.dhw2record)
          this->rec.dhw[node] = tdhw;
        this->hw[node] += tdhw;

        // ## water height cannot be 0
        if (this->hw[node] < 0)
          this->hw[node] = 0;
      }

      // std::cout << max_fac << "<--- max fac||min fac ---->" << min_fac <<
      // std::endl;

      if (erosion) {
        float_t sumvm = 0;
        // float_t vmottot_o = vmottot;
        for (int i = 0; i < this->graph->nnodes; ++i) {
          if (this->connector->boundaries.forcing_io(i))
            continue;
          this->topography[i] += vmot[i];
          if (this->rec.vmot2record)
            this->rec.vmot[i] = vmot[i];
          vmottot -= vmot[i];
          sumvm += vmot[i];
        }
        // std::cout << "vmottot theoretical = " << vmottot_o << " vmotot
        // applied = " << sumvm << std::endl;;
      }

      // this->softbound(true,true);
    }
  }

  void run_MFD_static_SPL(float_t Ker, float_t Kerl, float_t Kdepl,
                          float_t mexp, float_t nexp, float_t transport_length,
                          float_t dt_hydro, int n_erosion, float_t dt_erosion,
                          bool stochastic_slope, bool no_erosion_at_all) {
    std::cout << "JUST FOR TESTING PURPOSES" << std::endl;
    bool erosion = false;
    this->stochaslope = 0;

    float_t vmottot = 0;

    for (int nit = 1; nit <= n_erosion; ++nit) {

      if (nit == n_erosion && !no_erosion_at_all) {
        if (stochastic_slope)
          this->stochaslope = 1.;
        erosion = true;
      }

      // Reinit the sediment fluxes
      this->init_erosion();
      if (erosion)
        this->rec.init_geo();
      this->rec.init_water();

      // This vector records and apply the vertical motions at the end of the
      // time step THis reduce the stochasticity but ensure consistency in the
      // sediment flux
      std::vector<float_t> vmot(this->connector->nnodes, 0.);

      // Getting the filled topo and processing the graph in MFD (only_SS =
      // false)
      this->graph_automator(false);

      // reinitialising the Qw in and out
      for (auto &v : this->Qwin)
        v = 0;
      for (auto &v : this->Qwout)
        v = 0;

      // pre-caching the receivers links
      std::vector<int> reclink = this->connector->get_empty_neighbour();
      // pre-caching the receivers weights
      std::vector<float> recweight(reclink.size(), 0.);
      // pre-caching the receivers slopes
      std::vector<float> recslope(reclink.size(), 0.);

      // debugging variables recording the factors regulating sediment fluxes
      float_t max_fac = 0;
      float_t min_fac = 1;

      // MAIN LOOP
      for (int i = this->graph->nnodes - 1; i >= 0; --i) {
        // Getting ht enext upstreamest non processed node
        int node = this->graph->stack[i];

        // If the current node is not active: ignore
        if (this->connector->flow_out_or_pit(node))
          continue;

        // No water? no flow
        // if(this->hw[i] == 0) continue;

        // getting the receiving links
        int nr = this->connector->get_receivers_idx_links(node, reclink);

        // ## (Rare) no receivers?? then I skip
        if (nr == 0) {
          continue;
        }

        // ## Accumulating local discharge
        this->Qwin[node] += this->Qbase[node];

        // SECTION OF THE CODE MANAGING WHO GIVES TO THE BOUNDARY
        // bool is_to_boundary = false;
        // for(int j = 0; j< nr; ++j)
        // {
        // 	int to = this->graph->get_to_links(reclink[j]);
        // 	if(!this->graph->flow_out_or_pit(to,*this->connector) &&
        // this->connector->boundary[to] != 3) 		is_to_boundary = true;
        // }

        // if(is_to_boundary)
        // {
        // 	int nj = 0;
        // 	for(int j = 0; j< nr; ++j)
        // 	{
        // 		int to = this->graph->get_to_links(reclink[j]);

        // 		if(!this->graph->flow_out_model(this->graph->get_to_links(reclink[j]),*this->connector)
        // && this->connector->boundary[to] != 3)
        // 		{
        // 			reclink[nj] = reclink[j];
        // 			++nj;
        // 		}

        // 	}
        // 	nr = nj;
        // 	//## (Rare) no receivers?? then I skip
        // 	if(nr == 0)
        // 	{
        // 		continue;
        // 	}

        // }

        // First calculating the weights and slope
        // # sumslopes sums the real slopes and wsumslopes is the weighted
        // equivalent (the parting coeff is not always one)
        float_t wsumslopes = 0., sumslopes = 0., maxslope = 0.;
        float_t facgrad = 0;

        for (int j = 0; j < nr; ++j) {

          // #getting the link ID
          int lix = reclink[j];

          if (this->connector->is_link_valid(lix) == false)
            continue;

          // #getting the slope and casting it to a minimum numeraicla value
          //  float_t tsw = (is_to_boundary)? this->boundary_slope :
          //  std::max(this->get_Sw_at_link(lix), 1e-6);
          float_t tsw = std::max(this->get_Sw_at_link(lix), 1e-6);

          // float_t tsw = std::max(this->get_Sw_at_link_custom_dx(lix,
          // this->connector->dx), 1e-6);
          // if(!this->graph->flow_out_model(this->graph->get_to_links(lix),*this->connector)
          // == false) tsw = 1e-6;

          recslope[j] = tsw;

          if (tsw > maxslope)
            maxslope = tsw;

          // # Registering and summing the slopes
          float_t weisw =
              (this->parting_coeff == 1) ? tsw : std::pow(tsw, parting_coeff);

          if (this->stochaslope > 0)
            weisw *= this->randu.get() + 1e-6;

          wsumslopes += weisw;
          recweight[j] = weisw;
          sumslopes += tsw;
        }

        // normalising weights
        float_t sumsum = 0;
        for (int j = 0; j < nr; ++j) {
          if (this->connector->is_link_valid(reclink[j]) == false)
            continue;

          int to = this->connector->get_to_links(reclink[j]);
          float_t tdx = this->connector->get_dx_from_links_idx(reclink[j]);

          recweight[j] = recweight[j] / wsumslopes;

          sumsum += recweight[j];

          // ### accumulating the equation factor
          if (this->hflow) {
            facgrad +=
                std::pow(this->get_hflow_at_link(reclink[j]), FIVETHIRD) *
                recslope[j] / tdx;
          } else
            facgrad += recslope[j] / tdx;

          this->Qwin[to] += this->Qwin[node] * recweight[j];
        }
        float_t cA = this->connector->get_area_at_node(node);

        // ## Finally calculating Qout
        if (this->hflow)
          this->Qwout[node] = this->topological_number * facgrad /
                              std::sqrt(maxslope) * cA /
                              this->mannings; // * (1/(4 * this->connector->dxy
                                              // )+ 1/ (2*this->connector->dx) +
                                              // 1/(2 * this->connector->dy));
        else
          this->Qwout[node] = this->topological_number * facgrad /
                              std::sqrt(maxslope) *
                              std::pow(this->hw[node], FIVETHIRD) * cA /
                              this->mannings; // * (1/(4 * this->connector->dxy
                                              // )+ 1/ (2*this->connector->dx) +
                                              // 1/(2 * this->connector->dy));

        // ## Divergence of Q to apply changes to water height
        float_t tdhw = (this->Qwin[node] - this->Qwout[node]) * dt_hydro /
                       this->connector->get_area_at_node(node);
        if (this->rec.dhw2record)
          this->rec.dhw[node] = tdhw;
        this->hw[node] += tdhw;

        // ## water height cannot be 0
        if (this->hw[node] < 0)
          this->hw[node] = 0;

        // Dealing with the divergence of sediments
        this->Qsout[node] = this->Qsin[node];

        for (int j = 0; j < nr; ++j) {
          // ### Not valid 4 some reasons? skip
          if (this->connector->is_link_valid(reclink[j]) == false)
            continue;

          int to = this->connector->get_to_links(reclink[j]);

          // bool is_to_boundary =
          // !!this->graph->flow_out_model(to,*this->connector);

          float_t tdl =
              this->connector->get_traverse_dx_from_links_idx(reclink[j]);

          float_t tdx = this->connector->get_dx_from_links_idx(reclink[j]);
          // float_t tdx = (is_to_boundary)? this->connector->dx :
          // this->connector->get_dx_from_links_idx(reclink[j]); float_t tdx =
          // this->connector->dx;

          // float_t tdx = this->connector->get_dx_from_links_idx(reclink[j]);
          float_t tQin = this->Qsin[node] * recweight[j];

          // float_t thflow = this->get_hflow_at_link(reclink[j]);

          // // EXPERIMENTAL!!!!!!!!!!
          float_t regulator = 1.;
          // if(thflow < this->hw[node])
          // {
          // 	regulator=std::pow(thflow/this->hw[node],
          // this->sensibility_to_flowdepth);
          // 	// throw std::runtime_error("exp RAN");
          // 	tau = tau * regulator;
          // 	// if(std::pow(thflow/this->hw[node],
          // this->sensibility_to_flowdepth) > 1)
          // 	// 	throw std::runtime_error("!");
          // }

          // first deposition:
          // float_t d_dot = (regulator>0) ? recweight[j]/regulator * tQin/(tdl
          // * transport_length):tQin/cA;

          // if(thflow < this->hw[i]/2) tau * thflow/this

          if (erosion && this->hw[node] > 0 &&
              this->connector->boundaries.forcing_io(node) == false) {

            float_t d_dot = tQin / (transport_length);

            float_t edot = 0.;

            float_t latedotA = 0;
            float_t latddotA = 0;
            float_t latedotB = 0;
            float_t latddotB = 0;

            auto orthonodes = this->connector->get_orthogonal_nodes(node, to);

            edot = Ker * std::pow(recslope[j], nexp) *
                   std::pow(this->Qwout[node] * recweight[j], mexp);

            if (orthonodes.first >= 0) {
              if (this->connector->flow_out_model(orthonodes.first) == false) {
                if (this->topography[orthonodes.first] >
                    this->topography[node]) {
                  latedotA = regulator *
                             (this->topography[orthonodes.first] -
                              this->topography[node]) /
                             tdl * edot * Kerl;
                } else if (tQin > 0 && Kdepl > 0) {
                  latddotA = ((this->topography[node] -
                               this->topography[orthonodes.first]) /
                              tdl) *
                             Kdepl * d_dot;
                }
              }
            }

            if (orthonodes.second >= 0) {
              if (this->connector->flow_out_model(orthonodes.second) == false) {
                if (this->topography[orthonodes.second] >
                    this->topography[node]) {
                  latedotB = regulator *
                             (this->topography[orthonodes.second] -
                              this->topography[node]) /
                             tdl * edot * Kerl;
                } else if (tQin > 0 && Kdepl > 0) {
                  latddotB = ((this->topography[node] -
                               this->topography[orthonodes.second]) /
                              tdl) *
                             Kdepl * d_dot;
                  ;
                }
              }
            }

            float_t fac = 1.;
            // float endQ = tQin/cA + edot + latedotA + latedotB - d_dot -
            // latddotA - latddotB; if(endQ < 0 && (abs(endQ) + tQin/cA) > 0)
            // fac = tQin/(cA * abs(endQ) + tQin);
            float_t endQ = (d_dot * tdx + latddotA * tdl + latddotB * tdl);
            if (endQ > tQin)
              fac = tQin / endQ;

            d_dot *= fac;
            latddotA *= fac;
            latddotB *= fac;

            bool quit = false;

            if (std::isfinite(fac) == false) {
              std::cout << "fac" << std::endl;
              quit = true;
            }

            if (std::isfinite(edot) == false) {
              std::cout << "edot" << std::endl;
              quit = true;
            }

            if (std::isfinite(tQin) == false) {
              std::cout << "tQin" << std::endl;
              quit = true;
            }

            if (std::isfinite(d_dot) == false) {
              std::cout << "d_dot" << std::endl;
              quit = true;
            }

            if (std::isfinite(latedotA) == false) {
              std::cout << "latedotA" << std::endl;
              quit = true;
            }

            if (std::isfinite(latedotB) == false) {
              std::cout << "latedotB" << std::endl;
              quit = true;
            }
            if (quit) {
              std::cout << fac << "|" << edot << "|" << latedotA << "|"
                        << latedotB << "|" << d_dot << "|" << latddotA << "|"
                        << latddotB << "|" << cA << "|" << tQin << std::endl;
              throw std::runtime_error("non finite erosion");
            }

            if (max_fac < fac)
              max_fac = fac;
            if (min_fac > fac)
              min_fac = fac;

            vmot[node] -= (edot - d_dot) * dt_erosion;
            if (orthonodes.first >= 0)
              vmot[orthonodes.first] -= (latedotA - latddotA) * dt_erosion;
            if (orthonodes.second >= 0)
              vmot[orthonodes.second] -= (latedotB - latddotB) * dt_erosion;

            if (std::isfinite(
                    (edot + latedotB + latddotA - d_dot - latddotA - latddotB) *
                        cA +
                    tQin) == false) {
              std::cout << fac *
                                   (edot + latedotB + latddotA - d_dot -
                                    latddotA - latddotB) *
                                   cA +
                               tQin
                        << std::endl;
              std::cout << fac << "|" << edot << "|" << latedotA << "|"
                        << latedotB << "|" << d_dot << "|" << latddotA << "|"
                        << latddotB << "|" << cA << "|" << tQin << std::endl;
              throw std::runtime_error("nonfinitestuff");
            }

            this->Qsin[to] +=
                (edot + latedotB + latddotA - d_dot - latddotA - latddotB) *
                    tdx +
                tQin;
            if (this->Qsin[to] < 0) {
              std::cout << "happens" << std::endl;
              this->Qsin[to] = 0;
            }

            vmottot -=
                (edot + latedotB + latddotA - d_dot - latddotA - latddotB) *
                dt_erosion;

            if (this->rec.edot2record)
              this->rec.edot[node] += edot;
            if (this->rec.ddot2record)
              this->rec.ddot[node] += d_dot;
            if (this->rec.lateral_edot2record && orthonodes.first >= 0) {
              this->rec.lateral_edot[orthonodes.first] += latedotA;
            }
            if (this->rec.lateral_edot2record && orthonodes.second >= 0) {
              this->rec.lateral_edot[orthonodes.second] += latedotB;
            }
            if (this->rec.lateral_ddot2record && orthonodes.first >= 0) {
              this->rec.lateral_ddot[orthonodes.first] += latddotA;
            }
            if (this->rec.lateral_ddot2record && orthonodes.second >= 0) {
              this->rec.lateral_ddot[orthonodes.second] += latddotB;
            }
          }
        }
      }

      // std::cout << max_fac << "<--- max fac||min fac ---->" << min_fac <<
      // std::endl;

      if (erosion) {
        float_t sumvm = 0;
        // float_t vmottot_o = vmottot;
        for (int i = 0; i < this->graph->nnodes; ++i) {
          if (this->connector->boundaries.forcing_io(i))
            continue;
          this->topography[i] += vmot[i];
          if (this->rec.vmot2record)
            this->rec.vmot[i] = vmot[i];
          vmottot -= vmot[i];
          sumvm += vmot[i];
        }
        // std::cout << "vmottot theoretical = " << vmottot_o << " vmotot
        // applied = " << sumvm << std::endl;;
      }

      // this->softbound(true,true);
    }
  }

  void out_boundary_match_donors() {
    // std::vector<int> dodons = this->connector->get_empty_neighbour();
    // for(int i=0; i<this->graph->nnodes;++i)
    // {
    // 	if(!this->graph->flow_out_model(i,*this->connector) == false &&
    // this->connector->boundary[i] != 3)
    // 	{
    // 		int nn = this->connector->get_neighbour_idx(i, dodons);
    // 		float_t minsurf = std::numeric_limits<float_t>::max();
    // 		for(int j=0; j<nn; ++j)
    // 		{

    // 			int dono = dodons[j];

    // 			if(this->connector->boundary[dono] <= 0 ||
    // this->connector->boundary[dono] == 3) continue;

    // 			// if(topo_)
    // 			// 	if(this->topography[dono]<this->topography[i])
    // 			// 		this->topography[i] =
    // this->topography[dono]
    // - 1e-6;
    // 			// // if(hw_)
    // 			// // 	if(this->topography[dono] +
    // this->hw[dono]>this->topography[i]
    // + this->hw[i])
    // 			// // 		this->topography[i] =
    // this->topography[dono]
    // - 1e-6;

    // 				if(this->topography[dono] + this->hw[dono] <
    // minsurf)
    // 				{
    // 					// this->topography[i] =
    // this->topography[dono]
    // - 1e-6;
    // 					// this->hw[i] = this->hw[dono] - 1e-6;
    // 					minsurf = this->topography[dono] +
    // this->hw[dono];
    // 				}

    // 		}

    // 		if(std::numeric_limits<float_t>::max() != minsurf)
    // 		{
    // 			if(minsurf < this->topography[i])
    // 			{
    // 				this->hw[i] = 0;
    // 				this->topography[i] = minsurf - 1e-6;
    // 			}
    // 			else
    // 			{
    // 				this->hw[i] = minsurf - this->topography[i] -
    // 1e-6;
    // 			}

    // 		}

    // 	}
    // }
  }

  // void softbound(bool topo_, bool hw_)
  // {
  // 	std::vector<int> dodons = this->connector->get_empty_neighbour();
  // 	for(int i=0; i<this->graph->nnodes;++i)
  // 	{
  // 		if(!this->graph->flow_out_model(i,*this->connector) == false)
  // 		{
  // 			this->graph->get_donors_idx(i, *this->connector,
  // dodons); 			for(auto dono:dodons)
  // 			{
  // 				if(topo_)
  // 					if(this->topography[dono]<this->topography[i])
  // 						this->topography[i] =
  // this->topography[dono]
  // - 1e-6;
  // 				// if(hw_)
  // 				// 	if(this->topography[dono] +
  // this->hw[dono]>this->topography[i] + this->hw[i])
  // 				// 		this->topography[i] =
  // this->topography[dono]
  // - 1e-6;

  // 			}
  // 		}
  // 	}
  // }

  float_t get_Sw_at_link(int linkidx) {
    int n1 = this->connector->get_from_links(linkidx);
    int n2 = this->connector->get_to_links(linkidx);
    return (this->hw[n1] - this->hw[n2] + this->topography[n1] -
            this->topography[n2]) /
           this->connector->get_dx_from_links_idx(linkidx);
  }
  //
  float_t get_Sw_at_link_custom_dx(int linkidx, float_t dx) {
    int n1 = this->connector->get_from_links(linkidx);
    int n2 = this->connector->get_to_links(linkidx);
    return (this->hw[n1] - this->hw[n2] + this->topography[n1] -
            this->topography[n2]) /
           dx;
  }

  float_t get_hflow_at_link(int linkidx) {
    int n1 = this->connector->get_from_links(linkidx);
    int n2 = this->connector->get_to_links(linkidx);
    return this->get_hflow_from_nodes(n1, n2);
  }

  float_t get_hflow_from_nodes(int n1, int n2) {
    return (std::max(this->hw[n1] + this->topography[n1],
                     this->hw[n2] + this->topography[n2]) -
            std::max(this->topography[n1], this->topography[n2]));
  }

  float_t get_surface_at_node(int n1) {
    return this->hw[n1] + this->topography[n1];
  }

  void fill_topo() {
    std::vector<float_t> surfpp(this->topography);
    auto savep = this->graph->depression_resolver;

    this->graph->depression_resolver = DAGGER::DEPRES::cordonnier_fill;
    this->graph->_compute_graph(surfpp, false, true);
    this->graph->depression_resolver = savep;

    for (int i = 0; i < this->connector->nnodes; ++i) {

      if (surfpp[i] > this->topography[i]) {
        float_t delat = surfpp[i] - this->topography[i];

        this->hw[i] -= delat;
        this->topography[i] += delat;
        if (this->hw[i] < 0)
          this->hw[i] = 0;
      }
    }
  }

  // {

  // 	this->init_erosion();

  // 	// std::cout << "DEBUGDEBUGA" << std::endl;
  // 	std::vector<float_t> surf(this->connector->nnodes,0);

  // 	// float_t sumgrad_sqrt = 0, sumHW = 0;

  // 	for(int i=0; i< this->connector->nnodes; ++i)
  // 		surf[i] = this->topography[i] + this->hw[i];

  // 	// std::cout << "DEBUGDEBUGA1" << std::endl;

  // 	// std::vector<float_t> surfpp(surf);
  // 	std::vector<float_t> surfpp(surf);
  // 	// std::cout << "DEBUGDEBUGA11" << std::endl;
  // 	(*this->graph)._compute_graph(surfpp, *(this->connector), true, true);
  // 	// py::array temp  = this->graph->compute_graph<Connector_t,
  // std::vector<double> , py::array>(depression_solver, surfpp,
  // *(this->connector), true, true);
  // 	// std::cout << "DEBUGDEBUGA2" << std::endl;
  // 	// float_t totaddedbystuff = 0;
  // 	for(int i=0; i<this->connector->nnodes; ++i)
  // 	{

  // 		if(surfpp[i] > surf[i])
  // 		{
  // 			this->hw[i] += surfpp[i] - surf[i];
  // 			// totaddedbystuff += surfpp[i] - surf[i];
  // 		}

  // 	}
  // 	// std::cout << "Tot added by filling::" << totaddedbystuff <<
  // std::endl;
  // 	// std::cout << "DEBUGDEBUGA3" << std::endl;

  // 	this->Qwin =
  // this->graph->_accumulate_variable_downstream_SFD(*this->connector,
  // this->Qbase); 	std::vector<float_t> gradients =
  // this->graph->_get_SFD_gradient(surfpp);

  // 	// float_t sumgrad2 = 0;
  // 	// for(auto v:gradients)
  // 	// {

  // 	// 	sumgrad2+= std::sqrt(std::max(v,1e-8));
  // 	// }
  // 	// std::cout << "sumgradient = " << sumgrad2 << std::endl;

  // 	// std::cout << "DEBUGDEBUGA4" << std::endl;

  // 	for(auto&v:this->Qwout)
  // 		v=0;

  // 	// std::cout << "DEBUGDEBUGA5" << std::endl;
  // 	int nneg = 0;
  // 	std::vector<bool> done(this->connector->nnodes,false);
  // 	for(int ti = this->connector->nnodes -1; ti >= 0; --ti)
  // 	{
  // 		int i = int(this->graph->Sstack[ti]);

  // 		if(gradients[i] < 0)
  // 		{
  // 			++nneg;
  // 			gradients[i] = 0;
  // 		}
  // 		// sumgrad_sqrt += std::sqrt(gradients[i]);

  // 		if(!this->graph->flow_out_model(i,*this->connector))
  // 		{

  // 			this->Qwout[i] =
  // this->graph->Sdistance2receivers[i]/this->mannings *
  // std::pow(this->hw[i],5./3) * std::pow(gradients[i],1./2);
  // this->hw[i] += dt
  // * ( this->Qwin[i] - this->Qwout[i])/this->connector->get_area_at_node(i);

  // 			if(this->hw[i] < 0)
  // 				this->hw[i] = 0;

  // 			if(this->Qwout[i] > 0)
  // 			{
  // 				float_t E = K * std::pow((this->hw[i] +
  // this->topography[i])/this->graph->Sdistance2receivers[i], n ) *
  // std::pow(this->Qwout[i], m ); 				float_t D =
  // this->qs[i]/(L
  // *
  // this->Qwout[i]); 				this->qs[i] += E - D;
  // if (this->qs[i] < 0) this->qs[i] = 0;
  // this->topography[i] += (D - E) * dt;
  // 			}
  // 		}
  // 	}
  // 	// std::cout << "man::" << 1/mannings << " hwpow::" << sumHW << "
  // grad::" << sumgrad_sqrt << " sumgrad2 " << sumgrad2 << " nmneg = " <<
  // nneg<< std::endl;
  // 	// std::cout << "DEBUGDEBUGA5" << std::endl;

  // }

  inline int random_node_index() {
    return std::floor(this->randu.get() * this->graph->nnodes);
  }

  void basicFloodos(int n_prec, float_t base_dt, float_t dV, int n_noout) {
    // vector tracking the last time passage of a single precipiton
    std::vector<float_t> last_t_precipitons(this->graph->nnodes, 0.);

    // Debug thingy
    int n_breaks = 0;

    // neighbour/rec/link fectcher and gradient trackers
    auto neighbours_links = this->connector->get_empty_neighbour();
    std::vector<int> reclink(neighbours_links.size(), -1);
    std::vector<float_t> gradhw(reclink.size(), 0),
        sqrtgradhw(reclink.size(), 0);

    // Main loop "launching" precipitons
    for (int i_prec = 0; i_prec < n_prec; ++i_prec) {

      // starting index at a random position
      int i = this->random_node_index();
      // checking if the precipiton is valid
      if (!this->connector->flow_out_model(i) == false) {
        // "Cancelling" the iteration
        --i_prec;
        continue;
      }

      // getting precipiton's time
      float_t this_time = base_dt * (i_prec + 1);

      // decout, ignore
      // std::cout << i << "/" << this->graph->nnodes << std::endl;

      // Walking the precipitons while in the landscape
      // -> tracking the precipiton's path
      int nwalk = 0;
      // -> running
      while (!this->connector->flow_out_model(i)) {
        // Recording the number of pixel
        nwalk++;
        // std::cout << i << "/" << this->graph->nnodes << "|";
        // Geting the neighbouring links:
        int nl = this->connector->get_neighbour_idx_links(i, neighbours_links);

        // Tracking the receivers link
        int nr = 0; // number of recs

        // The following variables (with absolutely not instinctive names) are
        // storing the information of the selected link for the precipiton's
        // path xam is the gradient, xamselect is the stochastic squareroot
        // gradient to use for selecting the right gradient
        float_t xam = 0, xamselect = 0;
        // xi is the link index of the selected link
        int xi = -1;
        // node index of the selected receiver (other side of the link)
        int node_xi = -1;

        // Calculating the effective timestep, corresponding to the time
        // difference from the last precipiton path
        float_t dt = this_time - last_t_precipitons[i];

        // Tracking the last precipiton's time
        last_t_precipitons[i] = this_time;

        // checking every neighbours
        for (int tn = 0; tn < nl; ++tn) {
          // current link
          int li = neighbours_links[tn];
          // other node (ie the node that is not the current one yo)
          int ri = (this->connector->linknodes[li * 2] == i)
                       ? this->connector->linknodes[li * 2 + 1]
                       : this->connector->linknodes[li * 2];
          // if hydrologic gradient > 0
          if (this->topography[i] + this->hw[i] >
              this->topography[ri] + this->hw[ri]) {
            // saving the link index as receiver
            reclink[nr] = li;
            // calculating the gradient
            gradhw[nr] = (this->topography[i] + this->hw[i] -
                          this->topography[ri] - this->hw[ri]) /
                         this->connector->get_dx_from_links_idx(li);
            // calculating the stochastic gradient
            sqrtgradhw[nr] =
                std::max(0., std::sqrt(gradhw[nr]) * this->randu.get() * 8);
            // selecting the current neighbour as receiver to use for the walk
            // if the stochastic gradient is higher
            if (sqrtgradhw[nr] > xamselect) {
              // stochastic gradient
              xamselect = sqrtgradhw[nr];
              // normal gradient
              xam = gradhw[nr];
              // link index
              xi = li;
              // recnode index
              node_xi = ri;
            }
            // incrementing the number of recs
            ++nr;
          }
          // end of loop checking neighbours
        }

        // decout
        // std::cout << i << "|" << node_xi << "|" << dV << " ||| ";

        // if I do have some receivers
        if (nr > 0) {
          // calculating the discharge using mannings
          double tQout =
              this->connector->get_dx_from_links_idx(xi) * 1. / this->mannings *
              std::pow(this->hw[i], double(5. / 3.)) * std::sqrt(xam);
          // ignore, speed up technique
          if (i_prec < n_noout) {
            // if speed up activated, only increase and ignore discharge
            this->hw[i] += (dV) / this->connector->get_area_at_node(i);
            // this->hw[i] += 1;
          } else {
            // the new hw is equal to the divergence of the fluxes (expressed as
            // volume there)
            this->hw[i] +=
                (dV - tQout * dt) / this->connector->get_area_at_node(i);
          }
          // removing negative discharges
          this->hw[i] = std::max(0., this->hw[i]);
        } else {
          // No recs?
          // increment with precipitons and leave
          this->hw[i] += dV / this->connector->get_area_at_node(i);
          this->hw[i] = std::max(0., this->hw[i]);
          ++n_breaks;
          break;
        }
        // next node is the receiving node
        i = node_xi;
      }
      // std::cout   << std::endl;
    }
  }

  void basicFloodos_v2(int n_prec, float_t base_dt, float_t dV, int n_noout) {
    // vector tracking the last time passage of a single precipiton
    std::vector<float_t> last_t_precipitons(this->graph->nnodes, 0.);
    auto surfpp = this->get_surface();
    this->graph->_compute_graph(surfpp, true, true);
    std::vector<float_t> last_S = this->connector->_get_SFD_gradient(surfpp);
    // std::vector<float_t> last_S(surfpp.size(),0.);
    std::vector<float_t> last_dx = this->connector->Sdistance2receivers;
    surfpp.clear();

    // Debug thingy
    int n_breaks = 0;

    // neighbour/rec/link fectcher and gradient trackers
    auto neighbours_links = this->connector->get_empty_neighbour();
    std::vector<int> reclink(neighbours_links.size(), -1);
    std::vector<float_t> gradhw(reclink.size(), 0),
        sqrtgradhw(reclink.size(), 0);

    // Main loop "launching" precipitons
    for (int i_prec = 0; i_prec < n_prec; ++i_prec) {

      // starting index at a random position
      int i = this->random_node_index();
      // checking if the precipiton is valid
      if (!this->connector->flow_out_model(i) == false) {
        // "Cancelling" the iteration
        --i_prec;
        continue;
      }

      // getting precipiton's time
      float_t this_time = base_dt * (i_prec + 1);

      // decout, ignore
      // std::cout << i << "/" << this->graph->nnodes << std::endl;

      // Walking the precipitons while in the landscape
      // -> tracking the precipiton's path
      int nwalk = 0;
      // -> running
      float_t MANNING_COEF = 0.6666666666666666666;
      float_t MANNING_COEF_INV = 1.5;
      while (!this->connector->flow_out_model(i)) {
        // Recording the number of pixel
        nwalk++;
        // std::cout << i << "/" << this->graph->nnodes << "|";
        // Geting the neighbouring links:
        int nl = this->connector->get_neighbour_idx_links(i, neighbours_links);

        // first updating the hw on all these links
        for (int tn = 0; tn < nl; ++tn) {

          // current link
          int li = neighbours_links[tn];
          // other node (ie the node that is not the current one yo)
          int ri = (this->connector->linknodes[li * 2] == i)
                       ? this->connector->linknodes[li * 2 + 1]
                       : this->connector->linknodes[li * 2];
          // this->hw[ri] = std::max(std::sqrt(last_S[ri]) *
          // this->mannings/this->connector->get_dx_from_links_idx(li), 0.1);
          float_t dx = this->connector->get_dx_from_links_idx(li);

          float_t tdt = this_time - last_t_precipitons[ri];
          last_t_precipitons[ri] = this_time;
          float_t W_coef =
              std::max(MANNING_COEF * std::sqrt(std::abs(last_S[ri])) *
                           (1 / this->mannings) / dx,
                       1e-6);

          this->hw[ri] =
              this->hw[ri] *
              std::pow(1.0 +
                           W_coef * std::pow(this->hw[ri], MANNING_COEF) * tdt,
                       -MANNING_COEF_INV);
          if (std::isfinite(this->hw[ri]) == false) {

            std::cout << W_coef << "," << last_S[ri] << "," << tdt << ","
                      << std::endl;
            std::cout << "sdafgsff::" << dx << std::endl;
            throw std::runtime_error("Bawulf");
          }
        }

        float_t tdt = this_time - last_t_precipitons[i];
        last_t_precipitons[i] = this_time;
        float_t dx = last_dx[i];
        float_t W_coef =
            std::max(MANNING_COEF * std::sqrt(std::abs(last_S[i])) *
                         (1 / this->mannings) / dx,
                     1e-6);

        this->hw[i] =
            this->hw[i] *
            std::pow(1.0 + W_coef * std::pow(this->hw[i], MANNING_COEF) * tdt,
                     -MANNING_COEF_INV);

        if (std::isfinite(this->hw[i]) == false) {

          std::cout << W_coef << "," << tdt << "," << std::endl;
          std::cout << "sdaf::" << dx << std::endl;
          throw std::runtime_error("Bawulf");
        }

        // Tracking the receivers link
        int nr = 0; // number of recs

        // The following variables (with absolutely not instinctive names) are
        // storing the information of the selected link for the precipiton's
        // path xam is the gradient, xamselect is the stochastic squareroot
        // gradient to use for selecting the right gradient
        float_t xam = 0, xamselect = 0;
        // xi is the link index of the selected link
        // int xi = -1;

        float_t dxi = 0;
        // node index of the selected receiver (other side of the link)
        int node_xi = -1;

        // Calculating the effective timestep, corresponding to the time
        // difference from the last precipiton path

        // Tracking the last precipiton's time
        // last_t_precipitons[i] = this_time;

        // checking every neighbours
        for (int tn = 0; tn < nl; ++tn) {
          // current link
          int li = neighbours_links[tn];
          // other node (ie the node that is not the current one yo)
          int ri = (this->connector->linknodes[li * 2] == i)
                       ? this->connector->linknodes[li * 2 + 1]
                       : this->connector->linknodes[li * 2];

          float_t tgrad = std::abs((this->topography[i] + this->hw[i] -
                                    this->topography[ri] - this->hw[ri]) /
                                   this->connector->get_dx_from_links_idx(li));
          // float_t coeff = std::max(std::sqrt(xam) *
          // this->mannings/this->connector->get_dx_from_links_idx(xi), 0.1);
          // float_t exp = 0.66666;
          // float_t invexp = 1.5;
          // this->hw[ri] =
          // this->hw[ri]*std::pow(1+coeff*std::pow(this->hw[ri],exp) *
          // tdt,-invexp);

          // if hydrologic gradient > 0
          if (this->topography[i] + this->hw[i] >
              this->topography[ri] + this->hw[ri]) {
            // saving the link index as receiver
            reclink[nr] = li;
            // calculating the gradient
            gradhw[nr] = tgrad;
            // calculating the stochastic gradient
            sqrtgradhw[nr] =
                std::max(0., std::sqrt(gradhw[nr]) * this->randu.get() * 1e-3);
            // selecting the current neighbour as receiver to use for the walk
            // if the stochastic gradient is higher
            if (sqrtgradhw[nr] > xamselect) {
              // stochastic gradient
              xamselect = sqrtgradhw[nr];
              // normal gradient
              xam = gradhw[nr];
              // link index
              // xi = li;
              dxi = this->connector->get_dx_from_links_idx(li);
              // recnode index
              node_xi = ri;
            }
            // incrementing the number of recs
            ++nr;
          }
          // end of loop checking neighbours
        }

        // decout
        // std::cout << i << "|" << node_xi << "|" << dV << " ||| ";

        // if I do have some receivers
        if (nr > 0) {
          last_S[i] = xam;
          last_dx[i] = dxi;
          // calculating the discharge using mannings
          // double tQout = this->connector->get_dx_from_links_idx(xi)
          // * 1./this->mannings * std::pow(this->hw[i],double(5./3.)) *
          // std::sqrt(xam); ignore, speed up technique
          if (i_prec < n_noout) {
            // if speed up activated, only increase and ignore discharge
            this->hw[i] += (dV) / this->connector->get_area_at_node(i);
            // this->hw[i] += 1;
          } else {
            // the new hw is equal to the divergence of the fluxes (expressed as
            // volume there) this->hw[i] += (dV - tQout* dt)
            // /this->connector->get_area_at_node(i); float_t coeff =
            // std::max(std::sqrt(last_S[i]) *
            // this->mannings/this->connector->get_dx_from_links_idx(xi), 0.1);
            // // float_t C = 25;
            // float_t exp = 0.66666;
            // float_t invexp = 1.5;
            // float_t tdx = this->connector->get_dx_from_links_idx(xi);
            // if(dt > 0)
            // 	this->hw[i] =
            // this->hw[i]*std::pow(1+coeff*std::pow(this->hw[i],exp) *
            // dt,-invexp); else
            this->hw[i] += (dV) / this->connector->get_area_at_node(i);
          }
          // removing negative discharges
          this->hw[i] = std::max(0., this->hw[i]);
        } else {
          // No recs?
          // increment with precipitons and leave
          this->hw[i] += dV / this->connector->get_area_at_node(i);
          this->hw[i] = std::max(0., this->hw[i]);
          ++n_breaks;
          break;
        }
        // last_t_increase[i] = this_time;
        // next node is the receiving node
        i = node_xi;
      }
      // std::cout   << std::endl;
    }
  }

  void basicFloodos_v3(int n_prec, float_t base_dt,
                       float_t precipitation_rates) {

    // init the variables if it has not been done yet
    if (this->lastdl.size() == 0) {
      this->lastdl = std::vector<float_t>(this->graph->nnodes, 0.);
      this->lastS = std::vector<float_t>(this->graph->nnodes, 0.);
      this->lastdt = std::vector<float_t>(this->graph->nnodes, 0.);
    }

    // Calculating the Volume of precipitons function of uniform precipitation
    // rates
    float_t dV = precipitation_rates * (base_dt) * this->graph->nnodes *
                 this->connector->get_area_at_node(0);
    float_t C = 1. / this->mannings, mannings_coeff = 2. / 3., invmann = -1.5;

    auto neighbours_links = this->connector->get_empty_neighbour();
    auto neighbours_idx = this->connector->get_empty_neighbour();
    std::vector<int> reclink(neighbours_links.size(), -1);
    std::vector<float_t> gradhw(reclink.size(), 0),
        sqrtgradhw(reclink.size(), 0);

    // int n_stop_middle = 0;
    for (int iprec = 1; iprec <= n_prec; ++iprec) {
      float_t time = maxdt + base_dt;
      this->maxdt = time;
      int i = this->random_node_index();

      if (!this->connector->flow_out_model(i) == false) {
        --iprec;
        continue;
      }

      while (!this->connector->flow_out_model(i)) {

        float_t dt = time - this->lastdt[i];
        if (this->lastdt[i] > 0 && this->lastdt[i] > 0 &&
            std::abs(lastS[i]) > 0)
          this->hw[i] =
              this->hw[i] *
              std::pow(1 + mannings_coeff * C / this->lastdl[i] *
                               std::sqrt(std::abs(lastS[i])) *
                               std::pow(this->hw[i], mannings_coeff) * dt,
                       invmann);

        int nl = this->connector->get_neighbour_idx_links(i, neighbours_links);

        float_t xam = 0, xamselect = 0;
        // xi is the link index of the selected link
        // int xi = -1;

        float_t dxi = 0;
        // node index of the selected receiver (other side of the link)
        int node_xi = -1;
        int nr = 0;

        // first updating the hw on all these links
        for (int tn = 0; tn < nl; ++tn) {
          // current link
          int li = neighbours_links[tn];
          // other node (ie the node that is not the current one yo)
          int ri = (this->connector->linknodes[li * 2] == i)
                       ? this->connector->linknodes[li * 2 + 1]
                       : this->connector->linknodes[li * 2];
          // this->hw[ri] = std::max(std::sqrt(last_S[ri]) *
          // this->mannings/this->connector->get_dx_from_links_idx(li), 0.1);
          float_t dx = this->connector->get_dx_from_links_idx(li);
          float_t tgrad = std::abs((this->topography[i] + this->hw[i] -
                                    this->topography[ri] - this->hw[ri]) /
                                   dx);

          // if hydrologic gradient > 0
          if (this->topography[i] + this->hw[i] >
              this->topography[ri] + this->hw[ri]) {
            // saving the link index as receiver
            reclink[nr] = li;
            // calculating the gradient
            gradhw[nr] = tgrad;
            // calculating the stochastic gradient
            sqrtgradhw[nr] =
                std::max(0., std::sqrt(gradhw[nr]) * this->randu.get());
            // selecting the current neighbour as receiver to use for the walk
            // if the stochastic gradient is higher
            if (sqrtgradhw[nr] > xamselect) {
              // stochastic gradient
              xamselect = sqrtgradhw[nr];
              // normal gradient
              xam = gradhw[nr];
              // link index
              // xi = li;
              dxi = dx;
              // recnode index
              node_xi = ri;
            }
            ++nr;
          }
        }

        this->lastdt[i] = time;

        this->hw[i] += dV / this->connector->get_area_at_node(i);

        if (node_xi != -1) {
          this->lastdl[i] = dxi;
          this->lastS[i] = xam;
          i = node_xi;
        } else {
          this->hw[i] += dV / this->connector->get_area_at_node(i);
          // n_stop_middle++;
          while (true) {
            int nn = this->connector->get_neighbour_idx(i, neighbours_idx);
            int maxsurf_i = i;
            float_t maxsurf = this->hw[i] + this->topography[i];
            bool hasrec = false;
            for (int tn = 0; tn < nn; ++tn) {
              int oi = neighbours_idx[tn];
              if (this->hw[oi] + this->topography[oi] >= maxsurf) {
                maxsurf = this->hw[oi] + this->topography[oi];
                maxsurf_i = oi;
              } else {
                hasrec = true;
                break;
              }
            }
            if (hasrec || i == maxsurf_i)
              break;
            this->hw[i] = this->hw[maxsurf_i] + 1e-3 +
                          this->topography[maxsurf_i] - this->topography[i];
            i = maxsurf_i;
          }
        }
      }
    }

    // std::cout << "N stop middle:: " << n_stop_middle << std::endl;
  }

  void set_mannings(float_t new_mannings) { this->mannings = new_mannings; }

  void testDebugWalk(int n_prec, float_t base_dt, float_t K, float_t m,
                     float_t n, float_t tolerance) {
    // last time passage of a single precipiton
    std::vector<float_t> last_t_precipitons(this->graph->nnodes, 0.);

    std::vector<int> stack(this->graph->nnodes, -1),
        recs(this->graph->nnodes, -1);
    std::vector<float_t> Dxs(this->graph->nnodes, 0.);
    std::vector<float_t> Area(this->graph->nnodes, 0.);

    auto neighbours_links = this->connector->get_empty_neighbour();
    std::vector<int> reclink(neighbours_links.size(), -1);
    std::vector<float_t> gradhw(reclink.size(), 0),
        sqrtgradhw(reclink.size(), 0);
    for (int i_prec = 0; i_prec < n_prec; ++i_prec) {
      // getting precipiton's time
      float_t this_time = base_dt * (i_prec + 1);
      // starting index
      int i = this->random_node_index();

      int nn_proc = 0;

      // Walking the precipitons while in the landscape
      float_t tA = 0;
      while (!this->connector->flow_out_model(i)) {
        int nl = this->connector->get_neighbour_idx_links(i, neighbours_links);
        int nr = 0;
        float_t xam = 0, xamselect = 0;
        int node_xi = -1;
        tA += this->connector->get_area_at_node(i);

        for (int tn = 0; tn < nl; ++tn) {
          int li = neighbours_links[tn];
          int ri = (this->connector->linknodes[li * 2] == i)
                       ? this->connector->linknodes[li * 2 + 1]
                       : this->connector->linknodes[li * 2];
          if (this->topography[i] > this->topography[ri]) {
            reclink[nr] = li;
            gradhw[nr] = (this->topography[i] - this->topography[ri]) /
                         this->connector->get_dx_from_links_idx(li);
            sqrtgradhw[nr] = std::max(0., std::sqrt(gradhw[nr]));
            if (sqrtgradhw[nr] > xamselect) {
              xamselect = sqrtgradhw[nr];
              xam = this->connector->get_dx_from_links_idx(li);
              node_xi = ri;
            }
            ++nr;
          }
        }
        if (nr > 0) {
          // double E = 3e-5 * std::pow(xam,1) *
          // std::pow(this->connector->get_area_at_node(i),0.45) * dt;
          // this->topography[i] -= E;
          stack[nn_proc] = i;
          recs[nn_proc] = node_xi;
          Dxs[nn_proc] = xam;
          Area[nn_proc] = tA;
          ++nn_proc;
        } else {
          break;
        }

        i = node_xi;
      }

      for (int j = nn_proc - 1; j >= 0; --j) {
        int node = stack[j];
        int rec = recs[j];
        float_t dx = Dxs[j];

        float_t dt = this_time - last_t_precipitons[node];
        last_t_precipitons[node] = this_time;

        float_t factor = K * dt * std::pow(Area[node], m) / std::pow(dx, n);
        ;

        float_t node_elevation = this->topography[node];
        float_t rec_elevation = this->topography[rec];

        float_t elevation_k = node_elevation;
        float_t elevation_prev = std::numeric_limits<float_t>::max();

        while (std::abs(elevation_k - elevation_prev) > tolerance) {
          elevation_prev = elevation_k;
          float_t slope = std::max(elevation_k - rec_elevation, 1e-6);
          float_t diff =
              ((elevation_k - node_elevation + factor * std::pow((slope), n)) /
               (1. + factor * n * std::pow(slope, (n - 1))));
          elevation_k -= diff;
        }
        this->topography[node] = elevation_k;
      }
    }
  }

  void check_SD_val() {

    std::vector<std::pair<int, float_t>> neighbours(8);
    for (int i = 0; i < this->graph->nnodes; ++i) {
      if (!this->connector->flow_out_model(i)) {
        int NN = this->connector->get_neighbour_idx_distance(i, neighbours);
        int steepest_i = i;
        float_t steepest_slope = 0;
        for (int j = 0; j < NN; ++j) {
          int ri = neighbours[i].first;
          float_t dri = neighbours[i].second;
          float_t tS = (this->topography[i] + this->hw[i] - this->hw[ri] -
                        this->topography[ri]) /
                       dri;
          if (tS > steepest_slope) {
            steepest_slope = tS;
            steepest_i = ri;
          } else if (tS == steepest_slope)
            throw std::runtime_error("EQUALITYYYYYY");
        }

        if (this->connector->_Sreceivers[i] != steepest_i)
          throw std::runtime_error("QQQQWAGYNIARD");
      }
    }
    std::cout << "All good with SFD checks." << std::endl;
  }

#ifdef OPENMP_YOLO
  void check_devices() {
    int A[1] = {-1};
#pragma omp target
    { A[0] = omp_is_initial_device(); }

    if (!A[0]) {
      std::cout << "Able to use offloading!\n" << std::endl;
      ;
    } else
      std::cout << "No offloading :(" << std::endl;
    ;
  }

  // Implementing Caesar-Lisflood flooding algorithm in DAGGER's logic
  // For the physics see paper from Bates et al., 2010 and the implementation is
  // adapted and modified from HAIL-CAESAR implementation by Declan Valters
  float_t caesar_lisflood_OMP() {

    // FIrst checking if the Qlinks has been initialised
    if (this->Qlink.size() == 0) {
      this->Qlink = std::vector<float_t>(this->connector->links.size(), 0.);
    }

    // Initializing the dhw
    std::vector<float_t> dhw(this->graph->nnodes, 0.);

    std::vector<float_t> surf(this->connector->nnodes, 0);

    for (int i = 0; i < this->connector->nnodes; ++i)
      surf[i] = this->topography[i] + this->hw[i];

    // calculationg the maximum timestep respecting the CFD
    float_t maxhw =
        std::max(*std::max_element(this->hw.begin(), this->hw.end()), 1e-3);

    float_t dt =
        this->connector->dxmin * this->alpha / std::sqrt(GRAVITY * maxhw);

    if (dt == 0)
      dt = 1e-4;

// main loop
#pragma omp parallel for
    for (size_t i = 0; i < this->Qlink.size(); ++i) {
      // valid link??
      if (this->connector->linknodes[i * 2] == -1)
        continue;

      // by convention n1 and n2 are in the arbitrary order of the linknode
      // array (not really important)
      int n1 = this->connector->linknodes[i * 2],
          n2 = this->connector->linknodes[i * 2 + 1];

      // if both hw are below 0 I do not compute and qlink is null
      if (this->hw[n1] <= 0 && this->hw[n2] <= 0) {
        this->Qlink[i] = 0;
        continue;
      }

      // calculating gradient
      float_t Sw =
          (surf[n1] - surf[n2]) / this->connector->get_dx_from_links_idx(i);

      // Getting hflow, the flow height between the two cell:
      float_t hflow = std::max(surf[n1], surf[n2]) -
                      std::max(this->topography[n1], this->topography[n2]);

      // no flow, no algo
      if (hflow <= 0) {
        this->Qlink[i] = 0;
        continue;
      }

      // Get the trasverse dx
      float_t tdx = this->topological_number *
                    this->connector->get_traverse_dx_from_links_idx(i);

      // Get the discharge
      this->Qlink[i] =
          tdx * (this->Qlink[i] - GRAVITY * hflow * dt * Sw) /
          (1 + (GRAVITY * hflow * dt * std::pow(this->mannings, 2) *
                std::abs(this->Qlink[i]) / std::pow(hflow, 10. / 3.)));

      int n_neighbours = (this->topological_number == 1.) ? 4 : 8;
      // int n_neighbours = 4;

      // Discahrge checkers
      if (std::abs(this->Qlink[i]) / hflow / std::sqrt(GRAVITY * hflow) >
          this->froude_number)
        this->Qlink[i] = DAGGER::sgn(this->Qlink[i]) * hflow *
                         std::sqrt(GRAVITY * hflow) * this->froude_number;

      if (std::abs(this->Qlink[i]) * dt / tdx > this->hw[n2] / n_neighbours &&
          DAGGER::sgn(this->Qlink[i]) > 0)
        this->Qlink[i] = this->hw[n2] * this->connector->get_area_at_node(n2) /
                         (n_neighbours + 1) / dt;
      else if (std::abs(this->Qlink[i]) * dt / tdx >
                   this->hw[n1] / n_neighbours &&
               DAGGER::sgn(this->Qlink[i]) < 0)
        this->Qlink[i] = -1 * this->hw[n1] *
                         this->connector->get_area_at_node(n1) /
                         (n_neighbours + 1) / dt;
    }

    for (size_t i = 0; i < this->Qlink.size(); ++i) {
      // valid link??
      if (this->connector->linknodes[i * 2] == -1)
        continue;

      // by convention n1 and n2 are in the arbitrary order of the linknode
      // array (not really important)
      int n1 = this->connector->linknodes[i * 2],
          n2 = this->connector->linknodes[i * 2 + 1];
      dhw[n1] += dt * this->Qlink[i];
      dhw[n2] -= dt * this->Qlink[i];
    }

    for (int i = 0; i < this->graph->nnodes; ++i) {
      if (!this->connector->flow_out_model(i)) {
        this->hw[i] += (dhw[i]) / this->connector->get_area_at_node(i);
        if (this->hw[i] < 0)
          this->hw[i] = 0;
      }
    }

    for (int i = 0; i < this->graph->nnodes; ++i) {
      if (!this->connector->flow_out_model(i))
        this->hw[i] +=
            (this->Qbase[i] * dt) / this->connector->get_area_at_node(i);
    }

    return dt;
  }
#endif
};

// End of namespace fastflood
}; // namespace DAGGER

#endif
