//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef DEPHIERCH_HPP
#define DEPHIERCH_HPP

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
#include <unordered_map>
#include <vector>

// local includes
// -> General routines and data structures
#include "utils.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

namespace DAGGER {

// Useful namespace for my priority queue
using PQ_i_d = std::priority_queue<PQ_helper<int, double>,
                                   std::vector<PQ_helper<int, double>>,
                                   std::greater<PQ_helper<int, double>>>;
using PQH = PQ_helper<int, double>;

enum class DEPSTATE : std::uint8_t {
  NODEP,
  BUILDING,
  CHILD,
  OPEN,
};

template <class fT> class deptree {
public:
  // number of subdepressions in the tree
  int ndeps = 0;
  // Global index of the tree
  int index;

  // is a branch (subtree) or full
  bool branch = false;
  int motherbranch = -1;

  // level of the tree: increases dynamically during publication to reach its
  // max level
  fT level = std::numeric_limits<fT>::min();
};

// OLDAER TEST, UNFINISHED OR DO NOT WORK

// template<class float_t>
// class depression
// {
// public:
// 	int pitnode;
// 	int index;
// 	std::vector<int> inoutlets;
// 	int outlinks;
// 	float_t outscore = std::numeric_limits<float_t>::max();

// 	depression(){};
// 	depression(int pitnode, int index)
// 	{
// 		this->pitnode = pitnode;
// 		this->index = index;
// 	}

// 	// template<class topo_t>
// 	void internal_link(int outlet)
// 	{
// 		bool valid = true;
// 		// for(auto node:this->outlets)
// 		// {
// 		// 	if(topography[outlet] > topography[node])
// 		// 		valid = false;
// 		// 	else if (topography[outlet] < topography[node])
// 		// 		throw std::runtime_error("Depression Hierarchy
// exception unhandled #342vf8");
// 		// }
// 		// if(valid)
// 		// 	this->outlet.emplace_back(outlet);

// 		for(auto node:this->outlets)
// 		{
// 			if(outlet == node)
// 				valid = false;
// 		}
// 		if(valid)
// 			this->inoutlets.emplace_back(outlet);
// 	}

// 	// template<class topo_t>
// 	void external_link(int outlet, float_t score)
// 	{
// 		if(this->outscore > score)
// 		{
// 			this->outlink = outlet;
// 			this->outscore = score;
// 		}
// 	}

// };

// void depression_hierarchy_v0(Connector_t* connector, topo_t& topography)
// {
// 	// Setting up the original priority queue
// 	// for depression labelling
// 	PQ_i_d tpq;

// 	// Main depression index
// 	int lab = 0;

// 	// depression holder
// 	std::vector<depression> quik;
// 	std::vector<std::uint8_t> state;
// 	quik.reserve(200);
// 	state.reserve(200);

// 	//
// 	std::vector<std::uint8_t> vis(connector->nnodes,false);

// 	// keeping track of depression labels
// 	std::vector<int> depression_labels(connector->nnodes, -2);

// 	//Step 1: label ocean-connected cells in steepest descent
// 	for(int i=0; i<connector.nodes; ++i)
// 	{
// 		int node = connector.Sstack[i];
// 		if(connector.boundaries.no_data(node)) continue;
// 		int rec = connector._Sreceivers[node];
// 		if(node == rec)
// 		{
// 			if(connector.flow_out_model(node))
// 			{
// 				depression_labels[node] = -1;
// 			}
// 			else
// 			{
// 				tpq.emplace_back( PQH(node, topography[node]) );
// 			}
// 			vis[node] = true;
// 		}
// 		else
// 		{
// 			if(depression_labels[rec] == -1)
// 			{
// 				depression_labels[node] = -1;
// 			}
// 		}
// 	}

// 	// Step 2 add the MFD cells

// 	// Step 2.1: feeding the original queue with all the nodes to reconnect
// 	std::stack<int, std::vector<int> > mfdlab;
// 	int n1,n2;
// 	for(int i=0; i<connector->nlinks(); ++i)
// 	{
// 		if(connector->is_link_valid(i))
// 		{
// 			connector->node_idx_from_link_idx_nocheck(i,n1,n2);
// 			if(depression_labels[n1] != depression_labels[n2])
// 			{
// 				if(depression_labels[n1] == -1 && topography[n2]
// >= topography[n1])
// 				{
// 					mfdlab.emplace(n2);
// 					depression_labels[n2] = -1;
// 				}
// 				else if(depression_labels[n2] == -1 &&
// topography[n1]
// >= topography[n2])
// 				{
// 					mfdlab.emplace(n1);
// 					depression_labels[n1] = -1;
// 				}
// 			}
// 		}
// 	}

// 	// Step 2.2 label donors
// 	auto neighbourer = connector->get_empty_neighbour();
// 	while(mfdlab.empty() == false)
// 	{
// 		int node = mfdlab.top();
// 		mfdlab.pop();
// 		int nn = connector->get_donors_idx(node, neighbours);
// 		for(int j =0; j< nn; ++j )
// 		{
// 			if(depression_labels[neighbours[j]] != -1)
// 			{
// 				depression_labels[neighbours[j]] = -1;
// 				mfdlab.emplace(neighbours[j]);
// 			}
// 		}
// 	}

// 	// At tha stage all the ocean cells should be labelled
// 	// and the depressions in the priority Queue

// 	// Step 3: label all the depressions
// 	while(tpq.empty() == false)
// 	{
// 		// Getting the next node
// 		// in line
// 		auto next = tpw.top();
// 		tpq.pop();

// 		// getting the neighbours
// 		int nn = connector->get_neighbour_idx_links(next.node,
// neighbours);

// 		// Is the cell already part of a depression?
// 		if(depression_labels[next.node] == -2)
// 		{
// 			// if not I label it to new depressions (depression ID
// is the node id of the pit) quik.emplace_back(depression(next.node,lab));
// state.emplace_back(0);

// 			// Labelling the depression
// 			depression_labels[nextnode] = lab;

// 			// incrementing the global depression label
// 			++lab;
// 		}
// 		else if(depression_labels[next.node] == -1)
// 		{
// 			// case where cell is connected to the ocean
// 			for(int j =0; j<nn; ++j)
// 			{
// 				int lix = neighbours[j];
// 				int neigh =
// connector->get_other_node_from_links(lix,next.node);
// 				if(depression_labels[neigh] < 0) continue;
// 				// if a neighbour is a depression, this node is
// a potential outlet 				bool valid = false;
// if(topography[neigh] >=
// topography[next.node]) 					valid =
// quik[depression_labels[neigh]].external_link(lix, topography[neigh]);
// 				if(valid) state[depression_labels[neigh]] =
// true;
// 			}

// 			// Stopping the process here;
// 			continue;
// 		}
// 		else
// 		{
// 			// If the depression is already connected to the ocean,
// I abort the process 			if(state[depression_labels[next.node]])
// 			{
// 				depression_labels[next.node] = -1;
// 			}

// 			// I continue the loop (not pushing the neighbours)
// 			continue;
// 		}

// 		// Looping through cell's neighbours
// 		for(int j =0; j<nn; ++j)
// 		{
// 			int neigh = neighbours[j];

// 			// is the neighbour in a depression?
// 			if(depression_labels[neigh] < 0 && vis[neigh] == false)
// 			{
// 				vis[neigh] = true;
// 				tpq.emplace_back(PQH(neigh,topography[neigh]));
// 			}
// 		}
// 	}

// }

template <class Connector_t>
int reroute(Connector_t *connector,
            // std::vector<int>& baslab,
            std::vector<int> &_Sreceivers, int from, int pass) {
  // std::cout << "rerouting" << std::endl;
  int A = from;
  int B = _Sreceivers[A];
  bool goon = true;
  while (goon) {
    int C = _Sreceivers[B];
    if (C == B)
      goon = false;
    _Sreceivers[B] = A;
    A = B;
    B = C;
  }

  // int i=0;
  // do
  // {
  // 	// ++i;
  // 	int C = _Sreceivers[B];
  // 	_Sreceivers[B] = A;
  // 	connector->debug_print_row_col(B);
  // 	// std::cout << " now gives to "; connector->debug_print_row_col(A);
  // std::cout << "||"; 	A = B; 	B = C;
  // 	// if(i >1000)
  // 		// std::cout << A  << "|" << B << std::endl;
  // }while(_Sreceivers[B] != B);

  if (pass != -1)
    _Sreceivers[from] = pass;
  else
    _Sreceivers[pass] = pass;

  return B;

  // if(passout) _Sreceivers[pass] = pass;
}

template <class float_t, class topo_t, class Connector_t>
bool simple_depression_hierarchy(topo_t &topography, Connector_t *connector,
                                 std::vector<size_t> &Sstack,
                                 std::vector<int> &_Sreceivers) {

  // PROBLEM IDENTIFIED:: WHEN REROUTING THE SF, THE SFD DOES NOT NECESSARILY
  // DRAINS TO THE BASIN, ITSELF CALCULATED IN MFD LIKE

  // std::cout << "running?" << std::endl;

  // record the basin labels
  std::vector<int> baslab(connector->nnodes, -1);

  // the main priority queue
  PQ_i_d PQ;

  // std::unordered_map<int,int> pit2outlet;

  // Labelling all the nodes draining to the sea, pushing the local pits to the
  // priority queue
  int lab = 1;
  for (int i = 0; i < connector->nnodes; ++i) {
    int node = int(Sstack[i]);

    if (connector->boundaries.no_data(node))
      continue;

    int rec = _Sreceivers[node];

    if (node == rec) {

      if (connector->flow_out_model(node)) {
        // std::cout << "OUT!!! :: ";
        // connector->debug_print_row_col(node);
        baslab[node] = 0;
      } else {
        baslab[node] = lab;
        PQ.emplace(PQH(node, topography[node]));
        std::cout << "EMPLACING !!! :: ";
        connector->debug_print_row_col(node);
        ++lab;
      }
    }

    if (baslab[rec] == 0)
      baslab[node] = 0;
  }

  if (PQ.empty())
    return false;

  std::vector<std::uint8_t> basopen(lab, 0);
  basopen[0] = true;
  std::cout << "BASOPEN::" << std::to_string(basopen.size()) << std::endl;

  std::vector<std::vector<std::array<int, 2>>> bas2passes(
      lab, std::vector<std::array<int, 2>>());

  // std::cout << "PQ size = " << PQ.size() << std::endl;;

  auto neighbours = connector->get_empty_neighbour();
  // int yolo = 0;
  while (PQ.empty() == false) {
    // std::cout << PQ.size() << std::endl;
    // getting the next node
    auto next = PQ.top();
    PQ.pop();

    int tbas = baslab[next.node];

    // std::cout << std::endl << "processing ";
    // connector->debug_print_row_col(next.node);
    // std::cout << ":: basin " << tbas << " (out:" <<
    // std::to_string(basopen[tbas]) << ")"; if(tbas == 0) 	throw
    // std::runtime_error("should not happen");

    if (basopen[tbas] == true) {
      // std::cout << " pass!";
      continue;
    }

    if (connector->boundaries.can_out(next.node)) {
      // std::cout << "A" << std::endl;
      // std::cout << " reroute from edges";
      reroute(connector, _Sreceivers, next.node, -1);
      // std::cout << " done";
      // std::cout << "B" << std::endl;
      basopen[tbas] = true;
      continue;
    }

    int nn = connector->get_neighbour_idx(next.node, neighbours);

    bool does_it_outlets = false;
    PQH node_outlet(-1, topography[next.node]);
    int restor_Srec = -1;
    // if(tbas == 2)
    // 	std::cout << "Processing node ";
    // connector->debug_print_row_col(next.node);

    for (int j = 0; j < nn; ++j) {
      int nei = neighbours[j];
      int bnei = baslab[nei];

      // std::cout << tbas << "|" << bnei << std::endl;

      // if(basopen[bnei] == true || (topography[nei] < node_outlet.score &&
      // bnei != tbas))
      if (bnei >= 0) {
        if (bnei != tbas) {

          does_it_outlets = true;
          node_outlet.node = nei;
          node_outlet.score = topography[nei];
          restor_Srec = _Sreceivers[nei];
        }
      }

      if (bnei == -1) {
        PQ.emplace(PQH(nei, topography[nei]));
        baslab[nei] = tbas;
        bnei = tbas;
        _Sreceivers[nei] = next.node;
      }
    }

    // if(does_it_outlets && tbas ==2)
    // {
    // 	std::cout << "Processed " << next.node << " | outlet: " <<
    // node_outlet.node << std::endl;
    // }

    // found an outlet, rereouting
    if (does_it_outlets) {
      // std::cout << " potential outlet. ";

      int bnei = baslab[node_outlet.node];
      _Sreceivers[node_outlet.node] = restor_Srec;
      // if(bnei == 0)
      // 	throw

      if (basopen[bnei] == true) {
        // std::cout << "C::" << "|" <<  node_outlet.node << "|" << next.node <<
        // "|" << basopen[bnei] << "|" << std::to_string(bnei == tbas)<<
        // std::endl; std::cout << " rerouting: link is " ;
        // connector->debug_print_row_col(node_outlet.node); std::cout << " to
        // "; connector->debug_print_row_col(next.node);
        int lnode =
            reroute(connector, _Sreceivers, next.node, node_outlet.node);
        // std::cout << "D" << std::endl;
        // basopen[tbas] = true;
        std::queue<int> basb;

        basb.emplace(tbas);

        while (basb.empty() == false) {
          int nextbas = basb.front();
          basb.pop();
          basopen[nextbas] = true;

          for (auto &link : bas2passes[nextbas]) {
            // if(baslab[link[0]] == baslab[link[1]])
            // 	std::cout << "HAPPENS::" << baslab[link[0]]  << "|" <<
            // baslab[link[1]]  << std::endl;

            if (basopen[baslab[link[1]]] == true)
              continue;

            // std::cout << "C1::" << "|" <<  baslab[link[0]] << "|" <<
            // next.node << std::endl;
            lnode = reroute(connector, _Sreceivers, link[1], link[0]);
            // std::cout << "D" << std::endl;
            basb.emplace(baslab[lnode]);
          }
        }

      } else {
        // if(tbas == bnei)
        // 	throw std::runtime_error("???");

        // if(baslab[next.node] == baslab[node_outlet.node])
        // 	throw std::runtime_error("???");

        std::cout << "Link::" << tbas << "|" << bnei << std::endl;

        bas2passes[tbas].emplace_back(
            std::array<int, 2>{next.node, node_outlet.node});
        bas2passes[bnei].emplace_back(
            std::array<int, 2>{node_outlet.node, next.node});
      }
    }

    // for(int j = 0; j<nn; ++j)
    // {
    // 	int nei = neighbours[j];

    // }

    // double continue;
    // cnt:;
  }

  int nopen = 0, nclose = 0;
  for (int i = 0; i < basopen.size(); ++i) {
    if (basopen[i] == false) {
      ++nclose;
      std::cout << "UNREROUTED BASIN ::" << i << std::endl;
    } else {
      ++nopen;
    }
  }

  std::cout << "yolo::" << nclose << " vs " << nopen << std::endl;
  return true;
};

} // namespace DAGGER

#endif
