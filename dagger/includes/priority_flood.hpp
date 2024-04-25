//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef priority_Q_HPP
#define priority_Q_HPP

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

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

// #ifdef BOOST_AVAILABLE
// #include "boost/heap/priority_queue.hpp"
// #include "boost/heap/binomial_heap.hpp"
// #include "boost/heap/fibonacci_heap.hpp"
// #include "boost/heap/pairing_heap.hpp"
// #include "boost/heap/skew_heap.hpp"
// #endif

namespace DAGGER {

// Useful namespace for my priority queue
using PQFORPF = std::priority_queue<PQ_helper<int, double>,
                                    std::vector<PQ_helper<int, double>>,
                                    std::greater<PQ_helper<int, double>>>;
using PQH = PQ_helper<int, double>;
// using PQFORPF = boost::heap::priority_queue<PQH,
// boost::heap::compare<std::greater<PQH> > >; using PQFORPF =
// boost::heap::binomial_heap<PQH, boost::heap::compare<std::greater<PQH> > >;
// using PQFORPF = boost::heap::fibonacci_heap<PQH,
// boost::heap::compare<std::greater<PQH> > >; using PQFORPF =
// boost::heap::pairing_heap<PQH, boost::heap::compare<std::greater<PQH> > >;
// using PQFORPF = boost::heap::skew_heap<PQH,
// boost::heap::compare<std::greater<PQH> > >;

template <class Connector_t, class int_t, class out_t, class fT>
out_t standalone_priority_flood(int_t &topography, Connector_t &connector) {
  auto ttopo = format_input(topography);
  auto vtopo = DAGGER::to_vec(ttopo);
  PriorityFlood(vtopo, connector);
  return DAGGER::format_output<std::vector<fT>, out_t>(vtopo);
}

template <class fT, class Connector_t>
void PriorityFlood(std::vector<fT> &topography, Connector_t &connector) {
  PQFORPF open;
  std::queue<PQH> pit;

  std::vector<int8_t> closed(connector.nnodes, false);

  for (int i = 0; i < connector.nnodes; ++i) {
    if (connector.boundaries.can_out(i) == false)
      continue;

    open.emplace(PQH(i, topography[i]));
    closed[i] = true;
  }

  auto neighbours = connector.get_empty_neighbour();
  while (open.size() > 0 || pit.size() > 0) {
    PQH c;
    // if(pit.size()>0 && open.size()>0 && open.top().score ==
    // pit.front().score)
    // {
    // 	c=open.top();
    // 	open.pop();
    // 	PitTop=-9999;
    // }
    // else if(pit.size()>0)
    // {
    // 	c=pit.front();
    // 	pit.pop();
    // 	if(PitTop==-9999)
    // 		PitTop=topography[c.node];
    // // } else {
    c = open.top();
    open.pop();
    // fT PitTop = topography[c.node];
    // }

    int nn = connector.get_neighbour_idx(c.node, neighbours);
    fT ttopo = topography[c.node];

    for (int j = 0; j < nn; ++j) {
      int n = neighbours[j];

      if (connector.is_in_bound(n) == false ||
          connector.boundaries.can_create_link(n) == false)
        continue;

      if (closed[n])
        continue;

      closed[n] = true;
      ttopo += connector.randu->get() * 1e-6 + 1e-8;
      // fT ntopo = ;
      topography[n] = std::max(ttopo, topography[n]);
      open.emplace(n, topography[n]);
    }
  }
}

template <class fT, class Connector_t>
void _PriorityFool(std::vector<fT> &topography, Connector_t *connector,
                   std::vector<size_t> &stack, std::vector<size_t> &Sstack) {
  PQFORPF open;
  std::queue<PQH> pit;
  // std::cout << "DEBUGDEBUG1" << std::endl;
  // connector->update_links_MFD_only(topography);

  // std::cout << "DEBUGDEBUG2" << std::endl;

  // auto topography = format_input(ttopography);
  // uint64_t pitc            = 0;
  // auto     PitTop          = -9999;
  // int      false_pit_cells = 0;

  std::vector<int8_t> closed(connector->nnodes, 0);

  for (int i = 0; i < connector->nnodes; ++i) {
    if (connector->boundaries.can_out(i) == false) {
      // connector->Sreceivers[i] = i+1;
      continue;
    }

    open.emplace(PQH(i, topography[i]));
    closed[i] = 1;
  }

  // std::cout << "DEBUGDEBUG3" << std::endl;
  int stack_incrementor = 0;
  auto neighbours = connector->get_empty_neighbour();
  while (open.size() > 0 || pit.size() > 0) {
    // std::cout << "DEBUGDEBUGstack_incrementor:" << stack_incrementor << "|"
    // << connector->nnodes << std::endl;
    PQH c;
    // if(pit.size()>0 && open.size()>0 && open.top().score ==
    // pit.front().score)
    // {
    // 	c=open.top();
    // 	open.pop();
    // 	PitTop=-9999;
    // }
    // else if(pit.size()>0)
    // {
    // 	c=pit.front();
    // 	pit.pop();
    // 	if(PitTop==-9999)
    // 		PitTop=topography[c.node];
    // // } else {
    c = open.top();
    open.pop();
    // PitTop=topography[c.node];
    stack[stack_incrementor] = c.node;
    ++stack_incrementor;
    // }

    int nn = connector->get_neighbour_idx_links(c.node, neighbours);
    fT ttopo = topography[c.node];

    for (int j = 0; j < nn; ++j) {
      int lix = neighbours[j];
      if (connector->is_link_valid(lix) == false)
        continue;

      int n = connector->get_other_node_from_links(lix, c.node);

      if (closed[n] == 0) {
        // if(closed[n] == 2)
        // 	connector->update_local_link(lix,topography);
        ttopo += connector->randu->get() * 1e-6 + 1e-8;
        // fT ntopo = ;

        if (ttopo > topography[n]) {
          topography[n] = ttopo;
          // connector->update_local_link(lix,topography);
          closed[n] = 2;

        } else {
          closed[n] = 1;
        }

        open.emplace(n, topography[n]);
        connector->update_local_link(lix, topography);
        fT dx = connector->get_dx_from_links_idx(lix);
        fT ts = (topography[n] - topography[c.node]) / dx;
        if (ts >= connector->SS[n]) {
          connector->SS[n] = ts;
          connector->_Sreceivers[n] = c.node;
          connector->Sdistance2receivers[n] = dx;
        }
      }
    }
  }
  connector->compute_SF_donors_from_receivers();
}

// Grayscale morphological reconstruction algorithm
void morphologicalReconstruction(std::vector<double> &data, int width,
                                 int height) {
  // Initialize the marker image with the original data
  std::vector<double> marker(data);

  // Structuring element
  int se_width = 3;
  int se_height = 3;
  std::vector<double> se_data = {1, 1, 1, 1, 1, 1, 1, 1, 1};

  // Iteratively dilate and take intersection
  bool changed = true;
  while (changed) {
    changed = false;

    std::vector<double> dilated(data.size());

    // Perform dilation
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        double max_val = 0;

        // Apply structuring element
        for (int j = -se_height / 2; j <= se_height / 2; j++) {
          for (int i = -se_width / 2; i <= se_width / 2; i++) {
            int nx = x + i;
            int ny = y + j;

            // Skip out-of-bounds pixels
            if (nx < 0 || nx >= width || ny < 0 || ny >= height) {
              continue;
            }

            double val =
                marker[ny * width + nx] +
                se_data[(j + se_height / 2) * se_width + i + se_width / 2];
            max_val = std::max(max_val, val);
          }
        }

        dilated[y * width + x] = max_val;
      }
    }

    // Take intersection with original data
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        double old_val = marker[y * width + x];
        double new_val = std::min(dilated[y * width + x], data[y * width + x]);
        marker[y * width + x] = new_val;

        // Check if any values changed
        if (new_val != old_val) {
          changed = true;
        }
      }
    }
  }

  // Copy the marker image back into the original data
  for (size_t i = 0; i < data.size(); i++) {
    data[i] = marker[i];
  }
}

// DOES NOT WORK!
template <class Connector_t, class graph_t, class int_t, class out_t, class fT>
out_t standalone_priority_flood_opti(int_t &topography, Connector_t &connector,
                                     graph_t &GRAPH) {
  auto ttopo = format_input(topography);
  auto vtopo = DAGGER::to_vec(ttopo);
  PriorityFlood_opti(vtopo, connector, GRAPH.Sstack);
  return DAGGER::format_output<std::vector<fT>, out_t>(vtopo);
}

template <class fT, class Connector_t>
bool PriorityFlood_opti(std::vector<fT> &topography, Connector_t &connector,
                        std::vector<size_t> &Sstack) {
  PQFORPF open;
  std::queue<PQH> pit;

  // auto topography = format_input(ttopography);
  uint64_t pitc = 0;
  auto PitTop = -9999;
  int false_pit_cells = 0;

  std::vector<int8_t> closed(connector.nnodes, false);
  int checker = 0;
  for (int i = 0; i < connector.nnodes; ++i) {
    int node = int(Sstack[i]);
    if (connector.boundaries.no_data(node)) {
      closed[node] = true;
      ++checker;
      continue;
    }

    if (connector.boundaries.can_out(node) == false) {
      closed[node] = closed[connector._Sreceivers[node]];

      if (closed[node])
        ++checker;

      continue;
    }

    // open.emplace(PQH(node,topography[node]));
    closed[node] = true;
    ++checker;
  }
  // std::cout << "CHECKER::" << checker << "|" << connector.nnodes <<
  // std::endl;

  if (checker == connector.nnodes) {
    return false;
  }

  for (int i = 0; i < int(connector.links.size()); ++i) {
    if (connector.is_link_valid(i) == false)
      continue;

    int u, o;
    connector.node_idx_from_link_idx_nocheck(i, u, o);

    if (closed[u] == closed[o])
      continue;

    // open.emplace(PQH(i,topography[i]));

    if (closed[u]) {
      open.emplace(PQH(u, topography[u]));
    } else {
      open.emplace(PQH(o, topography[o]));
    }
  }

  // std::cout << "open size::" << open.size() << std::endl;

  auto neighbours = connector.get_empty_neighbour();
  while (open.size() > 0 || pit.size() > 0) {
    PQH c;
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

    int nn = connector.get_neighbour_idx(c.node, neighbours);

    for (int j = 0; j < nn; ++j) {
      int n = neighbours[j];

      if (connector.is_in_bound(n) == false ||
          connector.boundaries.can_create_link(n) == false)
        continue;

      if (closed[n])
        continue;

      closed[n] = true;

      fT ntopo = topography[c.node] + connector.randu->get() * 1e-6 + 1e-8;

      if (topography[n] == -9999)
        pit.emplace(PQH(n, -9999));

      else if (topography[n] <= topography[c.node]) {

        if (PitTop != -9999 && PitTop < topography[n] &&
            ntopo >= topography[n]) {
          ++false_pit_cells;
        }

        ++pitc;

        topography[n] = ntopo;

        pit.emplace(n, topography[n]);
      } else
        open.emplace(n, topography[n]);
    }
  }

  return true;
}

// Largely adapted from RichDEM
template <class topo_t, class Connector_t>
std::vector<double> PriorityFlood_old(topo_t &ttopography,
                                      Connector_t &connector) {

  PQFORPF open;
  std::queue<PQ_helper<int, double>> pit;

  double PitTop = -9999;

  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator
  std::uniform_real_distribution<> distr(1e-7, 1e-6); // define the range

  auto topography = format_input(ttopography);

  std::vector<int8_t> closed(connector.nnodes, false);

  for (int i = 0; i < connector.nnodes; ++i) {
    if (connector.boundaries.can_out(i)) {
      open.emplace(PQ_helper<int, double>(i, topography[i]));
      closed[i] = true;
    }
  }

  auto neighbours = connector.get_empty_neighbour();
  while (open.size() > 0 || pit.size() > 0) {
    PQ_helper<int, double> c;
    if (pit.size() > 0) {
      c = pit.front();
      pit.pop();
    } else {
      c = open.top();
      open.pop();
    }

    int nn = connector.get_neighbour_idx(c.node, neighbours);

    for (int j = 0; j < nn; ++j) {
      int n = neighbours[j];

      if (connector.boundaries.no_data(n) ||
          connector.boundaries.force_giving(n))
        continue;

      if (closed[n])
        continue;

      closed[n] = true;

      if (topography[n] <= c.score) {
        topography[n] = std::nextafter(c.score, c.score + 1);
        pit.emplace(PQ_helper<int, double>(n, topography[n]));
      } else
        open.emplace(PQ_helper<int, double>(n, topography[n]));
    }
  }
  return to_vec(topography);
}

// Largely adapted from RichDEM
template <class topo_t, class Connector_t, class Graph_t>
std::vector<double> PriorityFlood_broken(topo_t &ttopography,
                                         Connector_t &connector,
                                         Graph_t &graph) {
  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator
  std::uniform_real_distribution<> distr(1e-7, 1e-6); // define the range

  auto tttopography = format_input(ttopography);
  std::vector<double> topography = to_vec(tttopography);

  PQFORPF open;
  std::queue<PQ_helper<int, double>> pit;

  uint64_t processed_cells = 0;
  uint64_t pitc = 0;
  int false_pit_cells = 0;
  double PitTop = -9999;

  std::vector<bool> closed(connector.nnodes, false);

  // RDLOG_PROGRESS<<"Adding cells to the priority queue...";
  for (int i = 0; i < connector.nnodes; ++i) {
    // if(connector.can_flow_even_go_there(i) && i == graph.Sreceivers[i])
    if (connector.boundaries.can_out(i)) {
      closed[i] = true;
      open.emplace(i, topography[i]);
    }
  }

  // RDLOG_PROGRESS<<"Performing Priority-Flood+Epsilon...";
  auto neighbours = connector.get_empty_neighbour();

  // throw std::runtime_error("???A");
  // int ii = 0;
  while (open.size() > 0 || pit.size() > 0) {
    // ++ii;
    // if(ii % 100000 == 0) std::cout << open.size() << "||" << pit.size() <<
    // std::flush << std::endl;

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

    int nn = connector.get_neighbour_idx(c.node, neighbours);
    for (int tnn = 0; tnn < nn; ++nn) {
      int n = neighbours[tnn];
      if (connector.can_flow_even_go_there(n) == false)
        continue;

      if (closed[n])
        continue;

      closed[n] = true;

      if (topography[n] <=
          std::nextafter(c.score, std::numeric_limits<double>::infinity())) {
        if (PitTop != -9999 && PitTop < topography[n] &&
            std::nextafter(c.score, std::numeric_limits<double>::infinity()) >=
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

} // namespace DAGGER

#endif
