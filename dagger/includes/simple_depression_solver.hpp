/*
This header file extends the graph to provide routines to process depressions
using Cordonnier et al, 2019
*/

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef SIMPLE_DEP_SOLVER_HPP
#define SIMPLE_DEP_SOLVER_HPP

// STL imports
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "utils.hpp"
#include "veque.hpp"
// #include "graph.hpp"

// #ifdef BOOST_AVAILABLE
// #include <boost/circular_buffer.hpp>
// #endif

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

namespace DAGGER {
// Useful namespace for my priority queue
using PQ_i_d = std::priority_queue<PQ_helper<int, double>,
                                   std::vector<PQ_helper<int, double>>,
                                   std::greater<PQ_helper<int, double>>>;
using PQH = PQ_helper<int, double>;

// I KEEP THAT HERE IN CASE
// STOP TRYING BORIS THERE IS NO BETTER QUEUE THAT DEQUE FOR MOST CASES
// ONLY CIRCULAR RINGS FITTING THE MEMORY COULD ADD SOME SPEED (i.e < few
// hundreds elements)
template <class T> class vecQ {
public:
  vecQ(){};
  vecQ(int resa) { this->_data.reserve(resa); }
  size_t reader = 0;
  // int reader = 0;
  std::vector<T> _data;
  void emplace(T tt) { this->_data.emplace_back(tt); }
  T front() { return this->_data[this->reader]; }
  T back() { return this->_data.back(); }
  void pop() { ++this->reader; }
  size_t size() { return this->_data.size() - reader; }
  bool empty() { return this->_data.size() == reader; }
};

// Let's try the same with a vector of queues of fixed size

// enum class LDSTATE : std::uint8_t
// {
// 	UNPROC,
// 	DEP,
// 	BORD
// };

template <class in_t, class out_t, class Connector_t>
out_t label_depressions_PQ(Connector_t *connector, in_t &ttopography) {

  // ingesting the topography in an acceptable format
  auto topography = format_input(ttopography);

  // getting the labels
  int ndep = 0;
  DBBTree tree;
  std::vector<int> labels =
      _label_depressions_PQ(connector, topography, ndep, tree);

  // done
  return format_output<decltype(labels), out_t>(labels);
}

template <class in_t, class Connector_t>
std::vector<int> _label_depressions_PQ(Connector_t *connector, in_t &topography,
                                       int &ndep, DBBTree &tree) {
  // init label output
  // std::vector<std::uint8_t> out(connector->nnodes, 0);

  std::vector<int> baslab(connector->nnodes, 0);
  std::vector<std::uint8_t> inQ(connector->nnodes, 0);

  // DBBTree tree;
  // Init with label 0
  tree.add();

  // priority queue
  PQ_i_d sirius;

  // checking each and every nodes
  for (int i = 0; i < connector->nnodes; ++i) {
    // Ignore the ones who escape the model and/or aren't valid
    if (connector->boundaries.no_data(i) || connector->flow_out_model(i))
      continue;

    // then emplace if pit
    if (connector->_Sreceivers[i] == i) {
      sirius.emplace(PQH(i, topography[i]));
      int tlab = tree.add();
      baslab[i] = tlab;
      inQ[i] = 1;
    }
  }

  // Processing the PQ
  auto neighbours = connector->get_empty_neighbour();
  while (sirius.empty() == false) {
    // getting the next node
    auto next = sirius.top();
    sirius.pop();

    int bas = baslab[next.node];

    // First checking if the current depression has been merged in the meantime
    if (tree.top[bas] != bas) {
      // yes, correcting the depression
      bas = tree.top[bas];
      baslab[next.node] = bas;
    }

    if (bas == 0)
      throw std::runtime_error("Should not happen");

    // then checking if it is resolved
    if (tree.label[bas] == 2) {
      // if it is, it means this node is technically not in the depression
      // anymore if(next.node == 172627)
      // {
      // 	bool iswq = topography[next.node] == tree.outlet_Z[bas] ;
      // 	std::cout << "topo:" << topography[next.node] << " || outlet: "
      // <<  tree.outlet_Z[bas]  << " || " << iswq<< std::endl; 	throw
      // std::runtime_error("kjfsdalk;jfsda");
      // }
      baslab[next.node] = 0;
      continue;
    }

    // If I got there, it means the depression is not opened yet

    // getting the neighbours
    int nn = connector->get_neighbour_idx(next.node, neighbours);

    // for each of them
    // bool merge = false;
    std::set<int> tomerge;
    for (int j = 0; j < nn; ++j) {
      int nei = neighbours[j];
      bool lowerz = topography[nei] < topography[next.node];
      int nbas = tree.top[baslab[nei]];

      if (connector->boundaries.can_out(next.node) && tree.label[bas] != 2) {
        // LEAVES MODEL
        tree.label[bas] = 2;
        tree.outlet_node[bas] = next.node;
        tree.outlet_Z[bas] = topography[next.node];
        break;
      } else if ((nbas == 0 && lowerz == false &&
                  inQ[nei] == false)) // || connector->boundaries.can_out(nei))
      {

        // if(inQ[nei] == 1)
        // 	throw std::runtime_error("happens");

        // CASE UNLABELLED AND LOWER ELEVATION
        sirius.emplace(PQH(nei, topography[nei]));
        inQ[nei] = 1;
        baslab[nei] = bas;
      } else if (nbas != bas && nbas != 0) {
        // CAS MERGE OR OUTLET IN DIFFERENT DEP
        if (tree.label[nbas] == 0) {
          tomerge.insert(nbas);
        } else if (tree.label[nbas] == 2) {
          tree.label[bas] = 2;
          tree.outlet_node[bas] = tree.outlet_node[nbas]; // = next.node;
          tree.outlet_Z[bas] = tree.outlet_Z[nbas]; // topography[next.node];
        }
      } else if (lowerz && nbas == 0) {
        // CASE OUTLET
        tree.label[bas] = 2;
        tree.outlet_node[bas] = next.node;
        tree.outlet_Z[bas] = topography[next.node];
      }
    }

    // finally merging if needed
    if (tomerge.size() > 0) {
      tomerge.insert(bas);
      int lbas = tree.merge(tomerge);
      if (tree.label[bas] == 2) {
        tree.label[lbas] = 2;
        tree.outlet_node[lbas] = tree.outlet_node[bas];
        tree.outlet_Z[lbas] = tree.outlet_Z[bas];
      }
    }
  }

  // ndep = 0;

  for (int i = 0; i < connector->nnodes; ++i)
    baslab[i] = tree.top[baslab[i]];

  ndep = tree.nnodes;

  return baslab;
}

template <class in_t, class out_t, class Connector_t>
out_t label_ocean(Connector_t *connector, in_t &ttopography) {

  // ingesting the topography in an acceptable format
  auto topography = format_input(ttopography);

  // getting the labels
  // int ndep = 0;
  // DBBTree tree;
  std::vector<int> labels = _label_ocean(connector, topography);

  // done
  return format_output<decltype(labels), out_t>(labels);
}

template <class in_t, class Connector_t>
std::vector<int> _label_ocean(Connector_t *connector, in_t &topography) {
  // TEST LABELLING HERE
  std::vector<int> ocean(connector->nnodes, -1);

  connector->update_links_MFD_only(topography);
  // #ifdef BOOST_AVAILABLE
  // // std::queue<int, boost::circular_buffer<int> >
  // nexus(boost::circular_buffer<int>(10000)); std::queue<int, std::list<int> >
  // nexus(std::list<int>); #else std::queue<int> nexus; #endif

  // std::queue<int > _, OcQ, PiQ;

  std::stack<int, std::vector<int>> _, OcQ, PiQ;

  // std::queue<int > _, OcQ, PiQ;
  // std::queue<int, boost::circular_buffer<int> >
  // OcQ(boost::circular_buffer<int>(1000));

  PQ_i_d divide_PQ;
  int npits = 0;
  // vecQ<int> OcQ(5000);
  for (int i = 0; i < connector->nnodes; ++i) {
    if (connector->boundaries.can_out(i)) {
      _.emplace(i);
      ocean[i] = 0;
    } else if (connector->flow_out_or_pit(i)) {
      PiQ.emplace(i);
      ocean[i] = i;
      if (i == 150050)
        std::cout << "GASSS" << std::endl;

      ++npits;
    }
  }

  // std::cout << "NPITS::" << npits << std::endl;

  auto neighbours = connector->get_empty_neighbour();
  // int maxsize = 0;

  while (_.empty() == false) {
    OcQ.emplace(_.top());
    _.pop();
    while (OcQ.empty() == false) {
      int next = OcQ.top();
      OcQ.pop();
      int nn = connector->get_donors_idx(next, neighbours);
      for (int j = 0; j < nn; ++j) {
        int oi = neighbours[j];
        if (connector->boundaries.no_data(oi) || ocean[oi] == 0)
          continue;
        ocean[oi] = 0;
        OcQ.emplace(oi);
      }
      // if(OcQ.size() > maxsize) maxsize = OcQ.size();
    }
  }

  while (PiQ.empty() == false) {
    int pit = PiQ.top();
    PiQ.pop();
    OcQ.emplace(pit);
    while (OcQ.empty() == false) {
      int next = OcQ.top();
      OcQ.pop();
      int nn = connector->get_donors_idx(next, neighbours);
      for (int j = 0; j < nn; ++j) {
        int oi = neighbours[j];
        if (connector->boundaries.no_data(oi) || ocean[oi] > 0 ||
            ocean[oi] == -2)
          continue;
        if (ocean[oi] == 0) {
          ocean[oi] = -2;
          divide_PQ.emplace(PQH(oi, topography[oi]));
        } else {
          ocean[oi] = pit;
          OcQ.emplace(oi);
        }

        // if(OcQ.size() > maxsize) maxsize = OcQ.size();
      }
    }
  }

  // std::cout << "divide_PQ size is " << divide_PQ.size() << " " << std::endl;;
  // std::cout << "max size is " << maxsize << " " << std::endl;;

  while (divide_PQ.empty() == false) {
    int tirnext = divide_PQ.top().node;

    divide_PQ.pop();
    int rec = -1;
    while (rec != tirnext) {
      ocean[tirnext] = 0;
      double minZ = topography[tirnext];
      int nn = connector->get_receivers_idx_links(tirnext, neighbours);
      for (int j = 0; j < nn; ++j) {
        int lix = neighbours[j];
        int node = connector->get_to_links(lix);
        if (ocean[node] == 0)
          continue;
        if (topography[node] < minZ) {
          rec = node;
          minZ = topography[node];
          connector->reverse_link(lix);
        }
      }
      // if(tirnext == 150050) std::cout << "GBAA" << std::endl;

      if (rec == -1)
        break;
      tirnext = rec;
    }

    OcQ.emplace(tirnext);
    while (OcQ.empty() == false) {
      int next = OcQ.top();
      OcQ.pop();
      int nn = connector->get_donors_idx(next, neighbours);
      for (int j = 0; j < nn; ++j) {
        int oi = neighbours[j];
        if (connector->boundaries.no_data(oi) || ocean[oi] == 0)
          continue;
        ocean[oi] = 0;
        OcQ.emplace(oi);
      }
      // if(OcQ.size() > maxsize) maxsize = OcQ.size();
    }
  }

  return ocean;
}

template <class Connector_t, class topo_t, class fT>
std::vector<std::int8_t> _dagger_fill(Connector_t *connector,
                                      topo_t &topography,
                                      std::vector<size_t> &stack) {
  // Ocean contains label of nodes connected to the ocean (0), not connected to
  // anything yet (-1), temp, unresolved drainage divide (-2) or connected to an
  // unresolved depression (value = node index of the pit)
  std::vector<std::int8_t> ocean(connector->nnodes, -1);
  std::vector<std::int8_t> inQ(connector->nnodes, false);

  // Precomputing links
  // connector->update_links(topography);
  int stackindex = 0;

  // initialising the different stacks (LIFO data structures)
  // The code juggles between the three to enhance locality and maximise CPU
  // caching Not sure I understand the subbtleties, but it is definitly faster
  // than a single, bigger stack (I benchmarked)
  std::stack<int, std::vector<int>> oc_LIFO, subLIFO, pit_LIFO;

  // Readying a Queue for the filling
  std::queue<int> phil_collins;
  // PQ_i_d phil_collins;

  // Stores the drainage edges
  PQ_i_d divide_PQ;

  std::vector<std::int8_t> trackneighb(connector->nnodes, 0);

  auto neighbours = connector->get_empty_neighbour();
  auto neighbours_links = connector->get_empty_neighbour();
  auto neighbours_Z = connector->template get_empty_neighbour<fT>();

  // Step I) pushing to the Queue the ocean/pit cells and labelling them
  // accordingly

  int npits = 0;
  for (int i = 0; i < connector->nnodes; ++i) {
    int nn = connector->get_neighbour_idx_nodes_and_links(i, neighbours,
                                                          neighbours_links);
    trackneighb[i] += nn;
    // If flow can out the cell -> labels as ocean, emplace in the ocean LIFO
    if (connector->boundaries.can_out(i)) {
      oc_LIFO.emplace(i);
      ocean[i] = 0;
      stack[stackindex] = i;
      ++stackindex;
      for (int j = 0; j < nn; ++j) {
        --trackneighb[i];
        --trackneighb[i];
      }
    }
  }

  // Step II) DFS traversal from ocean nodes to the donors direction to label
  // ocean cells

  while (oc_LIFO.empty() == false) {
    // Getting the first main node from the parent LIFO queue
    subLIFO.emplace(oc_LIFO.top());
    oc_LIFO.pop();

    // Sub-DFS traversal to label ocean cells linked to dat one
    while (subLIFO.empty() == false) {
      int next = subLIFO.top();
      subLIFO.pop();
      int nn = connector->get_neighbour_idx_nodes_links_external_array(
          next, neighbours, neighbours_links, neighbours_Z, topography);
      for (int j = 0; j < nn; ++j) {
        int oi = neighbours[j];
        int lix = neighbours_links[j];
        connector->update_local_link(lix, topography);
        if (ocean[oi] == 0)
          continue;
        // if(ocean[oi] == -1)	connector->update_local_link(lix, topography);
        // basically if the donors is not labelled yet as ocean, I label it and
        // emplace it in the LIFO
        if (topography[oi] > topography[next]) {
          if (ocean[oi] == -1)
            subLIFO.emplace(oi);
          ocean[oi] = 0;
        } else if (ocean[oi] == -1)
          ocean[oi] = -2;
      }
    }
  }
  // Done, all the primary ocean cells are labelled
  for (int i = 0; i < connector->nnodes; ++i) {
    if (ocean[i] == -2) {
      divide_PQ.emplace(i, topography[i]);
      ocean[i] = -1;
      inQ[i] = true;
    }
  }

  // // Now I need to label the cells connected to their respective pits
  // while(pit_LIFO.empty() == false)
  // {
  // 	// getting next pit
  // 	int pit = pit_LIFO.top(); pit_LIFO.pop();

  // 	// reusing the subLIFO to enhance locality
  // 	subLIFO.emplace(pit);
  // 	while(subLIFO.empty() == false)
  // 	{
  // 		int next = subLIFO.top();
  // 		subLIFO.pop();
  // 		int nn = connector->get_donors_idx(next, neighbours);
  // 		for(int j=0; j<nn;++j)
  // 		{
  // 			int oi = neighbours[j];
  // 			// Ignoring nodes already labelled as divide or other
  // pit 			if(connector->boundaries.no_data(oi) ||
  // ocean[oi]
  // >
  // 0
  // || ocean[oi] == -2) continue;
  // 			// If the node is connected to the ocean I push it into
  // the divide PQ and label it as divide 			if(ocean[oi] ==
  // 0)
  // 			{
  // 				ocean[oi] = -2;
  // 				divide_PQ.emplace(PQH(oi,topography[oi]));
  // 			}
  // 			// Otherwise, is connected to dat pit AND EMPLACED TO
  // THE CURRENT SUBlifo whoops cap lock 			else
  // 			{
  // 				ocean[oi] = ocean[pit];
  // 				subLIFO.emplace(oi);
  // 			}
  // 		}
  // 	}
  // }

  // At that stage I have a fully labelled landscapes where each node is either
  // connected to the ocean or to a pit, and the divides are sorted by ascending
  // elevation in the PQ note that the labels are valid in MFD but only
  // represent ONE of the possible state of the nodes: an ocean node can be
  // connected to a pit, it's just that we jsut want to know if it is connected
  // to the ocean and hence can outlet this is also valid for the pit labels: in
  // MFD nodes can be (and will be) connected to multiple pits, but again as we
  // are interested by the FIRST nodes connected to another pit/ocean it's fine
  // return ocean;
  // std::cout << "DEBUGSS4" << ocean[649605] << std::endl;
  // Starting with the lowest divide
  while (divide_PQ.empty() == false) {
    // next node
    int tirnext = divide_PQ.top().node;
    divide_PQ.pop();
    // std::cout << topography[tirnext] << std::endl;

    // if the node is not a divide anymore: skip - it has already been processed
    if (ocean[tirnext] == 0)
      continue;
    phil_collins.emplace(tirnext);
    // phil_collins.emplace(PQH(tirnext, topography[tirnext]));
    while (phil_collins.empty() == false) {
      int next = phil_collins.front();
      phil_collins.pop();
      ocean[next] = 0;

      int nn = connector->get_neighbour_idx_nodes_and_links(next, neighbours,
                                                            neighbours_links);
      double ttopo = topography[next];

      for (int j = 0; j < nn; ++j) {
        int oi = neighbours[j];
        int lix = neighbours_links[j];
        if (ocean[oi] == 0) {
          connector->update_local_link(lix, topography);
          continue;
        }
        if (topography[oi] <= ttopo) {
          ttopo = ttopo + 1e-6 * connector->randu->get() + 1e-8;
          topography[oi] = ttopo;
          phil_collins.emplace(oi);
        } else {
          oc_LIFO.emplace(oi);
        }
        ocean[oi] = 0;
        connector->update_local_link(lix, topography);
      }
    }

    while (oc_LIFO.empty() == false) {
      int tirnext = oc_LIFO.top();
      oc_LIFO.pop();
      subLIFO.emplace(tirnext);
      // Sub-DFS traversal to label ocean cells linked to dat one
      while (subLIFO.empty() == false) {
        int next = subLIFO.top();
        subLIFO.pop();
        int nn = connector->get_neighbour_idx_nodes_and_links(next, neighbours,
                                                              neighbours_links);
        for (int j = 0; j < nn; ++j) {
          int oi = neighbours[j];
          int lix = neighbours_links[j];
          if (ocean[oi] != 0) {
            if (topography[oi] > topography[next]) {
              if (ocean[oi] == -1)
                subLIFO.emplace(oi);
              ocean[oi] = 0;
            } else if (ocean[oi] == -1)
              ocean[oi] = -2;
          }
          connector->update_local_link(lix, topography);
        }
      }
    }

    // Done, all the primary ocean cells are labelled
    for (int i = 0; i < connector->nnodes; ++i) {
      if (ocean[i] == -2 && inQ[i] == false) {
        divide_PQ.emplace(i, topography[i]);
        ocean[i] = -1;
        inQ[i] = true;
      }
    }
  }

  // connector->recompute_SF_donors_from_receivers();
  return ocean;
}

// TEST WITH PQ
// template<class in_t, class Connector_t>
// std::vector<int> _label_ocean(Connector_t* connector, in_t& topography )
// {
// 	// TEST LABELLING HERE
// 	std::vector<int> ocean(connector->nnodes,0);
// 	// #ifdef BOOST_AVAILABLE
// 	// // std::queue<int, boost::circular_buffer<int> >
// nexus(boost::circular_buffer<int>(10000));
// 	// std::queue<int, std::list<int> > nexus(std::list<int>);
// 	// #else
// 	// std::queue<int> nexus;
// 	// #endif

// 	PQ_i_d nexus;
// 	// vecQ<int> nexus(5000);
// 	for(int i=0; i<connector->nnodes; ++i)
// 	{
// 		if(connector->boundaries.can_out(i))
// 		{
// 			nexus.emplace(PQH(i,topography[i]));
// 			ocean[i] = 1;
// 		}
// 	}
// 	auto neighbours = connector->get_empty_neighbour();
// 	// int maxsize = 0;
// 	while(nexus.empty() == false)
// 	{
// 		int next = nexus.top().node;
// 		nexus.pop();
// 		int nn = connector->get_neighbour_idx(next, neighbours);
// 		for(int j=0; j<nn;++j)
// 		{
// 			int oi = neighbours[j];
// 			if(connector->boundaries.no_data(oi) || ocean[oi] == 1
// ||
// topography[oi] < topography[next]) continue; ocean[oi] = 1;
// nexus.emplace(PQH(oi,topography[oi]));
// 		}

// 		// if(nexus.size() > maxsize) maxsize = nexus.size();

// 	}

// 	// std::cout << "max size was " << maxsize << " " << std::endl;;

// 	return ocean;

// }

template <class in_t, class Connector_t>
std::vector<int> _yolo45(Connector_t *connector, in_t &topography) {

  std::vector<int> ocean(connector->nnodes, 0);
  std::queue<int> nexus;
  PQ_i_d inverter;

  // Step 1: label all the out points
  for (int i = 0; i < connector->nnodes; ++i) {
    if (connector->boundaries.can_out(i)) {
      nexus.emplace(i);
      ocean[i] = 1;
    }
  }

  auto neighbours = connector->get_empty_neighbour();
  // int maxsize = 0;
  while (nexus.empty() == false) {
    int next = nexus.front();
    nexus.pop();
    int nn = connector->get_neighbour_idx(next, neighbours);
    for (int j = 0; j < nn; ++j) {
      int oi = neighbours[j];
      if (connector->boundaries.no_data(oi) || ocean[oi] == 1 ||
          topography[oi] < topography[next])
        continue;
      ocean[oi] = 1;
      nexus.emplace(oi);
    }

    // if(nexus.size() > maxsize) maxsize = nexus.size();
  }

  // std::cout << "max size was " << maxsize << " " << std::endl;;

  return ocean;
}

template <class Connector_t, class fT>
bool simple_depression_solver(Connector_t *connector,
                              std::vector<fT> &topography,
                              std::vector<size_t> &Sstack) {
  int ndep = 0;
  DBBTree tree;
  std::vector<int> baslab =
      _label_depressions_PQ(connector, topography, ndep, tree);
  std::vector<std::uint8_t> inQ(connector->nnodes, false);

  std::queue<int> filler;
  auto neighbours = connector->get_empty_neighbour();
  for (int bas = 1; bas < ndep; ++bas) {
    if (tree.top[bas] == bas) {
      int outnode = tree.outlet_node[bas];
      // if(outnode == -1)
      // 	throw std::runtime_error("Bite");

      filler.emplace(outnode);
      inQ[outnode] = true;
      fT topo =
          std::nextafter(tree.outlet_Z[bas], std::numeric_limits<fT>::max());
      topography[outnode] = topo;
      // topo = topo + 1e-7 + connector->randu->get() * 1e-8;
      topo = std::nextafter(topo, std::numeric_limits<fT>::max());
      // topo = std::nextafter(topo, std::numeric_limits<fT>::max()) + 1e-7 +
      // connector->randu->get() * 1e-8;

      while (filler.empty() == false) {
        int next = filler.front();
        filler.pop();
        int nn = connector->get_neighbour_idx(next, neighbours);
        for (int j = 0; j < nn; ++j) {
          int nei = neighbours[j];
          if ((baslab[nei] != bas && topography[nei] != tree.outlet_Z[bas]) ||
              inQ[nei])
            continue;
          inQ[nei] = true;
          topography[nei] = topo;
          // topo = topo + 1e-7 + connector->randu->get() * 1e-8;
          topo = std::nextafter(topo, std::numeric_limits<fT>::max());
          filler.emplace(nei);
        }
      }
    }
  }

  if (ndep > 1)
    return true;
  return false;
}

// std::vector<fT> minscore(connector->nnodes, std::numeric_limits<fT>::max());
// std::vector<int> minnode(connector->nnodes, -1);
// std::vector<std::uint8_t> processed(connector->nnodes, false);
// for(auto node:Sstack)
// {
// 	if(connector->flow_out_model(node)) processed[node] = true;
// 	else processed[node] = processed[connector->_Sreceivers[node]];
// }

// PQ_i_d nexus;

// for(int i=0; i<connector->nlinks(); ++i)
// {
// 	if(connector->is_link_valid(i) == false) continue;
// 	int A,B;connector->node_idx_from_link_idx_nocheck(i,A,B);
// 	if(baslab[A] != baslab[B])
// 	{
// 		fT score = std::min(topography[A],topography[B]);
// 		if(baslab[A] > 0)
// 		{
// 			if(minscore[baslab[A]] > score )
// 			{
// 				minscore[baslab[A]] = score;
// 				minnode[baslab[A]] = A;
// 			}
// 		}
// 		if(baslab[B] > 0)
// 		{
// 			if(minscore[baslab[B]] > score )
// 			{
// 				minscore[baslab[B]] = score;
// 				minnode[baslab[B]] = B;
// 			}
// 		}
// 	}
// }

// for(int i=1;i<ndep; ++i)
// {
// 	if(minnode[i] != -1)
// 	{
// 		nexus.emplace(PQH(minnode[i],minscore[i]));
// 		processed[i] = true;
// 	}
// }

// std::queue<int> filler;
// auto neighbours = connector->get_empty_neighbour();
// auto neighbours2 = connector->get_empty_neighbour();
// while(nexus.empty() == false)
// {
// 	auto next = nexus.top();
// 	nexus.pop();

// 	filler.emplace(next.node);
// 	fT topo = std::nextafter(next.score,std::numeric_limits<fT>::max());
// 	topography[next.node] = topo;
// 	topo = std::nextafter(topo,std::numeric_limits<fT>::max());
// 	int targbas = baslab[next.node];
// 	// baslab[next.node] = 0;

// 	while(filler.empty() == false)
// 	{
// 		int n2fill = filler.front();
// 		filler.pop();
// 		int nn = connector->get_neighbour_idx(n2fill,neighbours);
// 		for(int j=0; j<nn; ++j)
// 		{
// 			int nei1 = neighbours[j];
// 			if(processed[nei1])
// 			{
// 				continue;
// 			}

// 			if((baslab[nei1] == targbas))// || (topography[nei1] <=
// minscore[targbas]))
// 			{

// 				// baslab[nei1] = 0;
// 				topography[nei1] = topo;
// 				topo =
// std::nextafter(topo,std::numeric_limits<fT>::max());
// filler.emplace(nei1); 				processed[nei1] = true;

// 			}
// 			// else if(topography[nei1] <= topo && baslab[nei1] ==
// 0)
// 			// {
// 			// 	bool newLM = true;
// 			// 	int mm = connector->get_neighbour_idx(nei1,
// neighbours2);
// 			// 	for(int k=0; k<mm && newLM; ++k)
// 			// 	{
// 			// 		int nei2 = neighbours2[k];
// 			// 		if(topography[nei2] < topography[nei2])
// 			// 			newLM = false;

// 			// 	}

// 			// 	if(newLM)
// 			// 	{
// 			// 	// baslab[neighbours[j]] = 0;
// 			// 		topography[nei1] = topo;
// 					// topo =
// std::nextafter(topo,std::numeric_limits<fT>::max());
// 			// 		filler.emplace(nei1);
// 			// 		processed[nei1] = true;
// 			// 	}
// 			// }
// 		}

// 	}

// }

// std::vector<std::uint8_t> isdone(ndep, false);

// PQ_i_d spqr;

// auto neighbours = connector->get_empty_neighbour();
// for(int i=0; i<connector->nnodes; ++i)
// {

// 	if(connector->boundaries.no_data(i) || baslab[i] == 0) continue;

// 	int nn = connector->get_neighbour_idx(i,neighbours);
// 	// bool isin = false;
// 	for(int j=0; j< nn ; ++j)
// 	{
// 		int oi = neighbours[j];
// 		if(baslab[oi] != baslab[i])
// 		{
// 			spqr.emplace(PQH(i,topography[i]));
// 			break;
// 		}
// 	}
// }

// std::queue<int> nexus;

// while(spqr.empty() == false)
// {
// 	auto next = spqr.top();
// 	spqr.pop();
// 	int bas = baslab[next.node];
// 	if(isdone[bas]) continue;
// 	isdone[bas] = true;
// 	nexus.emplace(next.node)
// }

// while(nexus.empty() == false)
// {

// 	int nn = connector->get_neighbour_idx(next.node,neighbours);
// 	fT mintopo = topography[next.node];
// 	for(int j=0; j< nn ; ++j)
// 	{
// 		int oi = neighbours[j];
// 		if(baslab[oi] != baslab[next.node] && topography[oi] < mintopo)
// 		{
// 			mintopo = topography[oi];
// 		}
// 		else if (baslab[oi] != 0)
// 			spqr.emplace(PQH(oi,topography[oi]));

// 	}
// 	baslab[next.node] = 0;
// 	topography[next.node] = std::nextafter(mintopo,
// std::numeric_limits<fT>::max());

// }

// return true;

// };
// {

// 	// STEP 1 LABEL ALL BASINS and fetch the one giving out
// 	std::vector<int> basin_key(this->connector->nnodes, -1);
// 	int lab = 0;
// 	for(auto node : Sstack)
// 	{
// 		int rec = connector->_Sreceivers[node];
// 		if(connector->flow_out_model(node))
// 		{
// 			basin_key[node] = 0;
// 		}
// 		else if (node == rec)
// 		{
// 			++lab;
// 			basin_key[node] = lab;
// 		}
// 		else
// 		{
// 			basin_key[node] = basin_key[rec];
// 		}
// 	}

// 	// if no outside basins, returning false (i.e. does not need
// reprocessing) 	if(lab == 0) 		return false;

// 	int n_basins = lab+1;
// 	std::vector<int> connode(n_basins, -1);
// 	std::vector<std::uint8_t> isopened(n_basins, 0);
// 	isopened[0] = 1;
// 	std::vector<fT> connode_z(n_basins, std::numeric_limits<fT>::max());

// 	while( lab > 0 )
// 	{
// 		for(int i = 0; i < int(connector->links.size()) ++i)
// 		{
// 			if(connector->is_link_valid(i) == false)
// 				continue;
// 			int a,b; node_idx_from_link_idx_nocheck(i,a,b);
// 			int ba = basin_key[a], bb = basin_key[b];
// 			if(isopened[ba] != isopened[bb])
// 			{
// 				int tclosed = (isopened[ba] == 1) ? bb : ba;
// 				fT tmitopo =
// std::min(topography[a],topography[b]);
// if(connode_z[tclosed] > tmitopo)
// 				{
// 					connode//...
// 				}
// 			}
// 		}
// 	}

// }

} // namespace DAGGER

#endif
