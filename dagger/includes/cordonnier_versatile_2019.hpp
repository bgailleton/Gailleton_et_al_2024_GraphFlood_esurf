/*
This header file extends the graph to provide routines to process depressions
using Cordonnier et al, 2019
*/

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef cordonnier_versatile_2019_HPP
#define cordonnier_versatile_2019_HPP

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
// #include "graph.hpp"

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

namespace DAGGER {

template <class float_t> class basinLink {
public:
  int b1, b2, from, to;
  float_t score;
  basinLink() { ; }
  basinLink(int b1, int b2, int from, int to, float_t score) {
    this->b1 = b1;
    this->b2 = b2;
    this->from = from;
    this->to = to;
    this->score = score;
  }
};

// Custom operator sorting the nodes by scores
template <class float_t>
inline bool operator>(const basinLink<float_t> &lhs,
                      const basinLink<float_t> &rhs) {
  return lhs.score > rhs.score;
}
template <class float_t>
inline bool operator<(const basinLink<float_t> &lhs,
                      const basinLink<float_t> &rhs) {
  return lhs.score < rhs.score;
}

template <class float_t> class SparseStorer {
public:
  int nel = 0;
  int nlinks = 0;
  bool optimize_memory = false;
  std::vector<bool> isin;
  std::unordered_map<std::string, basinLink<float_t>> blinks;
  SparseStorer() { ; }
  SparseStorer(int nel) {
    this->nel = nel;
    if (this->optimize_memory == false)
      this->isin = std::vector<bool>(this->nel * this->nel, false);
  }

  void insert(std::string idx, basinLink<float_t> &basl) {
    bool tisin = false;
    if (this->optimize_memory) {
      tisin = (this->blinks.count(idx) == 0) ? false : true;
    }

    if (tisin) {
      // basinLink<float_t>& cubasl = this->blinks[idx];
      if (basl.score < this->blinks[idx].score) {
        // std::cout << "happens::" << basl.score << " vs " <<
        // this->blinks[idx].score << std::endl;
        this->blinks[idx] = basl;
      }
    } else {
      ++this->nlinks;
      this->blinks[idx] = basl;
    }
  }
};

template <class n_t, class dist_t, class Connector_t, class topo_t, class LM_t>
class UnionFind {
public:
  UnionFind(int size, LM_t &cod) {
    this->cod = &cod;
    this->_parent = std::vector<int>(size);
    this->_open = std::vector<bool>(size);
    for (int i = 0; i < size; i++) {
      this->_parent[i] = i;
      this->_open[i] = this->cod->is_open_basin[i];
    }
    this->_rank = std::vector<int>(size, 0);
  };

  void Union(int &x, int &y) {
    int xroot = this->Find(x);
    int yroot = this->Find(y);

    if (xroot != yroot) {
      if (this->_rank[xroot] < this->_rank[yroot])
        this->_parent[xroot] = yroot;
      else {
        this->_parent[yroot] = xroot;
        if (this->_rank[xroot] == this->_rank[yroot])
          this->_rank[xroot]++;
      }

      if (this->_open[xroot] || this->_open[yroot]) {
        this->_open[xroot] = true;
        this->_open[yroot] = true;
      }
    }
  }

  int Find(int &x) {
    int xp = x, xc;
    while (true) {
      xc = xp;
      xp = this->_parent[xc];
      if (xp == xc)
        break;
    }
    this->_parent[x] = xc;
    return xc;
  }

  std::vector<int> _parent;
  std::vector<int> _rank;
  std::vector<bool> _open;
  LM_t *cod;
};

// Only for pairs of std::hash-able types for simplicity.
// You can of course template this struct to allow other hash functions
struct pair_hash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &p) const {
    auto h1 = std::hash<T1>{}(p.first);
    auto h2 = std::hash<T2>{}(p.second);

    // Mainly for demonstration purposes, i.e. works but is overly simple
    // In the real world, use sth. like boost.hash_combine
    return h1 ^ h2;
  }
};

template <class float_t> class LMRerouter {

public:
  int nbas;

  // node size
  std::vector<int> basins;

  // nbasins size
  std::vector<bool> is_open_basin;
  std::vector<int> receivers;
  std::vector<int> pitnode;
  std::vector<std::pair<int, int>> receivers_node;
  std::vector<int> stack;
  std::vector<std::vector<int>> donors;

  float_t minimum_slope = 1e-4;
  float_t slope_randomness = 1e-6;

  bool opti_sparse_border = false;

  //
  // std::unordered_map<int , float_t> edges;
  // std::unordered_map<int , std::pair<int,int> > edges_nodes;

  SparseStorer<float_t> stostor;

  LMRerouter() { ; };

  template <class topo_t, class Connector_t>
  bool run(DEPRES method, topo_t &topography, Connector_t *connector,
           std::vector<size_t> &Sstack) {
    // std::cout << "DEBUGLM_II::1" <<std::endl;
    // tracking the number of basins
    this->nbas = -1;
    // tracking to which basin each node belongs to
    this->basins = std::vector<int>(connector->nnodes, 0);
    // number of internal basins to reroute
    int nbas2solve = 0;

    // First comuting the basin array
    for (int i = 0; i < connector->nnodes; ++i) {
      // getting the current nodes
      size_t node = Sstack[i];
      if (connector->boundaries.no_data(node))
        continue;

      // If it is its own receiver, then
      if (connector->_Sreceivers[node] == int(node)) {
        // incrementing the basin label
        ++this->nbas;
        // is it a base level
        if (connector->boundaries.can_out(node)) {
          // then it is an open basin
          this->is_open_basin.emplace_back(true);
          // saving its pit node
          this->pitnode.emplace_back(node);
        } else {
          // otherwise it is a basin to solve
          ++nbas2solve;
          // not open
          this->is_open_basin.emplace_back(false);
          // saving its pit node
          this->pitnode.emplace_back(node);
        }
      }

      // labelling the node
      this->basins[node] = this->nbas;
    }

    // need a last increment as first label is 0
    ++this->nbas;

    // std::cout << "DEBUGLM_II::nbas2solve" << nbas2solve <<std::endl;

    // Relabelling 0 all the open basins to gain time
    if (this->opti_sparse_border) {
      for (int i = 0; i < connector->nnodes; ++i) {
        if (this->is_open_basin[this->basins[i]])
          this->basins[i] = 0;
      }
    }

    // std::cout << "DEBUGLM_II::3" <<std::endl;
    // if there is literally no basins to solve, then I am done
    if (nbas2solve == 0) {
      // std::cout << "DEBUG::NODEP" << std::endl;
      return false;
    }

    // tracking the number of connector->linknodes between a basin to another
    // int nconnector->linknodes = 0;

    // std::cout << "DEBUGLM_II::4" <<std::endl;
    // this->stostor =
    // SparseStorer<float_t>(std::round(std::pow(this->nbas,2)/2));
    this->stostor = SparseStorer<float_t>();
    this->stostor.optimize_memory = true;

    // going through each and every link
    for (int i = 0; i < int(connector->links.size()); ++i) {
      // if the link is not valid, I continue
      if (connector->is_link_valid(i) == false) {
        continue;
      }

      int o1, o2;
      connector->node_idx_from_link_idx_nocheck(i, o1, o2);

      if (!connector->boundaries.can_receive(o1) ||
          !connector->boundaries.can_receive(o2))
        continue;

      // translating the nodes to basin IDs
      int bj = this->basins[o1];
      int bk = this->basins[o2];

      // if in same basin or both open -> I skip
      if (bj == bk || (this->is_open_basin[bj] && this->is_open_basin[bk]))
        continue;

      // The score is the minimum elevation of the pass or the elevation of the
      // outlet if the flow can out a place
      float_t score;
      if (connector->boundaries.can_out(o1)) {
        score = topography[o1];
      } else if (connector->boundaries.can_out(o2)) {
        score = topography[o2];
      } else {
        score = std::min(topography[o1], topography[o2]);
      }

      // is bj < bk (the std::pair storing the pass always starts from the lowes
      // to the highest by convention to keep the std::pair map keys unique)
      bool bjmin = bj < bk;
      if (bjmin == false) {
        std::swap(bj, bk);
        std::swap(o1, o2);
      }

      // if (bj<bk) pair is {bj,bk} else {bk,bj} (I love ternary operators)
      // std::pair<int,int> tp = {(bjmin)?bj:bk, (bjmin)?bk:bj};
      std::string tp = this->uniqueBasid(bj, bk);

      basinLink<float_t> tbasl(bj, bk, o1, o2, score);
      this->stostor.insert(tp, tbasl);
    }
    // std::cout << "I had " << this->stostor.nlinks << std::endl;
    // Done with the link construction

    // Gathering all the links in a vector
    std::vector<basinLink<float_t>> thesebasinlinks;
    thesebasinlinks.reserve(this->stostor.nlinks);
    for (auto it : this->stostor.blinks) {
      thesebasinlinks.emplace_back(it.second);
    }

    // And sorting it
    std::sort(thesebasinlinks.begin(), thesebasinlinks.end());

    // This will track which links are active or not
    std::vector<bool> isactive(thesebasinlinks.size(), false);

    // std::cout << "DEBUGLM_II::6" <<std::endl;

    // This part is applying the kruskal algorithm (I think)
    UnionFind<int, float_t, Connector_t, topo_t, LMRerouter> uf(this->nbas,
                                                                (*this));

    // trackng the receiver of all basins
    this->receivers = std::vector<int>(this->nbas);
    // and the subsequent node pair
    this->receivers_node = std::vector<std::pair<int, int>>(this->nbas);
    // init recs to themselves (base level)
    for (int i = 0; i < this->nbas; ++i)
      this->receivers[i] = i;

    // ok going through all the links from lowest pass to the highest
    for (size_t i = 0; i < thesebasinlinks.size(); ++i) {
      // getting next link
      auto &next = thesebasinlinks[i];

      // getting basin IDs
      int b1 = next.b1;
      int b2 = next.b2;

      // getting basin IDs unionised ☭
      int fb1 = uf.Find(b1);
      int fb2 = uf.Find(b2);

      // If they are united, I skip (they already merged)
      if (fb1 != fb2) {

        // if both are open, I skip
        if (uf._open[fb1] && uf._open[fb2])
          continue;

        // Unification of both basin
        uf.Union(b1, b2);
        // this link is active then
        isactive[i] = true;

        // // If basin one is open
        // if(this->is_open_basin[b1])
        // {
        // 	std::cout << "gulg::1" << this->is_open_basin[b2]	<<
        // std::endl;
        // 	// rec of b2 is b1
        // 	this->receivers[b2] = b1;
        // 	// connecting node are node b2 to node b1
        // 	this->receivers_node[b2] =
        // std::pair<int,int>{this->edges_nodes[next.node].second
        // ,this->edges_nodes[next.node].first};
        // 	// b2 is now open
        // 	this->is_open_basin[b2] = true;
        // }
        // else if(this->is_open_basin[b2])
        // {
        // 	std::cout << "gulg::2" << this->is_open_basin[b1] << std::endl;
        // 	this->receivers[b1] = b2;
        // 	this->receivers_node[b1] =
        // std::pair<int,int>{this->edges_nodes[next.node].first
        // ,this->edges_nodes[next.node].second}; 	this->is_open_basin[b1]
        // = true;
        // }
      }
    }

    // std::cout << "DEBUGLM_II::7" <<std::endl;

    int p_nlinkignored = -1;
    while (true) {
      bool alltrue = true;
      int nlinkignored = 0;
      for (size_t i = 0; i < thesebasinlinks.size(); ++i) {
        if (isactive[i] == false) {
          nlinkignored++;
          continue;
        }

        int b1 = thesebasinlinks[i].b1, b2 = thesebasinlinks[i].b2;
        if (this->is_open_basin[b1] && this->is_open_basin[b2]) {
          nlinkignored++;
          continue;
        }

        auto &next = thesebasinlinks[i];

        // std::cout << "bulf";
        // if(this->basins[next.to] != b2)
        // 	std::cout << "IJKOPFDKLFDSL:DSFL:DFSL:JKDFDFLHJKDFLK" <<
        // std::endl;

        // if(this->basins[next.from] != b1)
        // 	std::cout << "IJKOPFDKLFDSL:DSFL:DFSL:JKDFDFLHJKDFLK2" <<
        // std::endl;

        if (this->is_open_basin[b1]) {
          // std::cout << "pluf" << std::endl;
          this->receivers[b2] = b1;
          this->receivers_node[b2] = std::pair<int, int>{next.to, next.from};
          this->is_open_basin[b2] = true;
        } else if (this->is_open_basin[b2]) {
          // std::cout << "pluf" << std::endl;
          this->receivers[b1] = b2;
          this->receivers_node[b1] = std::pair<int, int>{next.from, next.to};
          this->is_open_basin[b1] = true;
        } else
          alltrue = false;
      }
      // std::cout << nlinkignored << "/" << thesebasinlinks.size() <<
      // std::endl;

      if (alltrue || p_nlinkignored == nlinkignored)
        break;
      p_nlinkignored = nlinkignored;
    }

    // std::cout << "DEBUGLM_II::8" <<std::endl;

    this->donors =
        std::vector<std::vector<int>>(this->nbas, std::vector<int>());
    for (int i = 0; i < this->nbas; ++i) {
      if (this->receivers[i] != i) {
        this->donors[this->receivers[i]].emplace_back(i);
      }
    }

    // std::cout << "DEBUGLM_II::9" <<std::endl;
    this->compute_TO_SF_stack_version();

    // std::cout << "DEBUGLM_II::10::" << this->stack.size() <<std::endl;

    if (method == DEPRES::cordonnier_carve ||
        true) // artificially enabled (see cordonnier fill bellow, temporarily
              // broken)
    {
      for (int i = this->nbas - 1; i >= 0; --i) {
        int bas = this->stack[i];
        // std::cout << bas << "/" << this->nbas << std::endl;;
        // HOTFIX TO CHECK!!!!!: || this->is_open_basin[bas] == false, needed
        // when p_nlinkignored is triggered
        if (connector->boundaries.can_out(
                this->pitnode[bas])) // || this->is_open_basin[bas] == false)
          continue;

        // std::cout << "A" << std::endl;
        int from = this->receivers_node[bas].first;
        int to = this->receivers_node[bas].second;
        // std::cout << "B" << std::endl;
        // std::cout << connector->_Sreceivers[this->pitnode[bas]] << "|";

        int A = from;
        int B = connector->_Sreceivers[A];
        int C = B;
        // std::cout << "C" << std::endl;

        while (A != this->pitnode[bas]) {
          // std::cout << B << std::endl;
          C = connector->_Sreceivers[B];
          connector->_Sreceivers[B] = A;
          // connector->Sdistance2receivers[B] =
          // connector->Sdistance2receivers[A];
          connector->Sdistance2receivers[B] = connector->dx;

          A = B;
          B = C;
        }
        // std::cout << "D" << std::endl;

        connector->_Sreceivers[from] = to;
        connector->Sdistance2receivers[from] = connector->dx;

        // std::cout << connector->_Sreceivers[this->pitnode[bas]] << std::endl;
      }

      for (int i = 0; i < connector->nnodes; ++i) {
        if (i == connector->_Sreceivers[i]) {
          if (connector->boundaries.force_giving(i))
            throw std::runtime_error("HAPPENS::force_giving_is_outlet");
        }
      }

    } else if (method == DEPRES::cordonnier_fill &&
               false) // artificially disconnected, see bellow
    {
      // std::cout << "NEED FIXING" << std::endl;

      // THIS IS THE "normal" cordonnier_fill but is broken at the moment.
      // Switching back to a fake carve-fill
      std::vector<std::uint8_t> inQ(connector->nnodes, 0);
      auto neighbours = connector->get_empty_neighbour();
      for (int i = 0; i < this->nbas; ++i)
      // for(int i =	this->nbas-1; i>=0; --i)
      {
        int bas = this->stack[i];

        // HOTFIX TO CHECK!!!!!: || this->is_open_basin[bas] == false, needed
        // when p_nlinkignored is triggered
        if (connector->boundaries.can_out(this->pitnode[bas]) ||
            this->is_open_basin[bas] == false)
          continue;

        int from = this->receivers_node[bas].first;
        int to = this->receivers_node[bas].second;

        // float_t cost = topography[to];
        connector->_Sreceivers[from] = to;
        connector->Sdistance2receivers[from] = connector->dx;
        topography[from] = std::max(topography[to], topography[from]) +
                           connector->randu->get() * 1e-6 + 1e-8;
        ;

        std::queue<int> Q;
        Q.emplace(from);

        while (Q.empty() == false) {

          int n = Q.front();
          Q.pop();

          inQ[n] = 1;

          int nn = connector->get_neighbour_idx(n, neighbours);

          for (int j = 0; j < nn; ++j) {
            int tn = neighbours[j];
            int tbas = this->basins[tn];

            if (tbas != bas) {
              continue;
            }

            if (inQ[tn] == 1) {
              continue;
            }

            if (topography[tn] > topography[n]) {
              continue;
            }

            if (connector->boundaries.forcing_io(tn))
              continue;

            // if( tn == 0 )
            // 	std::cout << "before::" << topography[tn];
            // topography[tn] = topography[n] + connector->randu->get() * 1e-6 +
            // 1e-8; if( tn == 0 ) 	std::cout << "after::" << topography[tn]
            // << std::endl;;
            connector->_Sreceivers[tn] = n;
            connector->Sdistance2receivers[tn] =
                connector->Sdistance2receivers[n];
            Q.emplace(tn);
            inQ[tn] = 1;
          }
        }
      }

      // TO KEEP IN CASE
      // std::vector<bool> isinQ(connector->nnodes,false);
      // std::vector<bool> isfilled(connector->nnodes,false);
      // std::vector<bool> basinDone(this->nbas,false);
      // std::vector<int> basfam(this->nbas,-1);
      // for(int i = 0; i< this->nbas; ++i)
      // {
      // 	if(connector->boundaries.can_out(this->pitnode[i]))
      // 		basinDone[i] = true;

      // 	int node = this->stack[i];
      // 	if(this->receivers[node] != node)
      // 		basfam[node] = basfam[this->receivers[node]];
      // 	else
      // 		basfam[node] = node;
      // }

      // auto neighbours = connector->get_empty_neighbour();
      // for(int i = 0; i < this->nbas; ++i)
      // {
      // 	int bas = this->stack[i];
      // 	if(connector->boundaries.can_out(this->pitnode[bas]))
      // 		continue;

      // 	int from = this->receivers_node[bas].first;
      // 	int to = this->receivers_node[bas].second;
      // 	float_t zref = std::max(topography[from], topography[to]);
      // 	connector->_Sreceivers[from] = to;
      // 	connector->Sdistance2receivers[from] = connector->dx;
      // 	isinQ[from] = true;
      // 	std::queue<int> Q;Q.emplace(from);
      // 	while(Q.empty() == false)
      // 	{
      // 		int next = Q.front();Q.pop();
      // 		isfilled[next] = true;
      // 		int nn = connector->get_neighbour_idx(next, neighbours);

      // 		float_t lowest_z =
      // std::max(topography[connector->_Sreceivers[next]],topography[next]);
      // int nznodeext = connector->_Sreceivers[next];

      // 		for(int tnn = 0 ; tnn<nn; ++tnn )
      // 		{

      // 			int n = neighbours[tnn];

      // 			if (n<0 || n >= connector->nnodes)
      // 				std::cout <<
      // "JDFKHJLDSFHJKDSFHJKLDFSHJKLDFS" << std::endl;

      // 			int basn = this->basins[n];

      // 			if(isfilled[n] || basinDone[basn] ||
      // basfam[basn]
      // != basfam[bas])
      // 			{
      // 				if(lowest_z > topography[n])
      // 				{
      // 					lowest_z = topography[n];
      // 					nznodeext = n;
      // 				}
      // 			}

      // 			if(isinQ[n])
      // 				continue;
      // 			// if(basfam[basn] != basfam[bas])
      // 			if(basn != bas)
      // 				continue;
      // 			if(basinDone[basn])
      // 				continue;
      // 			if(topography[n] <= zref)
      // 			{
      // 				connector->_Sreceivers[n] = next;
      // 				isinQ[n] = true;
      // 				Q.emplace(n);
      // 			}
      // 		}

      // 		topography[next] = std::max(lowest_z + minimum_slope +
      // connector->randu->get() * slope_randomness, topography[next]);
      // zref = std::max(topography[next],zref);
      // connector->_Sreceivers[next] = nznodeext;
      // connector->Sdistance2receivers[next] = connector->dx;
      // 	}

      // 	basinDone[bas] = true;
      // }

    }

    else if (method == DEPRES::cordonnier_simple) {
      // std::cout << "KJSDFJHKDAFAKJHL" << std::endl;
      for (int i = 0; i < this->nbas; ++i) {
        int bas = this->stack[i];
        if (connector->boundaries.can_out(this->pitnode[bas]))
          continue;

        int to = this->receivers_node[bas].second;
        // std::cout << "rerouting " << this->pitnode[bas] << " to " << to <<
        // std::endl;
        connector->_Sreceivers[this->pitnode[bas]] = to;
        connector->Sdistance2receivers[this->pitnode[bas]] = connector->dx;
      }
    }

    // std::cout << "DEBUGLM_II::11" <<std::endl;

    return true;
  }

  void compute_TO_SF_stack_version() {
    // Initialising the stack
    this->stack.clear();
    // reserving the amount of stuff
    this->stack.reserve(this->nbas);

    // The stack container helper
    std::stack<int, std::vector<int>> stackhelper;
    // std::vector<bool> isdone(this->nbas,false);
    // going through all the nodes
    for (int i = 0; i < this->nbas; ++i) {
      // if they are base level I include them in the stack
      if (this->receivers[i] == i) {
        stackhelper.emplace(i);
      }

      // While I still have stuff in the stack helper
      while (stackhelper.empty() == false) {
        // I get the next node and pop it from the stack helper
        int nextnode = stackhelper.top();
        stackhelper.pop();
        // std::cout << stackhelper.size() << "->" << nextnode << std::endl;

        // // I emplace it in the stack
        // if(isdone[nextnode] == true)
        //	throw std::runtime_error("node-duplicate");

        // isdone[nextnode] = true;
        this->stack.emplace_back(nextnode);

        // as well as all its donors which will be processed next
        for (size_t j = 0; j < this->donors[nextnode].size(); ++j) {
          stackhelper.emplace(this->donors[nextnode][j]);
        }
      }
    }

    // if(this->nbas != this->stack.size())
    // 	throw std::runtime_error("stacksize issue in LMRerouter::" +
    // std::to_string(this->stack.size()) + " vs " +
    // std::to_string(this->stack.size()));

    // for(auto v:this->stack)
    // 	std::cout << v << "|";
    // std::cout << std::endl;
  }

  std::string uniqueBasid(int b1, int b2) {
    return std::to_string(b1) + "_" + std::to_string(b2);
  }
  // void basid2bas(int basid, int& b1, int& b2)
  // {
  // 	b1 = basid % this->nbas;
  // 	b2 = (int)std::floor(basid/this->nbas);
  // 	if(b2<b1)
  // 		std::swap(b2,b1);
  // }
};

// // SOMETHING TO TEST HERE
// template<class Graph_t, class Connector_t, class float_t>
// class dagger_LM
// {
// public:

// 	Graph_t* graph;
// 	Connector_t* con;
// 	std::vector<size_t> Sstack;
// 	std::vector<int> baslab;
// 	std::vector<int> neighbours;
// 	std::vector<bool> is_Sdonor;
// 	std::vector<float_t>* topography;
// 	std::vector<bool> basout;
// 	std::stack<size_t, std::vector<size_t> > stackhelper;
// 	int idx_stack_checker = 0;

// 	dagger_LM(){;}
// 	dagger_LM(Graph_t& tgraph, Connector_t& tcon)
// 	{
// 		this->graph = &tgraph;
// 		this->con = &tcon;
// 		this->Sstack.reserve(this->graph->nnodes);
// 		this->baslab = std::vector<int>(this->graph->nnodes, -1);
// 		this->neighbours = this->con->get_empty_neighbour();
// 		this->is_Sdonor =
// std::vector<bool>(this->neighbours.size(),false);
// 	}

// 	void run()
// 	{
// 		// first selecting the base nodes
// 		std::priority_queue< PQ_helper<std::pair<int,int>,float_t>,
// std::vector<PQ_helper<std::pair<int,int>,float_t> >,
// std::greater<PQ_helper<std::pair<int,int>,float_t> > > maPQ;

// 		std::vector<int> next_nodes;
// 		next_nodes.reserve(this->con->nx * 4);
// 		for(int i=0; i<this->graph->nnodes; ++i) if(i ==
// this->graph->connector->_Sreceivers[i]) next_nodes.emplace_back(i);
// bool keepon = true; 		do
// 		{
// 			this->toposort_section(next_nodes);
// 			while(this->idx_stack_checker <
// int(this->Sstack.size()))
// 			{

// 				++this->idx_stack_checker;
// 			}

// 		}while(keepon);

// 	}

// 	void check_pass_at_node(int i)
// 	{
// 		int nn = this->con->get_neighbour_idx(i,this->neighbours);
// 		for(int j=0; j< nn; ++j)
// 		{
// 			if()
// 		}
// 	}

// 	void toposort_section(std::vector<int>& startnodes)
// 	{

// 		// going through all the starting nodes
// 		for(auto i:startnodes)
// 		{

// 			stackhelper.emplace(i);

// 			// While I still have stuff in the stack helper
// 			while(stackhelper.empty() == false)
// 			{
// 				// I get the next node and pop it from the stack
// helper 				int nextnode =
// stackhelper.top();stackhelper.pop();
// this->Sstack.emplace_back(nextnode);
// this->basout[nextnode] = true;

// 				int nn = this->get_Sdonors_manual(nextnode);

// 				// as well as all its donors which will be
// processed next 				for(int j = 0; j < nn; ++j)
// stackhelper.emplace(this->neighbours[j]);

// 			}

// 		}
// 	}

// 	int get_Sdonors_manual(int i)
// 	{
// 		// this->rinit_is_Sdonor();
// 		int nn = this->con->get_neighbour_idx(i,this->neighbours);
// 		int jn = 0;
// 		for(int j=0; j< nn; ++j)
// 		{
// 			if(this->graph->connector->_Sreceivers[this->neighbours[j]]
// == i)
// 			{
// 				// this->is_Sdonor[j] = true;
// 				this->neighbours[jn] = this->neighbours[j];
// 				++jn;
// 			}
// 		}
// 		return jn;
// 	}

// 	void rinit_is_Sdonor()
// 	{
// 		for(size_t j =0; j < this->is_Sdonor.size(); ++j)
// 			this->is_Sdonor[j] = false;
// 	}

// };

// bool dagger_rerouter_v1(std::vector<float_t>& ttopo, Graph_t& graph,
// Connector_t& connector)
// {

// 	// Dynamic Stack
// 	std::vector<size_t> Sstack;Sstack.reserve(graph.nnodes);
// 	// Basin Labels
// 	std::vector<int> baslab;baslab.reserve(std::ceil(graph.nnodes/10));
// 	// is a basin connected to the outer world?
// 	std::vector<bool> basout;basout.reserve(std::ceil(graph.nnodes/10));
// 	// indices
// 	int idinsert = 0;

// 	// first, let's compute a partial toposort

// 	return true;
// }

// end of DAGGER namesapce
}; // namespace DAGGER

#endif
