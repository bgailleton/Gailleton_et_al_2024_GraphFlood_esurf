//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef graph_HPP
#define graph_HPP

/*
This file contains the graph class.
The graph manages anything linked to the DAG and is the main interface to the
code. B.G. 2022
*/

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
#include <iomanip>

// local includes
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
#include "depression_hierarchy.hpp"
// -> The connector classes
#include "D8connector.hpp"
#include "priority_flood.hpp"

#include "wrap_helper.hpp"

namespace DAGGER {

template<class fT,
				 class Connector_t,
				 class dummy_t =
					 int> // the class type dummy_t is there to bypass pybind11 issue
								// of not being able to bind the same object twice
class graph
{

	// Everything goes public, more straighforward
public:
	// Number of nodes in the graph
	int nnodes;

	Connector_t* connector;

	//======================================================================
	//======================================================================
	//======================================================================
	//======================================================================
	// MOVED TO THE CONNECTOR OBJECT
	// refactored 02/2023
	// Keeping that in for a few commits in case
	//=====================================================================
	// // uint8_t vector for each link indicating its directionality:
	// // - 0 is inverse (linknode 2 --> linknode 1)
	// // - 1 is normal (linknode 1 --> linknode 2)
	// // - 3 is invalid (link index may exist, but is either inexsitant (index is
	// conserved for speed reason when the link index is calculated from node
	// index) ) or temporarily invalid due to dynamic boundary conditions
	// // The index itself is calculated from a connector
	// std::vector<std::uint8_t> links;

	// // Single graph receivers
	// // -> Sreceivers: steepest recervers (nnodes size),
	// // -> number of donors (nnodes size),
	// // -> Steepest donors (nnodes * n_neighbours size)
	// // --> Sdonors of node i are located from index i*n_neighbours to index
	// i*n_neighbours + nSdonors[i] not included std::vector<int>
	// Sreceivers,nSdonors,Sdonors;

	// // Single graph distance to receivers
	// std::vector<fT> Sdistance2receivers;

	// // Steepest slope
	// std::vector<fT> SS;

	// // integer vector of 2*links size with the node indices of each link
	// // for example, the nodes of link #42 would be indices 84 and 85
	// std::vector<int> linknodes;
	//======================================================================
	//======================================================================
	//======================================================================
	//======================================================================

	// Topological order and Single graph topological order
	// While the MD stack can be used for single flow topology,
	// the SD graph sensu Braun and Willett 2013 is built in a very comprehensive
	// way making operations such as watershed labelling or connected component
	// gathering particularly efficient
	std::vector<size_t> stack, Sstack;

	bool debug_mask = false;
	std::vector<std::uint8_t> _debug_mask;
	std::vector<std::uint8_t> get_debug_mask() { return this->_debug_mask; }

	std::vector<int> _debug_int;
	std::vector<int> get_debug_int() { return this->_debug_int; }

	// What depression resolver to use when computing the graph
	DEPRES depression_resolver = DEPRES::cordonnier_carve;

	// hte minimum slope to impose on the fake topography
	fT minimum_slope = 1e-4;
	fT slope_randomness = 1e-6;

	std::shared_ptr<easyRand> randu;

	// Solving large number of local minima with cordonnier can be expensive and
	// this optimises But the Algorithm may become unstable - especially if
	// boundary conditions are weird
	bool opti_sparse_border = false;

	int n_pits = 0;
	int get_n_pits() { return this->n_pits; }

	fT _rdid = 0.;
	void gen_rdid() { this->_rdid = this->randu->get(); };
	fT get_rdid() { return this->_rdid; };

	// Other options

	// # if opt_stst_rerouting is true, only the multiple flow receivers are fully
	// recalculated after preprocessing the dem # The steepest descent are
	// approximated by the rerouting and may not represent the steepest slope in
	// areas where the topography has been modified # Defaulted to true because
	// the differences are minor and not necessarily wanted
	bool opt_stst_rerouting = true;
	void set_opt_stst_rerouting(bool nval) { this->opt_stst_rerouting = nval; }

	// default empty constructor
	graph(){};

	// Classic constructor, simply giving the number of nodes and the number of
	// neighbours per node
	graph(Connector_t& con) { this->init(con); }
	void init(Connector_t& con)
	{
		this->connector = &con;
		this->nnodes = this->connector->nnodes;
		this->randu = con.randu;
		this->init_graph();
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
	Graph functions	- everything about computing, updating and managing the
graph, links, nodes and local minima The core of the code in other words.
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	// Initialisation of the graph structure
	// Uses the nnodes, n_neighbours and connector to alloate the single vectors
	// and the irec/linknodes attribute template<class Connector_t>
	void init_graph()
	{
		// Allocate vectors to node size
		this->_allocate_vectors();
	}

	// used to reinitialise the the vectors
	// template<class Connector_t>
	void reinit_graph()
	{
		// reinitialise the vector without reallocating the full memory
		this->_reallocate_vectors();
	}

	// Compute the graph using the cordonnier metod to solve the depressions
	// template arguments are the connector type, the wrapper input type for
	// topography and the wrapper output type for topography
	template<class topo_t, class out_t>
	out_t compute_graph(
		topo_t& ttopography, // the input topography
		bool only_SD,				 // only computes the single flow graph if true
		bool quicksort			 // computes the MF toposort with a quicksort algo
												 // if true, else uses a BFS-based algorithm
												 // (which one is better depends on many things)
	)
	{
		// Formatting the input to match all the wrappers
		auto topography = format_input<topo_t>(ttopography);
		// Formatting the output topo
		std::vector<fT> faketopo(this->nnodes, 0);

		for (int i = 0; i < this->nnodes; ++i) {
			faketopo[i] = topography[i];
		}

		this->_compute_graph(faketopo, only_SD, quicksort);

		// Finally I format the output topography to the right wrapper
		return format_output<decltype(faketopo), out_t>(faketopo);
	}

	// template<class Connector_t>
	void _compute_graph(
		std::vector<fT>& faketopo,
		bool only_SD,	 // only computes the single flow graph if true
		bool quicksort // computes the MF toposort with a quicksort
									 // algo if true, else uses a BFS-based algorithm
									 // (which one is better depends on many things)
	)
	{

		// Checking if the depression method is cordonnier or node
		bool isCordonnier = this->is_method_cordonnier();
		// bool isDagger = this->is_method_dagger();

		if (this->depression_resolver == DEPRES::priority_full_MFD) {
			this->reinit_graph();
			this->stack = std::vector<size_t>(this->nnodes, 0);
			_PriorityFool(faketopo, this->connector, this->stack, this->Sstack);
			this->topological_sorting_SF();
			return;
		}

		// if(this->depression_resolver == DEPRES::dagger_fill)
		// {
		// 	this->reinit_graph();
		// 	// this->stack = std::vector<size_t>(this->nnodes,0);
		// 	this-> _debug_int = _dagger_fill(this->connector, faketopo);
		// 	this->topological_sorting_dag();
		// 	this->topological_sorting_SF();

		// 	return;
		// }

		// reinit to 0 n+pits
		this->n_pits = 0;

		if (this->depression_resolver == DEPRES::priority_flood) {
			PriorityFlood(faketopo, *this->connector);
		}

		// if the method is not Cordonnier -> apply the other first
		// if(isCordonnier == false && this->depression_resolver != DEPRES::none)
		// {
		// 	// filling the topography with a minimal slope using Wei et al., 2018
		// 	// if(this->depression_resolver == DEPRES::priority_flood)
		// 	faketopo = PriorityFlood(faketopo,  *(this));
		// 	// else
		// 	// 	faketopo = connector->PriorityFlood(faketopo);
		// }

		// Making sure the graph is not inheriting previous values
		this->reinit_graph();

		// Updates the links vector and the Srecs vector by checking each link new
		// elevation
		this->connector->update_links(faketopo);

		this->compute_npits();

		// Compute the topological sorting for single stack
		// Braun and willett 2014 (modified)
		this->topological_sorting_SF();

		if (this->depression_resolver == DEPRES::priority_flood_opti) {
			bool need_reco =
				PriorityFlood_opti(faketopo, *this->connector, this->Sstack);
			if (need_reco) {
				this->connector->update_links(faketopo);
				this->topological_sorting_SF();
			}
		}

		// manages the Cordonnier method if needed
		if (isCordonnier) {
			bool need_recompute;

			if (isCordonnier) {
				// LMRerouter is the class managing the different cordonnier's mthod
				LMRerouter<fT> depsolver;

				if (this->opti_sparse_border)
					depsolver.opti_sparse_border = true;

				depsolver.minimum_slope = this->minimum_slope;
				depsolver.slope_randomness = this->slope_randomness;
				// Execute the local minima solving, return true if rerouting was
				// necessary, meaning that some element needs to be recomputed Note that
				// faketopo are modified in place.
				need_recompute = depsolver.run(
					this->depression_resolver, faketopo, this->connector, this->Sstack);
			}
			// else
			// {
			// 	// need_recompute = simple_depression_hierarchy<fT,
			// std::vector<fT>, Connector_t >(faketopo, this->connector, this->Sstack,
			// this->connector->_Sreceivers); 	need_recompute =
			// simple_depression_solver
			// }

			// Right, if reomputed needs to be
			if (need_recompute) {

				// Re-inverting the Sreceivers into Sdonors
				this->connector->recompute_SF_donors_from_receivers();

				// Recomputing Braun and willett 2014 (modified)
				this->topological_sorting_SF();

				// This is a bit confusing and needs to be changed but filling in done
				// in one go while carving needs a second step here
				if (this->depression_resolver == DEPRES::cordonnier_carve)
					this->carve_topo_v2(faketopo);
				else if (this->depression_resolver == DEPRES::cordonnier_fill)
					this->fill_topo_v2(1e-5, faketopo);

				// And updating the receivers (Wether the Sreceivers are updated or not
				// depends on opt_stst_rerouting)
				if (this->opt_stst_rerouting == false)
					this->connector->update_links(
						faketopo); // up to 30% slower - slightly more accurate for Sgraph

				// My work here is done if only SD is needed
				if (only_SD)
					return;

				// if(this->opt_stst_rerouting)
				this->connector->update_links_MFD_only(
					faketopo); // up to 30% faster - slightly less accurate for Sgraph
			}
		}
		// else if(isDagger)
		// {
		// 	bool need_recompute = simple_depression_solver( this->connector,
		// faketopo, this->Sstack); 	if(need_recompute)
		// 	{
		// 		this->reinit_graph();
		// 		this->connector->update_links(faketopo);
		// 		this->topological_sorting_SF();
		// 	}
		// }

		// else if (this->depression_resolver == DEPRES::priority_flood)
		// {
		// 	PriorityFlood<fT, Connector_t>(faketopo, *(this->connector));
		// 	// faketopo = connector->PriorityFlood(faketopo);
		// 	this->connector->update_links(faketopo); // up to 30% slower - slightly
		// more accurate for Sgraph

		// 	this->connector->recompute_SF_donors_from_receivers();

		// 	this->topological_sorting_SF();

		// }

		// if there is no need to recompute neighbours, then I can only calculate
		// the topological sorting for multiple as the toposort for SF is already
		// done
		if (only_SD == false) {
			if (quicksort)
				this->topological_sorting_quicksort(faketopo);
			else
				this->topological_sorting_dag();
		}
	}

	void compute_npits()
	{
		if (this->debug_mask)
			this->_debug_mask = std::vector<std::uint8_t>(this->nnodes, 0);

		for (int i = 0; i < this->nnodes; ++i) {
			if (this->connector->boundaries.can_out(i) == false &&
					i == this->connector->_Sreceivers[i]) {
				if (this->debug_mask)
					this->_debug_mask[i] = 1;
				++this->n_pits;
			}
		}
	}

	// To adapt to the new graph structure
	// // Compute the graph using the cordonnier metod to solve the depressions
	// // template arguments are the connector type, the wrapper input type for
	// topography and the wrapper output type for topography template<class
	// topo_t, class out_t> out_t compute_graph_timer(
	//   topo_t& ttopography, // the input topography
	//   // the input connector to use (e.g. D8connector)
	//   bool only_SD, // only computes the single flow graph if true
	//   bool quicksort, // computes the MF toposort with a quicksort algo if
	//   true, else uses a BFS-based algorithm (which one is better depends on
	//   many things) int NN // N iteration of the process
	//   )
	// {
	// 	ocarina timer;
	// 	std::map<std::string, std::pair<int,double> > timer_by_stuff;
	// 	timer_by_stuff["initialising_data"] = {0,0.};
	// 	timer_by_stuff["reinit_graph"] = {0,0.};
	// 	timer_by_stuff["local_minima"] = {0,0.};
	// 	timer_by_stuff["update_links"] = {0,0.};
	// 	timer_by_stuff["TS_SF"] = {0,0.};
	// 	timer_by_stuff["TS_MF"] = {0,0.};
	// 	auto topography = format_input<topo_t>(ttopography);
	// 	std::vector<fT> faketopo(this->nnodes,0);

	// 	for (int tn = 0; tn<NN ; ++tn)
	// 	{
	// 		// updating the timer increment
	// 		for(auto& kv:timer_by_stuff)
	// 			++kv.second.first;

	// 		// Formatting the input to match all the wrappers
	// 		timer.tik();
	// 		topography = format_input<topo_t>(ttopography);
	// 		// Formatting the output topo
	// 		faketopo = std::vector<fT>(this->nnodes,0);

	// 		for(int i=0; i<this->nnodes; ++i)
	// 		{
	// 			faketopo[i] = topography[i];
	// 		}
	// 		timer_by_stuff["initialising_data"].second += timer.tok();

	// 		this->_compute_graph_timer(faketopo,only_SD,quicksort,timer,
	// timer_by_stuff);
	// 	}

	// 	fT total = 0;
	// 	for(auto& kv:timer_by_stuff)total+=kv.second.second;

	// 	std::cout << "Result for the timer - mean time on " << NN << "
	// iterations:" << std::endl; 	for(auto& kv:timer_by_stuff)
	// std::cout
	// << kv.first << ": " << kv.second.second/kv.second.first << " ms on average,
	// which is " << 100*kv.second.second/total << " %% of total time" <<
	// std::endl;

	// 	// Finally I format the output topography to the right wrapper
	// 	return format_output<decltype(faketopo), out_t >(faketopo);
	// }

	void activate_opti_sparse_border_cordonnier()
	{
		this->opti_sparse_border = true;
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
					Admin functions	managing attribute and variable initialisations
					=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
					=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	// Helper functions to allocate and reallocate vectors when
	// computing/recomputing the graph
	void _allocate_vectors()
	{
		this->Sstack = std::vector<size_t>(this->nnodes, 0);
	}

	void _reallocate_vectors()
	{
		this->connector->_reallocate_vectors();
		this->Sstack = std::vector<size_t>(this->nnodes, 0);
	}

	// Set the local minima resolver method
	void set_LMR_method(DEPRES method) { this->depression_resolver = method; }

	// Sets the minimum slope to be imposed by the local minima solver
	void set_minimum_slope_for_LMR(fT slope) { this->minimum_slope = slope; }

	// Sets the magnitude of the random noise within the local minima solver
	// (needs to be order of magnitude lower than the actual minimal slope to make
	// sense)
	void set_slope_randomness_for_LMR(fT slope)
	{
		if (slope >= this->minimum_slope)
			throw std::runtime_error(
				"slope randomness cannot be >= to the minimum slope and is even "
				"reccomended to be at least one or two order of magnitude lower");
		this->slope_randomness = slope;
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
	Functions performing topological sorting
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

	// This is my implementation of Braun and willett 2014
	// Slightly modified:
	// - First it is based on the fastscapelib version, which bypass the need of
	// the delta vectors and all by simply applying the recursion directly feeding
	// the stack
	// - Secondly, recursion is not the best practice for most languages so I am
	// using a stack data structure instead Results are identical and performances
	// similar (marginally better, but that is linked to c++ not being a heavy
	// recursion friendly language)
	void topological_sorting_SF()
	{
		// The stack container helper
		std::stack<size_t, std::vector<size_t>> stackhelper;
		// std::vector<bool> isdone(this->nnodes,false);
		// going through all the nodes
		int istack = 0;
		for (int i = 0; i < this->nnodes; ++i) {
			// if they are base level I include them in the stack
			if (this->connector->_Sreceivers[i] == i) {
				stackhelper.emplace(i);
				// ++istack;
			}

			// While I still have stuff in the stack helper
			while (stackhelper.empty() == false) {
				// I get the next node and pop it from the stack helper
				int nextnode = stackhelper.top();
				stackhelper.pop();
				this->Sstack[istack] = nextnode;
				++istack;

				// as well as all its donors which will be processed next
				for (int j = 0; j < this->connector->nSdonors[nextnode]; ++j) {
					stackhelper.emplace(
						this->connector->get_nth_steepest_donors_of_node(nextnode, j));
				}
			}
		}
	}

	// performs the topological sorting using a BFS algorithm
	// It starts from the receivers-less nodes, then uses a queue to add all the
	// donors of these nodes. Each time a donor pops our of the queue, it
	// increment a visited array recording the number of time a node is visited.
	// if the number of visits equals the number of receiver of the node, it is
	// saved in the stack/ The result is a stack of node from the most downstream
	// to the most upstream one template< class Connector_t>
	void topological_sorting_dag()
	{
		// nrecs tracks the number of receivers for each nodes
		std::vector<int> nrecs(this->nnodes, 0);
		// The queueß
		std::stack<int, std::vector<int>> toproc;

		// preparing the stack
		this->stack.clear();
		this->stack.reserve(this->nnodes);

		// Iterating through the links
		// int debug_cpt = 0;
		int node = 0;
		int incr = 0;
		for (int i = 0; i < int(this->connector->links.size()); ++i) {

			if (this->connector->is_link_valid(i)) {
				// Otherwise incrementing the number of receivers for the right link
				++nrecs[this->connector->get_from_links(i, node)];
			}

			++incr;
			if (incr == 4) {
				++node;
				incr = 0;
			}
		}

		// this->connector->debug_print_neighbours(12);

		// Now checking the receiverless nodes and emplacing them in the queueß
		int asdgfkjhl = 0;
		for (int i = 0; i < this->nnodes; ++i) {
			if (nrecs[i] == 0) {
				toproc.emplace(i);
				++asdgfkjhl;
			}
		}

		// then as lon g as there are nodes in the queue:
		auto dons = connector->get_empty_neighbour();
		while (toproc.empty() == false) {
			// getting the next node
			// int next = toproc.front();
			int next = toproc.top();
			// (and popping the node from the queue)
			toproc.pop();
			// if the node is poped out of the queue -> then it is ready to be in the
			// stack Because we are using a FIFO queue, they are sorted correctly in
			// the queue
			this->stack.emplace_back(next);
			// getting the idx of the dons
			int nn = this->connector->get_donors_idx(next, dons);
			for (int td = 0; td < nn; ++td) {
				int d = dons[td];
				// Decrementing the number of receivers (I use it as a visited vector)
				--nrecs[d];
				// if it has reached 0, the rec is ready to be put in the FIFO
				if (nrecs[d] == 0) {
					toproc.emplace(d);
				}
			}
		}

		if (int(this->stack.size()) != this->nnodes) {
			std::cout << "WARNING::Stack->" << this->stack.size() << "/"
								<< this->nnodes << std::endl;
			throw std::runtime_error("STACK ISSUE");
		}

		// std::cout << debug_cpt << " <-- Nlinks invalids" << std::endl;
	}

	// Multiple flow topological sorting using quicksort algorithm
	// Topography is simply sorted by absolute elevation keeping track of the
	// indices
	template<class topo_t>
	void topological_sorting_quicksort(topo_t& ttopography)
	{
		// Formatting hte topographic input from the wrapper
		auto topography = format_input<topo_t>(ttopography);

		// Dortng by index
		auto yolo = sort_indexes(topography);

		// the sorted index is now the stack
		this->stack = std::move(yolo);
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
	Functions affecting topography from corrected receivers
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

	/// this function enforces minimal slope
	/// It starts from the most upstream part of the landscapes and goes down
	/// following the Sreceiver route it carve on the go, making sure the
	/// topography of a receiver is lower
	template<class topo_t>
	void carve_topo_v2(topo_t& topography)
	{

		// Traversing the (SFD) stack on the reverse direction
		for (int i = this->nnodes - 1; i >= 0; --i) {
			// Getting the node
			int node = this->Sstack[i];
			// Checking its validiyt AND if it is not a base level
			if (this->connector->flow_out_model(node))
				continue;
			// Getting the single receiver info
			int rec = this->connector->_Sreceivers[node];
			// Checking the difference in elevation
			fT dz = topography[node] - topography[rec];
			// if the difference in elevation is bellow 0 I need to carve
			if (dz <= 0) {
				// And I do ! Note that I add some very low-grade randomness to avoid
				// flat links
				topography[rec] =
					topography[node] - this->minimum_slope +
					this->randu->get() * this->slope_randomness; // * d2rec;
			}
		}
	}

	/// Opposite of the above function
	/// It starts from the most dowstream nodes and climb its way up.
	/// when a node is bellow its receiver, we correct the slope
	template<class topo_t>
	void fill_topo_v2(fT slope, topo_t& topography)
	{
		for (int i = 0; i < this->nnodes; ++i) {
			int node = this->Sstack[i];
			if (this->connector->flow_out_model(node))
				continue;

			int rec = this->connector->_Sreceivers[node];
			fT dz = topography[node] - topography[rec];

			if (dz <= 0) {
				topography[node] = std::nextafter(
					topography[rec], std::numeric_limits<fT>::max()); // * d2rec;
			}
		}
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
					Accessing receivers/donors/... for a single node
					=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
					=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	// DEbugging function checking the validity of the single flow stack
	// It tests a number of checks like double values or if it has the right
	// number of nodes
	bool is_Sstack_full()
	{
		// IS the stack the right size
		if (int(this->Sstack.size()) != this->nnodes) {
			std::cout << "stack size (" << this->Sstack.size() << ") is invalid."
								<< std::endl;
			return false;
		}

		// Are there values appearing multiple times
		std::vector<int> ntimenodes(this->nnodes, 0);
		for (auto v : this->Sstack) {
			++ntimenodes[v];
		}
		int n_0 = 0, n_p1 = 0;
		for (int i = 0; i < this->nnodes; ++i) {
			if (ntimenodes[i] == 0)
				++n_0;
			else if (ntimenodes[i] > 1)
				++n_p1;
		}

		if (n_0 > 0 || n_p1 > 0) {
			std::cout << "Stack issue: " << n_p1
								<< " nodes appearing more than once and " << n_0
								<< " nodes not appearing" << std::endl;
			return false;
		}

		// IS the stack in the right order in accordance with the SStack array
		std::vector<bool> isdone(this->nnodes, false);
		for (int i = this->nnodes - 1; i >= 0; --i) {
			auto v = this->Sstack[i];
			isdone[v] = true;
			if (int(v) != this->connector->_Sreceivers[v]) {
				if (isdone[this->connector->_Sreceivers[v]]) {
					std::cout << "Receiver processed before node stack is fucked"
										<< std::endl;
					return false;
				}
			}
		}

		return true;
	}

	std::vector<bool> has_Srecs()
	{
		std::vector<bool> haSrecs(this->nnodes, true);
		for (int i = 0; i < this->nnodes; ++i) {
			if (this->connector->_Sreceivers[i] == i)
				haSrecs[i] = false;
		}
		return haSrecs;
	}

	bool is_method_cordonnier()
	{
		if (this->depression_resolver == DEPRES::cordonnier_carve ||
				this->depression_resolver == DEPRES::cordonnier_fill ||
				DEPRES::cordonnier_simple == this->depression_resolver)
			return true;
		else
			return false;
	}

	bool is_method_dagger()
	{
		if (this->depression_resolver == DEPRES::dagger_carve ||
				this->depression_resolver == DEPRES::dagger_fill ||
				DEPRES::dagger_simple == this->depression_resolver)
			return true;
		else
			return false;
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
					Functions to access sets of nodes draining to/from a single point
					=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
					=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	template<class out_t>
	out_t get_all_nodes_upstream_of(int node, bool only_SD = false)
	{
		std::vector<int> out;
		out = this->_get_all_nodes_upstream_of_using_graph(node, only_SD);
		return format_output<decltype(out), out_t>(out);
	}

	std::vector<int> _get_all_nodes_upstream_of_using_graph(int node,
																													bool only_SD)
	{

		// Formatting the output vector
		std::vector<int> out;
		out.reserve(round(this->nnodes / 4));
		// Creating a visited vector tracking which node has been visited
		std::vector<bool> vis(this->nnodes, false);
		// marking the initial node as true
		vis[node] = true;

		for (auto tnode : this->Sstack) {
			// ignoring the not ode and outlets
			if (this->connector->flow_out_or_pit(tnode) == false) {
				// Getting the receiver
				int rec = this->connector->_Sreceivers[tnode];
				// checkng if receiver is visited but not node
				if (vis[rec]) {
					// current noer is visited
					vis[tnode] = true;
					// and is draining to this node
					out.emplace_back(tnode);
				}
			}
		}

		// if only steepest descent is needed, we stop there
		if (only_SD)
			return out;

		// else, we have to use a queue to add all the donors
		std::queue<int> toproc;

		// first checking if all the steepest descent nodes I already have there
		// have a not-SD donor
		auto donors = connector->get_empty_neighbour();
		for (auto v : out) {
			// gettign the donors
			int nn = this->connector->get_donors_idx(v, donors);
			// for all donors of dat nod
			for (int td = 0; td < nn; ++td) {
				int d = donors[td];
				// if not visited
				if (vis[d] == false) {
					// becomes visited
					vis[d] = true;
					// and I emplace it in the queues
					toproc.emplace(d);
				}
			}
		}

		// once this is done I work until the queue is empty
		while (toproc.empty() == false) {
			// getting the next node in line
			int next = toproc.front();
			toproc.pop();
			// recording it as draining to the original node
			out.emplace_back(next);
			// getting all its donors
			int nn = this->connector->get_donors_idx(next, donors);
			// for all donors of dat nod
			for (int td = 0; td < nn; ++td) {
				int d = donors[td];
				// if not visited including it (see above)
				if (vis[d] == false) {
					vis[d] = true;
					toproc.emplace(d);
				}
			}
		}

		// LOK queue is empty and I have everything I need

		return out;
	}

	template<class out_t>
	out_t get_all_nodes_downstream_of(int node, bool only_SD = false)
	{
		std::vector<int> out;
		out = this->_get_all_nodes_downstream_of_using_graph(node, only_SD);

		return format_output<decltype(out), out_t>(out);
	}

	std::vector<int> _get_all_nodes_downstream_of_using_graph(int node,
																														bool only_SD)
	{

		// Formatting the output vector
		std::vector<int> out;
		out.reserve(round(this->nnodes / 4));
		// Creating a visited vector tracking which node has been visited
		std::vector<bool> vis(this->nnodes, false);
		// marking the initial node as true
		vis[node] = true;

		for (int i = this->nnodes - 1; i >= 0; --i) {
			int tnode = this->Sstack[i];
			// ignoring the not ode and outlets
			if (this->connector->flow_out_or_pit(tnode) == false) {
				// Getting the receiver
				int rec = this->connector->_Sreceivers[tnode];
				// checkng if receiver is visited but not node
				if (vis[node] && rec != node && tnode != node) {
					// current noer is visited
					vis[rec] = true;
					// and is draining to this node
					out.emplace_back(tnode);
				}
			}
		}

		// if only steepest descent is needed, we stop there
		if (only_SD)
			return out;

		// else, we have to use a queue to add all the receivers
		std::queue<int> toproc;

		auto receivers = this->connector->get_empty_neighbour();

		// first checking if all the steepest descent nodes I already have there
		// have a not-SD rec
		for (auto v : out) {
			// gettign the receivers
			int nn = this->connector->get_receivers_idx(v, receivers);
			// for all receivers of dat nod
			for (int tr = 0; tr < nn; ++tr) {
				int r = receivers[tr];
				// if not visited
				if (vis[r] == false) {
					// becomes visited
					vis[r] = true;
					// and I emplace it in the queues
					toproc.emplace(r);
				}
			}
		}

		// once this is done I work until the queue is empty
		while (toproc.empty() == false) {
			// getting the next node in line
			int next = toproc.front();
			toproc.pop();
			// recording it as draining to the original node
			out.emplace_back(next);
			// getting all its receivers
			int nn = this->connector->get_receivers_idx(next, receivers);
			for (int tr = 0; tr < nn; ++tr) {
				int r = receivers[tr];
				// if not visited including it (see above)
				if (vis[r] == false) {
					vis[r] = true;
					toproc.emplace(r);
				}
			}
		}

		// LOK queue is empty and I have everything I need

		return out;
	}

	std::vector<int> _get_flow_acc()
	{
		std::vector<int> flowacc(this->nnodes, 0);
		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->Sstack[i];
			if (this->connector->boundaries.no_data(node) == false)
				continue;

			int rec = this->connector->_Sreceivers[node];
			if (this->connector->flow_out_or_pit(node) == false) {
				flowacc[rec] += flowacc[node] + 1;
			}
		}

		return flowacc;
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
					Functions to access to bulk receivers/links/donors/...
					=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
					=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	template<class out_t>
	out_t get_SFD_stack()
	{
		return format_output<std::vector<size_t>, out_t>(this->Sstack);
	}

	template<class out_t>
	out_t get_MFD_stack()
	{
		return format_output<std::vector<size_t>, out_t>(this->stack);
	}

	template<class topo_t>
	std::vector<fT> _get_max_val_link_array(topo_t& array)
	{
		std::vector<fT> tmax(this->nnodes, 0);
		for (size_t i = 0; i < this->connector->links.size(); ++i) {
			if (this->connector->is_link_valid(i) == false)
				continue;
			int go = this->connector->get_from_links(i);
			if (tmax[go] < array[i])
				tmax[go] = array[i];
		}
		return tmax;
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
	Functions to propagate signal upstream and downstream
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	template<class out_t>
	out_t accumulate_constant_downstream_SFD(fT var)
	{
		std::vector<fT> out = this->_accumulate_constant_downstream_SFD(var);
		return format_output<decltype(out), out_t>(out);
	}

	std::vector<fT> _accumulate_constant_downstream_SFD(fT var)
	{
		fT bal = 0, tot = 0;
		std::vector<fT> out(this->nnodes, 0);
		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->Sstack[i];
			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += var;
			bal += var;
			tot += var;

			if (this->connector->flow_out_or_pit(node)) {
				// if(this->connector->is_on_dem_edge(node) == false)
				// {
				// 	std::cout << "WARNING_DEBUG_45b::FLOW OUTS - NOT ON DEM EDGE "
				// << std::endl;; 	std::cout << node << "||" <<
				// this->connector->_Sreceivers[node] << std::endl;
				// 	this->connector->debug_print_neighbours(node);

				// }
				if (this->connector->flow_out_model(node))
					bal -= out[node];

				continue;
			}

			out[this->connector->_Sreceivers[node]] += out[node];
		}

		return out;
	}

	std::vector<fT> _accumulate_constant_downstream_area_SFD(fT var)
	{

		std::vector<fT> out(this->nnodes, 0);

		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->Sstack[i];
			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += var + this->connector->get_area_at_node(node);

			if (this->connector->flow_out_or_pit(node)) {
				continue;
			}

			out[this->connector->_Sreceivers[node]] += out[node];
		}

		return out;
	}

	template<class out_t, class topo_t>
	out_t accumulate_variable_downstream_SFD(topo_t& tvar)
	{
		auto var = format_input<topo_t>(tvar);
		std::vector<fT> out = this->_accumulate_variable_downstream_SFD(var);
		return format_output<decltype(out), out_t>(out);
	}

	template<class topo_t>
	std::vector<fT> _accumulate_variable_downstream_SFD(topo_t& var)
	{
		std::vector<fT> out(this->nnodes, 0);
		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->Sstack[i];

			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += var[node];

			if (this->connector->flow_out_or_pit(node))
				continue;

			out[this->connector->_Sreceivers[node]] += out[node];
		}

		return out;
	}

	template<class topo_t>
	std::vector<fT> _accumulate_variable_downstream_area_SFD(topo_t& var)
	{
		std::vector<fT> out(this->nnodes, 0);
		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->Sstack[i];

			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += var[node] + this->connector->get_area_at_node(node);

			if (this->connector->flow_out_or_pit(node))
				continue;

			out[this->connector->_Sreceivers[node]] += out[node];
		}

		return out;
	}

	template<class topo_t, class out_t>
	out_t accumulate_constant_downstream_MFD(topo_t& tweights, fT var)
	{
		auto weights = format_input<topo_t>(tweights);
		std::vector<fT> out =
			this->_accumulate_constant_downstream_MFD(weights, var);
		return format_output<decltype(out), out_t>(out);
	}

	template<class topo_t>
	std::vector<fT> _accumulate_constant_downstream_MFD(topo_t& weights, fT var)
	{
		std::vector<fT> out(this->nnodes, 0);
		auto reclinks = this->connector->get_empty_neighbour();
		for (int i = this->nnodes - 1; i >= 0; --i) {

			int node = this->stack[i];
			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += var;

			if (this->connector->flow_out_or_pit(node))
				continue;

			int nn = this->connector->get_receivers_idx_links(node, reclinks);
			for (int ttl = 0; ttl < nn; ++ttl) {
				int ti = reclinks[ttl];
				if (this->connector->is_link_valid(ti)) {
					int rec = this->connector->get_to_links(ti);
					out[rec] += out[node] * weights[ti];
				}
			}
		}

		return out;
	}

	template<class topo_t>
	std::vector<fT> _accumulate_constant_downstream_area_MFD(topo_t& weights,
																													 fT var)
	{
		std::vector<fT> out(this->nnodes, 0);
		auto reclinks = this->connector->get_empty_neighbour();
		for (int i = this->nnodes - 1; i >= 0; --i) {

			int node = this->stack[i];
			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += var + this->connector->get_area_at_node(node);

			if (this->connector->flow_out_or_pit(node))
				continue;

			int nn = this->connector->get_receivers_idx_links(node, reclinks);
			for (int ttl = 0; ttl < nn; ++ttl) {
				int ti = reclinks[ttl];
				if (this->connector->is_link_valid(ti)) {
					int rec = this->connector->get_to_links(ti);
					out[rec] += out[node] * weights[ti];
				}
			}
		}

		return out;
	}

	template<class topo_t, class out_t>
	out_t accumulate_variable_downstream_MFD(topo_t& tweights, topo_t& tvar)
	{
		auto weights = format_input<topo_t>(tweights);
		auto var = format_input<topo_t>(tvar);
		std::vector<fT> out =
			this->_accumulate_variable_downstream_MFD(weights, var);
		return format_output<decltype(out), out_t>(out);
	}

	template<class topo_t>
	std::vector<fT> _accumulate_variable_downstream_MFD(topo_t& weights,
																											topo_t& var)
	{
		std::vector<fT> out(this->nnodes, 0);
		auto reclinks = this->connector->get_empty_neighbour();
		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->stack[i];
			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += var[node];

			if (this->connector->flow_out_or_pit(node))
				continue;

			int nn = this->connector->get_receivers_idx_links(node, reclinks);
			for (int tr = 0; tr < nn; ++tr) {
				int ti = reclinks[tr];
				int rec = this->connector->get_to_links(ti);
				if (connector->is_in_bound(rec))
					out[rec] += out[node] * weights[ti];
			}
		}

		return out;
	}

	template<class topo_t>
	std::vector<fT> _accumulate_variable_downstream_area_MFD(topo_t& weights,
																													 topo_t& var)
	{
		std::vector<fT> out(this->nnodes, 0);
		auto reclinks = this->connector->get_empty_neighbour();
		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->stack[i];
			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += var[node] + this->connector->get_area_at_node(node);

			if (this->connector->flow_out_or_pit(node))
				continue;

			int nn = this->connector->get_receivers_idx_links(node, reclinks);
			for (int tr = 0; tr < nn; ++tr) {
				int ti = reclinks[tr];
				int rec = this->connector->get_to_links(ti);
				if (connector->is_in_bound(rec))
					out[rec] += out[node] * weights[ti];
			}
		}

		return out;
	}

	template<class out_t>
	out_t get_drainage_area_SFD()
	{
		std::vector<fT> out = this->_get_drainage_area_SFD();
		return format_output<decltype(out), out_t>(out);
	}

	std::vector<fT> _get_drainage_area_SFD()
	{
		std::vector<fT> out(this->nnodes, 0.);
		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->Sstack[i];

			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += this->connector->get_area_at_node(node);
			int rec = this->connector->_Sreceivers[node];
			if (node != rec)
				out[rec] += out[node];
		}
		return out;
	}

	template<class out_t, class topo_t>
	out_t get_drainage_area_MFD(topo_t& ttopo, fT exp_slope)
	{
		auto topo = format_input(ttopo);
		std::vector<fT> out = this->_get_drainage_area_MFD(topo, exp_slope);
		return format_output<decltype(out), out_t>(out);
	}

	template<class topo_t>
	std::vector<fT> _get_drainage_area_MFD(topo_t& topo, fT exp_slope)
	{
		std::vector<fT> out(this->nnodes, 0.);

		auto lgrad = this->connector->_get_links_gradient(topo, 1e-6);
		std::vector<fT> weights(lgrad.size(), 0.);

		this->connector->_get_link_weights_exp(weights, lgrad, exp_slope);

		auto reclinks = this->connector->get_empty_neighbour();
		for (int i = this->nnodes - 1; i >= 0; --i) {
			int node = this->stack[i];
			if (this->connector->boundaries.no_data(node))
				continue;

			out[node] += this->connector->get_area_at_node(node);

			if (this->connector->flow_out_or_pit(node))
				continue;

			int nn = this->connector->get_receivers_idx_links(node, reclinks);
			for (int tr = 0; tr < nn; ++tr) {
				int ti = reclinks[tr];
				int rec = this->connector->get_to_links(ti);
				if (connector->is_in_bound(rec))
					out[rec] += out[node] * weights[ti];
			}
		}

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
	Links/node utility functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

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
	Distance utility functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

	/// this function computes the flow distance from model outlets using the
	/// siungle direction graph
	template<class out_t>
	out_t get_SFD_distance_from_outlets()
	{
		std::vector<fT> distfromoutlet(this->nnodes, 0.);
		this->_get_SFD_distance_from_outlets(distfromoutlet);
		return format_output<decltype(distfromoutlet), out_t>(distfromoutlet);
	}

	void _get_SFD_distance_from_outlets(std::vector<fT>& distfromoutlet)
	{
		// just iterating through the Sstack in the upstream direction adding dx to
		// the receiver
		for (int i = 0; i < this->nnodes; ++i) {
			// next node in the stack
			int node = this->Sstack[i];
			// checking if active
			if (this->connector->flow_out_or_pit(node))
				continue;
			// Getting the receiver
			int rec = this->connector->_Sreceivers[node];
			// And integrating the distance from outlets
			distfromoutlet[node] =
				distfromoutlet[rec] + this->connector->Sdistance2receivers[node];
		}
	}

	/// this function computes the flow distance from model outlets using the
	/// siungle direction graph
	template<class out_t>
	out_t get_SFD_min_distance_from_sources()
	{
		std::vector<fT> distfromsources(this->nnodes, 0.);
		this->_get_SFD_min_distance_from_sources(distfromsources);
		return format_output<decltype(distfromsources), out_t>(distfromsources);
	}

	void _get_SFD_min_distance_from_sources(std::vector<fT>& distfromsources)
	{
		// just iterating through the Sstack in the upstream direction adding dx to
		// the receiver
		for (int i = this->nnodes - 1; i >= 0; --i) {
			// next node in the stack
			int node = this->Sstack[i];
			// checking if active
			if (this->connector->flow_out_or_pit(node))
				continue;
			int rec = this->connector->_Sreceivers[node];
			if (distfromsources[rec] == 0 ||
					distfromsources[rec] >
						distfromsources[node] + this->connector->Sdistance2receivers[node])
				distfromsources[rec] =
					distfromsources[node] + this->connector->Sdistance2receivers[node];
		}
	}

	/// this function computes the flow distance from model outlets using the
	/// siungle direction graph
	template<class out_t>
	out_t get_SFD_max_distance_from_sources()
	{
		std::vector<fT> distfromsources(this->nnodes, 0.);
		this->_get_SFD_max_distance_from_sources(distfromsources);
		return format_output<decltype(distfromsources), out_t>(distfromsources);
	}

	void _get_SFD_max_distance_from_sources(std::vector<fT>& distfromsources)
	{
		// just iterating through the Sstack in the upstream direction adding dx to
		// the receiver
		for (int i = this->nnodes - 1; i >= 0; --i) {
			// next node in the stack
			int node = this->Sstack[i];
			// checking if active
			if (this->connector->flow_out_or_pit(node))
				continue;
			int rec = this->connector->_Sreceivers[node];
			if (distfromsources[rec] == 0 ||
					distfromsources[rec] <
						distfromsources[node] + this->connector->Sdistance2receivers[node])
				distfromsources[rec] =
					distfromsources[node] + this->connector->Sdistance2receivers[node];
		}
	}

	/// this function computes the flow distance from model outlets using the
	/// siungle direction graph
	template<class out_t>
	out_t get_MFD_max_distance_from_sources()
	{
		std::vector<fT> distfromsources(this->nnodes, 0.);
		this->_get_MFD_max_distance_from_sources(distfromsources);
		return format_output<decltype(distfromsources), out_t>(distfromsources);
	}

	void _get_MFD_max_distance_from_sources(std::vector<fT>& distfromsources)
	{
		// just iterating through the Sstack in the upstream direction adding dx to
		// the receiver
		auto receilink = this->connector->get_empty_neighbour();
		for (int i = this->nnodes - 1; i >= 0; --i) {
			// next node in the stack
			int node = this->stack[i];
			// checking if active
			if (this->connector->flow_out_or_pit(node))
				continue;

			int nl = this->connector->get_receivers_idx_links(node, receilink);

			for (int tl = 0; tl < nl; ++tl) {
				int rec = this->connector->get_to_links(receilink[tl]);
				fT dx = this->connector->get_dx_from_links_idx(receilink[tl]);
				if (distfromsources[rec] == 0 ||
						distfromsources[rec] < distfromsources[node] + dx)
					distfromsources[rec] = distfromsources[node] + dx;
			}
		}
	}

	/// this function computes the flow distance from model outlets using the
	/// siungle direction graph
	template<class out_t>
	out_t get_MFD_min_distance_from_sources()
	{
		std::vector<fT> distfromsources(this->nnodes, 0.);
		this->_get_MFD_min_distance_from_sources(distfromsources);
		return format_output<decltype(distfromsources), out_t>(distfromsources);
	}

	void _get_MFD_min_distance_from_sources(std::vector<fT>& distfromsources)
	{
		// just iterating through the Sstack in the upstream direction adding dx to
		// the receiver
		auto receilink = this->connector->get_empty_neighbour();
		for (int i = this->nnodes - 1; i >= 0; --i) {
			// next node in the stack
			int node = this->stack[i];
			// checking if active
			if (this->connector->flow_out_or_pit(node))
				continue;

			int nl = this->connector->get_receivers_idx_links(node, receilink);

			for (int tl = 0; tl < nl; ++tl) {
				int rec = this->connector->get_to_links(receilink[tl]);
				fT dx = this->connector->get_dx_from_links_idx(receilink[tl]);
				if (distfromsources[rec] == 0 ||
						distfromsources[rec] > distfromsources[node] + dx)
					distfromsources[rec] = distfromsources[node] + dx;
			}
		}
	}

	/// this function computes the flow distance from model outlets using the
	/// siungle direction graph
	template<class out_t>
	out_t get_MFD_max_distance_from_outlets()
	{
		std::vector<fT> distfromoutlets(this->nnodes, 0.);
		this->_get_MFD_max_distance_from_outlets(distfromoutlets);
		return format_output<decltype(distfromoutlets), out_t>(distfromoutlets);
	}

	void _get_MFD_max_distance_from_outlets(std::vector<fT>& distfromoutlets)
	{
		// just iterating through the Sstack in the upstream direction adding dx to
		// the receiver
		auto receilink = this->connector->get_empty_neighbour();
		for (int i = 0; i < this->nnodes; ++i) {

			// next node in the stack
			int node = this->stack[i];

			// checking if active
			if (this->connector->flow_out_or_pit(node))
				continue;

			int nl = this->connector->get_receivers_idx_links(node, receilink);

			for (int tl = 0; tl < nl; ++tl) {
				int rec = this->connector->get_to_links(receilink[tl]);
				fT dx = this->connector->get_dx_from_links_idx(receilink[tl]);
				if (distfromoutlets[node] == 0 ||
						distfromoutlets[node] < distfromoutlets[rec] + dx)
					distfromoutlets[node] = distfromoutlets[rec] + dx;
			}
		}
	}

	/// this function computes the flow distance from model outlets using the
	/// siungle direction graph
	template<class out_t>
	out_t get_MFD_min_distance_from_outlets()
	{
		std::vector<fT> distfromoutlets(this->nnodes, 0.);
		this->_get_MFD_min_distance_from_outlets(distfromoutlets);
		return format_output<decltype(distfromoutlets), out_t>(distfromoutlets);
	}

	void _get_MFD_min_distance_from_outlets(std::vector<fT>& distfromoutlets)
	{
		// just iterating through the Sstack in the upstream direction adding dx to
		// the receiver
		auto receilink = this->connector->get_empty_neighbour();
		for (int i = 0; i < this->nnodes; ++i) {

			// next node in the stack
			int node = this->stack[i];

			// checking if active
			if (this->connector->flow_out_or_pit(node))
				continue;

			int nl = this->connector->get_receivers_idx_links(node, receilink);

			for (int tl = 0; tl < nl; ++tl) {
				int rec = this->connector->get_to_links(receilink[tl]);
				fT dx = this->connector->get_dx_from_links_idx(receilink[tl]);
				if (distfromoutlets[node] == 0 ||
						distfromoutlets[node] > distfromoutlets[rec] + dx)
					distfromoutlets[node] = distfromoutlets[rec] + dx;
			}
		}
	}

	/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
									. - ~ ~ ~ - .
				_     .-~               ~- .
				\ `..~                       ` .
					}  }              /       \   \
(\   \\ \~^..'                 |       }  \
\`.-~  o      /       }       |        /  \
(__          |       /        |       /    `.
`- - ~ ~ -._|      /_ - ~ ~ ^|   ____/- _    `.
				|     /          |     /     ~-.     ~- _
				|_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Watershed labelling functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

	template<class out_t>
	out_t get_SFD_basin_labels()
	{

		std::vector<int> baslab = this->_get_SFD_basin_labels();
		return format_output<decltype(baslab), out_t>(baslab);
	}

	std::vector<int> _get_SFD_basin_labels()
	{
		int nobasin_value = -1;
		std::vector<int> baslab(this->nnodes, nobasin_value);
		int label = -1;
		for (int i = 0; i < this->nnodes; ++i) {
			int node = this->Sstack[i];
			if (this->connector->boundaries.no_data(node))
				continue;
			if (this->connector->flow_out_or_pit(node))
				++label;
			baslab[node] = label;
		}
		return baslab;
	}

	// EXPERIMENTAL SECTION. Trying stuffs. Without breaking the compiler.

	// // Compute the graph using the cordonnier metod to solve the depressions
	// // template arguments are the connector type, the wrapper input type for
	// topography and the wrapper output type for topography template<class
	// topo_t, class int_t, class out_t> out_t compute_graph_exp(
	//   topo_t& ttopography, // the input topography
	//   bool only_SD, // only computes the single flow graph if true
	//   in_t& tstaters
	//   )
	// {
	// 	// Formatting the input to match all the wrappers
	// 	auto topography = format_input<topo_t>(ttopography);
	// 	auto starters = format_input<in_t>(tstaters);

	// 	// Formatting the output topo
	// 	std::vector<fT> faketopo(this->nnodes,0);

	// 	for(int i=0; i<this->nnodes; ++i)
	// 	{
	// 		faketopo[i] = topography[i];
	// 	}

	// 	this->_compute_graph_exp(faketopo,only_SD,starters);

	// 	// Finally I format the output topography to the right wrapper
	// 	return format_output<decltype(faketopo), out_t >(faketopo);
	// }

	// // template<class Connector_t>
	// void _compute_graph_exp(
	// 	std::vector<fT>& faketopo,
	//   bool only_SD, // only computes the single flow graph if true
	//   std::vector<int>& faketopo, // computes the MF toposort with a quicksort
	//   algo if true, else uses a BFS-based algorithm (which one is better
	//   depends on many things)
	// 	)
	// {

	// 	// Checking if the depression method is cordonnier or node
	// 	bool isCordonnier = this->is_method_cordonnier();
	// 	bool isDagger = this->is_method_dagger();

	// 	//reinit to 0 n+pits
	// 	this->n_pits = 0;

	// 	if(this->depression_resolver == DEPRES::priority_flood)
	// 	{
	// 		PriorityFlood(faketopo,*this->connector);
	// 	}

	// 	// if the method is not Cordonnier -> apply the other first
	// 	// if(isCordonnier == false && this->depression_resolver !=
	// DEPRES::none)
	// 	// {
	// 	// 	// filling the topography with a minimal slope using Wei et al.,
	// 2018
	// 	// 	// if(this->depression_resolver == DEPRES::priority_flood)
	// 	// 	faketopo = PriorityFlood(faketopo,  *(this));
	// 	// 	// else
	// 	// 	// 	faketopo = connector->PriorityFlood(faketopo);
	// 	// }

	// 	// Making sure the graph is not inheriting previous values
	// 	this->reinit_graph();

	// 	// Updates the links vector and the Srecs vector by checking each link
	// new elevation 	this->connector->update_links(faketopo);

	// 	this->compute_npits();

	// 	// Compute the topological sorting for single stack
	// 	// Braun and willett 2014 (modified)
	// 	this->topological_sorting_SF();

	// 	if(this->depression_resolver == DEPRES::priority_flood_opti)
	// 	{
	// 		bool need_reco = PriorityFlood_opti(faketopo,*this->connector,
	// this->Sstack); 		if(need_reco)
	// 		{
	// 			this->connector->update_links(faketopo);
	// 			this->topological_sorting_SF();
	// 		}
	// 	}

	// 	// manages the Cordonnier method if needed
	// 	if(isCordonnier || isDagger)
	// 	{
	// 		bool need_recompute;

	// 		if(isCordonnier)
	// 		{
	// 			// LMRerouter is the class managing the different
	// cordonnier's mthod 			LMRerouter<fT> depsolver;

	// 			if(this->opti_sparse_border)
	// 				depsolver.opti_sparse_border = true;

	// 			depsolver.minimum_slope = this->minimum_slope;
	// 			depsolver.slope_randomness = this->slope_randomness;
	// 			// Execute the local minima solving, return true if
	// rerouting was necessary, meaning that some element needs to be recomputed
	// 			// Note that faketopo are modified in place.
	// 			need_recompute =
	// depsolver.run(this->depression_resolver, faketopo, this->connector,
	// this->Sstack);
	// 		}
	// 		else
	// 		{
	// 			need_recompute = simple_depression_hierarchy<fT,
	// std::vector<fT>, Connector_t >(faketopo, this->connector, this->Sstack,
	// this->connector->_Sreceivers);
	// 		}

	// 		// Right, if reomputed needs to be
	// 		if(need_recompute)
	// 		{

	// 			// Re-inverting the Sreceivers into Sdonors
	// 			this->connector->recompute_SF_donors_from_receivers();

	// 			// Recomputing Braun and willett 2014 (modified)
	// 			this->topological_sorting_SF();

	// 			// This is a bit confusing and needs to be changed but
	// filling in done in one go while carving needs a second step here
	// 			if(this->depression_resolver ==
	// DEPRES::cordonnier_carve)
	// this->carve_topo_v2(faketopo); else if (this->depression_resolver ==
	// DEPRES::cordonnier_fill) 				this->fill_topo_v2(1e-5,
	// faketopo);

	// 			// And updating the receivers (Wether the Sreceivers are
	// updated or not depends on opt_stst_rerouting)
	// if(this->opt_stst_rerouting == false)
	// 				this->connector->update_links(faketopo); // up
	// to 30% slower
	// - slightly more accurate for Sgraph

	// 			// My work here is done if only SD is needed
	// 			if(only_SD)
	// 				return;

	// 			// if(this->opt_stst_rerouting)
	// 			this->connector->update_links_MFD_only(faketopo); // up
	// to 30% faster - slightly less accurate for Sgraph

	// 		}
	// 	}
	// 	// else if (this->depression_resolver == DEPRES::priority_flood)
	// 	// {
	// 	// 	PriorityFlood<fT, Connector_t>(faketopo, *(this->connector));
	// 	// 	// faketopo = connector->PriorityFlood(faketopo);
	// 	// 	this->connector->update_links(faketopo); // up to 30% slower -
	// slightly more accurate for Sgraph

	// 	// 	this->connector->recompute_SF_donors_from_receivers();

	// 	// 	this->topological_sorting_SF();

	// 	// }

	// 	// if there is no need to recompute neighbours, then I can only
	// calculate the topological sorting
	// 	// for multiple as the toposort for SF is already done
	// 	if (only_SD == false)
	// 	{
	// 		this->topological_sorting_dag();
	// 	}

	// }

	// end of the graph class
};

// end of Dagger namespace
}; // namespace DAGGER

#endif
