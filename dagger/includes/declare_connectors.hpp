#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_D8connector(py::module& m, std::string typestr)
{

	py::class_<Connector8<int, double>>(m, "Connector8")

		.def(py::init<int, int, double, double, Hermes<int, double>&>(),
				 py::arg("nx"),
				 py::arg("ny"),
				 py::arg("dx"),
				 py::arg("dy"),
				 py::arg("dbag"))

		.def(
			"init",
			&Connector8<int, double>::init,
			R"pdoc(Initialises the connector inner structure. Will affects the data bag with neighbour codes.)pdoc")
		.def("reinit",
				 &Connector8<int, double>::reinit,
				 R"pdoc(Internal function resetting receivers/donors/...)pdoc")
		.def("set_flowtopo",
				 &Connector8<int, double>::set_flowtopo,
				 R"pdoc(Set the flow topology parameter (pre-init).)pdoc")
		.def("set_condou",
				 &Connector8<int, double>::set_condou,
				 R"pdoc(Set the boundary parameter (pre-init).)pdoc")
		.def(
			"compute",
			&Connector8<int, double>::compute,
			R"pdoc(Main function recomputing the graph to match the current surface array in the data bag.)pdoc")
		.def(
			"PFcompute_all",
			&Connector8<int, double>::PFcompute_all,
			R"pdoc(Temp function while I code the graph v2 computing the full graph using priority flood.)pdoc")

		;

	// D8Connector: linked-based graph
	py::class_<D8connector<FLOATING_POINT_DAGGER>>(m, typestr.c_str(), R"pbdoc(
	D8 regular grid connector D8N

	Description:
	------------

	D8N class: the built-in connector for regular grids. D8N regular grids are the most used type of grids, like most global DEMs for example.
	They are defined by a number of rows (ny) and a number of columns (nx) with associated spacing between nodes in both directions (dy,dx).
	The D8N constructor defaults to ``4edges`` boundary preset and precompiles neighbouring. Neighbouring and indexing is optimised by geometrical relationships:
	Node index is simply a flat index in the row major direction (idx = row_id * nx + col_id) and the link ID for the top-right, right, bottomright and bottom neighbours are respectiveley idx *4, idx*4 +1, idx * 4+2 and idx*4 + 3

	The constructor only needs basics geometrical information.

	Parameters:
	-----------

		* nx (int): number of nodes in the X direction (columns)
		* ny (int): number of nodes in the Y direction (rows)
		* dx (float64): distance between nodes in the x directions
		* dy (float64): distance between nodes in the y directions
		* x_min (float64): X coordinates of the bottom left corner
		* y_min (float64): Y coordinates of the top left corner


	Authors:
	--------
	B.G.

	)pbdoc")

		.def(py::init<int,
									int,
									FLOATING_POINT_DAGGER,
									FLOATING_POINT_DAGGER,
									FLOATING_POINT_DAGGER,
									FLOATING_POINT_DAGGER>(),
				 py::arg("nx"),
				 py::arg("ny"),
				 py::arg("dx"),
				 py::arg("dy"),
				 py::arg("x_min"),
				 py::arg("x_max"),
				 R"pbdoc(
	D8 regular grid connector D8N

	Description:
	------------

	D8N class: the built-in connector for regular grids. D8N regular grids are the most used type of grids, like most global DEMs for example.
	They are defined by a number of rows (ny) and a number of columns (nx) with associated spacing between nodes in both directions (dy,dx).
	The D8N constructor defaults to ``4edges`` boundary preset and precompiles neighbouring. Neighbouring and indexing is optimised by geometrical relationships:
	Node index is simply a flat index in the row major direction (idx = row_id * nx + col_id) and the link ID for the top-right, right, bottomright and bottom neighbours are respectiveley idx *4, idx*4 +1, idx * 4+2 and idx*4 + 3

	The constructor only needs basics geometrical information.

	Parameters:
	-----------

		* nx (int): number of nodes in the X direction (columns)
		* ny (int): number of nodes in the Y direction (rows)
		* dx (float64): distance between nodes in the x directions
		* dy (float64): distance between nodes in the y directions
		* x_min (float64): X coordinates of the bottom left corner
		* y_min (float64): Y coordinates of the top left corner


	Authors:
	--------
	B.G.

	)pbdoc")
		.def("set_default_boundaries",
				 &D8connector<FLOATING_POINT_DAGGER>::set_default_boundaries,
				 py::arg("boundary_preset"),
				 R"pbdoc(

	Automatically sets the default boundary system based on predefined presets.

	Description:
	------------

	Boundary conditions can be tricky to set. Luckily, most use cases needs
	classical opened bounary at the edge of the dem, or periodic (sometimes called
	cyclic) bounaries at EW or NS edges. This functions automate these (i.e. set the
	right boundary codes and recompute the linknode array)

	Parameters:
	-----------

		* boundary_preset (str): the preset. Can be "4edges", "periodic_EW" or
			"periodic_NS"


	Authors:
	--------
	B.G.

	)pbdoc")

		.def("set_custom_boundaries",
				 &D8connector<FLOATING_POINT_DAGGER>::set_custom_boundaries<
					 py::array_t<int, 1>>,
				 py::arg("boundary_codes"),
				 R"pbdoc(

	Manually sets the boundary codes as uint8 1D array (see doc for options).

	Description:
	------------

	Function to manually specify boundary conditions, for cases where the classic
	presets are not providing the desired option. For example if one needs specific
	outlets/inlets location or deactivate some nodes. This is a good function to
	automate the ingestion of nodata, isolating a watershed, or defining entry
	points to the landscape.

	Parameters:
	-----------

		* boundary_codes (1D Array): the uint8 array of node size containing the
			boundary codes (see section :ref:`boundary`)


	Authors:
	--------
	B.G.

	)pbdoc"

				 )
		.def("print_dim",
				 &D8connector<FLOATING_POINT_DAGGER>::print_dim,
				 R"pdoc(Debugging function)pdoc")
		.def(
			"get_HS",
			&D8connector<FLOATING_POINT_DAGGER>::
				get_HS<std::vector<FLOATING_POINT_DAGGER>, py::array>,
			R"pdoc(Deprecated, kept for legacy (see hillshade function outside of ``connector``))pdoc")
		.def("get_mask_array",
				 &D8connector<FLOATING_POINT_DAGGER>::get_mask_array,
				 R"pdoc(Returns a 1D array of bool where false are nodata)pdoc")

		.def("set_values_at_boundaries",
				 &D8connector<FLOATING_POINT_DAGGER>::set_values_at_boundaries<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("array_to_modify"),
				 py::arg("value"),
				 R"pdoc(
	Modify an array in place with given value where connector can out flux

	Description:
	------------

	Every nodes on the input array where the connector satisfy the coundary condition
	``can_out`` will be set to a given value. It modify the array in place (i.e.
	does not return anything but modify the input).

	Useful to impose boundary condition in LEMs in a versatile way.

	Parameters:
	-----------

		* array_to_modify (1D Array): flot64 array to be modified in place
		* value (float64): value to impose


	Authors:
	--------
	B.G.

	)pdoc")

		.def("set_out_boundaries_to_permissive",
				 &D8connector<FLOATING_POINT_DAGGER>::set_out_boundaries_to_permissive,
				 R"pdoc(
	Converts OUT boundaries to CAN_OUT.

	Description:
	------------

	Converts OUT boundaries, where flow entering this cell has to out the model to
	CAN_OUT, where flow leaves the model if the cell has no downslope receivers.

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_boundary_at_node",
				 &D8connector<FLOATING_POINT_DAGGER>::get_boundary_at_node,
				 py::arg("node_idx"),
				 R"pdoc(Returns the boundary code at node index)pdoc")

		.def(
			"get_rowcol_Sreceivers",
			&D8connector<FLOATING_POINT_DAGGER>::get_rowcol_Sreceivers,
			py::arg("row_index"),
			py::arg("col_index"),
			R"pdoc(Debug function to get the receiver (node) indices of a node from its row and column index)pdoc")

		.def(
			"print_receivers",
			&D8connector<FLOATING_POINT_DAGGER>::template print_receivers<
				std::vector<FLOATING_POINT_DAGGER>>,
			py::arg("node_index"),
			py::arg("topography"),
			R"pdoc(Debuggin function printing to the terminal the receivers of a node index and their topography (post graph computation! so the topographic field may not be the one used for the receivers/LM computations))pdoc")

		.def("get_rec_array_size",
				 &D8connector<FLOATING_POINT_DAGGER>::get_rec_array_size,
				 R"pdoc(Debug function - ignore)pdoc")

		.def("update_links_MFD_only",
				 &D8connector<FLOATING_POINT_DAGGER>::template update_links_MFD_only<
					 std::vector<FLOATING_POINT_DAGGER>>,
				 py::arg("topography"),
				 R"pdoc(
	Updates all the link directionalities - but not the SFD receiver/donors.

	Description:
	------------

	Updates all the link directionalities based on a given topography. Note that it
	does not process the SFD receiver/donors and is (paradoxally) faster. However,
	the use of this function is reserved to experienced user who seek tuning perform
	ances as many routines are speed up by SFD info behind the scene -  even in MFD.

	Parameters
	----------

	- topography (1D array): the flat topography


	Authors:
	--------
	B.G.

	)pdoc")

		.def("update_links_from_topo",
				 &D8connector<FLOATING_POINT_DAGGER>::template update_links_from_topo<
					 std::vector<FLOATING_POINT_DAGGER>>,
				 py::arg("topography"),
				 R"pdoc(
	Updates all the link directionalities - but not the SFD receiver/donors.

	Description:
	------------

	Updates all the link directionalities based on a given topography. This is the
	full updating function computing MFD/SFD links/receivers/donors/... .

	Parameters
	----------

	- topography (1D array): the flat topography


	Authors:
	--------
	B.G.

	)pdoc")

		.def("sum_at_outlets",
				 &D8connector<FLOATING_POINT_DAGGER>::template sum_at_outlets<
					 py::array_t<FLOATING_POINT_DAGGER, 1>,
					 FLOATING_POINT_DAGGER>,
				 py::arg("array"),
				 py::arg("include_pits"),
				 R"pdoc(
	Sum the values contains in the input array where flux out the model.

	Description:
	------------

	Sum the values contains in the input array where flux out the model. It can also
	include the internal pit if needed. Internal pits are only referring to nodes
	where local minimas have not been resolved

	DEPRECATION WARNING: Will be detached to a standalone algorithm in a future
	update.

	Parameters
	----------

	- array (1D array): the array to sum (of node size)
	- include_pits (bool): if true, internal pits (receiverless nodes) are sumed too

	Returns
	----------

	The sum.

	Authors:
	--------
	B.G.

	)pdoc")

		.def("keep_only_at_outlets",
				 &D8connector<FLOATING_POINT_DAGGER>::template keep_only_at_outlets<
					 py::array_t<FLOATING_POINT_DAGGER, 1>,
					 py::array>,
				 py::arg("array"),
				 py::arg("include_pits"),
				 R"pdoc(
	return a copy of the input array where all the non-outting nodes are set to 0.

	Description:
	------------

	returns a copy of the input array where all the outting nodes are set to 0.
	It can also include the internal pit if needed. Internal pits are only referring
	to nodes where local minimas have not been resolved.

	DEPRECATION WARNING: Will be detached to a standalone algorithm in a future
	update.

	Parameters
	----------

	- array (1D array): the array to sum (of node size)
	- include_pits (bool): if true, internal pits (receiverless nodes) are sumed too

	Returns
	----------

	The modified array

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_SFD_receivers",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_SFD_receivers<
					 py::array_t<int, 1>>,
				 R"pdoc(
	returns the array of SFD receivers.

	Description:
	------------

	returns the array of SFD receivers according to the current state of the \
	connector. These can be obtained after computing links from a topography, but
	can also be modified when a graph resolves local minimas.

	Returns
	----------

	array of node size containing the node index of the SFD receivers

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_SFD_dx",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_SFD_dx<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(
	returns the array of SFD distance to receivers.

	Description:
	------------

	returns the array of SFD dist. to receivers according to the current state of
	the connector. These can be obtained after computing links from a topography,
	but can also be modified when a graph resolves local minimas.

	Returns
	----------

	array of floating point distance to the SFD receiver

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_SFD_ndonors",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_SFD_ndonors<
					 py::array_t<int, 1>>,
				 R"pdoc(
	returns the array of SFD number of donors.

	Description:
	------------

	returns the array of number of SFD donors for every nodes.

	Returns
	----------

	array of integer of node size with the number of SFD donors for each nodes

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_SFD_donors_flat",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_SFD_donors_flat<
					 py::array_t<int, 1>>,
				 R"pdoc(
	returns a flat array of SFD donors(read description for indexing!).

	Description:
	------------

	returns a flat array of SFD donors. It is a sparse array for the sake of simplicity:
	the Donors for a node are every >=0 nodes from the index node_index * 8 to
	node_index * 8 + nSdonors[node_index].

	Returns
	----------

	The flat array of donor indices

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_SFD_donors_list",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_SFD_donors_list<
					 std::vector<std::vector<int>>>,
				 R"pdoc(
	returns a list (not an array!) of irregular size with donor indices.

	Description:
	------------

	returns a list (not an array!) of irregular size with donor indices. The first
	index is the node index and it points to a list of ndonors size with all the
	SFD donors of the given node.

	Returns
	----------

	THe 2D irregular list of node size

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_links",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_links<
					 std::vector<std::uint8_t>>,
				 R"pdoc(
	returns an array of link size with link type.

	Description:
	------------

	returns an array of link size with link type. See the documentation for the
	relationships between nodes and links. the node indices for a link index can be
	found in the linknode array, at link_index*2 and link_index*2 + 1. The link type
	is an uint8 number: 0 for inverse (node 2 give to node 1), 1 for normal (node 1
	fives to node 2) or 3 invalid/inactive link.

	Returns
	----------

	the array of link size with link code (uint8)

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_linknodes_flat",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_linknodes_flat<
					 py::array_t<int, 1>>,
				 R"pdoc(
	returns a flat array of linknodes (node indices pair for each links).

	Description:
	------------

	1D flat array of linknodes. The node corresponding to link index li can be found
	at linknodes[li*2] and linknodes[li*2 + 1]. Non-existing or invalid links have
	node indices of -1.

	Returns
	----------

	the array of 2 * link size with node indicies for each links

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_linknodes_list",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_linknodes_list<
					 std::vector<std::vector<int>>>,
				 R"pdoc(
	returns a list (not an array!) of link nodes.

	Description:
	------------

	2D link of link nodes. index is link index and it provides a list of nodes
	composing the link.

	Returns
	----------

	List of link size.

	Authors:
	--------
	B.G.

	)pdoc")

		.def(
			"get_linknodes_list_oriented",
			&D8connector<FLOATING_POINT_DAGGER>::template get_linknodes_list_oriented<
				std::vector<std::vector<int>>>,
			R"pdoc(
	returns a list (not an array!) of link nodes, donor first, rec second.

	Description:
	------------

	2D link of link nodes. index is link index and it provides a list of nodes
	composing the link. This version orients the list with the donor always first
	(and the receiver always second).

	Returns
	----------

	List of link size.

	Authors:
	--------
	B.G.

	)pdoc")

		.def(
			"get_SFD_receivers_at_node",
			&D8connector<FLOATING_POINT_DAGGER>::get_SFD_receivers_at_node,
			R"pdoc(Returns the node index of the SFD receivers for a given node index)pdoc")

		.def(
			"get_SFD_dx_at_node",
			&D8connector<FLOATING_POINT_DAGGER>::get_SFD_dx_at_node,
			R"pdoc(Returns the distance to the SFD receivers for a given node index)pdoc")

		.def("get_SFD_ndonors_at_node",
				 &D8connector<FLOATING_POINT_DAGGER>::get_SFD_ndonors_at_node

				 )

		.def("get_SFD_donors_at_node",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_SFD_donors_at_node<
					 std::vector<int>>,
				 R"pdoc(Returns a list of SFD donors for a given node index)pdoc")

		.def("get_SFD_gradient",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_SFD_gradient<
					 py::array_t<FLOATING_POINT_DAGGER, 1>,
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("topography"),
				 R"pdoc(
	returns an array of node size with the topographic gradient

	Description:
	------------

	Takes a topography and returns the steepest descent gradient using precomputed
	SFD receivers informations. Note that if you feed the function with a different
	topography than the one the graph has been computed with - or if local minimas
	have been resolved, you may obtain local negative values.

	Parameters
	-----------

	- Topography (array 1D): the topography to use for the gradient computation

	Returns
	----------

	Array of node size fo topographic gradient

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_links_gradient",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_links_gradient<
					 py::array_t<FLOATING_POINT_DAGGER, 1>,
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("topography"),
				 py::arg("minimum_slope"),
				 R"pdoc(
	returns an array of link size with the topographic gradient for each of them.

	Description:
	------------

	Takes a topography and returns the gradient for each topographic link. It can
	recast teh negative/small gradients to a minimum value if needed.

	Parameters
	-----------

	- Topography (array 1D): the topography to use for the gradient computation
	- min_slope (float): the minimum slope - set to very negative value if not
	needed.

	Returns
	----------

	Array of link size of topographic gradient

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_MFD_mean_gradient",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_MFD_mean_gradient<
					 py::array_t<FLOATING_POINT_DAGGER, 1>,
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("topography"),
				 R"pdoc(
	returns an array of node size with the mean topographic gradient for each nodes.

	Description:
	------------

	returns an array of node size with the mean topographic gradient for each nodes.
	The gradient is computed for all receivers and averaged.

	Parameters
	-----------

	- Topography (array 1D): the topography to use for the gradient computation


	Returns
	----------

	Array of node size of mean topographic gradient

	Authors:
	--------
	B.G.

	)pdoc"

				 )

		.def(
			"get_MFD_weighted_gradient",
			&D8connector<FLOATING_POINT_DAGGER>::template get_MFD_weighted_gradient<
				py::array_t<FLOATING_POINT_DAGGER, 1>,
				py::array_t<FLOATING_POINT_DAGGER, 1>>,
			py::arg("topography"),
			py::arg("weights"),
			R"pdoc(
	returns an array of node size with the weighted mean gradient for each nodes.

	Description:
	------------

	returns an array of node size with the weighted mean gradient for each nodes.
	It requires a weight for each links which can be obtained using the
	``get_link_weights`` function.

	Parameters
	-----------

	- Topography (array 1D): the topography to use for the gradient computation
	- weight (array 1D): link size array obtained with ``get_link_weights`` function


	Returns
	----------

	Array of node size of weighted topographic gradient

	Authors:
	--------
	B.G.

	)pdoc")

		.def("get_link_weights",
				 &D8connector<FLOATING_POINT_DAGGER>::template get_link_weights<
					 py::array_t<FLOATING_POINT_DAGGER, 1>,
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("gradients"),
				 py::arg("exponent"),
				 R"pdoc(
	Computes partition weights for each link function of rec slopes per node basis.

	Description:
	------------

	Assign a [0,1] weight to each link to partition flow in MFD from every nodes to
	their receivers. The partition uses an exponent to set the sensitivity to slope
	differences. in a general manner, the higher the value, the more weight would
	go toward the steepest receiver. There are three particular cases (numerically
	and conceptually):

	- exp = 0: equally parted towards all receivers regardless of their slopes
	- exp = 0.5: proportional to the squareroot of the slopes
	- exp = 1: perfectly proportional to the slopes

	Parameters
	-----------

	- link_gradients (1D array): gradients per links (``get_link_gradients``)
	- exponent (float): a positive exponent setting the sensitivity to slope

	Returns
	----------

	Array of link size of partition weights

	Authors:
	--------
	B.G.)pdoc")

		.def(
			"set_stochaticiy_for_SFD",
			&D8connector<FLOATING_POINT_DAGGER>::set_stochaticiy_for_SFD,
			py::arg("magnitude"),
			R"pdoc(EXPERIMENTAL: adds stochasticity to the SFD receivers calculation. Best to ignore.)pdoc");
}

// ARCHIVE
// py::class_<D4connector<FLOATING_POINT_DAGGER>>(
// 		m,
// 		"D4N",
// 		R"pdoc(DEPRECATED - will be back at some points, keeping for legacy)pdoc")
// 		.def(py::init<int,
// 									int,
// 									FLOATING_POINT_DAGGER,
// 									FLOATING_POINT_DAGGER,
// 									FLOATING_POINT_DAGGER,
// 									FLOATING_POINT_DAGGER>())
// 		.def("set_default_boundaries",
// 				 &D4connector<FLOATING_POINT_DAGGER>::set_default_boundaries)
// 		.def("set_custom_boundaries",
// 				 &D4connector<FLOATING_POINT_DAGGER>::set_custom_boundaries<
// 					 py::array_t<int, 1>>)
// 		.def("print_dim", &D4connector<FLOATING_POINT_DAGGER>::print_dim)
// 		.def("get_HS",
// 				 &D4connector<FLOATING_POINT_DAGGER>::
// 					 get_HS<std::vector<FLOATING_POINT_DAGGER>, py::array>)
// 		// .def("fill_barne_2014",
// 		//
// &D4connector<FLOATING_POINT_DAGGER>::fill_barne_2014<std::vector<FLOATING_POINT_DAGGER>
// 		// >)
// 		.def("get_mask_array",
// &D4connector<FLOATING_POINT_DAGGER>::get_mask_array)
// 		.def("set_values_at_boundaries",
// 				 &D4connector<FLOATING_POINT_DAGGER>::set_values_at_boundaries<
// 					 py::array_t<FLOATING_POINT_DAGGER, 1>>)
// 		.def("set_out_boundaries_to_permissive",
// 				 &D4connector<FLOATING_POINT_DAGGER>::set_out_boundaries_to_permissive)
// 		.def("get_boundary_at_node",
// 				 &D4connector<FLOATING_POINT_DAGGER>::get_boundary_at_node);

// 	py::class_<numvec<FLOATING_POINT_DAGGER>>(m, "numvecf64")
// 		.def(py::init<py::array_t<FLOATING_POINT_DAGGER, 1>&>())
// 		.def("get", &numvec<FLOATING_POINT_DAGGER>::get)
// 		.def("set", &numvec<FLOATING_POINT_DAGGER>::set);
