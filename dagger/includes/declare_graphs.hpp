/*
Fucntion declaring the different graphs

*/

#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

// graph v1
template<typename CONNECTOR_T>
void
declare_graph(py::module& m, std::string typestr)
{

	m.def("compute_SFD_basin_labels",
				&compute_SFD_basin_labels<int, double, Connector8<int, double>>);

	m.def("fill_lake_in_hw",
				&fill_lake_in_hw<int, double, Connector8<int, double>>);

	m.def("compute_SFD_DA",
				&compute_SFD_DA<int, double, Connector8<int, double>>);

	m.def(
		"compute_min_distance_from_outlets",
		&compute_min_distance_from_outlets<int, double, Connector8<int, double>>);

	m.def(
		"compute_max_distance_from_outlets",
		&compute_max_distance_from_outlets<int, double, Connector8<int, double>>);

	m.def("compute_relief",
				&compute_relief<int, double, Connector8<int, double>>);

	m.def("replace_with_slopped_surface",
				&replace_with_slopped_surface<int, double, Connector8<int, double>>);

	m.def("recast_BC_bellow_Z",
				&recast_BC_bellow_Z<int, double, Connector8<int, double>>);

	m.def("recast_BC_from_outlet",
				&recast_BC_from_outlet<int, double, Connector8<int, double>>);

	m.def("compute_elevation_above_nearest_mask",
				&compute_elevation_above_nearest_mask<int,
																							double,
																							Connector8<int, double>,
																							py::array_t<std::uint8_t, 1>,
																							py::array_t<double, 1>>);

	m.def("compute_sources_D8",
				&compute_sources_D8<int,
														double,
														Connector8<int, double>,
														py::array_t<int, 1>>);

	py::class_<graph<FLOATING_POINT_DAGGER, CONNECTOR_T>>(m,
																												typestr.c_str(),
																												R"pdoc(
Full Graph module, to plug on a connector to unlock non-local topological operations.

Description:
------------

The graph manages all non local connectivity and ensures the graph is acyclic.
It contains all the routines for topological ordering, local minima resolving
(multiple methods), or accumulating values in the upstream/downstream directions


Authors:
--------
B.G.)pdoc")

		.def(py::init<CONNECTOR_T&>())
		.def(
			"init_graph",
			&graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::init_graph,
			R"pdoc(Initialise the data structure and allocate memory (mostly used internally).)pdoc")
		.def(
			"set_opt_stst_rerouting",
			&graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::set_opt_stst_rerouting,
			py::arg("onoff"),
			R"pdoc(Activate (true) or deactivate (false) an optimiser. Most of the time does not make a difference but can _eventually_ approximate a few link a bit more precisely when rerouting local minimas.)pdoc")
		.def("compute_graph",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::template compute_graph<
					 py::array_t<FLOATING_POINT_DAGGER, 1>,
					 py::array>,
				 py::arg("topography"),
				 py::arg("no_MFD"),
				 py::arg("quicksort_on"),
				 R"pdoc(
Full computation of the graph (connector updates of links included).

Description
-------------

COmpute the graph by (i) updating receivers/donors, (ii) resolving local minimas
with the method and (iii) computing topological sortings. Can be set to SFD only
to save time.

Parameters
-----------

- topography (1D array): topographic field (of node size)
- no_MFD: if true, no MFD info are computed (especially the toposort MFD and the
recomputing after local minima solver that can be time consuming
if done repeteadly - e.g. for LEMs)
- quicksort: if true, the toposort MFD is done by sorting "filled" topography by
elevation using the std::sort algorithm (aka quicksort), otherwise it uses an
homemade  topological sorting algorithm. Quicksort is O(nlogn) and toposort
O(n+l) so the most efficient is case dependent.


returns:
--------

PP_topography (1D array): a preprocessed topography where LM have been filled or
carved or processed in the wanted way.

Authors:
--------
B.G.

)pdoc"

				 )
		.def("gen_rdid",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::gen_rdid,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("get_rdid",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::get_rdid,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("is_Sstack_full",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::is_Sstack_full,
				 R"pdoc(Debugging function, to ignore)pdoc")
		.def("activate_opti_sparse_border_cordonnier",
				 &graph<FLOATING_POINT_DAGGER,
								CONNECTOR_T>::activate_opti_sparse_border_cordonnier,
				 R"pdoc(Debugging function, to ignore)pdoc")
		.def("get_all_nodes_upstream_of",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_all_nodes_upstream_of<py::array_t<int, 1>>,
				 py::arg("node"),
				 py::arg("only_SFD"),
				 R"pdoc(
Fecth all the nodes upstream of a given one.

Description
-------------

Use the SFD stack structure from Braun and Willett (2013) to gather all the
nodes upstream of a given one. Effectively labels a watershed from a custom
starting point. Can be extended to MFD using FIFO queues (slower). This is not
the fastest way to label all the watersheds and should be reserved for one-off
operations.

Parameters
-----------

- node (int): the starting node index
- only_SFD (bool): only label Steepest Descent watersheds (faster, less nodes)

returns:
--------

nodes_upstream (1D array: int): returns a 1D array containing all the node
indices of upstream locations.

Authors:
--------
B.G.

)pdoc")
		.def("get_all_nodes_downstream_of",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_all_nodes_downstream_of<py::array_t<int, 1>>,
				 py::arg("node"),
				 py::arg("only_SFD"),
				 R"pdoc(
Fecth all the nodes downstream of a given one.

Description
-------------

Use the SFD stack structure from Braun and Willett (2013) to gather all the
nodes downstream of a given one. Effectively labels a water paths from a custom
starting point. Can be extended to MFD using FIFO queues (slower). This is not
the fastest way to label all the flow paths and should be reserved for one-off
operations.

Parameters
-----------

- node (int): the starting node index
- only_SFD (bool): only label Steepest Descent watersheds (faster, less nodes)

returns:
--------

nodes_downstream (1D array: int): returns a 1D array containing all the node
indices of downstream locations.

Authors:
--------
B.G.

)pdoc")
		.def("get_SFD_stack",
				 &graph<FLOATING_POINT_DAGGER,
								CONNECTOR_T>::template get_SFD_stack<py::array_t<size_t, 1>>,
				 R"pdoc(
Returns the single flow stack (Braun and Willett. 2013) in "stack order".

Description
-------------

Returns the single flow direction stack (sensu Braun and Willett (2013).
Warning: it corresponds to the last computation (``compute_graph``) and will
return an empty array if not computed at all yet.

The stack is in "stack order", **i.e.** from the most dowstream nodes to the
most upstream ones.

This can also be call the topologicaly sorted list of nodes. be careful though,
the array is of node size but contains ordered node index.

returns:
--------

SFD stack array.

Authors:
--------
B.G.

)pdoc")
		.def("get_MFD_stack",
				 &graph<FLOATING_POINT_DAGGER,
								CONNECTOR_T>::template get_MFD_stack<py::array_t<size_t, 1>>,
				 R"pdoc(
Returns the Multiple flow stack in "stack order".

Description
-------------

Returns the multiple flow direction stack.

Warning: it corresponds to the last computation (``compute_graph`` with
only_SFD = false) and will return an empty array if not computed at all yet.

The stack is in "stack order", **i.e.** from the most dowstream nodes to the
most upstream ones.

This can also be call the topologicaly sorted list of nodes. be careful though,
the array is of node size but contains ordered node index. It does not follow
the same topology than Braun and Willett (2013) and is more difficult to compute
. Only use if MFD is truly needed.

returns:
--------

MFD stack array.

Authors:
--------
B.G.

)pdoc")

		.def("accumulate_constant_downstream_SFD",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template accumulate_constant_downstream_SFD<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("constant_value"),
				 R"pdoc(
Accumulates (integrate) a constant downstream in the SFD direction.

Description
-------------

Sums a constant downstream in the single flow direction (each node adds the
constant value and give it to its receivers following the inverse stack order).
For example, can be the area of a regular grid cell to compute drainage area.
The operation is very efficient (needs the computed graph).

Parameters
-----------

- constant_value (float): the value to accumulate

returns:
--------

The 1D accumulated array (float, node size)

Authors:
--------
B.G.

)pdoc")
		.def("accumulate_variable_downstream_SFD",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template accumulate_variable_downstream_SFD<
						 py::array_t<FLOATING_POINT_DAGGER, 1>,
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("values"),
				 R"pdoc(
Accumulates (integrate) a variable downstream in the SFD direction.

Description
-------------

Sums a variable downstream in the single flow direction (each node adds the
variable value and give it to its receivers following the inverse stack order).
For example, can be the cell area times variable precipitation rates of a
regular grid cell to compute weighted drainage area (~runoff).

The operation is very efficient (needs the computed graph).

Parameters
-----------

- values (1D array of node size, float): the value to accumulate

returns:
--------

The 1D accumulated array (float, node size)

Authors:
--------
B.G.

)pdoc")
		.def("accumulate_constant_downstream_MFD",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template accumulate_constant_downstream_MFD<
						 py::array_t<FLOATING_POINT_DAGGER, 1>,
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("weights"),
				 py::arg("constant_value"),
				 R"pdoc(
Accumulates (integrate) a constant downstream in the MFD direction.

Description
-------------

Sums a constant downstream in the multiple flow direction (each node adds the
variable value and give it to its receivers following the inverse stack order).

It needs to ingest partitionning weights, an array of n links size telling each
nodes how to partition their value to their multiple receivers. The latter can
be calculated from a ``Connector`` and is most of the time expressed function of
topographic gradient.

For example, can be the area of a regular grid cell to compute drainage area.

Parameters
-----------

- weigths (1D array of link size, float): the weights
- constant_value (float): the value to accumulate

returns:
--------

The 1D accumulated array (float, node size)

Authors:
--------
B.G.

)pdoc")
		.def("accumulate_variable_downstream_MFD",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template accumulate_variable_downstream_MFD<
						 py::array_t<FLOATING_POINT_DAGGER, 1>,
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 py::arg("weights"),
				 py::arg("values"),
				 R"pdoc(
Accumulates (integrate) a variable downstream in the MFD direction.

Description
-------------

Sums a variable value of node size downstream in the multiple flow direction
(each node adds the variable value and give it to its receivers following the
inverse stack order).

It needs to ingest partitionning weights, an array of n links size telling each
nodes how to partition their value to their multiple receivers. The latter can
be calculated from a ``Connector`` and is most of the time expressed function of
topographic gradient.

For example, can be the area of a regular grid cell times a variable
precipitation rate to compute runoff.

Parameters
-----------

- weigths (1D array of link size, float): the weights
- values (1D array of link size, float): the values to accumulate

returns:
--------

The 1D accumulated array (float, node size)

Authors:
--------
B.G.

)pdoc")
		.def("set_LMR_method",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::set_LMR_method,
				 py::arg("LMR_method"),
				 R"pdoc(
Sets the Local Minima Resolver method.

Description
-------------

Select which LMR (Local Minima Resolver) to use. It needs to be a LMR enum value
and will dictacte if/how the local minima (internal pits) are rerouted (or not).

Can be one of the following:

- dagger.LMR.cordonnier_fill: approximate filling from Cordonnier et al. (2019)
- dagger.LMR.cordonnier_carve: approximate carving from Cordonnier et al. (2019)
- dagger.LMR.priority_flood: Fill with Barnes et al., 2014
- dagger.LMR.none: Ignore pits

There is no better solution than another, performances as well as accuracy are
extremely case dependent.


Parameters
-----------

- LMR_method (LMR enum)


Authors:
--------
B.G.

)pdoc")
		.def(
			"set_minimum_slope_for_LMR",
			&graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::set_minimum_slope_for_LMR,
			py::arg("slope"),
			R"pdoc(LMR solvers impose a numerical topographic gradient to avoid 0 slopes. Default is 1e-5.)pdoc")
		.def(
			"set_slope_randomness_for_LMR",
			&graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::set_slope_randomness_for_LMR,
			py::arg("magnitude"),
			R"pdoc(Avoid falt surfaces by imposing a very small randomness when processing local minimas. Must be an order of magitude smaller than the  minimal slope.)pdoc")

		// Distance functions
		.def("get_SFD_distance_from_outlets",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_SFD_distance_from_outlets<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(
Calculates the distance from the outlets following the SFD.

Description
-------------

Uses the SFD stack order to integrate the distance from the outlets to the
sources. Useful for long profiles.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc"

				 )
		.def("get_SFD_min_distance_from_sources",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_SFD_min_distance_from_sources<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(
Calculates the minimum distance from the sources following the SFD.

Description
-------------

Uses the SFD stack in inverse order to integrate the distance from the sources
toward the outlets. It keeps the minimum distance. Useful to identify the
closest sources.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc")
		.def("get_SFD_max_distance_from_sources",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_SFD_max_distance_from_sources<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(
Calculates the maximum distance from the sources following the SFD.

Description
-------------

Uses the SFD stack in inverse order to integrate the distance from the sources
toward the outlets. It keeps the maximum distance. Useful to identify the
furthest sources.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc")
		.def("get_MFD_max_distance_from_sources",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_MFD_max_distance_from_sources<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(
Calculates the maximum distance from the sources following the MFD.

Description
-------------

Uses the MFD stack in inverse order to integrate the distance from the sources
toward the outlets. It keeps the maximum distance. Useful to identify the
furthest sources.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc")
		.def("get_MFD_min_distance_from_sources",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_MFD_min_distance_from_sources<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(
Calculates the minimum distance from the sources following the MFD.

Description
-------------

Uses the MFD stack in inverse order to integrate the distance from the sources
toward the outlets. It keeps the minimum distance. Useful to identify the
closest sources.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc")
		.def("get_MFD_max_distance_from_outlets",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_MFD_max_distance_from_outlets<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(
Calculates the maximum distance from the outlet following the MFD.

Description
-------------

Uses the MFD stack in inverse order to integrate the distance from the outlet
toward the outlets. It keeps the maximum distance. Useful to identify the
furthest outlet.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc")
		.def("get_MFD_min_distance_from_outlets",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
					 template get_MFD_min_distance_from_outlets<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(
Calculates the minimum distance from the outlet following the MFD.

Description
-------------

Uses the MFD stack in inverse order to integrate the distance from the outlet
toward the outlets. It keeps the minimum distance. Useful to identify the
closest outlet.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc")

		// Watershed labelling
		.def(
			"get_SFD_basin_labels",
			&graph<FLOATING_POINT_DAGGER,
						 CONNECTOR_T>::template get_SFD_basin_labels<py::array_t<int, 1>>,
			R"pdoc(
Labels SFD watersheds with unique ID.

Description
-------------

Uses the SFD stack order to label very efficiently basins.

returns:
--------

1D array of integer labels

Authors:
--------
B.G.

)pdoc")

		.def(
			"get_drainage_area_SFD",
			&graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
				template get_drainage_area_SFD<py::array_t<FLOATING_POINT_DAGGER, 1>>,
			R"pdoc(
Labels SFD watersheds with unique ID.

Description
-------------

Uses the SFD stack order to label very efficiently basins.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc")

		.def(
			"get_drainage_area_MFD",
			&graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::
				template get_drainage_area_MFD<py::array_t<FLOATING_POINT_DAGGER, 1>,
																			 py::array_t<FLOATING_POINT_DAGGER, 1>>,
			R"pdoc(
Labels MFD watersheds with unique ID.

Description
-------------

Uses the MFD stack order to label very efficiently basins.

returns:
--------

1D array of flow distance

Authors:
--------
B.G.

)pdoc")

		.def("get_n_pits",
				 &graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::get_n_pits,
				 R"pdoc(Return the number of internal pits (prior solving).)pdoc")

		.def(
			"get_debug_mask",
			&graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::get_debug_mask,
			R"pdoc(Ignore. Internal debugging mask process. Changes purpose and is mostly deactivated.)pdoc")

		.def(
			"get_debug_int",
			&graph<FLOATING_POINT_DAGGER, CONNECTOR_T>::get_debug_int,
			R"pdoc(Ignore. Internal debugging int process. Changes purpose and is mostly deactivated.)pdoc")

		;
}
