#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_algos(py::module& m)
{
	m.def("hillshade",
				&hillshade<D8connector<FLOATING_POINT_DAGGER>,
									 py::array_t<FLOATING_POINT_DAGGER, 1>,
									 py::array_t<FLOATING_POINT_DAGGER, 1>,
									 FLOATING_POINT_DAGGER>,
				py::arg("connector"),
				py::arg("topography"),
				R"pbdoc(
Hillshading function for visualisation

Description:
------------

Returns a [0,1] hillshade for a regular grid. Negative/nodata are set to 0

Parameters:
-----------

	* D8connector
	* Flat topography of node size (1D array)

Returns:
--------

	* Flat hillshade of node size (1D array)

Authors:
--------
B.G.

)pbdoc");

	m.def("rayshade",
				&rayshade<DAGGER::graph<FLOATING_POINT_DAGGER,
																DAGGER::D8connector<FLOATING_POINT_DAGGER>>,
									D8connector<FLOATING_POINT_DAGGER>,
									py::array_t<FLOATING_POINT_DAGGER, 1>,
									py::array_t<FLOATING_POINT_DAGGER, 1>,
									FLOATING_POINT_DAGGER>);

	m.def("set_BC_to_remove_seas",
				&set_BC_to_remove_seas<D8connector<FLOATING_POINT_DAGGER>,
															 py::array_t<FLOATING_POINT_DAGGER, 1>,
															 FLOATING_POINT_DAGGER>);

	m.def("label_depressions_PQ",
				&label_depressions_PQ<py::array_t<FLOATING_POINT_DAGGER, 1>,
															py::array_t<int, 1>,
															D8connector<FLOATING_POINT_DAGGER>>);

	m.def("label_ocean",
				&label_ocean<py::array_t<FLOATING_POINT_DAGGER, 1>,
										 py::array_t<int, 1>,
										 D8connector<FLOATING_POINT_DAGGER>>);

	m.def("standalone_priority_flood",
				&standalone_priority_flood<D8connector<FLOATING_POINT_DAGGER>,
																	 py::array_t<FLOATING_POINT_DAGGER, 1>,
																	 py::array_t<FLOATING_POINT_DAGGER, 1>,
																	 FLOATING_POINT_DAGGER>,
				py::arg("topography"),
				py::arg("connector"));

	m.def("standalone_priority_flood_opti",
				&standalone_priority_flood_opti<
					D8connector<FLOATING_POINT_DAGGER>,
					DAGGER::graph<FLOATING_POINT_DAGGER,
												DAGGER::D8connector<FLOATING_POINT_DAGGER>>,
					py::array_t<FLOATING_POINT_DAGGER, 1>,
					py::array_t<FLOATING_POINT_DAGGER, 1>,
					FLOATING_POINT_DAGGER>,
				py::arg("topography"),
				py::arg("connector"),
				py::arg("graph"));

	m.def(
		"RiverNetwork",
		&RiverNetwork<FLOATING_POINT_DAGGER,
									DAGGER::D8connector<FLOATING_POINT_DAGGER>,
									DAGGER::graph<FLOATING_POINT_DAGGER,
																DAGGER::D8connector<FLOATING_POINT_DAGGER>>>);

	m.def("RiverNetworkC8",
				&RiverNetworkC8<FLOATING_POINT_DAGGER,
												DAGGER::Connector8<int, FLOATING_POINT_DAGGER>>);

	m.def("DrainageDivides",
				&DrainageDivides<
					FLOATING_POINT_DAGGER,
					DAGGER::D8connector<FLOATING_POINT_DAGGER>,
					DAGGER::graph<FLOATING_POINT_DAGGER,
												DAGGER::D8connector<FLOATING_POINT_DAGGER>>>);
}
