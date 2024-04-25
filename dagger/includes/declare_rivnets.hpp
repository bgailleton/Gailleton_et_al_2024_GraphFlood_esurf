#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_rivnet(py::module& m)
{

	m.def("add_river_network_from_threshold",
				&add_river_network_from_threshold<int,double>);

	py::class_<RivNet1D<int, double>>(m, "RivNet1D")
		.def(py::init<Connector8<int,double>&, Hermes<int,double>& > ())
		.def("get_nodes",
				 &RivNet1D<int, double>::template get_nodes<py::array_t<int, 1>>)
		.def("get_nodes",
				 &RivNet1D<int, double>::template get_nodes<py::array_t<int, 1>>)
		.def("get_recs",
				 &RivNet1D<int, double>::template get_recs<py::array_t<int, 1>>)
		.def("get_source_keys",
				 &RivNet1D<int, double>::template get_source_keys<py::array_t<int, 1>>)
		.def("get_dXs",
				 &RivNet1D<int, double>::template get_dXs<py::array_t<double, 1>>)
		.def("get_flow_distance",
				 &RivNet1D<int, double>::template get_flow_distance<py::array_t<double, 1>>)
		.def("get_sources",
				 &RivNet1D<int, double>::template get_sources<py::array_t<int, 1>>)
		.def("get_outlets",
				 &RivNet1D<int, double>::template get_outlets<py::array_t<int, 1>>)
		.def("build_from_sources",
				 &RivNet1D<int, double>::template build_from_sources<py::array_t<int, 1>>)
			
		;
}
