#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

template<typename CONNECTOR_T>
void
declare_popscape(py::module& m, std::string typestr)
{
	py::class_<popscape<FLOATING_POINT_DAGGER,
											DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											CONNECTOR_T>>(m, typestr.c_str())
		.def(py::init<DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>&,
									CONNECTOR_T&>())
		.def("StSt",
				 &popscape<FLOATING_POINT_DAGGER,
									 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									 CONNECTOR_T>::StSt)
		.def("restriction",
				 &popscape<FLOATING_POINT_DAGGER,
									 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									 CONNECTOR_T>::restriction)
		.def("interpolation",
				 &popscape<FLOATING_POINT_DAGGER,
									 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									 CONNECTOR_T>::interpolation)
		.def("smooth",
				 &popscape<FLOATING_POINT_DAGGER,
									 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									 CONNECTOR_T>::smooth)
		.def(
			"set_topo",
			&popscape<
				FLOATING_POINT_DAGGER,
				DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
				CONNECTOR_T>::template set_topo<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def(
			"get_topo",
			&popscape<
				FLOATING_POINT_DAGGER,
				DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
				CONNECTOR_T>::template get_topo<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("get_QA",
				 &popscape<
					 FLOATING_POINT_DAGGER,
					 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
					 CONNECTOR_T>::template get_QA<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("get_chistar",
				 &popscape<FLOATING_POINT_DAGGER,
									 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									 CONNECTOR_T>::
					 template get_chistar<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("simple_Kfchi",
				 &popscape<FLOATING_POINT_DAGGER,
									 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									 CONNECTOR_T>::simple_Kfchi)
		.def("simple_Kfz",
				 &popscape<FLOATING_POINT_DAGGER,
									 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									 CONNECTOR_T>::simple_Kfz);
}

template<typename CONNECTOR_T>
void
declare_popscape_old(py::module& m, std::string typestr)
{
	py::class_<popscape_old<FLOATING_POINT_DAGGER,
													DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
													CONNECTOR_T>>(m, typestr.c_str())
		.def(py::init<RANDNOISE,
									int,
									int,
									FLOATING_POINT_DAGGER,
									FLOATING_POINT_DAGGER>())
		// .def_readwrite("graph",  &popscape_old<FLOATING_POINT_DAGGER,
		// DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>, CONNECTOR_T
		// >::graph) .def_readwrite("connector",
		// &popscape_old<FLOATING_POINT_DAGGER,
		// DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>, CONNECTOR_T
		// >::connector)
		.def("solve_generic",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::solve_generic)
		.def("get_topo",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::template get_topo<py::array>)
		.def("get_QA",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::template get_QA<py::array>)
		.def("compute_graph",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::compute_graph)
		.def("compute_DA_SFD",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::compute_DA_SFD)
		.def("apply_uplift",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::apply_uplift)
		.def(
			"apply_variable_uplift",
			&popscape_old<FLOATING_POINT_DAGGER,
										DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										CONNECTOR_T>::
				template apply_variable_uplift<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("solve_SFD_SPL_imp",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::solve_SFD_SPL_imp)
		.def("hydraulic_erosion_v0",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::hydraulic_erosion_v0)
		.def("normalise_topography",
				 &popscape_old<FLOATING_POINT_DAGGER,
											 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
											 CONNECTOR_T>::normalise_topography)
		// .def("run_SFD_exp_latmag", &popscape_old<FLOATING_POINT_DAGGER,
		// DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>, CONNECTOR_T
		// >::run_SFD_exp_latmag)

		;
}
