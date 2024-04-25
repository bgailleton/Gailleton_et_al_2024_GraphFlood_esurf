#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_param(py::module& m)
{
	py::class_<ParamBag<int, double>>(m, "ParamBag")

		.def(py::init<>())
		.def("get_morphomode", &ParamBag<int, double>::get_morphomode)
		.def("set_morphomode", &ParamBag<int, double>::set_morphomode)
		.def("get_hydromode", &ParamBag<int, double>::get_hydromode)
		.def("set_hydromode", &ParamBag<int, double>::set_hydromode)
		.def("set_ke", &ParamBag<int, double>::set_ke)
		.def("set_variable_ke",
				 &ParamBag<int, double>::set_variable_ke<py::array_t<double, 1>>)
		.def("enable_gf2_diffuse_Qwin",
				 &ParamBag<int, double>::enable_gf2_diffuse_Qwin)
		.def("disable_gf2_diffuse_Qwin",
				 &ParamBag<int, double>::disable_gf2_diffuse_Qwin)
		.def("enable_gf2_morpho", &ParamBag<int, double>::enable_gf2_morpho)
		.def("disable_gf2_morpho", &ParamBag<int, double>::disable_gf2_morpho)
		.def("set_kd", &ParamBag<int, double>::set_kd)
		.def("get_kd", &ParamBag<int, double>::get_kd)
		.def("set_kdl", &ParamBag<int, double>::set_kdl)
		.def("get_kdl", &ParamBag<int, double>::get_kdl)
		.def("set_kel", &ParamBag<int, double>::set_kel)
		.def("get_kel", &ParamBag<int, double>::get_kel)
		.def("calculate_ke_tau_c_from_MPM",
				 &ParamBag<int, double>::calculate_ke_tau_c_from_MPM)

		.def("enable_TSG_dist", &ParamBag<int, double>::enable_TSG_dist)

		.def("disable_TSG_dist", &ParamBag<int, double>::disable_TSG_dist)
		.def("set_TSG_distmax", &ParamBag<int, double>::set_TSG_distmax)
		.def("set_time_dilatation_morpho",
				 &ParamBag<int, double>::set_time_dilatation_morpho)
		.def("get_time_dilatation_morpho",
				 &ParamBag<int, double>::get_time_dilatation_morpho)
		.def("enable_transient_flow", &ParamBag<int, double>::enable_transient_flow)
		.def("disable_transient_flow",
				 &ParamBag<int, double>::disable_transient_flow)
		.def("enable_bank_erosion", &ParamBag<int, double>::enable_bank_erosion)
		.def("disable_bank_erosion", &ParamBag<int, double>::disable_bank_erosion)
		.def("set_gf2Bbval", &ParamBag<int, double>::set_gf2Bbval)
		.def("get_gf2Bbval", &ParamBag<int, double>::get_gf2Bbval)
		.def("set_minWeightQw", &ParamBag<int, double>::set_minWeightQw)
		.def("get_minWeightQw", &ParamBag<int, double>::get_minWeightQw)
		.def("set_capacityFacQs", &ParamBag<int, double>::set_capacityFacQs)
		.def("get_capacityFacQs", &ParamBag<int, double>::get_capacityFacQs)
		.def("set_capacityFacQw", &ParamBag<int, double>::set_capacityFacQw)
		.def("get_capacityFacQw", &ParamBag<int, double>::get_capacityFacQw)

		;
}
