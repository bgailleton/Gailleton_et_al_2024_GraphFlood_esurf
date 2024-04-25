#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
declare_dbag(py::module& m)
{
	py::class_<Hermes<int, double>>(m, "Hermes")
		.def(py::init<>())
		.def("set_surface",
				 &Hermes<int, double>::set_surface<py::array_t<double, 1>>)
		.def("get_surface",
				 &Hermes<int, double>::get_surface<py::array_t<double, 1>>)
		.def("set_hw", &Hermes<int, double>::set_hw<py::array_t<double, 1>>)
		.def("get_hw", &Hermes<int, double>::get_hw<py::array_t<double, 1>>)
		.def("get_ibag", &Hermes<int, double>::get_ibag<py::array_t<int, 1>>)
		.def("get_fbag", &Hermes<int, double>::get_fbag<py::array_t<double, 1>>)
		.def("set_Qwin", &Hermes<int, double>::set_Qwin<py::array_t<double, 1>>)
		.def("get_Qwin", &Hermes<int, double>::get_Qwin<py::array_t<double, 1>>)
		.def("set_DA", &Hermes<int, double>::set_DA<py::array_t<double, 1>>)
		.def("get_DA", &Hermes<int, double>::get_DA<py::array_t<double, 1>>)
		.def("set_Qwout", &Hermes<int, double>::set_Qwout<py::array_t<double, 1>>)
		.def("get_Qwout", &Hermes<int, double>::get_Qwout<py::array_t<double, 1>>)
		.def("set_Qsin", &Hermes<int, double>::set_Qsin<py::array_t<double, 1>>)
		.def("get_Qsin", &Hermes<int, double>::get_Qsin<py::array_t<double, 1>>)
		.def("set_Qsout", &Hermes<int, double>::set_Qsout<py::array_t<double, 1>>)
		.def("get_Qsout", &Hermes<int, double>::get_Qsout<py::array_t<double, 1>>)

		.def("set_theta_flow_in",
				 &Hermes<int, double>::set_theta_flow_in<py::array_t<double, 1>>)
		.def("get_theta_flow_in",
				 &Hermes<int, double>::get_theta_flow_in<py::array_t<double, 1>>)
		.def("set_theta_flow_out",
				 &Hermes<int, double>::set_theta_flow_out<py::array_t<double, 1>>)
		.def("get_theta_flow_out",
				 &Hermes<int, double>::get_theta_flow_out<py::array_t<double, 1>>)

		.def("set_boundaries",
				 &Hermes<int, double>::set_boundaries<py::array_t<std::uint8_t, 1>>)
		.def("get_boundaries",
				 &Hermes<int, double>::get_boundaries<py::array_t<std::uint8_t, 1>>)

		;
}
