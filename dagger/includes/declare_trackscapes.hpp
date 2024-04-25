#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

template<typename CONNECTOR_T>
void
declare_trackscape(py::module& m, std::string typestr)
{
	py::class_<trackscape<FLOATING_POINT_DAGGER,
												DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
												CONNECTOR_T>>(m, typestr.c_str())
		.def(py::init<>())
		.def_readwrite(
			"graph",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::graph)
		.def_readwrite(
			"connector",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::connector)
		.def("init_random",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::init_random)
		.def("init_perlin",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::init_perlin)
		.def("get_topo",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::template get_topo<py::array>)
		.def("get_hillshade",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::template get_hillshade<py::array>)
		.def("get_h_sed",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::template get_h_sed<py::array>)
		.def("get_Qw",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::template get_Qw<py::array>)
		.def("get_Qs_fluvial",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::template get_Qs_fluvial<py::array>)
		.def("get_Qs_hillslopes",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::template get_Qs_hillslopes<py::array>)
		.def("get_precipitations",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::template get_precipitations<py::array>)
		.def("set_full_stochastic",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_full_stochastic)
		.def("run",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::run)
		.def("run_SFD",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::run_SFD)

		.def("standalone_implicit_SPL",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::standalone_implicit_SPL)
		.def("standalone_cidre_hillslopes",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::standalone_cidre_hillslopes)
		.def("standalone_cidre",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::standalone_cidre)
		.def("standalone_davy2009",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::standalone_davy2009)

		.def("exp_std_cidre_sepA",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::exp_std_cidre_sepA)

		.def("block_uplift",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::block_uplift)
		.def("external_uplift",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template external_uplift<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("init_TSP_module",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template init_TSP_module<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("update_TSP_source",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template update_TSP_source<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("sample_carrot_TSP",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template sample_carrot_TSP<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def(
			"sample_carrot_Ch_MTSI",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::
				template sample_carrot_Ch_MTSI<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("get_transect_Ch_MTSI",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template get_transect_Ch_MTSI<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("get_transect_TSP",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template get_transect_TSP<py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("set_dt",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_dt)
		.def("set_transfer_rate_Qs_hs2fl",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_transfer_rate_Qs_hs2fl)
		.def("set_single_Ks",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_Ks)
		.def("set_single_Kr",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_Kr)
		.def("set_single_Ke",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_Ke)
		.def("set_single_Kle",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_Kle)
		.def("set_single_Kld",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_Kld)
		.def("set_single_depcoeff",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_depcoeff)
		.def("set_single_precipitations",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_precipitations)
		.def("set_single_kappa_s",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_kappa_s)
		.def("set_single_kappa_r",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_kappa_r)
		.def("set_single_Sc",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_Sc)
		.def("set_single_longdep",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_longdep)

		.def("set_single_Sc_M",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_Sc_M)
		.def("set_single_lambda",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_lambda)
		.def("set_single_sea_level",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_sea_level)
		.def("set_single_internal_friction",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_internal_friction)
		.def("set_single_tls",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_single_tls)
		.def("set_hillslopes_mode",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_hillslopes_mode)
		.def("set_fluvial_mode",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_fluvial_mode)
		.def("set_secondary_fluvial_mode",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_secondary_fluvial_mode)
		.def("set_flowtopo_mode",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_flowtopo_mode)
		.def("set_marine_mode",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_marine_mode)
		.def("fill_up",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::fill_up)
		.def("init_Ch_MTSI",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::init_Ch_MTSI)
		.def("rise_boundary_by",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::rise_boundary_by)
		.def("get_TSP_surface_concentrations",
				 &trackscape<
					 FLOATING_POINT_DAGGER,
					 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
					 CONNECTOR_T>::template get_TSP_surface_concentrations<py::array>)
		.def("get_Ch_MTIS_surface_age",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::template get_Ch_MTIS_surface_age<py::array>)
		.def("set_variable_precipitations",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template set_variable_precipitations<
						 py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("set_variable_Kr",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template set_variable_Kr<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("set_variable_Ks",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template set_variable_Ks<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def(
			"set_variable_depcoeff",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::
				template set_variable_depcoeff<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def(
			"set_variable_kappa_s",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::
				template set_variable_kappa_s<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def(
			"set_variable_kappa_r",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::
				template set_variable_kappa_r<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("set_variable_Sc",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template set_variable_Sc<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("set_variable_Sc_M",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template set_variable_Sc_M<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("set_variable_Ke",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template set_variable_Ke<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("set_variable_lambda",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template set_variable_lambda<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def(
			"set_variable_sea_level",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::
				template set_variable_sea_level<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("feed_topo",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::
					 template feed_topo<py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def(
			"set_extra_Qs_fluvial",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::
				template set_extra_Qs_fluvial<py::array_t<int, 1>&,
																			py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def(
			"set_extra_Qw_fluvial",
			&trackscape<FLOATING_POINT_DAGGER,
									DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
									CONNECTOR_T>::
				template set_extra_Qw_fluvial<py::array_t<int, 1>&,
																			py::array_t<FLOATING_POINT_DAGGER, 1>&>)
		.def("disable_Qs_fluvial",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::disable_Qs_fluvial)
		.def("disable_Qw_fluvial",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::disable_Qw_fluvial)

		.def("set_m",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_m)
		.def("set_n",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_n)
		.def("set_N_boundary_to",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::set_N_boundary_to)

		.def("run_SFD_implicit",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::run_SFD_implicit)
		.def("lithify",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::lithify)
		.def("strip_sediment",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::strip_sediment)

		// SEt of standalone functions
		.def("Standalone_hyland_landslides",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::Standalone_hyland_landslides)
		.def("Standalone_hylands_single_landslide",
				 &trackscape<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER, CONNECTOR_T>,
										 CONNECTOR_T>::Standalone_hylands_single_landslide)

		;
}
