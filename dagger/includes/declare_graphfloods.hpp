#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

template<typename fT, typename GRAPH_T, typename CONNECTOR_T>
void
declare_graphflood(py::module& m, std::string typestr)
{
	py::class_<Graphflood2<int,
												 double,
												 Connector8<int, double>,
												 int,
												 Hermes<int, double>,
												 ParamBag<int, double>>>(m, "GF2")
		.def(py::init<Connector8<int, double>&,
									int&,
									Hermes<int, double>&,
									ParamBag<int, double>&>())
		.def("init",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::init)
		.def("compute_entry_points_from_P",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::compute_entry_points_from_P)
		.def("compute_entry_points_sources",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::compute_entry_points_sources)

		.def("run",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::run)
		.def("run_subgraphflood",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::run_subgraphflood)

		.def("morphoton",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::morphoton)
		.def("smooth_bedrock",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::smooth_bedrock)

		.def("run_dynamic",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::run_dynamic)
		.def("feed_inputQs_with_out",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::feed_inputQs_with_out)

		.def("fillrun_subgraphflood",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::fillrun_subgraphflood)
		.def("anarun_subgraphflood",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::anarun_subgraphflood)
		.def("diffuse_topo",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::diffuse_topo)
		.def("evaluate_convergence",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::evaluate_convergence)
		.def("multiply_Qw_input_points",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::multiply_Qw_input_points)
		.def("get_entry_node_PQ",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::get_entry_node_PQ)
		.def("get_input_node_Qw",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::get_input_node_Qw)
		.def("get_input_Qw",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::get_input_Qw)
		.def("get_input_Qs",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::get_input_Qs)
		.def("standalone_Qwin_D8",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::
					 template standalone_Qwin_D8<py::array_t<double, 1>>)

		.def("get_active_nodes",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::get_active_nodes)

		.def("compute_qr",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::compute_qr)

		.def("chunk_by_distance_to_outlet",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::chunk_by_distance_to_outlet)

		.def("solve_analytically_if",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::solve_analytically_if)
		.def("get_MB_Qwin_out",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::get_MB_Qwin_out)
		.def("get_meandhstar",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::get_meandhstar)

		.def("sum_Qw_inputs",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::sum_Qw_inputs)

		.def("prepare_tsg",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::prepare_tsg)

		.def("standalone_Qwin",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::standalone_Qwin)

		.def("run_tinysubgraph",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::run_tinysubgraph)
		.def("run_tinysubgraph_dyn",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::run_tinysubgraph_dyn)

		.def("run_subgraphflood_expA",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::run_subgraphflood_expA)
		.def("set_uniform_P",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::set_uniform_P)
		.def("set_Qw_input_points",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::
					 set_Qw_input_points<py::array_t<int, 1>, py::array_t<double, 1>>)
		.def("set_Qw_Qsinput_points",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::
					 set_QwQs_input_points<py::array_t<int, 1>, py::array_t<double, 1>>)
		.def("set_dt",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::set_dt)
		.def("initial_fill",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::initial_fill)
		.def("get_dt",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::get_dt)
		.def("bluplift",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::bluplift)

		.def("_quick_slipos_from_point",
				 &Graphflood2<int,
											double,
											Connector8<int, double>,
											int,
											Hermes<int, double>,
											ParamBag<int, double>>::_quick_slipos_from_point)

		;

	py::class_<graphflood<fT, GRAPH_T, CONNECTOR_T>>(m, typestr.c_str())
		.def(py::init<GRAPH_T&, CONNECTOR_T&>())
		.def("run",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::run,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("dynarun",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::dynarun,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("set_maxdHw",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_maxdHw,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("set_mindHw",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_mindHw,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("set_out_QS_as_input_QS",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_out_QS_as_input_QS)

		.def("set_out_QS_as_input_QS_min",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_out_QS_as_input_QS_min)

		.def("set_glob_dynatime",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_glob_dynatime,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("get_glob_dynatime",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_glob_dynatime,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def(
			"setup_entry_points_dynagraph_farea",
			&graphflood<fT, GRAPH_T, CONNECTOR_T>::setup_entry_points_dynagraph_farea,
			R"pdoc(Main function running the model from all the input params)pdoc")

		.def("set_topo",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template set_topo<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("set_hw",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template set_hw<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("get_hw",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_hw<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("get_surface_topo",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_surface_topo<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("get_bedrock_topo",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_bedrock_topo<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("get_Qwin",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_Qwin<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("get_Qs",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_Qs<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("get_SSTACKDEBUG",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_SSTACKDEBUG<
					 py::array_t<size_t, 1>>,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("gen_rdid",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::gen_rdid,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("get_rdid",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_rdid,
				 R"pdoc(Main function running the model from all the input params)pdoc")

		.def("enable_MFD",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_MFD,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("enable_SFD",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_SFD,
				 R"pdoc(Main function running the model from all the input params)pdoc")
		.def("set_dt_hydro", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_dt_hydro)
		.def("fill_minima", &graphflood<fT, GRAPH_T, CONNECTOR_T>::fill_minima)
		.def("reroute_minima",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::reroute_minima)
		.def("ignore_minima", &graphflood<fT, GRAPH_T, CONNECTOR_T>::ignore_minima)
		.def("enable_morpho", &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_morpho)
		.def("disable_morpho",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_morpho)
		.def("set_dt_morpho_multiplier",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_dt_morpho_multiplier)

		.def("set_dt_morpho", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_dt_morpho)
		.def("set_single_aexp",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_single_aexp)
		.def("set_single_ke", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_single_ke)
		.def("set_single_ke_lateral",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_single_ke_lateral)
		.def("set_single_kd", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_single_kd)
		.def("set_single_kd_lateral",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_single_kd_lateral)
		.def("set_single_tau_c",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_single_kd_lateral)
		.def("set_variable_ke",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template set_variable_ke<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("set_mannings", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_mannings)

		.def("init_convergence_checker",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::init_convergence_checker)
		.def("get_conv_nodes",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_conv_nodes<
					 py::array_t<int, 1>>)
		.def("get_conv_ini_Qw",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_conv_ini_Qw<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("get_conv_mean_Qr",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_conv_mean_Qr<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)
		.def("get_conv_mean_dhw",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_conv_mean_dhw<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("compute_tuqQ",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template compute_tuqQ<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("compute_AD8",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template compute_AD8<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("compute_AD8_maxQw",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template compute_AD8_maxQw<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def(
			"compute_AD8_stochastic_Qw",
			&graphflood<fT, GRAPH_T, CONNECTOR_T>::template compute_AD8_stochastic_Qw<
				py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def(
			"compute_QW8_stochastic_Qw",
			&graphflood<fT, GRAPH_T, CONNECTOR_T>::template compute_QW8_stochastic_Qw<
				py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def(
			"compute_AD8_stochastic_Sw",
			&graphflood<fT, GRAPH_T, CONNECTOR_T>::template compute_AD8_stochastic_Sw<
				py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("compute_elemental_transfer",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::
					 template compute_elemental_transfer<
						 py::array_t<FLOATING_POINT_DAGGER, 1>,
						 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("set_water_input_by_entry_points",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::
					 template set_water_input_by_entry_points<
						 py::array_t<FLOATING_POINT_DAGGER, 1>,
						 py::array_t<int, 1>>)
		.def("set_water_input_by_constant_precipitation_rate",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::
					 set_water_input_by_constant_precipitation_rate)
		.def("set_water_input_by_variable_precipitation_rate",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::
					 template set_water_input_by_variable_precipitation_rate<
						 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("set_sed_input_by_entry_points",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::
					 template set_sed_input_by_entry_points<
						 py::array_t<FLOATING_POINT_DAGGER, 1>,
						 py::array_t<int, 1>>)

		.def("enable_Qwout_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_Qwout_recording)
		.def("disable_Qwout_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_Qwout_recording)
		.def("get_Qwout_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_Qwout_recording<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("enable_Sw_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_Sw_recording)
		.def("disable_Sw_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_Sw_recording)
		.def("get_Sw_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_Sw_recording<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("enable_dhw_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_dhw_recording)
		.def("disable_dhw_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_dhw_recording)
		.def("get_dhw_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_dhw_recording<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("enable_filling_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_filling_recording)
		.def("disable_filling_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_filling_recording)
		.def("get_filling_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_filling_recording<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("enable_edot_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_edot_recording)
		.def("disable_edot_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_edot_recording)
		.def("get_edot_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_edot_recording<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("enable_flowvec_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_flowvec_recording)
		.def("disable_flowvec_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_flowvec_recording)
		.def("get_flowvec_recording",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template get_flowvec_recording<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("get_tot_Qw_input",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_tot_Qw_input)
		.def("get_tot_Qw_output",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_tot_Qw_output)
		.def("get_tot_Qwin_output",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_tot_Qwin_output)
		.def("get_tot_Qs_output",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_tot_Qs_output)

		.def("set_stochaslope",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_stochaslope)
		.def("disable_stochaslope",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_stochaslope)
		.def("set_fixed_hw_at_boundaries",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_fixed_hw_at_boundaries)
		.def("set_fixed_slope_at_boundaries",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_fixed_slope_at_boundaries)
		.def("get_dt_hydro", &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_dt_hydro)

		.def("set_partition_method",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_partition_method)
		.def("set_topological_number",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_topological_number)
		.def("get_topological_number",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_topological_number)

		.def("get_courant_dt_hydro",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_courant_dt_hydro)
		.def("set_courant_number",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_courant_number)
		.def("set_max_courant_dt_hydro",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_max_courant_dt_hydro)
		.def("set_min_courant_dt_hydro",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_min_courant_dt_hydro)
		.def("enable_courant_dt_hydro",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_courant_dt_hydro)
		.def("disable_courant_dt_hydro",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_courant_dt_hydro)
		.def("set_Qwin_crit", &graphflood<fT, GRAPH_T, CONNECTOR_T>::set_Qwin_crit)
		.def("get_nT", &graphflood<fT, GRAPH_T, CONNECTOR_T>::get_nT)

		.def("enable_hydrostationary",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::enable_hydrostationary)
		.def("disable_hydrostationary",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::disable_hydrostationary)

		.def("block_uplift", &graphflood<fT, GRAPH_T, CONNECTOR_T>::block_uplift)
		.def("variable_uplift",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::template variable_uplift<
					 py::array_t<FLOATING_POINT_DAGGER, 1>>)

		.def("run_precipitions",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::run_precipitions)
		.def("run_precipitions_exp",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::run_precipitions_exp)
		.def("run_graphipiton",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::run_graphipiton)
		.def("run_exp", &graphflood<fT, GRAPH_T, CONNECTOR_T>::run_exp)
		.def("run_hydro_only",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::run_hydro_only)

		.def("define_precipitations_Ath",
				 &graphflood<fT, GRAPH_T, CONNECTOR_T>::define_precipitations_Ath)

		;
}
