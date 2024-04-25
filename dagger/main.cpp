#include "declare_algo.hpp"
#include "declare_connectors.hpp"
#include "declare_dbag.hpp"
#include "declare_enums.hpp"
#include "declare_graphfloods.hpp"
#include "declare_graphs.hpp"
#include "declare_includes.hpp"
#include "declare_parambag.hpp"
#include "declare_popscapes.hpp"
#include "declare_rivnets.hpp"
#include "declare_trackscapes.hpp"

using namespace DAGGER;

PYBIND11_MODULE(dagger, m)
{
	m.doc() = R"pbdoc(
		DAGGER - python API
		===================

		Quick API
		---------

		.. autosummary::

			graph
			graph.init_graph
			graph.set_opt_stst_rerouting
			graph.compute_graph
			graph.is_Sstack_full
			graph.activate_opti_sparse_border_cordonnier
			graph.get_all_nodes_upstream_of
			graph.get_all_nodes_downstream_of
			graph.get_SFD_stack
			graph.get_MFD_stack
			graph.accumulate_constant_downstream_SFD
			graph.accumulate_variable_downstream_SFD
			graph.accumulate_constant_downstream_MFD
			graph.accumulate_variable_downstream_MFD
			graph.set_LMR_method
			graph.set_minimum_slope_for_LMR
			graph.set_slope_randomness_for_LMR
			graph.get_SFD_distance_from_outlets
			graph.get_SFD_min_distance_from_sources
			graph.get_SFD_max_distance_from_sources
			graph.get_MFD_max_distance_from_sources
			graph.get_MFD_min_distance_from_sources
			graph.get_MFD_max_distance_from_outlets
			graph.get_MFD_min_distance_from_outlets
			graph.get_SFD_basin_labels

			D8N
			D8N.__init__
			D8N.set_default_boundaries
			D8N.set_custom_boundaries
			D8N.print_dim
			D8N.get_HS
			D8N.get_mask_array
			D8N.set_values_at_boundaries
			D8N.set_out_boundaries_to_permissive
			D8N.get_boundary_at_node
			D8N.get_rowcol_Sreceivers
			D8N.print_receivers
			D8N.get_rec_array_size
			D8N.update_links_MFD_only
			D8N.update_links
			D8N.update_links_from_topo
			D8N.sum_at_outlets
			D8N.keep_only_at_outlets
			D8N.get_SFD_receivers
			D8N.get_SFD_dx
			D8N.get_SFD_ndonors
			D8N.get_SFD_donors_flat
			D8N.get_SFD_donors_list
			D8N.get_links
			D8N.get_linknodes_flat
			D8N.get_linknodes_list
			D8N.get_linknodes_list_oriented
			D8N.get_SFD_receivers_at_node
			D8N.get_SFD_dx_at_node
			D8N.get_SFD_ndonors_at_node
			D8N.get_SFD_donors_at_node
			D8N.get_SFD_gradient
			D8N.get_links_gradient
			D8N.get_MFD_mean_gradient
			D8N.get_MFD_weighted_gradient
			D8N.get_link_weights
			D8N.set_stochaticiy_for_SFD


		Full API
		---------

		.. autoclass:: D8N
		:members:

		.. autoclass:: graph
		:members:


	)pbdoc";

	delclare_enums(m);
	declare_dbag(m);
	declare_param(m);
	declare_D8connector(m, "D8N");
	declare_graph<D8connector<FLOATING_POINT_DAGGER>>(m, "graph");

	//=============================================================================================
	//=============================================================================================
	//===================== Standalone Algorithms
	//=================================================
	//=============================================================================================
	//=============================================================================================

	declare_algos(m);

	// m.def(
	//   "check_connector_template",
	//   &check_connector_template< D8connector<FLOATING_POINT_DAGGER>,
	//   FLOATING_POINT_DAGGER >
	// );

	declare_popscape_old<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(
		m, "popscape_old");
	declare_popscape<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m, "popscape");
	declare_trackscape<DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m,
																																 "trackscape");

	declare_graphflood<FLOATING_POINT_DAGGER,
										 DAGGER::graph<FLOATING_POINT_DAGGER,
																	 DAGGER::D8connector<FLOATING_POINT_DAGGER>>,
										 DAGGER::D8connector<FLOATING_POINT_DAGGER>>(m,
																																 "graphflood");
	declare_rivnet(m);
};
;

// end of file
