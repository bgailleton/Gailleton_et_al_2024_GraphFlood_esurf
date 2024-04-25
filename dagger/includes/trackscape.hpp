//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Trackscape main code
//
// Trackscape is a static CHONK implementation using DAGGER as backend
// I tuned it for tracking provenance and sediment timing  at relatively long
// term. Fluvial and hillslope processes are using laws adapted from CIDRE
// (Carretier et al. 2016)
//
// B. Gailleton - 2022
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

#ifndef TRACKSCAPE_HPP
#define TRACKSCAPE_HPP

// STL imports
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <stack>
#include <stdlib.h>
#include <string>
#include <thread>
#include <vector>

// local includes
// -> General routines and data structures
#include "utils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> The connector classes
#include "D4connector.hpp"
#include "D8connector.hpp"
#include "hillshading.hpp"
#include "popscape_utils.hpp"

#include "trackscape_enum.hpp"
#include "trackscape_utils.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

// #include "sedcells.hpp"

namespace DAGGER {
// {
// 	typedef void (*tscfunc)();
// 	typedef int (*tscfuncstack)(int);

// templates are: floating point type, graph type (DAGGER graph) and connector
// type (DAGGER connector)
template<class fT, class Graph_t, class Connector_t>
class trackscape
{
public:
	// ###############################################
	// ###############################################
	// Fluxes and quantities
	// ###############################################
	// ###############################################

	// Surface topography (top of bedrock + sediment)
	std::vector<fT> z_surf;
	// Sediment height
	std::vector<fT> h_sed;
	// Water Discharge (as P * drainage area)
	std::vector<fT> Qw;
	// Hillslope sediment flux
	std::vector<fT> Qs_hs;
	// Fluvial sediment flux
	std::vector<fT> Qs_fluvial;

	// ###############################################
	// ###############################################
	// Intermediate vectors storing vertical changes to apply to the model later
	// on
	// ###############################################
	// ###############################################
	std::vector<fT> dbedrockdt, dhsdt;

	// ###############################################
	// ###############################################
	// parameters for the different functions
	// ###############################################
	// ###############################################

	// TSD_FLUVIAL::DAVY2019
	// Spatial erodibility (sediments)
	std::vector<fT> _Ks;
	// Spatial erodibility (bedrock)
	std::vector<fT> _Kr;
	// Spatial deposition coefficient
	std::vector<fT> _depcoeff;

	// TSD_FLUVIAL::LATERALDAVY
	std::vector<fT> _Kle;
	std::vector<fT> _Kld;

	// TSD_HILLSLOPES::CIDRE
	// Spatial diffusivity (sediments)
	std::vector<fT> _kappa_s;
	// Spatial diffusivity (bedrock)
	std::vector<fT> _kappa_r;
	// Spatial critical slope
	std::vector<fT> _Sc;
	std::vector<fT> _longdep;

	// TSD_HILLSLOPES::HYLANDS
	// soil cohesion
	std::vector<fT> _C;
	std::vector<fT> _rho;
	std::vector<fT> _tls;
	std::vector<fT> _internal_friction;
	fT min_dep_hillslopes = 0.01;
	fT max_dep_hillslopes = 0.8;
	fT gravity = 9.81;

	bool transfer_Qs_hs2fl = false;
	fT transfer_rate_Qs_hs2fl = 0.1;
	void set_transfer_rate_Qs_hs2fl(fT stuff)
	{
		this->transfer_rate_Qs_hs2fl = stuff;
	}

	// Spatial precipitations
	std::vector<fT> _precipitations;

	// Spatial marine _Sc
	std::vector<fT> _Sc_M;
	// Spatial marine Ke
	std::vector<fT> _Ke;
	// Spatial marine lamdba
	std::vector<fT> _lambda;
	// Spatial marine sea_level
	std::vector<fT> _sea_level;

	// SPL exponents
	fT mexp = 0.45;
	fT nexp = 1;

	// Switching on and/or off the hillslope and/or fluvial
	TSC_FLOW_TOPOLOGY flowtopo = TSC_FLOW_TOPOLOGY::SFD;
	TSC_HILLSLOPE hillslopes = TSC_HILLSLOPE::NONE;
	TSC_FLUVIAL fluvial = TSC_FLUVIAL::DAVY2009;
	TSC_FLUVIAL secondary_fluvial = TSC_FLUVIAL::NONE;
	TSC_MARINE marine = TSC_MARINE::NONE;

	// bool stochaslope = false;
	// void enable_stochaslope(){this->stochaslope = true;}
	// void disable_stochaslope(){this->stochaslope = false;}

	// Switching on and/or off the spatiality of the different parameters
	bool variable_Ks = false;
	bool variable_Kr = false;
	bool variable_Kle = false;
	bool variable_Kld = false;
	bool variable_kappa_s = false;
	bool variable_kappa_r = false;
	bool variable_precipitations = false;
	bool variable_depcoeff = false;
	bool variable_Sc = false;
	bool variable_longdep = false;
	bool variable_Sc_M = false;
	bool variable_Ke = false;
	bool variable_lambda = false;
	bool variable_sea_level = false;
	bool variable_C = false;
	bool variable_rho = false;
	bool variable_tls = false;
	bool variable_internal_friction = false;

	// at least one tracking module is actiated if this is true
	// and it can conflict with some functions
	bool at_least_one_tracking_module_is_activated = false;

	// Hylands-specific parameters
	std::vector<std::uint32_t> landslidesid;

	// Module TSP:
	// TRACKING SINGLE SOURCE
	// # IS it activated
	bool TSP_module = false;
	// # source concentration (bedrock)
	std::vector<fT> TSP_concentrations;
	// # tracking the suspended concentrations
	std::vector<BasePropStorer<fT>> TSP_Qsf, TSP_Qsh;
	// # Stratigraphic info
	VerticalStorer<fT, BasePropStorer<fT>> TSP_store;

	// Module Ch_MTSI
	// CHRONO MEAN TIME SINCE INCISION
	// # Activated??
	bool Ch_MTSI = false;
	// # tracking the mobile timer (hs vs fl)
	std::vector<BasePropStorer<fT>> Ch_MTSI_Qsf, Ch_MTSI_Qsh;
	// # tracking the strati timer (hs vs fl)
	VerticalStorer<fT, BasePropStorer<fT>> Ch_MTSI_store;

	// The graph (DAGGER)
	Graph_t* graph;

	// The connector (DAGGER)
	Connector_t* connector;

	// Other params
	fT noise_magnitude = 1.;

	fT tdt = 1e2;
	void set_dt(fT dt) { this->tdt = dt; };

	bool add_extra_Qs_fluvial = false;
	std::vector<int> extra_Qsf_nodes;
	std::vector<fT> extra_Qsf;

	bool add_extra_Qw_fluvial = false;
	std::vector<int> extra_Qwf_nodes;
	std::vector<fT> extra_Qwf;

	bool full_stochastic = false;
	void set_full_stochastic(bool val) { this->full_stochastic = val; }

	// Private parameters for local processing
	// Basically they simplify my life when juggling between function so ignore
	// them
private:
	fT tdx, tdy, tZ, tSS, ths, tqs, tqw;

	fT thEs, thEr, thDs, tfEs, tfEr, tfDs;

	int tnode, tSrec;

	int tnn, trn, tdb;

	int orthonodA, orthonodB;

	std::array<int, 8> tneighbours_node, tneighbours_links, treceivers_nodes,
		treceivers_links, tdonors_nodes, tdonors_links;

	std::array<fT, 8> tslopes;

	// int ( trackscape<fT,Graph_t, Connector_t>::* ) (int) stackgetter;
	std::vector<void (trackscape<fT, Graph_t, Connector_t>::*)()> prefuncs;
	std::vector<void (trackscape<fT, Graph_t, Connector_t>::*)()> downstreamfuncs;
	std::vector<void (trackscape<fT, Graph_t, Connector_t>::*)()>
		downstreamfuncs_marine;
	std::vector<void (trackscape<fT, Graph_t, Connector_t>::*)()> upstreamfuncs;
	// std::array<fT,8> tfneighbours, tfneighboursA,tfneighboursB;

public:
	// empty constructor with default values for the different param (default =
	// non spatial)
	trackscape()
	{
		this->_Ks = { 2e-5 };
		this->_Kr = { 1e-5 };
		this->_Kle = { 0.1 };
		this->_Kld = { 0.1 };
		this->_kappa_s = { 2e-3 };
		this->_kappa_r = { 1e-3 };
		this->_precipitations = { 1. };
		this->_depcoeff = { 10. };
		this->_Sc = { 0.6 };
		this->_longdep = { 10. };
		this->_Sc_M = { 0.6 };
		this->_Ke = { 2e-2 };
		this->_lambda = { 1000 };
		this->_sea_level = { 0. };
		this->_C = { 15 };
		this->_rho = { 2700 };
		this->_tls = { 5000 };
		this->_internal_friction = { 0.55 };
	};

	// constructor 1:
	// -> DEPRECATED: NOW THIS FUNCTION ONLY DOES WHITE NOISE noise type is from
	// the RANDOISE enum (white, red, perlin, ...)
	// -> nx/ny are the number of nodes in the x/y dir
	// -> dx dy are the related spacing in [L]
	void init_random(int nx,
									 int ny,
									 fT dx,
									 fT dy,
									 std::string boundary_conditions)
	{

		// total number of nodes (so far assuming only regular grids)
		int nxy = nx * ny;

		// init the topo to 0
		this->z_surf = std::vector<fT>(nxy, 0.);

		// init connector
		this->connector = _create_connector_fromptr<fT>(nx, ny, dx, dy, 0., 0.);

		// boundary conditions:
		this->connector->set_default_boundaries(boundary_conditions);

		// init graph
		// _create_graph(nxy, *(this->connector), *(this->graph));
		_create_graph_fromptr(nxy, *(this->connector), this->graph);
		// this->graph->init_graph(this->connector);

		// init random noise
		DAGGER::add_noise_to_vector(
			this->z_surf, static_cast<fT>(0.), noise_magnitude);

		this->h_sed = std::vector<fT>(this->graph->nnodes, static_cast<fT>(0.));

		for (int i = 0; i < this->connector->nnodes; ++i) {
			if (this->connector->boundaries.can_out(i))
				this->z_surf[i] = 0;
		}
	}

	// constructor 1:
	// -> DEPRECATED: NOW THIS FUNCTION ONLY DOES WHITE NOISE noise type is from
	// the RANDOISE enum (white, red, perlin, ...)
	// -> nx/ny are the number of nodes in the x/y dir
	// -> dx dy are the related spacing in [L]
	void init_perlin(int nx,
									 int ny,
									 fT dx,
									 fT dy,
									 std::string boundary_conditions,
									 fT frequency,
									 int octaves,
									 fT Zmax,
									 std::uint32_t seed,
									 bool noise_on_top)
	{

		// total number of nodes (so far assuming only regular grids)
		int nxy = nx * ny;

		// init the topo to 0

		// init connector
		this->connector = _create_connector_fromptr<fT>(nx, ny, dx, dy, 0., 0.);

		// boundary conditions:
		this->connector->set_default_boundaries(boundary_conditions);

		this->z_surf =
			_generate_perlin_noise_2D(*(this->connector), frequency, octaves, seed);
		for (auto& v : this->z_surf)
			v *= Zmax;

		// init graph
		_create_graph_fromptr(nxy, *(this->connector), this->graph);
		// std::cout << this->graph->nnodes << std::endl;

		// init random noise
		if (noise_on_top)
			DAGGER::add_noise_to_vector<std::vector<fT>, fT>(
				this->z_surf, 0., noise_magnitude);

		this->h_sed = std::vector<fT>(this->graph->nnodes, 0.);

		for (int i = 0; i < this->connector->nnodes; ++i) {
			if (this->connector->boundaries.can_out(i))
				this->z_surf[i] = 0;
		}
	}

	// constructor 1:
	// -> DEPRECATED: NOW THIS FUNCTION ONLY DOES WHITE NOISE noise type is from
	// the RANDOISE enum (white, red, perlin, ...)
	// -> nx/ny are the number of nodes in the x/y dir
	// -> dx dy are the related spacing in [L]
	void _init_from_con_and_Z(Connector_t& connector, std::vector<fT>& ttopo)
	{

		// init the topo to 0
		this->z_surf = ttopo;

		// init connector
		this->connector = &connector;

		// init graph
		_create_graph_fromptr(
			this->connector->nnodes, *(this->connector), this->graph);
		// this->graph->init_graph(this->connector);

		this->h_sed = std::vector<fT>(this->connector->nnodes, static_cast<fT>(0.));
	}

	// switch on and off the hillslopes and fluvial processes
	void set_hillslopes_mode(TSC_HILLSLOPE mode) { this->hillslopes = mode; };
	void set_fluvial_mode(TSC_FLUVIAL mode) { this->fluvial = mode; };
	void set_secondary_fluvial_mode(TSC_FLUVIAL mode)
	{
		this->secondary_fluvial = mode;
	};
	void set_marine_mode(TSC_MARINE mode) { this->marine = mode; }
	void set_flowtopo_mode(TSC_FLOW_TOPOLOGY mode) { this->flowtopo = mode; }

	// Set the single parameters
	void set_single_Ks(fT tks)
	{
		this->_Ks = { tks };
		this->variable_Ks = false;
	}
	void set_single_Kr(fT tkr)
	{
		this->_Kr = { tkr };
		this->variable_Kr = false;
	}
	void set_single_Kle(fT tKle)
	{
		this->_Kle = { tKle };
		this->variable_Kle = false;
	}
	void set_single_Kld(fT tKld)
	{
		this->_Kld = { tKld };
		this->variable_Kld = false;
	}
	void set_single_depcoeff(fT tdepcoeff)
	{
		this->_depcoeff = { tdepcoeff };
		this->variable_depcoeff = false;
	}
	void set_single_precipitations(fT tprecipitations)
	{
		this->_precipitations = { tprecipitations };
		this->variable_precipitations = false;
	}
	void set_single_kappa_s(fT tkappa_s)
	{
		this->_kappa_s = { tkappa_s };
		this->variable_kappa_s = false;
	}
	void set_single_kappa_r(fT tkappa_r)
	{
		this->_kappa_r = { tkappa_r };
		this->variable_kappa_r = false;
	}
	void set_single_Sc_M(fT tSc_M)
	{
		this->_Sc_M = { tSc_M };
		this->variable_Sc_M = false;
	}
	void set_single_Sc(fT tSc)
	{
		this->_Sc = { tSc };
		this->variable_Sc = false;
	}
	void set_single_longdep(fT tlongdep)
	{
		this->_longdep = { tlongdep };
		this->variable_longdep = false;
	}
	void set_single_Ke(fT TKe)
	{
		this->_Ke = { TKe };
		this->variable_Ke = false;
	}
	void set_single_lambda(fT tlambda)
	{
		this->_lambda = { tlambda };
		this->variable_lambda = false;
	}
	void set_single_sea_level(fT tsea_level)
	{
		this->_sea_level = { tsea_level };
		this->variable_sea_level = false;
	}
	void set_single_internal_friction(fT tinternal_friction)
	{
		this->_internal_friction = { tinternal_friction };
		this->variable_internal_friction = false;
	}
	void set_single_tls(fT ttls)
	{
		this->_tls = { ttls };
		this->variable_tls = false;
	}

	void set_m(fT m) { this->mexp = m; }
	void set_n(fT n) { this->nexp = n; }

	template<class in_t>
	void set_variable_Kr(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Kr = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_Kr = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_Kr[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_Ks(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Ks = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_Ks = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_Ks[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_precipitations(in_t& prec)
	{
		auto tprec = DAGGER::format_input(prec);
		this->_precipitations = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_precipitations = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_precipitations[i] = tprec[i];
		}
	}

	template<class in_t>
	void set_variable_depcoeff(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_depcoeff = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_depcoeff = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_depcoeff[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_kappa_s(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_kappa_s = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_kappa_s = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_kappa_s[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_kappa_r(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_kappa_r = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_kappa_r = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_kappa_r[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_Sc(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Sc = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_Sc = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_Sc[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_Sc_M(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Sc_M = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_Sc_M = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_Sc_M[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_Ke(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_Ke = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_Ke = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_Ke[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_lambda(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_lambda = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_lambda = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_lambda[i] = tarr[i];
		}
	}

	template<class in_t>
	void set_variable_sea_level(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->_sea_level = std::vector<fT>(this->connector->nnodes, 0.);
		this->variable_sea_level = true;
		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->_sea_level[i] = tarr[i];
		}
	}

	template<class in_t>
	void feed_topo(in_t& arr)
	{
		auto tarr = DAGGER::format_input(arr);
		this->z_surf = std::vector<fT>(this->connector->nnodes, 0.);

		for (int i = 0; i < this->graph->nnodes; ++i) {
			this->z_surf[i] = tarr[i];
		}
	}

	template<class in_t, class ft_t>
	void set_extra_Qs_fluvial(in_t& arrnodes, ft_t& arrval)
	{
		this->add_extra_Qs_fluvial = true;
		auto tarrnodes = DAGGER::format_input(arrnodes);
		auto tarrval = DAGGER::format_input(arrval);
		this->extra_Qsf_nodes = DAGGER::to_vec(tarrnodes);
		this->extra_Qsf = DAGGER::to_vec(tarrval);
	}

	void disable_Qs_fluvial() { this->add_extra_Qs_fluvial = false; }

	template<class in_t, class ft_t>
	void set_extra_Qw_fluvial(in_t& arrnodes, ft_t& arrval)
	{
		this->add_extra_Qw_fluvial = true;
		auto tarrnodes = DAGGER::format_input(arrnodes);
		auto tarrval = DAGGER::format_input(arrval);
		this->extra_Qwf_nodes = DAGGER::to_vec(tarrnodes);
		this->extra_Qwf = DAGGER::to_vec(tarrval);
	}

	void disable_Qw_fluvial() { this->add_extra_Qw_fluvial = false; }

	// ------------------------------------------------
	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   Access to param value           /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|

	// All the methods related to accessing parameter values in a generic way
	// No matter which parametrisation method is used, these functions are called
	// to get the parameter from its node indice If spatial is not activated it
	// returns the single value; if spatial it returns the value at teh location
	// other option would be also called from here (e.g.K f(Qs))
	// ------------------------------------------------

	fT Ks(int i)
	{
		if (variable_Ks == false)
			return this->_Ks[0];
		else
			return this->_Ks[i];
	}

	fT Kr(int i)
	{
		if (variable_Kr == false)
			return this->_Kr[0];
		else
			return this->_Kr[i];
	}

	fT Kle(int i)
	{
		if (variable_Kle == false)
			return this->_Kle[0];
		else
			return this->_Kle[i];
	}

	fT Kld(int i)
	{
		if (variable_Kld == false)
			return this->_Kld[0];
		else
			return this->_Kld[i];
	}

	fT kappa_s(int i)
	{
		if (variable_kappa_s == false)
			return this->_kappa_s[0];
		else
			return this->_kappa_s[i];
	}

	fT kappa_r(int i)
	{
		if (variable_kappa_r == false)
			return this->_kappa_r[0];
		else
			return this->_kappa_r[i];
	}

	fT precipitations(int i)
	{
		if (this->variable_precipitations == false)
			return this->_precipitations[0];
		else
			return this->_precipitations[i];
	}

	fT depcoeff(int i)
	{
		if (variable_depcoeff == false)
			return this->_depcoeff[0];
		else
			return this->_depcoeff[i];
	}

	fT Sc(int i)
	{
		if (variable_Sc == false)
			return this->_Sc[0];
		else
			return this->_Sc[i];
	}

	fT longdep(int i)
	{
		if (variable_longdep == false)
			return this->_longdep[0];
		else
			return this->_longdep[i];
	}

	fT Sc_M(int i)
	{
		if (variable_Sc_M == false)
			return this->_Sc_M[0];
		else
			return this->_Sc_M[i];
	}

	fT Ke(int i)
	{
		if (variable_Ke == false)
			return this->_Ke[0];
		else
			return this->_Ke[i];
	}

	fT lambda(int i)
	{
		if (variable_lambda == false)
			return this->_lambda[0];
		else
			return this->_lambda[i];
	}

	fT sea_level(int i)
	{
		if (variable_sea_level == false)
			return this->_sea_level[0];
		else
			return this->_sea_level[i];
	}

	fT C(int i)
	{
		if (variable_C == false)
			return this->_C[0];
		else
			return this->_C[i];
	}

	fT rho(int i)
	{
		if (variable_rho == false)
			return this->_rho[0];
		else
			return this->_rho[i];
	}

	fT tls(int i)
	{
		if (variable_tls == false)
			return this->_tls[0];
		else
			return this->_tls[i];
	}

	fT internal_friction(int i)
	{
		if (variable_internal_friction == false)
			return this->_internal_friction[0];
		else
			return this->_internal_friction[i];
	}

	// Degree version
	fT internal_friction_d(int i)
	{
		if (variable_internal_friction == false)
			return std::atan(this->_internal_friction[0]);
		else
			return std::atan(this->_internal_friction[i]);
	}

	fT Sc_hylands(int i, fT& slope)
	{
		return (this->internal_friction(i) + slope) / 2;
	}

	// ------------------------------------------------
	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   Initialisation METHODS          /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|
	//
	// All the methods related to initialising things
	// ------------------------------------------------

	void fill_up()
	{
		this->graph->depression_resolver = DAGGER::DEPRES::cordonnier_fill;
		// std::vector<bool> testlinkchange(this->graph->links),
		// nnodes2change(this->graph->nnodes,false);
		this->graph->_compute_graph(this->z_surf, true, false);
	}

	void init_vectors()
	{
		this->Qw = std::vector<fT>(this->graph->nnodes, 0.);

		if (this->fluvial != TSC_FLUVIAL::NONE) {
			this->Qs_fluvial = std::vector<fT>(this->graph->nnodes, 0.);
		}

		if (this->hillslopes != TSC_HILLSLOPE::NONE ||
				this->marine != TSC_MARINE::NONE) {
			this->Qs_hs = std::vector<fT>(this->graph->nnodes, 0.);
		}

		if (this->TSP_module)
			this->init_Qs_TSP();

		if (this->Ch_MTSI)
			this->init_Qs_Ch_MTSI();

		this->downstreamfuncs.clear();
		this->upstreamfuncs.clear();

		this->dbedrockdt = std::vector<fT>(this->graph->nnodes, 0.);
		this->dhsdt = std::vector<fT>(this->graph->nnodes, 0.);
	}

	// ------------------------------------------------
	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   Run METHODS          /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|
	//
	// All the methods related to running a timestep
	// ------------------------------------------------

	void run()
	{
		// Initialising the data vectors
		this->init_vectors();

		bool need_mfrecs = this->do_I_need_MFD();

		// Setting up the function pointers
		this->determine_functors();
		auto stackgetter =
			(this->flowtopo == TSC_FLOW_TOPOLOGY::SFD)
				? &trackscape<fT, Graph_t, Connector_t>::get_istack_node_SFD
				: &trackscape<fT, Graph_t, Connector_t>::get_istack_node_MFD;

		std::vector<fT> faketopo(this->z_surf);
		this->graph->_compute_graph(faketopo, !need_mfrecs, false);

		if (this->marine != TSC_MARINE::NONE)
			this->fix_small_marine_patches();

		// this->graph->_compute_graph(this->z_surf, !need_mfrecs, false);

		if (this->prefuncs.size() > 0)
			for (auto& v : this->prefuncs)
				(this->*v)();

		if (this->downstreamfuncs.size() > 0) {
			for (int i = this->graph->nnodes - 1; i >= 0; --i) {
				// std::cout << std::endl << "RUN" << std::endl;

				// Getting the next node in line
				this->tnode = (this->*stackgetter)(i);

				// Check if no data
				if (this->connector->boundaries.no_data(this->tnode))
					continue;

				// Check if base level
				if (this->connector->flow_out_or_pit(this->tnode)) {
					// manage the base level evolution here
					continue;
				}

				this->_ready_node_state();

				if (this->marine == TSC_MARINE::NONE ||
						this->z_surf[this->tnode] >= this->sea_level(this->tnode)) {
					for (auto& v : this->downstreamfuncs)
						(this->*v)();
				} else {
					for (auto& v : this->downstreamfuncs_marine)
						(this->*v)();
				}

				// Applying the fluxes
				// this->h_sed[this->tnode] += (this->thDs + this->tfDs - this->thEs -
				// this->tfEs) * this->tdt; this->z_surf[this->tnode] += (this->thDs +
				// this->tfDs - this->thEs - this->thEr - this->tfEs - this->tfEr) *
				// this->tdt; std::cout << (this->thDs + this->tfDs - this->thEs -
				// this->thEr - this->tfEs - this->tfEr) * this->tdt << std::endl;
				// if((this->thDs + this->tfDs - this->thEs - this->thEr - this->tfEs -
				// this->tfEr) * this->tdt > 1e-1) std::cout <<
				// this->Qs_fluvial[this->tnode]/this->Qw[this->tnode] << "|" <<
				// this->tfDs << "|" << this->tfEs << "|" << tfEr << std::endl;;
				// this->add_to_dbedrockdt(this->tnode, this->tdt * (this->tfDs +
				// this->thDs - this->tfEr - this->tfEs - this->thEr - this->thEs) );
			}
		}

		this->apply_dzhsdt();

		// fT tmax = 0.;
		// for(auto v:this->h_sed)
		// {
		// 	if(std::abs(v) >tmax) tmax = v;
		// }
		// std::cout << "maxelev:" << tmax << std::endl;
		// std::cout  << std::endl;

		if (this->upstreamfuncs.size() > 0) {
			for (int i = 0; i < this->graph->nnodes; ++i) {
				// Getting the next node in line
				this->tnode = (this->*stackgetter)(i);

				// Check if no data
				if (this->connector->boundaries.no_data(this->tnode))
					continue;

				// Check if base level
				if (this->connector->flow_out_or_pit(this->tnode)) {
					// manage the base level evolution here
					continue;
				}

				this->tSrec = this->connector->_Sreceivers[this->tnode];

				// feeding the local private variables
				this->tdx = this->connector->Sdistance2receivers[this->tnode];
				this->tdy = this->connector->get_travers_dy_from_dx(this->tdx);
				this->tZ = this->z_surf[this->tnode];
				// this->tSS = this->connector->SS[this->tnode];
				this->tSS = std::max(static_cast<fT>(1e-6),
														 (tZ - this->z_surf[this->tSrec]) / this->tdx);
				this->ths = this->h_sed[this->tnode];

				if (need_mfrecs)
					this->trn = this->connector->get_receivers_idx_nodes_and_links(
						this->tnode, this->treceivers_nodes, this->treceivers_links);

				this->reset_fhED();

				for (auto& v : this->upstreamfuncs)
					(this->*v)();
			}
		}
	}

	void determine_functors()
	{

		this->downstreamfuncs.clear();
		this->downstreamfuncs_marine.clear();
		this->upstreamfuncs.clear();
		this->prefuncs.clear();

		if (this->hillslopes == TSC_HILLSLOPE::HYLANDS)
			this->prefuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::hillslopes_hylands_trigger);

		// if(this->marine != TSC_MARINE::NONE)
		// 	this->prefuncs.emplace_back(&trackscape<fT,Graph_t,
		// Connector_t>::fix_small_marine_patches);

		if (this->add_extra_Qs_fluvial)
			this->prefuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::prec_extra_Qs_fluvial);

		if (this->add_extra_Qw_fluvial)
			this->prefuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::prec_extra_Qw_fluvial);

		const bool need_Qw = this->fluvial != TSC_FLUVIAL::NONE;
		const bool need_Qs = this->fluvial == TSC_FLUVIAL::DAVY2009 ||
												 this->hillslopes != TSC_HILLSLOPE::NONE;
		const bool SFD = this->flowtopo == TSC_FLOW_TOPOLOGY::SFD;

		if (need_Qw) {
			if (SFD)
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::prec_Qw_SFD);
			else
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::prec_Qw_MFD);
		}

		// Switch statements when using enums
		if (this->hillslopes == TSC_HILLSLOPE::CIDRE)
			this->downstreamfuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::hillslopes_cidre);
		else if (this->hillslopes == TSC_HILLSLOPE::HYLANDS)
			this->downstreamfuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::
					hillslopes_cidre_dep_only_for_highlands);
		else if (this->hillslopes == TSC_HILLSLOPE::CIDRE_NOCRIT)
			this->downstreamfuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::hillslopes_cidre_no_crit);
		// end the switch with exception

		if (this->marine == TSC_MARINE::CHARLIE) {
			this->downstreamfuncs_marine.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::marine_charlie_III);
			if (SFD)
				this->downstreamfuncs_marine.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qshs_SFD);
			else
				this->downstreamfuncs_marine.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qshs_MFD);
		}

		if (this->secondary_fluvial == TSC_FLUVIAL::LATERALSPL)
			this->downstreamfuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::fluvial_lateral_erosion_SPL);

		if (this->fluvial == TSC_FLUVIAL::FASTSCAPE)
			this->upstreamfuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::fluvial_fastscape_SFD);
		else if (this->fluvial == TSC_FLUVIAL::DAVY2009)
			this->downstreamfuncs.emplace_back(
				&trackscape<fT, Graph_t, Connector_t>::fluvial_davy2009);

		if (this->secondary_fluvial == TSC_FLUVIAL::LATERALDAVY) {
			if (this->full_stochastic == false)
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::
						fluvial_lateral_erosion_deposition_davy);
			else
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::
						fluvial_lateral_erosion_deposition_davy_stochastic);
		}

		if (this->hillslopes != TSC_HILLSLOPE::NONE &&
				this->fluvial != TSC_FLUVIAL::NONE) {
			if (this->transfer_Qs_hs2fl)
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::capture_Qs_hs_from_fluv);
		}

		if (need_Qw && !need_Qs) {
			if (SFD)
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qw_SFD);
			else
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qw_MFD);
		} else if (need_Qw && need_Qs) {
			if (SFD)
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qw_Qs_SFD);
			else if (this->full_stochastic == false)
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qw_Qs_MFD);
			else
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qw_Qs_MFD_stochastic);
		} else if (need_Qs) {
			if (SFD)
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qshs_SFD);
			else if (this->full_stochastic == false)
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qshs_MFD);
			else
				this->downstreamfuncs.emplace_back(
					&trackscape<fT, Graph_t, Connector_t>::trans_Qshs_MFD_stochastic);
		}
	}

	// this function will get more sophisticated as I add processes. It returns
	// true if some processes needs to cache multiple neighbours.
	bool do_I_need_MFD()
	{
		if (this->flowtopo == TSC_FLOW_TOPOLOGY::MFD)
			return true;
		else
			return false;
	}

	// this function will get more sophisticated as I add processes. It returns
	// true if some processes needs to cache multiple neighbours.
	bool do_I_need_preprocess_Qw()
	{
		if (this->fluvial == TSC_FLUVIAL::FASTSCAPE)
			return true;
		else
			return false;
	}

	bool do_I_need_upstream_traversal()
	{
		if (this->fluvial == TSC_FLUVIAL::FASTSCAPE)
			return true;
		return false;
	}

	bool do_I_need_downstream_traversal()
	{
		if (this->hillslopes != TSC_HILLSLOPE::NONE ||
				this->fluvial == TSC_FLUVIAL::DAVY2009)
			return true;
		return false;
	}

	void reset_fhED()
	{
		this->thEs = 0.;
		this->thEr = 0.;
		this->thDs = 0.;
		this->tfEs = 0.;
		this->tfEr = 0.;
		this->tfDs = 0.;
	}

	/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																																													. - ~ ~ ~ - .
..     _      .-~               ~-.
//|     \ `..~                      `.
|| |      }  }              /       \  \
(\   \\ \~^..'                 |         }  \
\`.-~  o      /       }       |        /    \
(__          |       /        |       /      `.
`- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
													|     /          |     /     ~-.     ~- _
													|_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Process functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

	// no proc function, used when a process is deactivated
	void noproc() { return; }

	void fluvial_davy2009()
	{
		if (this->Qw[this->tnode] <= 0)
			return;

		// Erosion
		fT stream_power =
			std::pow(this->Qw[this->tnode], this->mexp) *
			std::pow(std::max(static_cast<fT>(1e-6), this->tSS), this->nexp);

		fT propused = 0;

		this->tfEs = Ks(this->tnode) * stream_power;
		fT remains = 0;

		if (this->add_to_dhsdt(this->tnode, -this->tfEs * this->tdt, remains)) {
			remains = std::abs(remains);
			propused = remains / (this->tfEs * this->tdt);
			this->tfEr = propused * Kr(this->tnode) * stream_power;

			if (tfEr < 0 || tfEs < 0)
				throw std::runtime_error("DMLKNFDFSLKJ");

			this->add_to_dbedrockdt(this->tnode, -this->tfEr * this->tdt);

			this->tfEs -= remains / this->tdt;
		}

		// this->tfDs = std::min(this->depcoeff(this->tnode) *
		// this->Qs_fluvial[this->tnode]/(this->Qw[this->tnode]),
		// this->Qs_fluvial[this->tnode]/(this->connector->get_area_at_node(this->tnode)));
		this->tfDs = this->Qs_fluvial[this->tnode] /
								 (this->depcoeff(this->tnode) * this->Qw[this->tnode]);
		// this->tfDs = 0.;

		if (this->connector->get_area_at_node(this->tnode) * this->tfDs >
				this->Qs_fluvial[this->tnode]) {
			this->tfDs = this->Qs_fluvial[this->tnode] /
									 this->connector->get_area_at_node(this->tnode);
			// std::cout << "JIJI" << std::endl;
		}

		this->add_to_dhsdt_nocheck(this->tnode, this->tfDs * this->tdt);

		this->Qs_fluvial[this->tnode] +=
			(this->tfEs + this->tfEr - this->tfDs) *
			this->connector->get_area_at_node(this->tnode);

		if (this->Qs_fluvial[this->tnode] < 0) {
			this->Qs_fluvial[this->tnode] = 0;
		}
		return;
	}

	void fluvial_fastscape_SFD()
	{
		fT tK = this->Kr(this->tnode);
		fT factor =
			tK * this->tdt * std::pow(this->Qw[this->tnode], this->mexp) /
			std::pow(this->connector->Sdistance2receivers[this->tnode], this->nexp);
		fT ielevation = this->z_surf[this->tnode];
		fT irec_elevation = this->z_surf[this->tSrec];
		fT elevation_k = ielevation;
		fT elevation_prev = std::numeric_limits<fT>::max();
		fT tolerance = 1e-6;

		while (abs(elevation_k - elevation_prev) > tolerance) {
			elevation_prev = elevation_k;
			fT slope = std::max(elevation_k - irec_elevation, static_cast<fT>(1e-6));
			fT diff =
				(elevation_k - ielevation + factor * std::pow(slope, this->nexp)) /
				(1. + factor * this->nexp * std::pow(slope, this->nexp - 1));
			elevation_k -= diff;
		}
		this->z_surf[this->tnode] = elevation_k;
	}

	void hillslopes_cidre()
	{
		// Storing the proportion of bedrock-power used
		fT propused = 0;

		// If the slope is bellow the critical values
		if (this->tSS <= Sc(this->tnode) - 1e-6) {
			// If I have sediments
			this->thEs = kappa_s(this->tnode) * this->tSS;
			fT remains = 0.;
			if (this->add_to_dhsdt(this->tnode, -this->tdt * this->thEs, remains)) {
				// erosion power
				remains = std::abs(remains);
				propused = remains / (this->thEs * this->tdt);
				this->thEs -= remains / this->tdt;
				// end of the check for excessing sed height
			} else {
				// I don't have any sed? -> no propused
				propused = 0.;
			}

			// Remaining applied to bedrock
			this->thEr = propused * kappa_r(this->tnode) * this->tSS;
			this->add_to_dbedrockdt(this->tnode, -this->thEr * this->tdt);

			fT L = this->connector->get_area_at_node(this->tnode) /
						 (1 - std::pow(this->tSS / Sc(this->tnode), 2));
			this->thDs = this->Qs_hs[this->tnode] / L;

			if (this->connector->get_area_at_node(this->tnode) * this->thDs >
					this->Qs_hs[this->tnode]) {
				this->thDs = this->Qs_hs[this->tnode] /
										 this->connector->get_area_at_node(this->tnode);
				// std::cout << "JIJI" << std::endl;
			}

			this->add_to_dhsdt_nocheck(this->tnode, this->thDs * this->tdt);

			if (L < 1)
				throw std::runtime_error("hillslopes::cidre::exception1994");

		} else {
			fT tothE = (this->z_surf[this->tnode] -
									(this->tdx * Sc(this->tnode) + z_surf[this->tSrec])) /
								 this->tdt;
			fT remains = 0;
			if (this->add_to_dhsdt(this->tnode, -1 * tothE * this->tdt, remains)) {
				remains = std::abs(remains);
				this->thEs = (tothE + remains) / this->tdt;
				this->thEr = remains / this->tdt;
				this->add_to_dbedrockdt(this->tnode, -this->thEr * this->tdt);
			} else
				this->thEs = tothE;
		}

		this->Qs_hs[this->tnode] += (this->thEs + this->thEr - this->thDs) *
																this->connector->get_area_at_node(this->tnode);

		if (this->Qs_hs[this->tnode] < 0)
			this->Qs_hs[this->tnode] = 0;
	}

	void hillslopes_cidre_no_crit()
	{
		// Storing the proportion of bedrock-power used
		fT propused = 0;

		fT min_dep = 0.01;

		fT ttss = std::max(this->tSS, static_cast<fT>(1e-2));

		// If I have sediments
		this->thEs = kappa_s(this->tnode) * ttss;
		fT remains = 0.;
		if (this->add_to_dhsdt(this->tnode, -this->tdt * this->thEs, remains)) {
			// erosion power
			remains = std::abs(remains);
			propused = remains / (this->thEs * this->tdt);
			this->thEs -= remains / this->tdt;
			// end of the check for excessing sed height
		} else {
			// I don't have any sed? -> no propused
			propused = 0.;
		}

		// Remaining applied to bedrock
		this->thEr = propused * kappa_r(this->tnode) * ttss;
		this->add_to_dbedrockdt(this->tnode, -this->thEr * this->tdt);

		// If the slope is bellow the critical values
		if (ttss <= Sc(this->tnode) - 1e-6) {
			fT L = this->connector->get_area_at_node(this->tnode) /
						 (1 - std::pow(ttss / Sc(this->tnode), 2));
			this->thDs = std::max(this->Qs_hs[this->tnode] / L,
														min_dep * this->Qs_hs[this->tnode] /
															this->connector->get_area_at_node(this->tnode));

			if (this->connector->get_area_at_node(this->tnode) * this->thDs >
					this->Qs_hs[this->tnode]) {
				this->thDs = this->Qs_hs[this->tnode] /
										 this->connector->get_area_at_node(this->tnode);
				// std::cout << "JIJI" << std::endl;
			}

			this->add_to_dhsdt_nocheck(this->tnode, this->thDs * this->tdt);

		} else
			this->thDs = min_dep * this->Qs_hs[this->tnode] /
									 this->connector->get_area_at_node(this->tnode);

		this->Qs_hs[this->tnode] += (this->thEs + this->thEr - this->thDs) *
																this->connector->get_area_at_node(this->tnode);

		if (this->Qs_hs[this->tnode] < 0)
			this->Qs_hs[this->tnode] = 0;
	}

	void hillslopes_cidre_dep_only_for_highlands()
	{

		if (this->Qs_hs[this->tnode] <= 0)
			return;
		// If the slope is bellow the critical values

		if (this->tSS <= internal_friction(this->tnode) -
											 1e-6) // using the internal friction as this section is
														 // used in hylands
		{

			fT L = this->connector->get_area_at_node(this->tnode) /
						 (1 - std::pow(this->tSS / internal_friction(this->tnode), 2));
			// this->thDs =
			// std::min(this->max_dep_hillslopes/this->connector->get_area_at_node(this->tnode),
			// std::max(1/L,
			// this->min_dep_hillslopes/this->connector->get_area_at_node(this->tnode))
			// )* this->Qs_hs[this->tnode];
			this->thDs = 1 / L * this->Qs_hs[this->tnode];
		}
		// else
		// 	this->thDs =
		// std::min(this->max_dep_hillslopes/this->connector->get_area_at_node(this->tnode),
		// this->min_dep_hillslopes/this->connector->get_area_at_node(this->tnode))
		// * this->Qs_hs[this->tnode];

		if (this->connector->get_area_at_node(this->tnode) * this->thDs >
				this->Qs_hs[this->tnode]) {
			this->thDs = this->Qs_hs[this->tnode] /
									 this->connector->get_area_at_node(this->tnode);
		}

		this->add_to_dhsdt_nocheck(this->tnode, this->thDs * this->tdt);

		this->Qs_hs[this->tnode] -=
			(this->thDs) * this->connector->get_area_at_node(this->tnode);

		if (this->Qs_hs[this->tnode] < 0)
			this->Qs_hs[this->tnode] = 0;
	}

	void hillslopes_cidre_dep_only()
	{

		if (this->Qs_hs[this->tnode] <= 0)
			return;
		// If the slope is bellow the critical values

		if (this->tSS <= Sc(this->tnode) - 1e-6) // using the internal friction as
																						 // this section is used in hylands
		{

			fT L = this->connector->get_area_at_node(this->tnode) /
						 (1 - std::pow(this->tSS / Sc(this->tnode), 2));
			// this->thDs =
			// std::min(this->max_dep_hillslopes/this->connector->get_area_at_node(this->tnode),
			// std::max(1/L,
			// this->min_dep_hillslopes/this->connector->get_area_at_node(this->tnode))
			// )* this->Qs_hs[this->tnode];
			this->thDs = 1 / L * this->Qs_hs[this->tnode];
		}
		// else
		//  this->thDs =
		// std::min(this->max_dep_hillslopes/this->connector->get_area_at_node(this->tnode),
		// this->min_dep_hillslopes/this->connector->get_area_at_node(this->tnode))
		// * this->Qs_hs[this->tnode];

		if (this->connector->get_area_at_node(this->tnode) * this->thDs >
				this->Qs_hs[this->tnode]) {
			this->thDs = this->Qs_hs[this->tnode] /
									 this->connector->get_area_at_node(this->tnode);
		}

		this->add_to_dhsdt_nocheck(this->tnode, this->thDs * this->tdt);

		this->Qs_hs[this->tnode] -=
			(this->thDs) * this->connector->get_area_at_node(this->tnode);

		if (this->Qs_hs[this->tnode] < 0)
			this->Qs_hs[this->tnode] = 0;
	}

	void hillslopes_cidre_longdep()
	{
		// Storing the proportion of bedrock-power used
		fT propused = 0;

		// If the slope is bellow the critical values
		if (this->tSS <= Sc(this->tnode) - 1e-6) {
			// If I have sediments
			this->thEs = kappa_s(this->tnode) * this->tSS;
			fT remains = 0.;
			if (this->add_to_dhsdt(this->tnode, -this->tdt * this->thEs, remains)) {
				// erosion power
				remains = std::abs(remains);
				propused = remains / (this->thEs * this->tdt);
				this->thEs -= remains / this->tdt;
				// end of the check for excessing sed height
			} else {
				// I don't have any sed? -> no propused
				propused = 0.;
			}

			// Remaining applied to bedrock
			this->thEr = propused * kappa_r(this->tnode) * this->tSS;
			this->add_to_dbedrockdt(this->tnode, -this->thEr * this->tdt);

			fT L = this->connector->get_area_at_node(this->tnode) /
						 (1 - std::pow(this->tSS / Sc(this->tnode), 2));
			this->thDs = this->Qs_hs[this->tnode] / (L * this->longdep(this->tnode));

			if (this->connector->get_area_at_node(this->tnode) * this->thDs >
					this->Qs_hs[this->tnode]) {
				this->thDs = this->Qs_hs[this->tnode] /
										 this->connector->get_area_at_node(this->tnode);
				// std::cout << "JIJI" << std::endl;
			}

			this->add_to_dhsdt_nocheck(this->tnode, this->thDs * this->tdt);

			if (L < 1)
				throw std::runtime_error("hillslopes::cidre::exception1994");

		} else {
			fT tothE = (this->z_surf[this->tnode] -
									(this->tdx * Sc(this->tnode) + z_surf[this->tSrec])) /
								 this->tdt;
			fT remains = 0;
			if (this->add_to_dhsdt(this->tnode, -1 * tothE * this->tdt, remains)) {
				remains = std::abs(remains);
				this->thEs = (tothE + remains) / this->tdt;
				this->thEr = remains / this->tdt;
				this->add_to_dbedrockdt(this->tnode, -this->thEr * this->tdt);
			} else
				this->thEs = tothE;
		}

		this->Qs_hs[this->tnode] += (this->thEs + this->thEr - this->thDs) *
																this->connector->get_area_at_node(this->tnode);

		if (this->Qs_hs[this->tnode] < 0)
			this->Qs_hs[this->tnode] = 0;
	}

	void hillslopes_cidre_longdep_dep_only()
	{

		if (this->Qs_hs[this->tnode] <= 0)
			return;
		// If the slope is bellow the critical values

		if (this->tSS <= Sc(this->tnode) - 1e-6) // using the internal friction as
																						 // this section is used in hylands
		{

			fT L = this->connector->get_area_at_node(this->tnode) /
						 (1 - std::pow(this->tSS / Sc(this->tnode), 2));
			// this->thDs =
			// std::min(this->max_dep_hillslopes/this->connector->get_area_at_node(this->tnode),
			// std::max(1/L,
			// this->min_dep_hillslopes/this->connector->get_area_at_node(this->tnode))
			// )* this->Qs_hs[this->tnode];
			this->thDs = this->Qs_hs[this->tnode] / (L * this->longdep(this->tnode));
		}
		// else
		//  this->thDs =
		// std::min(this->max_dep_hillslopes/this->connector->get_area_at_node(this->tnode),
		// this->min_dep_hillslopes/this->connector->get_area_at_node(this->tnode))
		// * this->Qs_hs[this->tnode];

		if (this->connector->get_area_at_node(this->tnode) * this->thDs >
				this->Qs_hs[this->tnode]) {
			this->thDs = this->Qs_hs[this->tnode] /
									 this->connector->get_area_at_node(this->tnode);
		}

		this->add_to_dhsdt_nocheck(this->tnode, this->thDs * this->tdt);

		this->Qs_hs[this->tnode] -=
			(this->thDs) * this->connector->get_area_at_node(this->tnode);

		if (this->Qs_hs[this->tnode] < 0)
			this->Qs_hs[this->tnode] = 0;
	}

	void capture_Qs_hs_from_fluv()
	{
		fT tqc = this->Qs_hs[this->tnode] * this->transfer_rate_Qs_hs2fl;
		this->Qs_fluvial[this->tnode] += tqc;
		this->Qs_hs[this->tnode] -= tqc;
	}

	void fluvial_lateral_erosion_deposition_davy()
	{
		this->_fluvial_lateral_erosion_deposition_davy(this->orthonodA);
		this->_fluvial_lateral_erosion_deposition_davy(this->orthonodB);
	}

	void _fluvial_lateral_erosion_deposition_davy(int onode)
	{
		if (onode <= -1 || onode >= this->connector->nnodes)
			return;

		fT dZ = this->z_surf[this->tnode] - this->z_surf[onode];

		// case deposition
		if (dZ > 0) {
			fT tlatdep =
				dZ / this->tdy * this->Kld(onode) * this->tfDs; // * this->tdt;
			if (tlatdep * this->connector->get_area_at_node(this->tnode) >
					this->Qs_fluvial[this->tnode])
				tlatdep = this->Qs_fluvial[this->tnode] /
									this->connector->get_area_at_node(this->tnode);

			this->add_to_dhsdt_nocheck(onode, tlatdep * this->tdt);

			this->Qs_fluvial[this->tnode] -=
				tlatdep * this->connector->get_area_at_node(this->tnode);

			if (this->Qs_fluvial[this->tnode] < 0)
				this->Qs_fluvial[this->tnode] = 0;

		} else if (dZ < 0) {
			fT remains = 0;
			fT tlater =
				std::abs(dZ) / this->tdy * this->Kle(onode) * (this->tfEs + this->tfEr);
			this->add_to_dhsdt(onode,
												 -1 * tlater * this->tdt,
												 remains); // the -1 for the vertical motion is included
																	 // in the negative remain
			// the remains is already negative if not null
			this->add_to_dbedrockdt(onode, remains);

			this->Qs_fluvial[this->tnode] +=
				tlater * this->connector->get_area_at_node(this->tnode);
			if (this->Qs_fluvial[this->tnode] < 0)
				this->Qs_fluvial[this->tnode] = 0;
		}
	}

	void fluvial_lateral_erosion_deposition_davy_stochastic()
	{
		this->_fluvial_lateral_erosion_deposition_davy_stochastic(this->orthonodA);
		this->_fluvial_lateral_erosion_deposition_davy_stochastic(this->orthonodB);
	}

	void _fluvial_lateral_erosion_deposition_davy_stochastic(int onode)
	{
		if (onode <= -1 || onode >= this->connector->nnodes)
			return;

		fT dZ = this->z_surf[this->tnode] - this->z_surf[onode];

		// case deposition
		if (dZ > 0) {
			fT tlatdep = dZ / this->tdy * this->Kld(this->tnode) * this->tfDs *
									 (this->connector->randu->get() + 0.5); // * this->tdt;
			if (tlatdep * this->connector->get_area_at_node(this->tnode) >
					this->Qs_fluvial[this->tnode])
				tlatdep = this->Qs_fluvial[this->tnode] /
									this->connector->get_area_at_node(this->tnode);

			this->add_to_dhsdt_nocheck(onode, tlatdep * this->tdt);

			this->Qs_fluvial[this->tnode] -=
				tlatdep * this->connector->get_area_at_node(this->tnode);

			if (this->Qs_fluvial[this->tnode] < 0)
				this->Qs_fluvial[this->tnode] = 0;
		} else if (dZ < 0) {
			fT remains = 0;
			fT tlater = std::abs(dZ) / this->tdy * this->Kle(this->tnode) *
									(this->tfEs + this->tfEr) *
									(this->connector->randu->get() + 0.5);
			this->add_to_dhsdt(onode,
												 -1 * tlater * this->tdt,
												 remains); // the -1 for the vertical motion is included
																	 // in the negative remain
			// the remains is already negative if not null
			this->add_to_dbedrockdt(onode, remains);

			this->Qs_fluvial[this->tnode] +=
				tlater * this->connector->get_area_at_node(this->tnode);
			if (this->Qs_fluvial[this->tnode] < 0)
				this->Qs_fluvial[this->tnode] = 0;
		}
	}

	void fluvial_lateral_erosion_SPL()
	{
		this->_fluvial_lateral_erosion_SPL(this->orthonodA);
		this->_fluvial_lateral_erosion_SPL(this->orthonodB);
	}

	void _fluvial_lateral_erosion_SPL(int onode)
	{
		if (onode <= -1 || onode >= this->connector->nnodes)
			return;

		fT DeltaA = this->Qw[onode] - this->Qw[this->tnode];

		if (DeltaA < 0)
			return;

		fT dZ = this->z_surf[this->tnode] - this->z_surf[onode];

		if (dZ < 0) {
			fT remains = 0;
			fT tlater = std::pow(std::abs(dZ) / this->tdy, this->nexp) *
									std::pow(DeltaA, this->mexp) * this->Kle(this->tnode);

			this->add_to_dhsdt(onode,
												 -1 * tlater * this->tdt,
												 remains); // the -1 for the vertical motion is included
																	 // in the negative remain
			// the remains is already negative if not null
			this->add_to_dbedrockdt(onode, remains);

			this->Qs_fluvial[this->tnode] +=
				tlater * this->connector->get_area_at_node(this->tnode);
			if (this->Qs_fluvial[this->tnode] < 0)
				this->Qs_fluvial[this->tnode] = 0;
		}
	}

	void prec_extra_Qs_fluvial()
	{
		for (size_t i = 0; i < this->extra_Qsf.size(); ++i) {
			this->Qs_fluvial[this->extra_Qsf_nodes[i]] += this->extra_Qsf[i];
		}
	}

	void prec_extra_Qw_fluvial()
	{
		for (size_t i = 0; i < this->extra_Qsf.size(); ++i) {
			this->Qw[this->extra_Qwf_nodes[i]] += this->extra_Qwf[i];
		}
	}

	void prec_Qw_SFD()
	{
		// std::cout << "prec_Qw_SFD" << std::endl;
		this->Qw[this->tnode] += this->precipitations(this->tnode) *
														 this->connector->get_area_at_node(this->tnode);
	}

	void trans_Qw_SFD()
	{
		// std::cout << "trans_Qw_SFD" << std::endl;
		this->Qw[this->tSrec] += this->Qw[this->tnode];
	}

	void prec_Qw_MFD()
	{
		// std::cout << "prec_Qw_MFD" << std::endl;
		this->Qw[this->tnode] += this->precipitations(this->tnode) *
														 this->connector->get_area_at_node(this->tnode);
	}

	void trans_Qw_MFD()
	{
		// std::cout << "trans_Qw_MFD" << std::endl;
		// this->tslopes
		fT sumsum = 0.;
		for (int j = 0; j < this->trn; ++j) {
			this->tslopes[j] = std::max(
				static_cast<fT>(1e-6),
				(this->tZ - this->z_surf[this->treceivers_nodes[j]]) /
					this->connector->get_dx_from_links_idx(this->treceivers_links[j]));

			sumsum += this->tslopes[j];
		}

		for (int j = 0; j < this->trn; ++j) {
			if (sumsum > 0)
				this->Qw[this->treceivers_nodes[j]] +=
					this->Qw[this->tnode] * this->tslopes[j] / sumsum;
			else
				this->Qw[this->treceivers_nodes[j]] +=
					this->Qw[this->tnode] * 1. / this->trn;
		}
	}

	void trans_Qw_Qs_SFD()
	{
		// std::cout << "trans_Qw_Qs_SFD" << std::endl;
		this->Qw[this->tSrec] += this->Qw[this->tnode];
		this->Qs_fluvial[this->tSrec] += this->Qs_fluvial[this->tnode];
		if (this->hillslopes != TSC_HILLSLOPE::NONE)
			this->Qs_hs[this->tSrec] += this->Qs_hs[this->tnode];
	}

	void trans_Qshs_SFD()
	{
		// std::cout << "trans_Qs_SFD" << std::endl;
		this->Qs_hs[this->tSrec] += this->Qs_hs[this->tnode];
	}

	void trans_Qshs_MFD()
	{
		// std::cout << "trans_Qs_MFD" << std::endl;
		// this->tslopes
		fT sumsum = 0.;
		for (int j = 0; j < this->trn; ++j) {
			this->tslopes[j] = std::max(
				static_cast<fT>(1e-6),
				(this->tZ - this->z_surf[this->treceivers_nodes[j]]) /
					this->connector->get_dx_from_links_idx(this->treceivers_links[j]));

			sumsum += this->tslopes[j];
		}

		for (int j = 0; j < this->trn; ++j) {
			if (sumsum > 0)
				this->Qs_hs[this->treceivers_nodes[j]] +=
					this->Qs_hs[this->tnode] * this->tslopes[j] / sumsum;
			else
				this->Qs_hs[this->treceivers_nodes[j]] +=
					this->Qs_hs[this->tnode] * 1. / this->trn;
		}
	}

	void trans_Qw_Qs_MFD()
	{
		// this->tslopes
		// std::cout << "trans_Qw_Qs_MFD" << std::endl;
		fT sumsum = 0.;
		for (int j = 0; j < this->trn; ++j) {
			this->tslopes[j] = std::max(
				static_cast<fT>(1e-6),
				(this->tZ - this->z_surf[this->treceivers_nodes[j]]) /
					this->connector->get_dx_from_links_idx(this->treceivers_links[j]));

			sumsum += this->tslopes[j];
		}

		for (int j = 0; j < this->trn; ++j) {
			if (sumsum > 0) {
				this->Qw[this->treceivers_nodes[j]] +=
					this->Qw[this->tnode] * this->tslopes[j] / sumsum;
				this->Qs_fluvial[this->treceivers_nodes[j]] +=
					this->Qs_fluvial[this->tnode] * this->tslopes[j] / sumsum;
				if (this->hillslopes != TSC_HILLSLOPE::NONE)
					this->Qs_hs[this->treceivers_nodes[j]] +=
						this->Qs_hs[this->tnode] * this->tslopes[j] / sumsum;
			} else {
				this->Qw[this->treceivers_nodes[j]] +=
					this->Qw[this->tnode] * 1. / this->trn;
				this->Qs_fluvial[this->treceivers_nodes[j]] +=
					this->Qs_fluvial[this->tnode] * 1. / this->trn;
				if (this->hillslopes != TSC_HILLSLOPE::NONE)
					this->Qs_hs[this->treceivers_nodes[j]] +=
						this->Qs_hs[this->tnode] * 1. / this->trn;
			}
		}
	}

	void trans_Qshs_MFD_stochastic()
	{
		// std::cout << "trans_Qs_MFD" << std::endl;
		// this->tslopes
		fT sumsum = 0.;
		for (int j = 0; j < this->trn; ++j) {
			this->tslopes[j] =
				std::max(
					static_cast<fT>(1e-6),
					(this->tZ - this->z_surf[this->treceivers_nodes[j]]) /
						this->connector->get_dx_from_links_idx(this->treceivers_links[j])) *
				this->connector->randu->get();

			sumsum += this->tslopes[j];
		}

		for (int j = 0; j < this->trn; ++j) {
			if (sumsum > 0)
				this->Qs_hs[this->treceivers_nodes[j]] +=
					this->Qs_hs[this->tnode] * this->tslopes[j] / sumsum;
			else
				this->Qs_hs[this->treceivers_nodes[j]] +=
					this->Qs_hs[this->tnode] * 1. / this->trn;
		}
	}

	void trans_Qw_Qs_MFD_stochastic()
	{
		// this->tslopes
		// std::cout << "trans_Qw_Qs_MFD" << std::endl;
		fT sumsum = 0.;
		for (int j = 0; j < this->trn; ++j) {
			this->tslopes[j] =
				std::max(
					static_cast<fT>(1e-6),
					(this->tZ - this->z_surf[this->treceivers_nodes[j]]) /
						this->connector->get_dx_from_links_idx(this->treceivers_links[j])) *
				this->connector->randu->get();

			sumsum += this->tslopes[j];
		}

		for (int j = 0; j < this->trn; ++j) {
			if (sumsum > 0) {
				this->Qw[this->treceivers_nodes[j]] +=
					this->Qw[this->tnode] * this->tslopes[j] / sumsum;
				this->Qs_fluvial[this->treceivers_nodes[j]] +=
					this->Qs_fluvial[this->tnode] * this->tslopes[j] / sumsum;
				if (this->hillslopes != TSC_HILLSLOPE::NONE)
					this->Qs_hs[this->treceivers_nodes[j]] +=
						this->Qs_hs[this->tnode] * this->tslopes[j] / sumsum;
			} else {
				this->Qw[this->treceivers_nodes[j]] +=
					this->Qw[this->tnode] * 1. / this->trn;
				this->Qs_fluvial[this->treceivers_nodes[j]] +=
					this->Qs_fluvial[this->tnode] * 1. / this->trn;
				if (this->hillslopes != TSC_HILLSLOPE::NONE)
					this->Qs_hs[this->treceivers_nodes[j]] +=
						this->Qs_hs[this->tnode] * 1. / this->trn;
			}
		}
	}

	void marine_charlie_III()
	{

		// Storing the proportion of bedrock-power used
		// fT propused = 0;

		// updating the node sediment flux with anything that came (river +
		// hillslope)
		this->Qs_hs[this->tnode] += this->Qs_fluvial[this->tnode];

		// Storing the proportion of bedrock-power used
		fT propused = 0;

		// If the slope is bellow the critical values
		if (this->tSS <= Sc_M(this->tnode) - 1e-6) {
			// If I have sediments
			this->thEs = this->Ke(this->tnode) * this->tSS;
			fT remains = 0.;
			if (this->add_to_dhsdt(this->tnode, -this->tdt * this->thEs, remains)) {
				// erosion power
				remains = std::abs(remains);
				propused = remains / (this->thEs * this->tdt);
				this->thEs -= remains / this->tdt;
				// end of the check for excessing sed height
			} else {
				// I don't have any sed? -> no propused
				propused = 0.;
			}

			// Remaining applied to bedrock
			// this->thEr = propused * kappa_r(this->tnode) * this->tSS;
			// this->add_to_dbedrockdt(this->tnode, - this->thEr * this->tdt);

			fT L = this->connector->get_area_at_node(this->tnode) *
						 this->lambda(this->tnode) /
						 (1 - std::pow(this->tSS / Sc_M(this->tnode), 2));
			this->thDs = this->Qs_hs[this->tnode] / L;

			if (this->connector->get_area_at_node(this->tnode) * this->thDs >
					this->Qs_hs[this->tnode]) {
				this->thDs = this->Qs_hs[this->tnode] /
										 this->connector->get_area_at_node(this->tnode);
				// std::cout << "JIJI" << std::endl;
			}

			if (this->z_surf[this->tnode] + this->thDs * this->tdt >
					this->sea_level(this->tnode))
				this->thDs = (this->sea_level(this->tnode) + 1e-3 * this->tdt -
											this->z_surf[this->tnode]) /
										 this->tdt;

			this->add_to_dhsdt_nocheck(this->tnode, this->thDs * this->tdt);

			if (L < 1)
				throw std::runtime_error("hillslopes::cidre::exception1994");

		} else {
			fT tothE = (this->z_surf[this->tnode] -
									(this->tdx * Sc_M(this->tnode) + z_surf[this->tSrec])) /
								 this->tdt;
			fT remains = 0;
			if (this->add_to_dhsdt(this->tnode, -1 * tothE * this->tdt, remains)) {
				remains = std::abs(remains);
				this->thEs = (tothE + remains) / this->tdt;
				this->thEr = remains / this->tdt;
				this->add_to_dbedrockdt(this->tnode, -this->thEr * this->tdt);
			} else
				this->thEs = tothE;
		}

		this->Qs_hs[this->tnode] += (this->thEs + this->thEr - this->thDs) *
																this->connector->get_area_at_node(this->tnode);

		if (this->Qs_hs[this->tnode] < 0)
			this->Qs_hs[this->tnode] = 0;

		// lasjdfljasdl;fja;slkjdfka;sjdlfjasldkjf

		// // If the slope is bellow the critical values
		// if(S <= Sc_M(this->tnode) - 1e-6)
		// {
		// 	// If I have sediments
		// 	if(this->h_sed[this->tnode] > 0)
		// 	{

		// 		// erosion power
		// 		mEs = Ke(this->tnode) * S;

		// 		// Checking I am not strapping more than the sediment pile
		// 		if(mEs * dt > this->h_sed[this->tnode])
		// 		{
		// 			// Correcting to the max sed height removal possible and
		// calculating the proportion of sediment used
		// 			// propused = this->h_sed[this->tnode]/dt /mEs;
		// 			mEs = this->h_sed[this->tnode]/dt;
		// 		}
		// 		else
		// 		{
		// 			// I don't have any sed? -> no propused
		// 			// propused = 1.;
		// 		}

		// 	// end of the check for excessing sed height
		// 	}

		// 	// Remaining applied to bedrock
		// 	// mEr = (1. - propused) * Kr(this->tnode) * S;

		// 	fT L = (this->connector->get_travers_dy_from_dx(dx) *
		// this->lambda(this->tnode))/(1 - std::pow(S/Sc_M(this->tnode),2));
		// mDs = this->Qs_hs[this->tnode]/std::max(L,cellarea);

		// }
		// else
		// {
		// 	fT tothE = (z_surf[this->tnode] - (dx * Sc(this->tnode) + z_surf[rec])
		// )/dt; 	if(tothE > this->h_sed[this->tnode])
		// 	{
		// 		mEs = this->h_sed[this->tnode]/dt;
		// 		// mEr = tothE - mEs;
		// 	}
		// 	else
		// 		mEs = tothE;
		// }

		// if(this->TSP_module)
		// 	this->apply_TSP(this->tnode, rec, mEs, mEr, mDs, dt, true);

		// if(this->Ch_MTSI)
		// 	this->apply_Ch_MTSI_SFD(this->tnode, rec, mEs, mEr, mDs, dt, true);

		// this->Qs_hs[this->tnode] += (mEs + mEr - mDs) * cellarea;

		// // Applying the fluxes
		// this->h_sed[this->tnode] += (mDs - mEs) * dt;
		// this->z_surf[this->tnode] += (mDs - mEs - mEr) * dt;
		// if(this->Qs_hs[this->tnode] < 0) this->Qs_hs[this->tnode] = 0;
		// this->Qs_hs[rec] += this->Qs_hs[this->tnode];
	}

	void apply_dzhsdt()
	{
		for (int i = 0; i < this->connector->nnodes; ++i) {
			this->z_surf[i] += this->dhsdt[i];
			this->h_sed[i] += this->dhsdt[i];
			this->z_surf[i] += this->dbedrockdt[i];
		}
	}

	/*
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																																																													. - ~ ~ ~ - .
													..     _      .-~               ~-.
									 //|     \ `..~                      `.
									|| |      }  }              /       \  \
	(\   \\ \~^..'                 |         }  \
	 \`.-~  o      /       }       |        /    \
	 (__          |       /        |       /      `.
					`- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
																													|     /          | /
	~-.     ~- _
																													|_____| |_____| ~ - .
	_ _~_-_

																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	prefuncs functions
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	void hillslopes_hylands_trigger()
	{
		if (this->flowtopo != TSC_FLOW_TOPOLOGY::MFD)
			throw std::runtime_error(
				"Not compatible SFD at the moment (WIP, there is a fatal bug)");

		// reinit the landslides ID
		// this->landslidesid = std::vector<std::uint32_t>(this->connector->nnodes,
		// 0); std::vector<std::uint8> islandslide(this->connector->nnodes, 0);
		// std::uint32_t id_ldsl = 1;

		// Readying stack helper
		std::stack<int, std::vector<int>> stackhelper;

		auto neighbours_nodes = this->connector->get_empty_neighbour();
		auto neighbours_links = this->connector->get_empty_neighbour();
		for (int i = 0; i < this->connector->nnodes; ++i) {

			int node = this->get_istack_node_MFD(i);
			if (this->connector->flow_out_or_pit(node) ||
					this->connector->boundaries.no_data(node))
				continue;
			// if(islandslide[node]) continue;

			int nn = this->connector->get_donors_idx_nodes_and_links(
				node, neighbours_nodes, neighbours_links);

			bool crit_exceeded = false;
			fT Hs = 0.;
			fT tslope = 0.;
			for (int j = 0; j < nn; ++j) {
				fT tths = (this->z_surf[neighbours_nodes[j]] - this->z_surf[node]);
				fT ttslope =
					tths / this->connector->get_dx_from_links_idx(neighbours_links[j]);
				if (ttslope > this->internal_friction(node)) {
					crit_exceeded = true;
					if (ttslope > tslope)
						tslope = ttslope;
				}
				if (tths > Hs)
					Hs = tths;
			}

			if (!crit_exceeded)
				continue;

			fT tslope_d = std::atan(tslope); // * (180.0/3.141592653589793238463);

			fT Hc = 4 * this->C(node) / (this->gravity * this->rho(node)) *
							(std::sin(tslope_d) * std::cos(this->internal_friction_d(node))) /
							(1 - std::cos(tslope_d - this->internal_friction_d(node)));
			fT proba = this->connector->randu->get() * this->tls(node) / this->tdt;

			// std::cout <<  Hs/Hc << "|" << proba << std::endl;;

			if (proba < Hs / Hc) {
				this->_hylands_trigger_single_landslide(
					node, tslope, stackhelper, neighbours_nodes, neighbours_links);
			}
		}
	}

	void _hylands_trigger_single_landslide(
		int node,
		fT tslope,
		std::stack<int, std::vector<int>>& stackhelper,
		std::array<int, 8>& neighbours_nodes,
		std::array<int, 8>& neighbours_links)
	{
		fT ldsl_angle = (this->internal_friction(node) + tslope) / 2;
		stackhelper.emplace(node);
		int onode = node;

		while (stackhelper.empty() == false) {
			node = stackhelper.top();
			stackhelper.pop();
			int nn = this->connector->get_neighbour_idx_nodes_and_links(
				node, neighbours_nodes, neighbours_links);
			for (int j = 0; j < nn; ++j) {
				int donode = neighbours_nodes[j];
				fT ddx = this->connector->get_dx_from_links_idx(neighbours_links[j]);
				fT target_z = this->z_surf[node] + ldsl_angle * ddx;
				if (target_z < this->z_surf[donode]) {
					fT tdZ = (this->z_surf[donode] - target_z);
					this->Qs_hs[onode] +=
						tdZ * this->connector->get_area_at_node(donode) / this->tdt;
					this->remove_shear_height_at_once(donode, tdZ);
					stackhelper.emplace(donode);
				}
			}
		}
	}

	void fix_small_marine_patches()
	{
		return;

		std::vector<std::uint8_t> needs(this->connector->nnodes, false);
		auto neighbours = this->connector->get_empty_neighbour();
		for (int i = 0; i < this->connector->nnodes; ++i) {
			if (this->z_surf[i] < this->sea_level(i)) {
				int nsea = 0;
				int nn = this->connector->get_neighbour_idx(i, neighbours);
				if (nn == 0)
					continue;
				for (int j = 0; j < nn; ++j) {
					int tn = neighbours[j];
					if (this->z_surf[tn] < this->sea_level(i))
						++nsea;
				}

				if (nsea / nn < 0.4)
					needs[i] = true;
			}
		}
		for (int i = 0; i < this->connector->nnodes; ++i) {
			if (needs[i])
				this->z_surf[i] =
					this->sea_level(i) + this->connector->randu->get() * 1e-5;
		}
	}

	// Some processes need Qw to be processed beforehand
	void preprocess_Qw_SFD()
	{
		if (this->variable_precipitations) {
			std::vector<fT> tQA(this->connector->nnodes, 0.);

			for (int i = 0; i < this->connector->nnodes; ++i) {
				tQA[i] =
					this->_precipitations[i] * this->connector->get_area_at_node(i);
			}
			this->Qw = this->graph->_accumulate_variable_downstream_SFD(tQA);
		} else {
			this->Qw = this->graph->_accumulate_constant_downstream_SFD(
				this->connector->get_area_at_node(0) * this->precipitations(0));
		}
	}

	/*
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																																																													. - ~ ~ ~ - .
													..     _      .-~               ~-.
									 //|     \ `..~                      `.
									|| |      }  }              /       \  \
	(\   \\ \~^..'                 |         }  \
	 \`.-~  o      /       }       |        /    \
	 (__          |       /        |       /      `.
					`- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
																													|     /          | /
	~-.     ~- _
																													|_____| |_____| ~ - .
	_ _~_-_

																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	Helper functions
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	// add to the vertical changes in sediment height and return true if something
	// remains
	bool add_to_dhsdt(int i, fT ths, fT& remain)
	{
		fT remaining_sed_height = this->h_sed[i] + this->dhsdt[i];

		// checking if the changes over the sediment height is more than it can
		// ingest
		if (remaining_sed_height + ths > 0) {
			this->dhsdt[i] += ths;
			return false;
		}

		// then I need to calculate what remains
		// which should be my vertical change + my actual sed height
		remain = ths + remaining_sed_height;
		this->dhsdt[i] = -this->h_sed[i];
		return true;
	}

	// add to the vertical changes in sediment height and return
	void add_to_dhsdt_nocheck(int i, fT ths)
	{
		// first adding the changes in sediment heights
		this->dhsdt[i] += ths;
	}

	// add to the vertical changes in sediment height and return
	void add_to_dbedrockdt(int i, fT tdbedrock)
	{
		// first adding the changes in sediment heights
		this->dbedrockdt[i] += tdbedrock;
	}

	void remove_shear_height_at_once(int i, fT dZ)
	{
		this->h_sed[i] = std::max(static_cast<fT>(0.), this->h_sed[i] - dZ);
		this->z_surf[i] -= dZ;
	}

	void _ready_node_state()
	{
		this->tSrec = this->connector->_Sreceivers[this->tnode];

		// feeding the local private variables
		this->tdx = this->connector->Sdistance2receivers[this->tnode];
		this->tdy = this->connector->get_travers_dy_from_dx(this->tdx);
		this->tZ = this->z_surf[this->tnode];
		// this->tSS = this->connector->SS[this->tnode];
		if (this->marine == TSC_MARINE::NONE)
			this->tSS = std::max(static_cast<fT>(1e-6),
													 (tZ - this->z_surf[this->tSrec]) / this->tdx);
		else if (this->z_surf[this->tnode] > this->sea_level(this->tnode) &&
						 this->z_surf[this->tSrec] > this->sea_level(this->tSrec))
			this->tSS = std::max(static_cast<fT>(1e-6),
													 (tZ - this->z_surf[this->tSrec]) / this->tdx);
		else if (this->z_surf[this->tnode] < this->sea_level(this->tnode) &&
						 this->z_surf[this->tSrec] < this->sea_level(this->tSrec))
			this->tSS = std::max(static_cast<fT>(1e-6),
													 (tZ - this->z_surf[this->tSrec]) / this->tdx);
		else if (this->z_surf[this->tnode] > this->sea_level(this->tnode) &&
						 this->z_surf[this->tSrec] < this->sea_level(this->tSrec))
			this->tSS = std::max(static_cast<fT>(1e-6),
													 (tZ - this->sea_level(this->tnode)) / this->tdx);
		else
			this->tSS = static_cast<fT>(1e-6);

		this->ths = this->h_sed[this->tnode];

		if (this->secondary_fluvial != TSC_FLUVIAL::NONE) {
			auto onodes =
				this->connector->get_orthogonal_nodes(this->tnode, this->tSrec);
			this->orthonodA = onodes.first;
			this->orthonodB = onodes.second;
		}

		if (this->flowtopo == TSC_FLOW_TOPOLOGY::MFD)
			this->trn = this->connector->get_receivers_idx_nodes_and_links(
				this->tnode, this->treceivers_nodes, this->treceivers_links);

		// this->reset_fhED();
	}

	void _exp_ready_node_state()
	{

		if (this->flowtopo == TSC_FLOW_TOPOLOGY::MFD)
			this->trn = this->connector->get_receivers_idx_nodes_and_links(
				this->tnode, this->treceivers_nodes, this->treceivers_links);

		this->tZ = this->z_surf[this->tnode];

		int tempSrec = this->connector->_Sreceivers[this->tnode];
		fT tempS = 0.;

		for (int j = 0; j < this->trn; ++j) {
			fT tSt =
				this->connector->randu
					->get(); // * std::max(static_cast<fT>(1e-6), (tZ -
									 // this->z_surf[this->treceivers_nodes[j]]) /
									 // this->connector->get_dx_from_links_idx(this->treceivers_links[j])
									 // );
			if (tSt > tempS) {
				tempS = tSt;
				tempSrec = this->treceivers_nodes[j];
			}
		}

		this->tSrec = tempSrec;

		// feeding the local private variables
		this->tdx = this->connector->Sdistance2receivers[this->tnode];
		this->tdy = this->connector->get_travers_dy_from_dx(this->tdx);
		// this->tSS = this->connector->SS[this->tnode];
		if (this->marine == TSC_MARINE::NONE)
			this->tSS = std::max(static_cast<fT>(1e-6),
													 (tZ - this->z_surf[this->tSrec]) / this->tdx);
		else if (this->z_surf[this->tnode] > this->sea_level(this->tnode) &&
						 this->z_surf[this->tSrec] > this->sea_level(this->tSrec))
			this->tSS = std::max(static_cast<fT>(1e-6),
													 (tZ - this->z_surf[this->tSrec]) / this->tdx);
		else if (this->z_surf[this->tnode] < this->sea_level(this->tnode) &&
						 this->z_surf[this->tSrec] < this->sea_level(this->tSrec))
			this->tSS = std::max(static_cast<fT>(1e-6),
													 (tZ - this->z_surf[this->tSrec]) / this->tdx);
		else if (this->z_surf[this->tnode] > this->sea_level(this->tnode) &&
						 this->z_surf[this->tSrec] < this->sea_level(this->tSrec))
			this->tSS = std::max(static_cast<fT>(1e-6),
													 (tZ - this->sea_level(this->tnode)) / this->tdx);
		else
			this->tSS = static_cast<fT>(1e-6);

		this->ths = this->h_sed[this->tnode];

		if (this->secondary_fluvial != TSC_FLUVIAL::NONE) {
			auto onodes =
				this->connector->get_orthogonal_nodes(this->tnode, this->tSrec);
			this->orthonodA = onodes.first;
			this->orthonodB = onodes.second;
		}

		// this->reset_fhED();
	}

	void set_N_boundary_to(fT val)
	{
		for (int i = 0; i < this->connector->nx; ++i) {
			this->z_surf[i] = val;
			this->connector->boundaries.codes[i] = BC::OUT;
		}
		this->connector->precompute_links();
	}

	/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																																													. - ~ ~ ~ - .
..     _      .-~               ~-.
//|     \ `..~                      `.
|| |      }  }              /       \  \
(\   \\ \~^..'                 |         }  \
\`.-~  o      /       }       |        /    \
(__          |       /        |       /      `.
`- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
													|     /          |     /     ~-.     ~- _
													|_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Standalone functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	void Standalone_hyland_landslides()
	{
		// save state
		auto hillslopes_state = this->hillslopes;
		this->hillslopes = TSC_HILLSLOPE::HYLANDS;

		this->init_vectors();
		// first computing graph
		std::vector<fT> faketopo(this->z_surf);
		this->graph->_compute_graph(faketopo, false, false);
		// running the trigger
		this->hillslopes_hylands_trigger();
		for (int i = this->connector->nnodes - 1; i >= 0; --i) {

			this->tnode = this->get_istack_node_MFD(i);

			if (this->connector->boundaries.no_data(this->tnode) ||
					this->connector->flow_out_or_pit(this->tnode))
				continue;

			this->_ready_node_state();
			this->hillslopes_cidre_dep_only_for_highlands();
			this->trans_Qshs_MFD();
		}

		this->apply_dzhsdt();

		// stops here temporarily
	}

	void Standalone_hylands_single_landslide(int node)
	{
		// save state
		auto hillslopes_state = this->hillslopes;
		this->hillslopes = TSC_HILLSLOPE::HYLANDS;

		auto neighbours_links = this->connector->get_empty_neighbour();
		auto neighbours_nodes = this->connector->get_empty_neighbour();
		std::stack<int, std::vector<int>> stackhelper;

		this->init_vectors();
		// first computing graph
		std::vector<fT> faketopo(this->z_surf);
		this->graph->_compute_graph(faketopo, false, false);

		fT slope = 1e-6;
		int nn = this->connector->get_neighbour_idx_nodes_and_links(
			node, neighbours_nodes, neighbours_links);
		for (int j = 0; j < nn; ++j) {
			fT tss = (this->z_surf[neighbours_nodes[j]] - this->z_surf[node]) /
							 this->connector->get_dx_from_links_idx(neighbours_links[j]);
			if (tss > slope)
				slope = tss;
		}

		this->_hylands_trigger_single_landslide(
			node, slope, stackhelper, neighbours_nodes, neighbours_links);

		for (int i = this->connector->nnodes - 1; i >= 0; --i) {

			this->tnode = this->get_istack_node_MFD(i);

			if (this->connector->boundaries.no_data(this->tnode) ||
					this->connector->flow_out_or_pit(this->tnode))
				continue;

			this->_ready_node_state();
			this->hillslopes_cidre_dep_only_for_highlands();
			this->trans_Qshs_MFD();
		}

		this->apply_dzhsdt();

		this->hillslopes = hillslopes_state;
	}

	/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																																													. - ~ ~ ~ - .
..     _      .-~               ~-.
//|     \ `..~                      `.
|| |      }  }              /       \  \
(\   \\ \~^..'                 |         }  \
\`.-~  o      /       }       |        /    \
(__          |       /        |       /      `.
`- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
													|     /          |     /     ~-.     ~- _
													|_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Legacy functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	void run_SFD(fT dt)
	{

		// std::cout << "NEED TO APPLY FLUXES AT THE END OF TIMESTEP OR AT LEAST V.
		// MOTIONS OR CHECK" << std::endl;;

		if (this->Ch_MTSI)
			this->Ch_MTSI_age(dt);

		this->init_vectors();
		this->graph->depression_resolver = DAGGER::DEPRES::cordonnier_carve;
		// std::vector<bool> testlinkchange(this->graph->links),
		// nnodes2change(this->graph->nnodes,false);
		std::vector<fT> faketopo(this->z_surf);
		this->graph->_compute_graph(faketopo, true, false);

		for (int i = this->graph->nnodes - 1; i >= 0; --i) {
			// Getting geometrical info
			// # location
			int node = this->graph->Sstack[i];
			// # Aborting if outnode
			if (!this->connector->flow_out_or_pit(node) == false)
				continue;
			// # receiver
			int rec = this->connector->_Sreceivers[node];
			// # distance to receiver
			fT dx = this->connector->Sdistance2receivers[node];
			// # cell area
			fT cellarea = this->connector->get_area_at_node(node);
			// # local gradient
			fT S = std::max((this->z_surf[node] - this->z_surf[rec]) / dx,
											static_cast<fT>(1e-6));

			// Hydrology
			// # local addition
			this->Qw[node] += cellarea * precipitations(node);

			// preparing processes
			fT fEs = 0., fEr = 0., fDs = 0., hEs = 0., hEr = 0., hDs = 0., mEs = 0.,
				 mEr = 0., mDs = 0.;

			bool isNodeContinental = true;
			if (this->marine != TSC_MARINE::NONE) {
				if (this->z_surf[node] + this->h_sed[node] < this->sea_level(node))
					isNodeContinental = false;
			}

			if (isNodeContinental) {
				// Running the hillslope processes
				if (this->hillslopes == TSC_HILLSLOPE::CIDRE) {
					this->_compute_SFD_hillslopes(
						node, rec, S, hEs, hEr, hDs, dt, dx, cellarea);
				}

				// Fluvial processes
				if (this->fluvial != TSC_FLUVIAL::NONE) {
					this->_compute_SFD_fluvial(
						node, rec, S, fEs, fEr, fDs, dt, dx, cellarea);
				}
			} else {
				this->_compute_SFD_marine(
					node, rec, S, mEs, mEr, mDs, dt, dx, cellarea);
			}
		}
	}

	// Legacy function of the hillslope processes
	void _compute_SFD_hillslopes(int node,
															 int rec,
															 fT& S,
															 fT& hEs,
															 fT& hEr,
															 fT& hDs,
															 fT& dt,
															 fT& dx,
															 fT& cellarea)
	{
		// Storing the proportion of bedrock-power used
		fT propused = 0;

		// If the slope is bellow the critical values
		if (S <= Sc(node) - 1e-6) {
			// If I have sediments
			if (this->h_sed[node] > 0) {

				// erosion power
				hEs = Ks(node) * S;

				// Checking I am not strapping more than the sediment pile
				if (hEs * dt > this->h_sed[node]) {
					// Correcting to the max sed height removal possible and calculating
					// the proportion of sediment used
					propused = this->h_sed[node] / dt / hEs;
					hEs = this->h_sed[node] / dt;
				} else {
					// I don't have any sed? -> no propused
					propused = 1.;
				}

				// end of the check for excessing sed height
			}

			// Remaining applied to bedrock
			hEr = (1. - propused) * Kr(node) * S;
			fT L = cellarea / (1 - std::pow(S / Sc(node), 2));
			hDs = this->Qs_hs[node] / L;

		} else {
			fT tothE = (z_surf[node] - (dx * Sc(node) + z_surf[rec])) / dt;
			if (tothE > this->h_sed[node]) {
				hEs = this->h_sed[node] / dt;
				hEr = tothE - hEs;
			} else
				hEs = tothE;
		}

		if (this->TSP_module)
			this->apply_TSP(node, rec, hEs, hEr, hDs, dt, true);

		if (this->Ch_MTSI)
			this->apply_Ch_MTSI_SFD(node, rec, hEs, hEr, hDs, dt, true);

		this->Qs_hs[node] += (hEs + hEr - hDs) * cellarea;

		// Applying the fluxes
		this->h_sed[node] += (hDs - hEs) * dt;
		this->z_surf[node] += (hDs - hEs - hEr) * dt;
		if (this->Qs_hs[node] < 0)
			this->Qs_hs[node] = 0;
		this->Qs_hs[rec] += this->Qs_hs[node];
	}

	// Legacy function for the fluvial processes
	void _compute_SFD_fluvial(int node,
														int rec,
														fT& S,
														fT& fEs,
														fT& fEr,
														fT& fDs,
														fT& dt,
														fT& dx,
														fT& cellarea)
	{
		// Erosion
		fT stream_power =
			std::pow(this->Qw[node], this->mexp) * std::pow(S, this->nexp);
		fT propused = 0;
		if (this->h_sed[node] > 0) {
			fEs = Ks(node) * stream_power;
			if (this->h_sed[node] < fEs * dt) {
				propused = this->h_sed[node] / dt / fEs;
				fEs = this->h_sed[node] / dt;
			} else {
				propused = 1.;
			}
		}
		fEr = (1. - propused) * Kr(node) * stream_power;

		// Deposition
		fT L = std::max(this->depcoeff(node) * this->Qw[node], cellarea);

		fDs = this->Qs_fluvial[node] / L;

		if (this->TSP_module)
			this->apply_TSP(node, rec, fEs, fEr, fDs, dt, false);

		if (this->Ch_MTSI)
			this->apply_Ch_MTSI_SFD(node, rec, fEs, fEr, fDs, dt, false);

		// Applying the fluxes
		this->h_sed[node] += (fDs - fEs) * dt;
		this->z_surf[node] += (fDs - fEs - fEr) * dt;
		this->Qs_fluvial[node] += (fEs + fEr - fDs) * cellarea;

		if (this->Qs_fluvial[node] < 0)
			this->Qs_fluvial[node] = 0;

		// Transferring to receivers
		this->Qw[rec] += this->Qw[node];
		this->Qs_fluvial[rec] += this->Qs_fluvial[node];
	}

	void run_SFD_implicit(fT dt)
	{

		this->init_vectors();
		this->graph->depression_resolver = DAGGER::DEPRES::cordonnier_carve;
		// std::vector<bool> testlinkchange(this->graph->links),
		// nnodes2change(this->graph->nnodes,false);
		this->graph->_compute_graph(this->z_surf, true, false);

		if (this->variable_precipitations) {
			std::vector<fT> tQA(this->connector->nnodes, 0.);

			for (int i = 0; i < this->connector->nnodes; ++i) {
				tQA[i] =
					this->_precipitations[i] * this->connector->get_area_at_node(i);
			}

			this->Qw = this->graph->_accumulate_variable_downstream_SFD(tQA);

		} else
			this->Qw = this->graph->_accumulate_constant_downstream_SFD(
				this->connector->get_area_at_node(0));

		for (int i = 0; i < this->graph->nnodes; ++i) {

			int node = this->graph->Sstack[i];
			int rec = this->connector->_Sreceivers[node];

			if (this->connector->flow_out_or_pit(node))
				continue;

			fT tK = this->Kr(node);

			fT factor =
				tK * dt * std::pow(this->Qw[node], this->mexp) /
				std::pow(this->connector->Sdistance2receivers[node], this->nexp);

			fT ielevation = this->z_surf[node];
			fT irec_elevation = this->z_surf[rec];

			fT elevation_k = ielevation;
			fT elevation_prev = std::numeric_limits<fT>::max();
			fT tolerance = 1e-4;

			while (abs(elevation_k - elevation_prev) > tolerance) {
				elevation_prev = elevation_k;
				fT slope =
					std::max(elevation_k - irec_elevation, static_cast<fT>(1e-6));
				fT diff =
					(elevation_k - ielevation + factor * std::pow(slope, this->nexp)) /
					(1. + factor * this->nexp * std::pow(slope, this->nexp - 1));
				elevation_k -= diff;
			}

			this->z_surf[node] = elevation_k;
		}
	}

	void _compute_SFD_marine(int node,
													 int rec,
													 fT& S,
													 fT& mEs,
													 fT& mEr,
													 fT& mDs,
													 fT& dt,
													 fT& dx,
													 fT& cellarea)
	{
		// Storing the proportion of bedrock-power used
		// fT propused = 0;

		// updating the node sediment flux with anything that came (river +
		// hillslope)
		this->Qs_hs[node] += this->Qs_fluvial[node];

		// If the slope is bellow the critical values
		if (S <= Sc_M(node) - 1e-6) {
			// If I have sediments
			if (this->h_sed[node] > 0) {

				// erosion power
				mEs = Ke(node) * S;

				// Checking I am not strapping more than the sediment pile
				if (mEs * dt > this->h_sed[node]) {
					// Correcting to the max sed height removal possible and calculating
					// the proportion of sediment used propused = this->h_sed[node]/dt
					// /mEs;
					mEs = this->h_sed[node] / dt;
				} else {
					// I don't have any sed? -> no propused
					// propused = 1.;
				}

				// end of the check for excessing sed height
			}

			// Remaining applied to bedrock
			// mEr = (1. - propused) * Kr(node) * S;

			fT L =
				(this->connector->get_travers_dy_from_dx(dx) * this->lambda(node)) /
				(1 - std::pow(S / Sc_M(node), 2));
			mDs = this->Qs_hs[node] / std::max(L, cellarea);

		} else {
			fT tothE = (z_surf[node] - (dx * Sc(node) + z_surf[rec])) / dt;
			if (tothE > this->h_sed[node]) {
				mEs = this->h_sed[node] / dt;
				// mEr = tothE - mEs;
			} else
				mEs = tothE;
		}

		if (this->TSP_module)
			this->apply_TSP(node, rec, mEs, mEr, mDs, dt, true);

		if (this->Ch_MTSI)
			this->apply_Ch_MTSI_SFD(node, rec, mEs, mEr, mDs, dt, true);

		this->Qs_hs[node] += (mEs + mEr - mDs) * cellarea;

		// Applying the fluxes
		this->h_sed[node] += (mDs - mEs) * dt;
		this->z_surf[node] += (mDs - mEs - mEr) * dt;
		if (this->Qs_hs[node] < 0)
			this->Qs_hs[node] = 0;
		this->Qs_hs[rec] += this->Qs_hs[node];
	}

	template<class out_t>
	out_t get_topo()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->z_surf);
	}

	template<class out_t>
	out_t get_hillshade()
	{
		return DAGGER::hillshade<Connector_t, std::vector<fT>&, out_t, fT>(
			*(this->connector), this->z_surf);
	}

	template<class out_t>
	out_t get_Qw()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->Qw);
	}

	template<class out_t>
	out_t get_Qs_fluvial()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->Qs_fluvial);
	}

	template<class out_t>
	out_t get_Qs_hillslopes()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->Qs_hs);
	}

	template<class out_t>
	out_t get_precipitations()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_precipitations);
	}

	template<class out_t>
	out_t get_h_sed()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->h_sed);
	}

	template<class in_t>
	void external_uplift(in_t& tU, fT dt, bool apply_to_edges = false)
	{
		auto U = DAGGER::format_input(tU);
		for (int i = 0; i < this->graph->nnodes; ++i) {
			if (!this->connector->boundaries.can_out(i) || apply_to_edges) {
				this->z_surf[i] += U[i] * dt;
			}
		}
	}

	void block_uplift(fT rate, fT dt)
	{
		for (int i = 0; i < this->graph->nnodes; ++i) {
			if (this->connector->boundaries.can_out(i) == false)
				this->z_surf[i] += rate * dt;
		}
	}

	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	// TSP module: Tracking Single Provenance
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

	template<class in_t>
	void init_TSP_module(fT dz, in_t& tconcent)
	{
		this->TSP_module = true;
		this->at_least_one_tracking_module_is_activated = true;
		auto concent = DAGGER::format_input(tconcent);
		this->TSP_concentrations = DAGGER::to_vec(concent);
		this->TSP_store =
			VerticalStorer<fT, BasePropStorer<fT>>(dz, this->graph->nnodes);
		this->init_Qs_TSP();
	}

	template<class in_t>
	void update_TSP_source(in_t& tconcent)
	{
		if (this->TSP_module == false)
			throw std::runtime_error(
				"DAGGER::trackscape::update_TSP_source -> cannot update the source "
				"area if TSP is not set. Run DAGGER::TRACKSCAPE::init_TSP_module "
				"function first.");
		// Updating the concentration field
		auto concent = DAGGER::format_input(tconcent);
		this->TSP_concentrations = DAGGER::to_vec(concent);
	}

	// Function reinitialising the trackers for the TSP module
	void init_Qs_TSP()
	{
		if (this->fluvial != TSC_FLUVIAL::NONE)
			this->TSP_Qsf = std::vector<BasePropStorer<fT>>(this->graph->nnodes,
																											BasePropStorer<fT>());
		if (this->hillslopes != TSC_HILLSLOPE::NONE)
			this->TSP_Qsh = std::vector<BasePropStorer<fT>>(this->graph->nnodes,
																											BasePropStorer<fT>());
	}

	void apply_TSP(int i, int ir, fT Es, fT Er, fT Ds, fT dt, bool thillslopes)
	{

		Es *= dt;
		Er *= dt;
		Ds *= dt;

		// Removing from the sediment pile
		auto rem1 = this->TSP_store.remove(i, Es);

		// Depositing on the sediment pile
		BasePropStorer<fT>& tpropr =
			(thillslopes) ? this->TSP_Qsh[i] : this->TSP_Qsf[i];
		this->TSP_store.pile_up(i, Ds, tpropr);

		// Modifying the fluxes
		fT zfQ = (thillslopes)
							 ? this->Qs_hs[i] / this->connector->get_area_at_node(i)
							 : this->Qs_fluvial[i] / this->connector->get_area_at_node(i);
		zfQ *= dt;
		BasePropStorer<fT>::mix(zfQ, tpropr, Es, rem1);
		zfQ += Es;
		BasePropStorer<fT> bdr(this->TSP_concentrations[i]);
		BasePropStorer<fT>::mix(zfQ, tpropr, Er, bdr);
		zfQ += Er;
		zfQ -= Ds;
		if (zfQ < 0)
			zfQ = 0;

		if (!this->connector->flow_out_or_pit(ir)) {
			if (thillslopes)
				BasePropStorer<fT>::mix(this->Qs_hs[ir] /
																	this->connector->get_area_at_node(i) * dt,
																this->TSP_Qsh[ir],
																zfQ,
																tpropr);
			else
				BasePropStorer<fT>::mix(this->Qs_fluvial[ir] /
																	this->connector->get_area_at_node(i) * dt,
																this->TSP_Qsf[ir],
																zfQ,
																tpropr);
		}
	}

	template<class out_t>
	out_t get_TSP_surface_concentrations()
	{
		if (this->TSP_module == false)
			throw std::runtime_error("Cannot return surface TSP if there is no TSP "
															 "module activated (yo!)");

		std::vector<fT> props(this->graph->nnodes, 0.);
		for (int i = 0; i < this->graph->nnodes; ++i) {
			if (!this->connector->flow_out_or_pit(i) == false)
				continue;

			if (this->TSP_store.pile[i].size() > 0)
				props[i] = this->TSP_store.pile[i].back().prop;
		}

		return DAGGER::format_output<std::vector<fT>, out_t>(props);
	}

	template<class out_t>
	out_t sample_carrot_TSP(int row, int col)
	{
		int i = this->connector->nodeid_from_row_col(row, col);

		std::vector<fT> carrot(this->TSP_store.pile[i].size());
		for (size_t j = 0; j < carrot.size(); ++j) {
			carrot[j] = this->TSP_store.pile[i][j].prop;
		}

		return DAGGER::format_output<std::vector<fT>, out_t>(carrot);
	}

	// Returns a cross section in the X (rwo = true) or Y (row = False) drection
	// for nz number of depth cells
	template<class out_t>
	out_t get_transect_TSP(int rowcol, int nz, bool row = true)
	{
		// number of nodes in the desired dimentsion
		int nxy = (row) ? this->connector->nx : this->connector->ny;
		// total number of nodes
		int nn = nxy * nz;

		std::vector<fT> transect(nn, 0.);
		int trow = rowcol;
		int tcol = rowcol;

		fT maxz = std::numeric_limits<fT>::min();
		for (int r = 0; r < nxy; ++r) {
			trow = (row) ? rowcol : r;
			tcol = (row) ? r : rowcol;
			int i = this->connector->nodeid_from_row_col(trow, tcol);
			if (this->z_surf[i] > maxz)
				maxz = this->z_surf[i];
		}

		for (int r = 0; r < nxy; ++r) {
			trow = (row) ? rowcol : r;
			tcol = (row) ? r : rowcol;
			int i = this->connector->nodeid_from_row_col(trow, tcol);

			fT tz = maxz;
			int j2 = int(this->TSP_store.pile[i].size() - 1);
			;
			for (int j = 0; j < nz; ++j) {
				int idxj = r * nz + j;
				if (tz > this->z_surf[i]) {
					transect[idxj] = -1.;
				} else {
					transect[idxj] = (j2 > 0) ? this->TSP_store.pile[i][j2].prop : -1.;
					--j2;
				}

				tz -= this->TSP_store.dz;
			}
		}

		return DAGGER::format_output<std::vector<fT>, out_t>(transect);
	}

	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	// Ch_MTSI module: timing the time since incision
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

	void init_Ch_MTSI(fT dz)
	{
		this->Ch_MTSI = true;
		this->at_least_one_tracking_module_is_activated = true;
		this->Ch_MTSI_store =
			VerticalStorer<fT, BasePropStorer<fT>>(dz, this->graph->nnodes);
		this->init_Qs_TSP();
	}

	// Function reinitialising the trackers for the Ch_MTSI module
	void init_Qs_Ch_MTSI()
	{
		if (this->fluvial != TSC_FLUVIAL::NONE)
			this->Ch_MTSI_Qsf = std::vector<BasePropStorer<fT>>(this->graph->nnodes,
																													BasePropStorer<fT>());
		if (this->hillslopes != TSC_HILLSLOPE::NONE)
			this->Ch_MTSI_Qsh = std::vector<BasePropStorer<fT>>(this->graph->nnodes,
																													BasePropStorer<fT>());
	}

	void
	apply_Ch_MTSI_SFD(int i, int ir, fT Es, fT Er, fT Ds, fT dt, bool thillslopes)
	{

		Es *= dt;
		Er *= dt;
		Ds *= dt;

		// Removing from the sediment pile
		auto rem1 = this->Ch_MTSI_store.remove(i, Es);

		// Depositing on the sediment pile
		BasePropStorer<fT>& tpropr =
			(thillslopes) ? this->Ch_MTSI_Qsh[i] : this->Ch_MTSI_Qsf[i];
		this->Ch_MTSI_store.pile_up(i, Ds, tpropr);

		// Modifying the fluxes
		fT zfQ = (thillslopes)
							 ? this->Qs_hs[i] / this->connector->get_area_at_node(i)
							 : this->Qs_fluvial[i] / this->connector->get_area_at_node(i);
		zfQ *= dt;
		BasePropStorer<fT>::mix(zfQ, tpropr, Es, rem1);
		zfQ += Es;
		BasePropStorer<fT> bdr(0.);
		BasePropStorer<fT>::mix(zfQ, tpropr, Er, bdr);
		zfQ += Er;
		zfQ -= Ds;
		if (zfQ < 0)
			zfQ = 0;

		if (!this->connector->flow_out_or_pit(ir)) {
			if (thillslopes)
				BasePropStorer<fT>::mix(this->Qs_hs[ir] /
																	this->connector->get_area_at_node(i) * dt,
																this->Ch_MTSI_Qsh[ir],
																zfQ,
																tpropr);
			else
				BasePropStorer<fT>::mix(this->Qs_fluvial[ir] /
																	this->connector->get_area_at_node(i) * dt,
																this->Ch_MTSI_Qsf[ir],
																zfQ,
																tpropr);
		}
	}

	void Ch_MTSI_age(fT dt)
	{
		for (int i = 0; i < this->graph->nnodes; ++i) {
			for (size_t j = 0; j < this->Ch_MTSI_store.pile[i].size(); ++j)
				this->Ch_MTSI_store.pile[i][j].prop += dt;
		}
	}

	template<class out_t>
	out_t get_Ch_MTIS_surface_age()
	{
		if (this->Ch_MTSI == false)
			throw std::runtime_error("Cannot return surface Ch_MTSI if there is no "
															 "Ch_MTSI module activated (yo!)");

		std::vector<fT> timeryolo(this->graph->nnodes, 0.);
		for (int i = 0; i < this->graph->nnodes; ++i) {
			if (!this->connector->flow_out_or_pit(i) == false)
				continue;

			if (this->Ch_MTSI_store.pile[i].size() > 0)
				timeryolo[i] = this->Ch_MTSI_store.pile[i].back().prop;
		}

		return DAGGER::format_output<std::vector<fT>, out_t>(timeryolo);
	}

	template<class out_t>
	out_t sample_carrot_Ch_MTSI(int row, int col)
	{
		int i = this->connector->nodeid_from_row_col(row, col);

		std::vector<fT> carrot(this->Ch_MTSI_store.pile[i].size());
		for (size_t j = 0; j < carrot.size(); ++j) {
			carrot[j] = this->Ch_MTSI_store.pile[i][j].prop;
		}

		return DAGGER::format_output<std::vector<fT>, out_t>(carrot);
	}

	template<class out_t>
	out_t get_transect_Ch_MTSI(int rowcol, int nz, bool row = true)
	{
		// number of nodes in the desired dimentsion
		int nxy = (row) ? this->connector->nx : this->connector->ny;
		// total number of nodes
		int nn = nxy * nz;

		std::vector<fT> transect(nn, 0.);
		int trow = rowcol;
		int tcol = rowcol;

		fT maxz = std::numeric_limits<fT>::min();
		for (int r = 0; r < nxy; ++r) {
			trow = (row) ? rowcol : r;
			tcol = (row) ? r : rowcol;
			int i = this->connector->nodeid_from_row_col(trow, tcol);
			if (this->z_surf[i] > maxz)
				maxz = this->z_surf[i];
		}

		for (int r = 0; r < nxy; ++r) {
			trow = (row) ? rowcol : r;
			tcol = (row) ? r : rowcol;
			int i = this->connector->nodeid_from_row_col(trow, tcol);

			fT tz = maxz;
			int j2 = int(this->Ch_MTSI_store.pile[i].size() - 1);
			;
			for (int j = 0; j < nz; ++j) {
				int idxj = r * nz + j;
				if (tz > this->z_surf[i]) {
					transect[idxj] = -1.;
				} else {
					transect[idxj] =
						(j2 > 0) ? this->Ch_MTSI_store.pile[i][j2].prop : -1.;
					--j2;
				}

				tz -= this->Ch_MTSI_store.dz;
			}
		}

		return DAGGER::format_output<std::vector<fT>, out_t>(transect);
	}

	/*
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
																																																	. - ~ ~ ~ - .
	..     _      .-~               ~-.
	//|     \ `..~                      `.
	|| |      }  }              /       \  \
	(\   \\ \~^..'                 |         }  \
	\`.-~  o      /       }       |        /    \
	(__          |       /        |       /      `.
	`- - ~ ~ -._|      /_ - ~ ~ ^|      /- _      `.
																	|     /          |     /     ~-.     ~- _
																	|_____|          |_____|         ~ - . _ _~_-_

	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	Standalone functions
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	*/

	void standalone_implicit_SPL()
	{
		// Do I need to calculate MFD recs at all?
		bool need_mfrecs = this->do_I_need_MFD();

		// pointer to the function getting the next node
		auto stackgetter =
			(this->flowtopo == TSC_FLOW_TOPOLOGY::SFD)
				? &trackscape<fT, Graph_t, Connector_t>::get_istack_node_SFD
				: &trackscape<fT, Graph_t, Connector_t>::get_istack_node_MFD;

		// Computing the graph
		std::vector<fT> faketopo(this->z_surf);

		fT Smax = 0.;

		this->graph->_compute_graph(faketopo, !need_mfrecs, false);
		// throw std::runtime_error("blug");

		if (this->flowtopo == TSC_FLOW_TOPOLOGY::SFD) {
			if (this->variable_precipitations) {
				this->Qw = this->graph->_accumulate_variable_downstream_area_SFD(
					this->_precipitations);
			} else {
				this->Qw = this->graph->_accumulate_constant_downstream_area_SFD(
					this->_precipitations[0]);
			}
		} else {
			auto gradient = this->connector->_get_links_gradient(this->z_surf, 1e-6);
			auto weights = this->connector->_get_link_weights(gradient, 1);
			if (this->variable_precipitations) {
				this->Qw = this->graph->_accumulate_variable_downstream_area_MFD(
					weights, this->_precipitations);
			} else {
				this->Qw = this->graph->_accumulate_constant_downstream_area_MFD(
					weights, this->_precipitations[0]);
			}
		}

		for (int i = 0; i < this->connector->nxy(); ++i) {

			// Getting the next node in line
			this->tnode = (this->*stackgetter)(i);

			this->_ready_node_state();

			// Check if no data
			if (this->connector->boundaries.no_data(this->tnode))
				continue;

			// Check if base level
			if (this->connector->flow_out_or_pit(this->tnode)) {
				// manage the base level evolution here
				continue;
			}

			this->fluvial_fastscape_SFD();
		}
	}

	void standalone_cidre_hillslopes()
	{

		this->fluvial = TSC_FLUVIAL::NONE;
		this->hillslopes = TSC_HILLSLOPE::CIDRE;
		this->marine = TSC_MARINE::NONE;

		// Initialising the data vectors
		this->init_vectors();

		bool need_mfrecs = this->do_I_need_MFD();

		// Setting up the function pointers
		auto stackgetter =
			(this->flowtopo == TSC_FLOW_TOPOLOGY::SFD)
				? &trackscape<fT, Graph_t, Connector_t>::get_istack_node_SFD
				: &trackscape<fT, Graph_t, Connector_t>::get_istack_node_MFD;

		std::vector<fT> faketopo(this->z_surf);
		this->graph->_compute_graph(faketopo, !need_mfrecs, false);

		const bool need_Qs = this->fluvial == TSC_FLUVIAL::DAVY2009 ||
												 this->hillslopes != TSC_HILLSLOPE::NONE;
		const bool SFD = this->flowtopo == TSC_FLOW_TOPOLOGY::SFD;

		for (int i = this->graph->nnodes - 1; i >= 0; --i) {

			// Getting the next node in line
			this->tnode = (this->*stackgetter)(i);

			// Check if no data
			if (this->connector->boundaries.no_data(this->tnode))
				continue;

			// Check if base level
			if (this->connector->flow_out_or_pit(this->tnode)) {
				// manage the base level evolution here
				continue;
			}

			this->_ready_node_state();

			this->hillslopes_cidre();
			if (SFD)
				this->trans_Qshs_SFD();
			else if (this->full_stochastic == false)
				this->trans_Qshs_MFD();
			else
				this->trans_Qshs_MFD_stochastic();
		}

		this->apply_dzhsdt();
	}

	void standalone_davy2009()
	{
		this->fluvial = TSC_FLUVIAL::DAVY2009;
		this->hillslopes = TSC_HILLSLOPE::NONE;
		this->marine = TSC_MARINE::NONE;

		// Initialising the data vectors
		this->init_vectors();

		bool need_mfrecs = this->do_I_need_MFD();

		// Setting up the function pointers
		auto stackgetter =
			(this->flowtopo == TSC_FLOW_TOPOLOGY::SFD)
				? &trackscape<fT, Graph_t, Connector_t>::get_istack_node_SFD
				: &trackscape<fT, Graph_t, Connector_t>::get_istack_node_MFD;

		std::vector<fT> faketopo(this->z_surf);
		this->graph->_compute_graph(faketopo, !need_mfrecs, false);

		if (this->add_extra_Qs_fluvial)
			this->prec_extra_Qs_fluvial();

		if (this->add_extra_Qw_fluvial)
			this->prec_extra_Qw_fluvial();

		const bool need_Qs = this->fluvial == TSC_FLUVIAL::DAVY2009 ||
												 this->hillslopes != TSC_HILLSLOPE::NONE;
		const bool SFD = this->flowtopo == TSC_FLOW_TOPOLOGY::SFD;

		const bool need_Qw = this->fluvial != TSC_FLUVIAL::NONE;

		for (int i = this->graph->nnodes - 1; i >= 0; --i) {

			// Getting the next node in line
			this->tnode = (this->*stackgetter)(i);

			// Check if no data
			if (this->connector->boundaries.no_data(this->tnode))
				continue;

			// Check if base level
			if (this->connector->flow_out_or_pit(this->tnode)) {
				// manage the base level evolution here
				continue;
			}

			this->_ready_node_state();

			if (need_Qw) {
				if (SFD)
					this->prec_Qw_SFD();
				else
					this->prec_Qw_MFD();
			}

			this->fluvial_davy2009();

			if (SFD)
				this->trans_Qw_Qs_SFD();
			else
				this->trans_Qw_Qs_MFD();
		}

		this->apply_dzhsdt();
	}

	void standalone_cidre()
	{
		this->fluvial = TSC_FLUVIAL::DAVY2009;
		this->hillslopes = TSC_HILLSLOPE::CIDRE;
		this->marine = TSC_MARINE::NONE;

		// Initialising the data vectors
		this->init_vectors();

		bool need_mfrecs = this->do_I_need_MFD();

		// Setting up the function pointers
		auto stackgetter =
			(this->flowtopo == TSC_FLOW_TOPOLOGY::SFD)
				? &trackscape<fT, Graph_t, Connector_t>::get_istack_node_SFD
				: &trackscape<fT, Graph_t, Connector_t>::get_istack_node_MFD;

		std::vector<fT> faketopo(this->z_surf);
		this->graph->_compute_graph(faketopo, !need_mfrecs, false);

		if (this->add_extra_Qs_fluvial)
			this->prec_extra_Qs_fluvial();

		if (this->add_extra_Qw_fluvial)
			this->prec_extra_Qw_fluvial();

		const bool need_Qs = this->fluvial == TSC_FLUVIAL::DAVY2009 ||
												 this->hillslopes != TSC_HILLSLOPE::NONE;
		const bool SFD = this->flowtopo == TSC_FLOW_TOPOLOGY::SFD;

		const bool need_Qw = this->fluvial != TSC_FLUVIAL::NONE;

		for (int i = this->graph->nnodes - 1; i >= 0; --i) {

			// Getting the next node in line
			this->tnode = (this->*stackgetter)(i);

			// Check if no data
			if (this->connector->boundaries.no_data(this->tnode))
				continue;

			// Check if base level
			if (this->connector->flow_out_or_pit(this->tnode)) {
				// manage the base level evolution here
				continue;
			}

			this->_ready_node_state();

			if (need_Qw) {
				if (SFD)
					this->prec_Qw_SFD();
				else
					this->prec_Qw_MFD();
			}

			this->hillslopes_cidre();
			this->fluvial_davy2009();

			if (SFD)
				this->trans_Qw_Qs_SFD();
			else
				this->trans_Qw_Qs_MFD();
		}

		this->apply_dzhsdt();
	}

	void exp_std_cidre_sepA(fT th_Qw)
	{
		this->fluvial = TSC_FLUVIAL::DAVY2009;
		this->hillslopes = TSC_HILLSLOPE::CIDRE;
		this->marine = TSC_MARINE::NONE;

		// Initialising the data vectors
		this->init_vectors();

		bool need_mfrecs = this->do_I_need_MFD();

		// Setting up the function pointers
		auto stackgetter =
			(this->flowtopo == TSC_FLOW_TOPOLOGY::SFD)
				? &trackscape<fT, Graph_t, Connector_t>::get_istack_node_SFD
				: &trackscape<fT, Graph_t, Connector_t>::get_istack_node_MFD;

		std::vector<fT> faketopo(this->z_surf);
		this->graph->_compute_graph(faketopo, !need_mfrecs, false);

		if (this->add_extra_Qs_fluvial)
			this->prec_extra_Qs_fluvial();

		if (this->add_extra_Qw_fluvial)
			this->prec_extra_Qw_fluvial();

		const bool need_Qs = this->fluvial == TSC_FLUVIAL::DAVY2009 ||
												 this->hillslopes != TSC_HILLSLOPE::NONE;
		const bool SFD = this->flowtopo == TSC_FLOW_TOPOLOGY::SFD;

		const bool need_Qw = this->fluvial != TSC_FLUVIAL::NONE;

		for (int i = this->graph->nnodes - 1; i >= 0; --i) {

			// Getting the next node in line
			this->tnode = (this->*stackgetter)(i);

			// Check if no data
			if (this->connector->boundaries.no_data(this->tnode))
				continue;

			// Check if base level
			if (this->connector->flow_out_or_pit(this->tnode)) {
				// manage the base level evolution here
				continue;
			}

			this->_exp_ready_node_state();

			if (need_Qw) {
				if (SFD)
					this->prec_Qw_SFD();
				else
					this->prec_Qw_MFD();
			}

			if (this->Qw[this->tnode] >= th_Qw) {

				this->fluvial_davy2009();
				// this->hillslopes_cidre_longdep_dep_only();
				this->hillslopes_cidre_longdep();

			} else {
				this->hillslopes_cidre_longdep();
			}

			if (SFD)
				this->trans_Qw_Qs_SFD();
			else
				this->trans_Qw_Qs_MFD();
		}

		this->apply_dzhsdt();
	}

	// ##################################################
	// ##################################################
	// ##################################################
	// ##################################################
	//
	//  	         ___..._
	//      _,--'       "`-.
	//    ,'.  .            |
	//  ,/:. .     .       .'
	//  |;..  .      _..--'
	//  `--:...-,-'""|
	//          |:.  `.
	//          l;.   l
	//          `|:.   |
	//           |:.   `.,
	//          .l;.    j, ,
	//       `. |`;:.   //,/
	//        .||)`;,||'/(
	//
	//
	//  Helper functions
	// ##################################################
	// ##################################################
	// ##################################################
	// ##################################################

	int get_istack_node_SFD(int i) { return int(this->graph->Sstack[i]); }

	int get_istack_node_MFD(int i) { return int(this->graph->stack[i]); }

	// this function "transforme" instantly all the sediments into rocks
	// Not possible with tracking activated
	void lithify()
	{
		if (this->at_least_one_tracking_module_is_activated)
			throw std::runtime_error("Cannot lothify if tracking is activated");

		for (int i = 0; i < this->connector->nnodes; ++i) {
			this->z_surf[i] += this->h_sed[i];
			this->h_sed[i] = 0;
		}
	}

	// this function remove instantly all the sediments
	// Not possible with tracking activated
	void strip_sediment()
	{
		if (this->at_least_one_tracking_module_is_activated)
			throw std::runtime_error(
				"Cannot remove all the seds if tracking is activated");

		for (int i = 0; i < this->connector->nnodes; ++i)
			this->h_sed[i] = 0;
	}

	void rise_boundary_by(std::string wish, fT woosh)
	{
		for (int i = 0; i < this->graph->nnodes; ++i) {
			if (wish == "N" && this->connector->is_on_top_row(i))
				this->z_surf[i] += woosh;
			else if (wish == "E" && this->connector->is_on_rightest_col(i))
				this->z_surf[i] += woosh;
			else if (wish == "W" && this->connector->is_on_leftest_col(i))
				this->z_surf[i] += woosh;
			else if (wish == "S" && this->connector->is_on_bottom_row(i))
				this->z_surf[i] += woosh;
		}
	}
};

} // namespace DAGGER

#endif
