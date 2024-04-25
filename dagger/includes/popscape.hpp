#ifndef POPSCAPE_HPP
#define POPSCAPE_HPP

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
#include "popscape_utils.hpp"
// #include "D8connector.hpp"
#include "graph.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "wrap_helper.hpp"

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
// #include <pybind11/numpy.h>

// namespace py = pybind11;

namespace DAGGER {

enum class PARAM_MODE
{
	CONSTANT,
	VARIABLE,
};

template<class fT, class Graph_t, class Connector_t>
class popscape
{

public:
	// stored as direct children: they will be changed anyway
	Graph_t graph;
	Connector_t connector;

	std::vector<fT> QA;
	std::vector<fT> topography;

	PARAM_MODE Kbase_mode = PARAM_MODE::CONSTANT;
	std::vector<fT> _Kbase = { 1e-3 };
	PARAM_MODE Kmod_mode = PARAM_MODE::CONSTANT;
	std::vector<fT> _Kmod = { 1. };
	PARAM_MODE m_mode = PARAM_MODE::CONSTANT;
	std::vector<fT> _m = { 0.45 };
	PARAM_MODE n_mode = PARAM_MODE::CONSTANT;
	std::vector<fT> _n = { 1. };
	PARAM_MODE precip_mode = PARAM_MODE::CONSTANT;
	std::vector<fT> _precip = { 1. };

	PARAM_MODE UE_mode = PARAM_MODE::CONSTANT;
	std::vector<fT> _UE = { 1e-3 };

	std::string boundary_string = "periodic_EW";

	popscape() { ; }
	popscape(Graph_t& graph, Connector_t& con)
	{

		this->connector =
			_create_connector(con.nx, con.ny, con.dx, con.dy, fT(0.), fT(0.));
		_create_graph(graph.nnodes, this->connector, this->graph);

		// this->graph = graph;
	}

	template<class in_t>
	void set_topo(in_t& itopo)
	{
		auto ftopo = DAGGER::format_input(itopo);
		this->topography = to_vec(ftopo);
	}

	template<class out_t>
	out_t get_topo()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->topography);
	}
	template<class out_t>
	out_t get_QA()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->QA);
	}

	template<class out_t>
	out_t get_chistar()
	{
		std::vector<fT> chistar = this->_chi_star();
		return DAGGER::format_output<std::vector<fT>, out_t>(chistar);
	}

	// Parameters:
	fT Kbase(int i)
	{
		if (this->Kbase_mode == PARAM_MODE::CONSTANT)
			return this->_Kbase[0];
		else
			return this->_Kbase[i];
	}

	fT Kmod(int i)
	{
		if (this->Kmod_mode == PARAM_MODE::CONSTANT)
			return this->_Kmod[0];
		else
			return this->_Kmod[i];
	}

	fT m(int i)
	{
		if (this->m_mode == PARAM_MODE::CONSTANT)
			return this->_m[0];
		else
			return this->_m[i];
	}

	fT n(int i)
	{
		if (this->n_mode == PARAM_MODE::CONSTANT)
			return this->_n[0];
		else
			return this->_n[i];
	}

	// Parameters:
	fT precip(int i)
	{
		if (this->precip_mode == PARAM_MODE::CONSTANT)
			return this->_precip[0];
		else
			return this->_precip[i];
	}
	fT UE(int i)
	{
		if (this->UE_mode == PARAM_MODE::CONSTANT)
			return this->_UE[0];
		else
			return this->_UE[i];
	}

	void _init_vecs() { this->QA = std::vector<fT>(this->connector.nnodes, 0.); }

	void StSt(int n_iterations)
	{
		// running the code for n_iterations
		for (int nit = 0; nit < n_iterations; ++nit) {
			this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
			this->graph._compute_graph(this->topography, true, false);
			this->_init_vecs();
			this->QA = this->graph._accumulate_constant_downstream_SFD(
				this->connector.get_area_at_node(0));
			for (int i = 0; i < this->graph.nnodes; ++i) {
				int node = this->graph.Sstack[i];

				if (this->connector.flow_out_or_pit(node))
					continue;

				int rec = this->connector._Sreceivers[node];
				this->topography[node] =
					this->topography[rec] +
					this->connector.Sdistance2receivers[node] *
						std::pow(this->UE(node) / (this->Kbase(node) * this->Kmod(node)),
										 1. / this->n(node)) *
						std::pow(this->QA[node], -this->m(node) / this->n(node));
			}
		}
	}

	void restriction(fT noise_intensity)
	{

		int nxy = 4 * this->graph.nnodes;

		int nx = this->connector.nx * 2;
		int ny = this->connector.ny * 2;

		std::vector<fT> ntopo(nxy, 0.);

		for (int i = 0; i < this->graph.nnodes; ++i) {
			int row, col;
			this->connector.rowcol_from_node_id(i, row, col);
			int nid = row * 2 * nx + col * 2;
			ntopo[nid] =
				this->topography[i] + this->connector.randu->get() * noise_intensity;
			nid = row * 2 * nx + col * 2 + 1;
			ntopo[nid] =
				this->topography[i] + this->connector.randu->get() * noise_intensity;
			nid = (row * 2 + 1) * nx + col * 2 + 1;
			ntopo[nid] =
				this->topography[i] + this->connector.randu->get() * noise_intensity;
			nid = (row * 2 + 1) * nx + col * 2;
			ntopo[nid] =
				this->topography[i] + this->connector.randu->get() * noise_intensity;
		}

		fT dx = this->connector.dx / 2;
		fT dy = this->connector.dy / 2;
		// init connector
		this->connector = _create_connector(nx, ny, dx, dy, fT(0.), fT(0.));
		this->connector.set_default_boundaries(this->boundary_string);

		// init graph
		_create_graph(nxy, this->connector, this->graph);
		// graph.
		// this->graph.init_graph(this->connector);
		this->topography = std::move(ntopo);
		this->_init_vecs();
		this->border20();
	}

	// TODO
	void interpolation()
	{
		// int nxy = floor(this->graph.nnodes/4);

		int nx = floor(this->connector.nx / 2);
		int ny = floor(this->connector.ny / 2);
		int nxy = nx * ny;

		std::vector<fT> ntopo(nxy, 0.);

		D8connector<fT> ncon = _create_connector(
			nx, ny, this->connector.dx * 2, this->connector.dy * 2, fT(0.), fT(0.));
		ncon.set_default_boundaries(this->boundary_string);

		for (int i = 0; i < nxy; ++i) {
			int row, col;
			ncon.rowcol_from_node_id(i, row, col);
			int nid = row * 2 * this->connector.nx + col * 2;
			fT tntopo = 0.;
			tntopo += this->topography[nid];
			nid = row * 2 * this->connector.nx + col * 2 + 1;
			tntopo += this->topography[nid];
			nid = (row * 2 + 1) * this->connector.nx + col * 2 + 1;
			tntopo += this->topography[nid];
			nid = (row * 2 + 1) * this->connector.nx + col * 2;
			tntopo += this->topography[nid];

			ntopo[i] = tntopo / 4;
		}
		// init connector
		this->connector = std::move(ncon);
		// init graph
		_create_graph(nxy, this->connector, this->graph);
		// this->graph.init_graph(this->connector);
		this->topography = std::move(ntopo);
		this->_init_vecs();
		this->border20();
	}

	void smooth(fT rr)
	{
		// this->topography = On_gaussian_blur<fT>(rr, this->topography,
		// this->connector.nx, this->connector.ny);
		this->topography = On_gaussian_blur(
			rr, this->topography, this->connector.nx, this->connector.ny);
	}

	void border20()
	{
		for (int i = 0; i < this->connector.nnodes; ++i) {
			if (this->connector.flow_out_model(i))
				this->topography[i] = 0;
		}
	}

	std::vector<fT> _chi_star()
	{
		std::vector<fT> A = this->graph._accumulate_constant_downstream_SFD(
			this->connector.get_area_at_node(0));
		std::vector<fT> chistar(this->graph.nnodes, 0.);
		fT chimax = 0;
		for (int i = 0; i < this->graph.nnodes; ++i) {
			int node = this->graph.Sstack[i];
			int rec = this->connector._Sreceivers[node];
			if (node == rec)
				continue;

			// chistar[node] = chistar[rec] +
			// this->connector.Sdistance2receivers[node] * (std::pow(1/A[node],
			// this->m(node)/this->n(node)));
			chistar[node] = chistar[rec] + this->connector.Sdistance2receivers[node] *
																			 (std::pow(1 / A[node], 0.2));
			if (chistar[node] > chimax)
				chimax = chistar[node];
		}

		for (auto& v : chistar)
			v /= chimax;

		return chistar;
	}

	std::vector<fT> _z_star()
	{
		std::vector<fT> zstar(this->topography);
		fT zmax = 0;
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (zstar[i] > zmax)
				zmax = zstar[i];
		}

		for (auto& v : zstar)
			v /= zmax;

		return zstar;
	}

	void simple_Kfchi(fT tkmod, fT chimin, fT chimax)
	{

		auto chistar = this->_chi_star();

		this->Kmod_mode = PARAM_MODE::VARIABLE;
		this->_Kmod = std::vector<fT>(this->graph.nnodes, 1.);
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (chistar[i] > chimin && chistar[i] < chimax)
				this->_Kmod[i] = tkmod;
		}

		this->StSt(1);

		this->Kmod_mode = PARAM_MODE::CONSTANT;
		this->_Kmod = { 1 };
	}

	void simple_Kfz(fT tkmod, fT zmin, fT zmax)
	{

		auto zstar = this->_z_star();

		this->Kmod_mode = PARAM_MODE::VARIABLE;
		this->_Kmod = std::vector<fT>(this->graph.nnodes, 1.);
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (zstar[i] > zmin && zstar[i] < zmax)
				this->_Kmod[i] = tkmod;
		}

		this->StSt(1);

		this->Kmod_mode = PARAM_MODE::CONSTANT;
		this->_Kmod = { 1 };
	}

	void normalise_topography()
	{
		fT mini = std::numeric_limits<fT>::max();
		fT maxi = std::numeric_limits<fT>::min();
		for (auto v : this->topography) {
			if (v < mini)
				mini = v;
			if (v > maxi)
				maxi = v;
		}
		maxi -= mini;
		for (auto& v : this->topography)
			v = (v - mini) / maxi;
	}

	// std::vector<fT> chiculations()
	// {

	// }
};

// Former Popscape, to be deprecated soon
template<class fT, class Graph_t, class Connector_t>
class popscape_old
{
public:
	// the graph is directly integrated into popscape
	Graph_t graph;
	// And so is the connector as popscape is optimised to pop a final landscape.
	Connector_t connector;

	// Topography
	std::vector<fT> topography, QA;

	// randomiser helper
	std::shared_ptr<DAGGER::easyRand> randu =
		std::make_shared<DAGGER::easyRand>();

	// empty constructor
	popscape_old(){};

	// constructor 1:
	// -> noise type is from the RANDOISE enum (white, red, perlin, ...)
	// -> nx/ny are the number of nodes in the x/y dir
	// -> dx dy are the related spacing in [L]
	popscape_old(RANDNOISE noisetype, int start_nx, int nit, fT dx, fT dy)
	{

		if (start_nx % 8 != 0)
			throw std::runtime_error(
				"target nx and start nx needs to be a multiple of 8");

		// total number of nodes (so far assuming only regular grids)
		int nx = start_nx, ny = start_nx;
		int nxy = nx * ny;

		// init the topo to 0
		this->topography = std::vector<fT>(nxy, 0.);

		// init connector
		this->connector = _create_connector(nx, ny, dx, dy, fT(0.), fT(0.));

		// init graph
		_create_graph(nxy, this->connector, this->graph);

		// init random noise
		if (noisetype == RANDNOISE::WHITE)
			DAGGER::add_noise_to_vector(this->topography, 0., 1.);

		this->connector.set_default_boundaries("periodic_EW");

		for (int i = 0; i < nit; ++i) {
			// std::cout << "REFINING::" << i << std::endl;
			int nnit = (i == 0) ? 50 : 5;
			// std::cout << "REFINING::A" << i << std::endl;
			this->solve_generic(0.45, 1.11, nnit);
			// std::cout << "REFINING::B" << i << std::endl;
			if (i < nit - 1) {
				// std::cout << "REFINING::C" << i << std::endl;

				this->double_res(10, noisetype);
			}
			// std::cout << "done" << std::endl;
		}
	}

	void _init_vecs() { this->QA = std::vector<fT>(this->graph.nnodes, 0.); }

	void solve_generic(fT m, fT n, int n_iterations)
	{
		// std::cout << "REFINING::|||||" <<n_iterations << std::endl;
		// running the code for n_iterations
		for (int nit = 0; nit < n_iterations; ++nit) {
			// std::cout << nit << "|A" << std::endl;

			this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
			this->graph._compute_graph(this->topography, true, false);
			this->_init_vecs();

			this->QA = this->graph._accumulate_constant_downstream_SFD(
				this->connector.get_area_at_node(0));

			// std::cout << nit << "|C" << std::endl;

			for (int i = 0; i < this->graph.nnodes; ++i) {
				int node = this->graph.Sstack[i];
				if (this->connector.flow_out_or_pit(node))
					continue;
				int rec = this->connector._Sreceivers[node];
				this->topography[node] =
					this->topography[rec] + this->connector.Sdistance2receivers[node] *
																		std::pow(1e2, 1. / n) /
																		std::pow(this->QA[node], m / n);
				// this->topography[node] = this->topography[rec] + 1;
			}
			// std::cout << nit << "|D" << std::endl;
		}
		// std::cout << std::endl;
	}

	void double_res(fT noise_intensity, RANDNOISE noisetype)
	{

		int nxy = 4 * this->graph.nnodes;

		int nx = this->connector.nx * 2;
		int ny = this->connector.ny * 2;

		std::vector<fT> ntopo(nxy, 0.);

		if (noisetype == RANDNOISE::WHITE)
			for (int i = 0; i < this->graph.nnodes; ++i) {
				int row, col;
				this->connector.rowcol_from_node_id(i, row, col);
				int nid = row * 2 * nx + col * 2;
				ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
				nid = row * 2 * nx + col * 2 + 1;
				ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
				nid = (row * 2 + 1) * nx + col * 2 + 1;
				ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
				nid = (row * 2 + 1) * nx + col * 2;
				ntopo[nid] = this->topography[i] + this->randu->get() * noise_intensity;
			}

		fT dx = this->connector.dx / 2;
		fT dy = this->connector.dy / 2;
		// init connector
		this->connector = _create_connector(nx, ny, dx, dy, fT(0.), fT(0.));

		// init graph
		_create_graph(nxy, this->connector, this->graph);
		// this->graph.init_graph(this->connector);
		this->topography = std::move(ntopo);
		this->_init_vecs();
	}

	void compute_graph(bool SFD_only)
	{
		this->graph.depression_resolver = DAGGER::DEPRES::cordonnier_carve;
		this->graph._compute_graph(this->topography, SFD_only, false);
	}

	void compute_DA_SFD()
	{
		this->QA = this->graph._accumulate_constant_downstream_SFD(
			this->connector.get_area_at_node(0));
	}

	void solve_SFD_SPL_imp(fT m, fT n, fT K, fT dt)
	{

		for (int i = 0; i < this->graph.nnodes; ++i) {

			int node = this->graph.Sstack[i];
			int rec = this->connector._Sreceivers[node];

			if (!this->connector.flow_out_or_pit(node) == false)
				continue;

			fT factor = K * dt * std::pow(this->QA[node], m) /
									std::pow(this->connector.Sdistance2receivers[node], n);

			fT ielevation = this->topography[node];
			fT irec_elevation = this->topography[rec];

			fT elevation_k = ielevation;
			fT elevation_prev = std::numeric_limits<fT>::max();
			fT tolerance = 1e-4;
			int nit = 0;
			while (abs(elevation_k - elevation_prev) > tolerance && nit < 10) {
				elevation_prev = elevation_k;
				fT slope =
					std::max(elevation_k - irec_elevation, static_cast<fT>(1e-6));
				fT diff = (elevation_k - ielevation + factor * std::pow(slope, n)) /
									(1. + factor * n * std::pow(slope, n - 1));
				elevation_k -= diff;
				++nit;
			}

			this->topography[node] = elevation_k;
		}
	}

	template<class out_t>
	out_t get_topo()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->topography);
	}
	template<class out_t>
	out_t get_QA()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->QA);
	}

	void apply_uplift(fT dt, fT U)
	{
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (!this->connector.flow_out_or_pit(i))
				this->topography[i] += U * dt;
		}
	}

	template<class out_t>
	void apply_variable_uplift(fT dt, out_t tU)
	{
		auto U = DAGGER::format_input(tU);
		for (int i = 0; i < this->graph.nnodes; ++i) {
			if (!this->connector.flow_out_or_pit(i))
				this->topography[i] += U[i] * dt;
		}
	}

	void hydraulic_erosion_v0(int n_particules, fT erosor)
	{
		// old test
	}

	void normalise_topography()
	{
		fT mini = std::numeric_limits<fT>::max();
		fT maxi = std::numeric_limits<fT>::min();
		for (auto v : this->topography) {
			if (v < mini)
				mini = v;
			if (v > maxi)
				maxi = v;
		}
		maxi -= mini;
		for (auto& v : this->topography)
			v = (v - mini) / maxi;
	}
};

// generates a quick and dirty fluvial SFD topo for testing purposes
// uses a multigrid analytical solution to the SPL (see saleve algorithm from
// Steer (2022) on which I applied multigrid methods using random noise and
// cordonnier carving for the projections) final size is 16 * 2^ncycles
// boundaries are "4edges" "periodic_EW" or "periodic_NW"
template<class fT>
std::vector<fT>
_quick_fluvial_topo(int ncycles, std::string boundaries, int nrefine = 10)
{
	// init dx to get a final one = 50
	fT dx = std::pow(2, ncycles) * 50;
	// init connector and boundary conditions
	D8connector<fT> con(16, 16, dx, dx, 0, 0);
	con.set_default_boundaries(boundaries);
	// init graph
	graph<fT, D8connector<fT>> gf(con);
	// init Popscape
	popscape<fT, graph<fT, D8connector<fT>>, D8connector<fT>> psc(gf, con);
	// init randomnoise
	std::vector<fT> topo(16 * 16, 0.);
	add_noise_to_vector(topo, 0, 1);
	psc.set_topo(topo);

	for (int i = 0; i < ncycles + 1; ++i) {
		psc.StSt(nrefine);
		if (i < ncycles)
			psc.restriction(5);
	}

	return psc.topography;
}

template<class fT, class out_t>
out_t
quick_fluvial_topo(int ncycles, std::string boundaries)
{
	std::vector<fT> out = _quick_fluvial_topo<fT>(ncycles, boundaries);
	return DAGGER::format_output<std::vector<fT>, out_t>(out);
}

// End of namespace fastflood
}; // namespace DAGGER

#endif
