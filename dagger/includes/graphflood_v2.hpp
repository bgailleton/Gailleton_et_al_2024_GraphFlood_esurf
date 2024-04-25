#pragma once
#include "dodcontexts.hpp"
#include "enumutils.hpp"
#include "graphflood_enums.hpp"
#include "graphflood_parts.hpp"
#include "graphuncs.hpp"
#include "parambag.hpp"
#include "tinysubgraph.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, // type for integers (e.g. node index)
				 class f_t, // type for floating points number: float (saves memory,
										// more approximations for small numbers) or double
										// (reccomended)
				 class CONNECTOR_T, // type of connector (grid)
				 class GRAPH_T,			// Type of graph (might be deprecated)
				 class DBAG_T,			// type of data holder
				 class PARAM_T>			// type of param holder
/// Graphflood2 is the new version of graphflood calculator
/// It uses induced dynamic subgraph to offer another significant speed up
/// Getting stable, still afew mass balance issue for the segmented tsg
/// calculation as of 17/02/2024 B.G. 2023 - 2024
class Graphflood2
{

	// Everything is public. I know, not best practice. But watch me.
public:
	/*
	=+=+=+=+=+=+=+=+=+=+=+=+=+=
	=+=+=+=+=+=+=+=+=+=+=+=+=+=
	=+	ATTRIBUTES +=+=+=+=+=+=
	=+=+=+=+=+=+=+=+=+=+=+=+=+=
	=+=+=+=+=+=+=+=+=+=+=+=+=+=
	*/

	/// Ptr to connector (grid)
	CONNECTOR_T* con;
	/// Ptr to graph
	GRAPH_T* gra;
	/// Ptr to data holder
	DBAG_T* data;
	/// Ptr to parameter holder
	PARAM_T* param;

	// Should be moved to parameter bag
	/// Controls how water enters the model:
	/// --> WATER_INPUT::PRECIPITATIONS_CONSTANT:
	WATER_INPUT water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;
	f_t Prate = 1e-5; // mm.s^-1. Yarr.
	void set_uniform_P(f_t val)
	{
		this->water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;
		this->Prate = val;
	}

	// Attribute managing the inputs and starting points of the subgraphs
	/// initial nodes in the PQ
	std::vector<i_t> entry_node_PQ;
	/// External Qw/Qs input location
	std::vector<i_t> input_node_Qw;
	/// External Qw input value
	std::vector<f_t> input_Qw;
	/// External Qs input value
	std::vector<f_t> input_Qs;
	/// Tracks which nodes are in the PQ or not
	std::vector<std::uint8_t> isInQ;

	std::vector<f_t> xtraQwin;
	std::vector<f_t> xtraQsout;

	// compute options
	RUN_GF2 computor = RUN_GF2::NORMAL;

	f_t MB_Qwin_out = 0.;
	f_t get_MB_Qwin_out() { return this->MB_Qwin_out; }

	f_t meandhstar = 0.;
	f_t get_meandhstar() { return this->meandhstar; }

	std::vector<i_t> input_node_Qw_tsg;
	std::vector<f_t> input_Qw_tsg;

	TinySubGraph<i_t, f_t, CONNECTOR_T, DBAG_T, PARAM_T> tsg;

	std::vector<i_t> active_nodes;

	bool fillonly = false;
	SUBGRAPHMETHOD sbg_method = SUBGRAPHMETHOD::V1;

	// Water transfer params
	f_t min_part_Qw = 0.1;

	// time management thingies
	// std::vector<f_t> data->_timetracker; // MOVED TO THE DATA BAG
	f_t time = 0.;
	f_t dt = 1e-3;
	void set_dt(f_t val) { this->dt = val; }
	f_t get_dt() { return this->dt; }

	// Graphflood params
	f_t mannings = 0.033;
	f_t hw_increment_LM = 1e-3;

	bool prefill_lakes = true;

	int debugyolo = 0;

	// std::vector<i_t> active_nodes;

	/// Empty constructor
	Graphflood2(){};

	/// Classic constructor: links the graphflood module to existing constructor,
	/// graph ,... Arguments:
	/// ###--> Connector
	/// ###--> graph
	/// ###--> data bag
	/// ###--> param bag
	Graphflood2(CONNECTOR_T& con, GRAPH_T& gra, DBAG_T& data, PARAM_T& param)
	{
		this->con = &con;
		this->gra = &gra;
		this->data = &data;
		this->param = &param;
	}

	void init()
	{

		if (this->data == nullptr)
			throw std::runtime_error(
				"Cannot init Graphflood2 -> no data bag connected");
		if (this->con == nullptr)
			throw std::runtime_error(
				"Cannot init Graphflood2 -> no connector connected");
		// if(this->gra == nullptr) throw std::runtime_error("Cannot init
		// Graphflood2 -> no graph connected");

		if (this->data->_surface.size() == 0)
			throw std::runtime_error(
				"databag has no surface data, cannot init GraphFlood2");
		if (this->data->_hw.size() == 0)
			this->data->_hw = std::vector<f_t>(this->con->nxy(), 0);
		if (this->data->_Qwin.size() == 0)
			this->data->_Qwin = std::vector<f_t>(this->con->nxy(), 0);
		if (this->data->_Qwout.size() == 0)
			this->data->_Qwout = std::vector<f_t>(this->con->nxy(), 0);
		this->data->_vmot_hw = std::vector<f_t>(this->con->nxy(), 0);
		this->isInQ = std::vector<std::uint8_t>(this->con->nxy(), false);

		this->tsg = TinySubGraph<i_t, f_t, CONNECTOR_T, DBAG_T, PARAM_T>(
			*this->con, *this->data, *this->param);

		this->xtraQwin = std::vector<f_t>(this->con->nxy(), 0);
		this->xtraQsout = std::vector<f_t>(this->con->nxy(), 0);
	}

	void initial_fill()
	{
		// copying the bedrock elevation
		std::vector<f_t> bedrock(this->data->_surface);

		this->con->reinit();

		// filling and computing the graph from it
		this->con->PFcompute_all();
		if (this->prefill_lakes) {
			for (int i = 0; i < this->con->nxy(); ++i)
				this->data->_hw[i] += this->data->_surface[i] - bedrock[i];
		}

		// backpropagating the bedrock to topo
		this->data->_surface = bedrock;
	}

	void compute_entry_points_from_P(f_t Qw_threshold)
	{
		this->_computeEntryPoints_prec(Qw_threshold);
	}

	void _computeEntryPoints_prec(f_t Qw_threshold)
	{

		this->entry_node_PQ.clear();
		this->input_node_Qw.clear();
		this->input_Qs.clear();
		this->input_Qw.clear();

		// copying the bedrock elevation
		std::vector<f_t> bedrock(this->data->_surface);

		this->con->reinit();

		// filling and computing the graph from it
		this->con->PFcompute_all();
		if (this->prefill_lakes) {
			for (int i = 0; i < this->con->nxy(); ++i)
				this->data->_hw[i] += this->data->_surface[i] - bedrock[i];
		}

		// backpropagating the bedrock to topo
		this->data->_surface = bedrock;

		// temporary QWin
		std::vector<f_t> QW(this->con->_nxy, 0.);

		// Calculating QWin
		for (int i = this->con->nxy() - 1; i >= 0; --i) {
			int node = this->data->_Sstack[i];
			if (nodata(this->data->_boundaries[node]))
				;
			int rec = this->con->Sreceivers(node);
			QW[node] += Prate * this->con->area(node);
			if (node != rec) {
				QW[rec] += QW[node];
			}
		}

		std::vector<std::uint8_t> isdone(this->con->nxy(), false);

		fillvec(this->data->_Qwin, 0.);

		for (int i = this->con->nxy() - 1; i >= 0; --i) {
			int node = this->data->_Sstack[i];
			if (nodata(this->data->_boundaries[node]) || isdone[node])
				continue;

			int rec = this->con->Sreceivers(node);
			if (QW[rec] >= Qw_threshold && QW[node] < Qw_threshold &&
					isdone[node] == false) {
				this->input_node_Qw.emplace_back(node);
				this->input_Qw.emplace_back(QW[node]);
				this->input_Qs.emplace_back(0);
				this->entry_node_PQ.emplace_back(node);

				isdone[node] = true;

				// while (rec != node) {
				// 	isdone[rec] = true;
				// 	node = rec;
				// 	rec = this->con->Sreceivers(node);
				// }
				// } else if (QW[node] >= Qw_threshold) {
				// 	this->data->_Qwin[node] += Prate * this->con->area(node);

				// } else if (QW[rec] >= Qw_threshold) {
				// 	this->data->_Qwin[rec] += QW[node];
			}

			int nd = this->con->nSdonors(node);
			if (nd == 0 && QW[node] >= Qw_threshold && isdone[node] == false) {
				this->input_node_Qw.emplace_back(node);
				this->input_Qw.emplace_back(QW[node]);
				this->input_Qs.emplace_back(0);
				this->entry_node_PQ.emplace_back(node);

				isdone[node] = true;
			}
		}

		// for (int i = this->con->nxy() - 1; i >= 0; --i) {
		// 	if (this->data->_Qwin[i] > 0) {
		// 		this->input_node_Qw.emplace_back(i);
		// 		this->input_Qw.emplace_back(this->data->_Qwin[i]);
		// 	}
		// }

		// this->water_input_mode = WATER_INPUT::ENTRY_POINTS_QW;
	}

	void compute_entry_points_sources(f_t val)
	{
		this->entry_node_PQ.clear();
		this->input_node_Qw.clear();
		this->input_Qw.clear();

		this->con->PFcompute_all(false);

		for (int i = 0; i < this->con->nxy(); ++i) {
			if (this->is_node_active(i) == false)
				continue;
			int nd = this->con->nSdonors(i);
			if (nd == 0) {
				this->entry_node_PQ.emplace_back(i);
				this->input_node_Qw.emplace_back(i);
				this->input_Qw.emplace_back(val);
			}
		}
	}

	template<class arrin_i_t, class arrin_f_t>
	void set_Qw_input_points(arrin_i_t& tarri, arrin_f_t& tarrf)
	{
		auto arri = format_input<arrin_i_t>(tarri);
		this->input_node_Qw = to_vec(arri);
		this->entry_node_PQ = to_vec(arri);
		auto arrf = format_input<arrin_f_t>(tarrf);
		this->input_Qw = to_vec(arrf);
		this->input_Qs = std::vector<f_t>(this->input_Qw.size(), 0.);
		this->water_input_mode = WATER_INPUT::ENTRY_POINTS_QW;
	}

	void multiply_Qw_input_points(f_t factor)
	{
		for (auto& v : this->input_Qw) {
			v *= factor;
		}
	}
	void multiply_Qs_input_points(f_t factor)
	{
		for (auto& v : this->input_Qs) {
			v *= factor;
		}
	}
	f_t sum_Qw_inputs()
	{
		f_t susum = 0.;
		for (auto& v : this->input_Qw) {
			susum += v;
		}
		return susum;
	}

	// Small getters
	std::vector<i_t> get_entry_node_PQ() { return entry_node_PQ; };
	std::vector<i_t> get_active_nodes() { return active_nodes; };
	std::vector<i_t> get_input_node_Qw() { return input_node_Qw; };
	std::vector<f_t> get_input_Qw() { return input_Qw; };
	std::vector<f_t> get_input_Qs() { return input_Qs; };

	template<class arrin_i_t, class arrin_f_t>
	void set_QwQs_input_points(arrin_i_t& tarri,
														 arrin_f_t& tarrf,
														 arrin_f_t& tarrfqs)
	{
		auto arri = format_input<arrin_i_t>(tarri);
		this->input_node_Qw = to_vec(arri);
		this->entry_node_PQ = to_vec(arri);
		auto arrf = format_input<arrin_f_t>(tarrf);
		this->input_Qw = to_vec(arrf);
		auto arrfqs = format_input<arrin_f_t>(tarrfqs);
		this->input_Qs = to_vec(arrfqs);
	}

	template<class CELL, class DSTACK>
	void init_dstack(DSTACK& dstack)
	{

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			this->xtraQwin[this->input_node_Qw[i]] += this->input_Qw[i];
		}

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			if (this->param->gf2_morpho)
				this->data->_Qsin[this->input_node_Qw[i]] += this->input_Qs[i];
		}

		for (size_t i = 0; i < this->entry_node_PQ.size(); ++i) {
			dstack.emplace(CELL(this->entry_node_PQ[i],
													this->data->_surface[this->entry_node_PQ[i]]));
			this->isInQ[this->entry_node_PQ[i]] = true;
		}

		if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT ||
				this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE) {
			for (int i = 0; i < this->con->nxy(); ++i) {
				if (this->is_node_active(i)) {
					this->xtraQwin[i] += this->Prate * this->con->area(i);
				}
			}
		}
	}

	void add_init_flux()
	{

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			this->xtraQwin[this->input_node_Qw[i]] += this->input_Qw[i];
		}

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			if (this->param->gf2_morpho)
				this->data->_Qsin[this->input_node_Qw[i]] += this->input_Qs[i];
		}
	}

	void add_init_flux(std::vector<f_t>& tqwqw, std::vector<f_t>& tqsqs)
	{

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			tqwqw[this->input_node_Qw[i]] += this->input_Qw[i];
		}

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			if (this->param->gf2_morpho)
				tqsqs[this->input_node_Qw[i]] += this->input_Qs[i];
		}
	}

	template<class CELL, class DSTACK>
	void init_dstack_tsbdyn(DSTACK& dstack)
	{

		for (size_t i = 0; i < this->input_node_Qw_tsg.size(); ++i) {
			this->data->_Qwin[this->input_node_Qw_tsg[i]] = this->input_Qw_tsg[i];
		}

		// for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
		// 	if (this->param->gf2_morpho)
		// 		this->data->_Qsin[this->input_node_Qw[i]] = this->input_Qs[i];
		// }

		for (size_t i = 0; i < this->input_node_Qw_tsg.size(); ++i) {
			dstack.emplace(CELL(this->input_node_Qw_tsg[i],
													this->data->_surface[this->input_node_Qw_tsg[i]]));
			this->isInQ[this->input_node_Qw_tsg[i]] = true;
		}
	}

	template<class CELL, class DSTACK>
	void init_dstack_WaterSed(DSTACK& dstack)
	{

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			dstack.emplace(CELL(this->input_node_Qw[i],
													this->data->_surface[this->input_node_Qw[i]],
													this->input_Qw[i],
													this->input_Qs[i]));
		}
	}

	void fillrun_subgraphflood()
	{
		this->sbg_method == SUBGRAPHMETHOD::FILLONLY;
		this->run_subgraphflood();
		this->sbg_method = SUBGRAPHMETHOD::V1;
	}

	void compute_qr()
	{

		this->computor = RUN_GF2::COMPUTEqr;
		this->run_subgraphflood();
		this->computor = RUN_GF2::NORMAL;
	}

	void run()
	{
		if (this->param->gf2_morpho == false ||
				this->param->morphomode == MORPHOMODE::NONE)
		{	
			if(this->param->hydromode == HYDROMODE::MFD)
				this->run_hydro_mfd();
			else if (this->param->hydromode == HYDROMODE::SFD)
				this->run_hydro_sfd();
		}
		else if(this->param->hydromode == HYDROMODE::SFD && this->param->morphomode == MORPHOMODE::MPM){
			this->run_morpho_sfd();
		}

		else {

			this->run_subgraphflood();
		}
	}

	void run_subgraphflood()
	{

		std::vector<f_t> qr_out;
		if (this->computor == RUN_GF2::COMPUTEqr)
			qr_out = std::vector<f_t>(this->con->nxy(), 0.);

		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_debug = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_theta_flow_in = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_theta_flow_out = std::vector<f_t>(this->con->nxy(), 0.);

			this->data->fbag["rc"] = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->fbag["latel"] = std::vector<f_t>(this->con->nxy(), 0.);
		}

		if (this->param->gf2_morpho && this->data->_Qsin.size() == 0) {
			this->data->_Qsin = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_Qsout = std::vector<f_t>(this->con->nxy(), 0.);
		}

		std::priority_queue<WaCell<i_t, f_t>,
												std::vector<WaCell<i_t, f_t>>,
												std::less<WaCell<i_t, f_t>>>
			dynastack;

		// ocarina nook;

		// nook.tik();

		if (this->active_nodes.size() == 0) {

			fillvec(this->data->_Qwin, 0.);
			fillvec(this->data->_theta_flow_in, 0.);
			fillvec(this->data->_theta_flow_out, 0.);
			fillvec(this->xtraQwin, 0.);
			fillvec(this->isInQ, false);
			if (this->param->gf2_morpho) {
				fillvec(this->data->_Qsin, 0.);
				fillvec(this->data->_Qsout, 0.);
				fillvec(this->xtraQsout, 0.);
			}
		} else {

			fillvec(this->data->_Qwin, 0., this->active_nodes);
			fillvec(this->xtraQwin, 0., this->active_nodes);
			fillvec(this->data->_theta_flow_in, 0., this->active_nodes);
			fillvec(this->data->_theta_flow_out, 0., this->active_nodes);
			fillvec(this->isInQ, false, this->active_nodes);
			if (this->param->gf2_morpho) {
				fillvec(this->data->_Qsin, 0., this->active_nodes);
				fillvec(this->data->_Qsout, 0., this->active_nodes);
				fillvec(this->xtraQsout, 0., this->active_nodes);
			}
		}

		std::vector<int> labels(this->con->nxy(), 0);

		this->MB_Qwin_out = 0.;
		this->meandhstar = 0.;

		this->init_dstack<WaCell<i_t, f_t>, decltype(dynastack)>(dynastack);

		CT_neighbours<i_t, f_t> ctx;

		this->time += this->dt;

		std::vector<i_t> tempnodes, temprec;
		tempnodes.reserve(this->con->_nx);
		temprec.reserve(this->con->_nx);

		std::array<i_t, 8> receivers;
		std::array<f_t, 8> receiversWeights;
		std::array<f_t, 8> receiversWeightsQs;
		std::array<f_t, 8> receiversSlopes;

		f_t sumout = 0.;

		// nook.tok("init took ");
		// nook.tik();

		while (dynastack.empty() == false) {

			// Getting the next node
			auto next = this->_dstack_next<WaCell<i_t, f_t>>(dynastack);

			// deregistering it as being in the PQ stack
			this->isInQ[next.node] = false;

			// ispast is true if the node has not been processed yet
			bool ispast = this->data->_timetracker[next.node] != this->time;

			// Updating the timer
			this->data->_timetracker[next.node] = this->time;

			// std::cout << next.node << std::endl;
			ctx.update(next.node, *this->con);
			;

			++labels[next.node];

			// std::cout << BC2str(ctx.boundary) << std::endl;

			// if (ispast == false) {
			// 	// this->data->_Qwin[next.node] =
			// std::max(this->data->_Qwin[next.node],
			// 	// next.Qw);
			// 	// if (this->data->_Qwin[next.node] > next.Qw) {
			// 		// this->data->_Qwin[next.node] += next.Qw;
			// 	// } else {
			// 		this->data->_Qwin[next.node] = next.Qw;
			// 	// }
			// }

			if (can_out(ctx.boundary)) {
				// if (ispast == false)
				// 	this->MB_Qwin_out += this->data->_Qwin[ctx.node];
				continue;
			}

			if (nodata(ctx.boundary)) {
				// std::cout << "nodata reached" << std::endl;
				continue;
			}

			int nr = 0;
			int SSi = 0;
			f_t SS = 0;
			f_t SSdx = 1.;
			f_t SSdy = 1.;

			this->update_receivers(ctx,
														 receivers,
														 receiversWeights,
														 nr,
														 SS,
														 SSdx,
														 SSdy,
														 SSi,
														 receiversSlopes);

			if (ispast) {
				temprec.emplace_back(nr >= 0 ? SSi : next.node);
				tempnodes.emplace_back(next.node);
			}

			f_t& thw = this->data->_hw[next.node];
			f_t& tsurf = this->data->_surface[next.node];

			for (int j = 0; j < ctx.nn; ++j) {
				int tn = ctx.neighbours[j];
				if (this->data->_surface[tn] > tsurf && this->data->_hw[tn] > 0 &&
						this->data->_timetracker[tn] != this->time) {
					dynastack.emplace(
						WaCell<i_t, f_t>(tn, this->data->_surface[tn], 0., 0.));
				}
			}

			// Actual flux calculation

			f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
								std::sqrt(std::max(1e-6, SS));
			f_t tQwout = thw * u_w * SSdy;

			this->data->_Qwout[ctx.node] = tQwout;

			// Flux in coprrection
			this->data->_Qwin[ctx.node] = next.Qw + this->xtraQwin[ctx.node];
			this->xtraQwin[ctx.node] = 0.;

			// Calculation of theta out
			f_t vx = 0., vy = 0.;
			for (int j = 0; j < nr; ++j) {
				int ij = receivers[j];
				int rec = ctx.neighbours[ij];
				f_t tvx, tvy;
				tvx = ctx.neighboursDx[ij] * ctx.neighboursDxSign[ij];
				tvy = ctx.neighboursDy[ij] * ctx.neighboursDySign[ij];
				vx += receiversWeights[j] * tvx;
				vy += receiversWeights[j] * tvy;

				f_t ttheta;
				VEC2D::cart2pol(tvx, tvy, ttheta);
				if (this->data->_Qwin[rec] + this->xtraQwin[rec] != 0)
					this->data->_theta_flow_in[rec] =
						VEC2D::mean_theta(this->data->_Qwin[rec] + this->xtraQwin[rec],
															this->data->_theta_flow_in[rec],
															receiversWeights[j] * this->data->_Qwin[ctx.node],
															ttheta);
				else {
					this->data->_theta_flow_in[rec] = ttheta;
				}
				// if(nr>1)
				// 	std::cout << "vx=" << tvx << "|" << " vy=" << tvy << "|" << ttheta
				// << "|" << this->data->_theta_flow_in[rec] << std::endl;
			}

			VEC2D::cart2pol(vx, vy, this->data->_theta_flow_out[ctx.node]);

			// if(nr>1)
			// 	std::cout << "vx=" << vx << "|" << " vy=" << vy << "|" <<  <<
			// std::endl;

			f_t radius =
				VEC2D::radius_of_curvature_theta(this->data->_theta_flow_in[ctx.node],
																				 this->data->_theta_flow_out[ctx.node]);
			// std::cout << this->data->_theta_flow_in[ctx.node] << "|" <<
			// this->data->_theta_flow_out[ctx.node]<<"|"<<rc << std::endl;
			// this->data->fbag["rc"][ctx.node] = rc;

			// TEST IN PROGRESS TO ASSESS VECTOR REP
			// if(true && nr > 1 && ctx.canout == false && ispast){
			// 	i_t n1, n2;
			// 	f_t w1, w2;
			// 	con->NeighboursTheta2(
			// 		ctx.node, this->data->_theta_flow_out[ctx.node] * (1 +
			// this->data->randu->get() * 0.4 - 0.2), n1, n2, w1, w2);
			// 	if(can_receive(this->data->_boundaries[n1]) &&
			// can_receive(this->data->_boundaries[n2]))
			// 	{
			// 		nr = 2;
			// 		receivers[0] = 0;
			// 		receivers[1] = 1;
			// 		ctx.neighbours[0] = n1;
			// 		ctx.neighbours[1] = n2;
			// 		receiversWeights[0] = w1;
			// 		receiversWeights[1] = w2;
			// 	}
			// }

			if (this->computor == RUN_GF2::COMPUTEqr)
				qr_out[ctx.node] = thw * u_w;

			// Manage the morphodynamics of the model
			if (this->param->gf2_morpho) {

				// Calculating shear stress
				f_t tau =
					0.; // = this->param->rho_water * this->param->GRAVITY * SS * thw;

				if (this->param->morphomode == MORPHOMODE::EROS) {
					tau = this->param->rho_water * this->param->GRAVITY * SS * thw;

					// rates
					// # basal erosion
					f_t edot = 0.;
					// # lateral erosion
					f_t eldot = 0;
					// # basal dep
					f_t ddot = 0.;
					// # lateral dep
					f_t dldot = 0.;

					// erosion only happens if critical shear stress increases
					if (tau > this->param->tau_c) {

						// calculating basal erosion: ke x (shear - critshear) ^ alpha
						edot = this->param->ke(ctx.node) *
									 std::pow(tau - this->param->tau_c, this->param->alpha);

						// lateral erosion
						if (this->param->bank_erosion) {
							// Getting the highest neighbour to erode
							int hni = ctx.idxHighestNeighbour(
								this->data->_surface, this->data->_hw, SSi);

							// hni is -1 if there is no higher neighbours
							if (hni >= 0) {

								// neihbour ID
								int hn = ctx.neighbours[hni];

								// difference in substrate elevation defining the max lateral
								// erosion rate
								f_t delta_Z =
									(this->data->_surface[hn] - this->data->_hw[hn] -
									 this->data->_surface[ctx.node] + this->data->_hw[ctx.node]);

								// Lateral gradient
								f_t latslope = delta_Z / ctx.neighboursDx[hni];

								// Actual rate
								eldot = this->param->kel * latslope * edot;

								// checking we do not erode too much
								if (eldot * this->dt > delta_Z)
									eldot = delta_Z / this->dt;

								// applyin' stuff
								this->xtraQsout[hn] += eldot * this->con->area(hn);
							}
						}
					}

					// TEMP DEACTIVATION
					// ddot =
					// 	this->data->_Qsin[ctx.node] / (this->param->kd * SSdy); //;

					f_t latslope = 0.;
					if (this->param->kd > 0) {
						int hni = ctx.idxLowestNeighbour(
							this->data->_surface, this->data->_hw, SSi);

						// hni is -1 if there is no higher neighbours
						if (hni >= 0) {
							// neihbour ID
							int hn = ctx.neighbours[hni];

							// difference in substrate elevation defining the max lateral
							// erosion rate
							f_t delta_Z =
								(this->data->_surface[ctx.node] - this->data->_hw[ctx.node] -
								 this->data->_surface[hn] + this->data->_hw[hn]);

							f_t latslope = delta_Z / ctx.neighboursDx[hni];

							dldot = latslope * this->param->kdl *
											this->data->_Qsin[ctx.node] / SSdy;

							this->xtraQsout[hn] += dldot * this->con->area(hn);
						}
					}

					double K = (1. / this->param->kd + this->param->kdl * latslope);

					double edotpsy = (edot + eldot) / K;
					double C1 = this->data->_Qsin[ctx.node] / SSdy - edotpsy;

					this->data->_Qsout[ctx.node] =
						SSdy * (edotpsy + C1 * std::exp(-SSdx * K));

					if (this->data->_Qsout[ctx.node] < 0)
						throw std::runtime_error("Qso < 0");
				} else if (this->param->morphomode == MORPHOMODE::MPM) {

					f_t correction_factor = 0., sumtaus = 0.;

					for (int j = 0; j < nr; ++j) {
						receiversWeightsQs[j] = 0.;
						f_t ttau = this->param->rho_water * this->param->GRAVITY *
											 receiversSlopes[j] * thw;
						ttau = std::pow(ttau - this->param->tau_c, 1.5) *
									 ctx.neighboursDy[receivers[j]];

						if (ttau > this->param->tau_c) {
							tau += ttau;
							ttau *= this->data->randu->get();
							sumtaus += ttau;
							receiversWeightsQs[j] = ttau;
						}
					}

					if (tau > 0) {
						for (int j = 0; j < nr; ++j) {
							receiversWeightsQs[j] /= sumtaus;
						}

						correction_factor = std::pow(this->param->rho_water *
																						 this->param->GRAVITY * SS * thw -
																					 this->param->tau_c,
																				 1.5) *
																SSdy / tau;
					}

					f_t capacity = 0;

					// Alrady tau - tau_c ed!
					capacity = this->param->E * tau;
					this->data->_Qsout[ctx.node] = capacity;

				} else if (this->param->morphomode == MORPHOMODE::MPMVEC) {

					tau = this->param->rho_water * this->param->GRAVITY * SS * thw;
					f_t capacity = 0.;
					if (tau > this->param->tau_c) {
						capacity =
							std::pow(tau - this->param->tau_c, 1.5) * this->param->E * SSdy;
						i_t n1, n2;
						f_t w1, w2;
						con->NeighboursTheta2(ctx.node,
																	this->data->_theta_flow_out[ctx.node] *
																		(1 + this->data->randu->get() * 0.4 - 0.2),
																	n1,
																	n2,
																	w1,
																	w2);
						this->data->_Qsin[n1] += capacity * w1;
						this->data->_Qsin[n2] += capacity * w2;
						this->data->_Qsout[ctx.node] = capacity;

						// if(this->data->_surface[n1] > this->data->_surface[ctx.node])
						// std::cout << "ping" << std::endl;
					}
				} else if (this->param->morphomode == MORPHOMODE::EROSVEC) {
					tau = this->param->rho_water * this->param->GRAVITY * SS * thw;

					// rates
					// # basal erosion
					f_t edot = 0.;
					// # lateral erosion
					f_t eldot = 0;
					// # basal dep
					f_t ddot = 0.;
					// # lateral dep
					f_t dldot = 0.;

					// erosion only happens if critical shear stress increases
					if (tau > this->param->tau_c) {

						// calculating basal erosion: ke x (shear - critshear) ^ alpha
						edot = this->param->ke(ctx.node) *
									 std::pow(tau - this->param->tau_c, this->param->alpha);
					}

					// Lateral wrosion from Langston and Tucker
					if (this->data->_theta_flow_in[ctx.node] != 0 &&
							this->param->kel > 0) {
						f_t theta_lat, rlat, wl1, wl2, x1, x2, y1, y2, xr, yr, xn, yn;
						i_t nl1, nl2;

						VEC2D::pol2cart(this->data->_theta_flow_in[ctx.node], x1, y1);
						VEC2D::pol2cart(this->data->_theta_flow_out[ctx.node], x2, y2);

						// VEC2D::add_pol(1.,this->data->_theta_flow_in[ctx.node],1.,this->data->_theta_flow_out[ctx.node],
						// rlat, theta_lat);
						VEC2D::add(x1, y1, x2, y2, xr, yr);

						f_t cross = VEC2D::cross(x1, y1, x2, y2);

						VEC2D::normal(xr, yr, xn, yn, cross > 0);

						VEC2D::cart2pol(xn, yn, rlat, theta_lat);

						con->NeighboursTheta2(ctx.node, theta_lat, nl1, nl2, wl1, wl2);

						f_t omega_c;
						if (radius > 1e-6)
							this->param->rho_water* SSdy* std::pow(u_w, 3) / radius;
						f_t teldot = this->param->kel * this->param->ke(ctx.node) * omega_c;

						f_t tbed =
							this->data->_surface[ctx.node] - this->data->_hw[ctx.node];

						if (this->data->_surface[nl1] - this->data->_hw[nl1] > tbed) {
							f_t tteldot = std::min(
								wl1 * teldot,
								(this->data->_surface[nl1] - this->data->_hw[nl1] - tbed) /
									(this->dt * this->param->time_dilatation_morpho));
							this->xtraQsout[nl1] += tteldot * this->con->area(ctx.node);
							eldot += tteldot;
							this->data->fbag["latel"][nl1] += tteldot;
						}
						if (this->data->_surface[nl2] - this->data->_hw[nl2] > tbed) {
							f_t tteldot = std::min(
								wl2 * teldot,
								(this->data->_surface[nl2] - this->data->_hw[nl2] - tbed) /
									(this->dt * this->param->time_dilatation_morpho));
							this->xtraQsout[nl2] += tteldot * this->con->area(ctx.node);
							eldot += tteldot;
							this->data->fbag["latel"][nl2] += tteldot;
						}
					}

					this->data->fbag["rc"][ctx.node] = radius;

					f_t latslope = 0.;

					double K = (1. / this->param->kd + this->param->kdl * latslope);

					double edotpsy = (edot + eldot) / K;
					double C1 = this->data->_Qsin[ctx.node] / SSdy - edotpsy;

					this->data->_Qsout[ctx.node] =
						SSdy * (edotpsy + C1 * std::exp(-SSdx * K));

					i_t n1, n2;
					f_t w1, w2;
					con->NeighboursTheta2(ctx.node,
																this->data->_theta_flow_out[ctx.node] *
																	(1 + this->data->randu->get() * 0.4 - 0.2),
																n1,
																n2,
																w1,
																w2);
					this->data->_Qsin[n1] += this->data->_Qsout[ctx.node] * w1;
					this->data->_Qsin[n2] += this->data->_Qsout[ctx.node] * w2;

					if (this->data->_Qsout[ctx.node] < 0)
						throw std::runtime_error("Qso < 0");
				}
			}

			f_t baseQw = this->data->_Qwin[ctx.node];
			f_t baseQs;
			if (this->param->gf2_morpho)
				baseQs = (ispast) ? this->data->_Qsout[next.node] : next.Qs;
			// baseQs += this->extra

			if (this->param->morphomode == MORPHOMODE::EROS)
				receiversWeightsQs = receiversWeights;
			// if(this->param->morphomode != MORPHOMODE::NONE) receiversWeightsQs =
			// receiversWeights;

			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				int rec = ctx.neighbours[i];
				if (receiversWeights[j] == 0 && receiversWeightsQs[j] == 0)
					continue;

				if (nodata(this->data->_boundaries[rec]))
					throw std::runtime_error("fasdgfdgkhfgkdhsfg");

				bool tizdone = this->data->_timetracker[rec] == this->time;

				if (tizdone) {
					if (this->param->gf2_morpho &&
							this->param->morphomode != MORPHOMODE::MPMVEC &&
							this->param->morphomode != MORPHOMODE::EROSVEC) {
						if (this->param->TSG_dist == false || this->tsg.xtraMask[rec]) {
							if (this->isInQ[rec])
								this->xtraQwin[rec] += baseQw * receiversWeights[j];
							else
								dynastack.emplace(
									WaCell<i_t, f_t>(rec,
																	 this->data->_surface[rec],
																	 baseQw * receiversWeights[j],
																	 baseQs * receiversWeightsQs[j]));

							this->isInQ[rec] = true;
						} else
							this->MB_Qwin_out += baseQw * receiversWeights[j];
					} else {
						if (this->param->TSG_dist == false || this->tsg.xtraMask[rec]) {
							if (this->isInQ[rec])
								this->xtraQwin[rec] += baseQw * receiversWeights[j];
							else
								dynastack.emplace(
									WaCell<i_t, f_t>(rec,
																	 this->data->_surface[rec],
																	 baseQw * receiversWeights[j]));

							this->isInQ[rec] = true;
						} else
							this->MB_Qwin_out += baseQw * receiversWeights[j];
					}

				} else {

					if (this->isInQ[rec] == false) {

						if (this->param->TSG_dist == false || this->tsg.xtraMask[rec]) {
							dynastack.emplace(
								WaCell<i_t, f_t>(rec, this->data->_surface[rec]));
							this->isInQ[rec] = true;
						} else
							this->MB_Qwin_out += baseQw * receiversWeights[j];
					}

					this->xtraQwin[rec] += baseQw * receiversWeights[j];
					if (this->param->gf2_morpho &&
							this->param->morphomode != MORPHOMODE::MPMVEC &&
							this->param->morphomode != MORPHOMODE::EROSVEC)
						this->data->_Qsin[rec] += baseQs * receiversWeightsQs[j];
				}
			}
		}

		if (this->computor == RUN_GF2::COMPUTEqr)
			this->data->fbag["qr_out"] = std::move(qr_out);

		if (this->computor != RUN_GF2::NORMAL)
			return;

		// nook.tok("PQ took ");
		// nook.tik();

		// if(this->param->gf2_morpho && this->param->morphomode ==
		// MORPHOMODE::MPM){

		// 	for(size_t i=0;i<tempnodes.size(); ++i){
		// 		int node = tempnodes[i];
		// 		this->data->_Qsin[node] += this->data->_Qsout[node];
		// 	}
		// 	for(size_t i=0;i<tempnodes.size(); ++i){
		// 		int node = tempnodes[i];
		// 		int rec  = temprec[i];

		// 		this->data->_Qsout[node] = this->data->_Qsin[rec];
		// 	}

		// }

		// OLD WAY OF CHECKING DISCONNECTED NODES
		// CT_neighbourer_WaCell<i_t, f_t> ctx2;

		// for (int i = 0; i < this->con->nxy(); ++i) {

		// 	if (this->active_nodes.size() > 0 && i >= this->active_nodes.size())
		// 		break;

		// 	int node = (this->active_nodes.size() > 0) ? this->active_nodes[i] : i;

		// 	this->data->_Qwin[node] += this->xtraQwin[node];
		// 	this->xtraQwin[node] = 0.;

		// 	if (this->is_node_active(node) == false )
		// 		continue;

		// 	if (this->data->_timetracker[node] < this->time && this->data->_hw[node]
		// > 0) {

		// 		this->data->_timetracker[node] = this->time;

		// 		this->data->_Qwin[node] = 0.;

		// 		ctx.update(node, *this->con);

		// 		this->_calculate_Qwout_for_disconnected_nodes<decltype(ctx)>(ctx);

		// 		++labels[node];

		// 	}

		// }

		int NN = 0;

		for (int i = 0; i < tempnodes.size(); ++i) {

			// if (this->active_nodes.size() > 0 && i >= this->active_nodes.size())
			// 	break;

			// int node = (this->active_nodes.size() > 0) ? this->active_nodes[i] : i;

			int node = tempnodes[i];

			if (can_out(this->data->_boundaries[node])) {
				this->MB_Qwin_out += this->data->_Qwin[node];
			}

			if (this->is_node_active(node) == false)
				continue;

			f_t& thw = this->data->_hw[node];
			f_t& tsurf = this->data->_surface[node];

			f_t dhw = 0.;
			if (this->sbg_method == SUBGRAPHMETHOD::V1) {
				dhw = this->data->_Qwin[node] -
							(this->data->_Qwout[node] * this->param->capacityFacQw);
			} else if (this->sbg_method == SUBGRAPHMETHOD::FILLONLY) {
				dhw = this->data->_Qwin[node];
			}

			dhw *= this->dt;
			dhw /= this->con->area(i);

			this->meandhstar += dhw / this->dt;
			++NN;

			thw += dhw;
			tsurf += dhw;

			if (this->param->gf2_morpho) {
				this->data->_Qsout[node] += xtraQsout[node];
				f_t dz = this->data->_Qsin[node] -
								 (this->data->_Qsout[node] * this->param->capacityFacQs);
				dz *= this->dt * this->param->time_dilatation_morpho;
				dz /= this->con->area(node);
				tsurf += dz;
			}

			if (std::isfinite(tsurf) == false) {
				std::cout << this->data->_Qwout[node] << "|" << this->data->_Qwin[node]
									<< "|" << this->data->_Qsin[node] << "|"
									<< this->data->_Qsout[node] << std::endl;
				;
				throw std::runtime_error("blug");
			}

			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}
		}

		// nook.tok("post took ");

		if (NN > 0)
			this->meandhstar /= NN;

		this->data->ibag["labels"] = labels;

		// std::cout << "Sumout: " << sumout << std::endl; ;
	}

	void run_hydro_mfd()
	{

		// ###################################
		// Initialising the data structures if needed
		// ###################################

		//Timer stuff ignore
		ocarina nook;
		// nook.tik();

		// In case I am computing the discharge per unit width
		std::vector<f_t> qr_out;
		if (this->computor == RUN_GF2::COMPUTEqr)
			qr_out = std::vector<f_t>(this->con->nxy(), 0.);

		// In case it is the first run
		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
		}

		// The priority queue is the dynamic stack popping/storing nodes in the
		// right way
		std::priority_queue<WaCell<i_t, f_t>,
												std::vector<WaCell<i_t, f_t>>,
												std::less<WaCell<i_t, f_t>>>
			dynastack;

		// ctx is the context neighbourer (helps nabigating through neighbours as
		// the stack is not computed)
		CT_neighbours<i_t, f_t> ctx;

		// Creating local stacks to speed up some operations
		std::vector<i_t> tempnodes, temprec;
		tempnodes.reserve(this->con->_nx);
		temprec.reserve(this->con->_nx);

		// placeholder for receiver data
		std::array<i_t, 8> receivers; // id of rec in the contextual neighbours
		std::array<f_t, 8> receiversWeights;
		std::array<f_t, 8> receiversWeightsQs;
		std::array<f_t, 8> receiversSlopes;

		// ###################################
		// Reinitialising the vectors to 0 if not furst run
		// ###################################

		// If active_nodes is activated, I only work on a subset of nodes
		if (this->active_nodes.size() == 0) {

			fillvec(this->data->_Qwin, 0.);
			fillvec(this->xtraQwin, 0.);
			fillvec(this->isInQ, false);
			// else I work on all the nodes
		} else {
			fillvec(this->data->_Qwin, 0., this->active_nodes);
			fillvec(this->xtraQwin, 0., this->active_nodes);
			fillvec(this->isInQ, false, this->active_nodes);
		}

		// ###################################
		// State variables
		// ###################################

		// Mass balance checker
		this->MB_Qwin_out = 0.;
		// increment water checker
		this->meandhstar = 0.;
		f_t sumout = 0.;

		// ###################################
		// Step 1: Preparing the inputs
		// ###################################

		// initialising the dynamic stack and the input points of water
		this->init_dstack<WaCell<i_t, f_t>, decltype(dynastack)>(dynastack);

		// Incrementing the timer
		this->time += this->dt;

		// ###################################
		// Step 2: Traversing the landscape
		// ###################################

		// Timer stuff
		// nook.tok("init took ");
		// nook.tik();

		// the main loop is running as long as there are still nodes in the priority
		// queue
		while (dynastack.empty() == false) {

			// nook.tik();

			// Getting the next node
			auto next = this->_dstack_next<WaCell<i_t, f_t>>(dynastack);

			// deregistering it as being in the PQ stack
			this->isInQ[next.node] = false;

			// ispast is true if the node has not been processed yet
			bool ispast = this->data->_timetracker[next.node] != this->time;

			// Updating the timer
			this->data->_timetracker[next.node] = this->time;

			// Updating the context (fetching local neighbours index and otehr info)
			ctx.update(next.node, *this->con);

			// If the node is outletting the model I skip it
			if (can_out(ctx.boundary)) {
				continue;
			}

			// If the node is no_data I skip it
			if (nodata(ctx.boundary)) {
				continue;
			}

			// ###################################
			// Local state variables
			// ###################################

			// Number of receivers
			int nr = 0;
			// Index of the steepest receiver
			int SSi = 0;
			// Steepest slope
			f_t SS = 0;
			// dx in the steepest direction
			f_t SSdx = 1.;
			// dy in the steepest direction
			f_t SSdy = 1.;

			// Fectching local hw and hydraulic surface
			f_t& thw = this->data->_hw[next.node];
			f_t& tsurf = this->data->_surface[next.node];

			// nook.tok("init took ");
			// nook.tik();


			// ###################################
			// Computing the local graph
			// ###################################

			// This function takes care of filtering, selecting and calculating the
			// receivers' characteristics
			this->update_receivers(ctx,
														 receivers,
														 receiversWeights,
														 nr,
														 SS,
														 SSdx,
														 SSdy,
														 SSi,
														 receiversSlopes);

			// nook.tok("update took ");
			// nook.tik();

			// Local stack
			if (ispast) {
				temprec.emplace_back(nr >= 0 ? SSi : next.node);
				tempnodes.emplace_back(next.node);
			}

			// nook.tok("locstack took ");
			// nook.tik();

			// Adding unprocessed upstream neighbours to the Queue to connect
			// disconnected node
			for (int j = 0; j < ctx.nn; ++j) {
				int tn = ctx.neighbours[j];
				if (this->data->_surface[tn] > tsurf && this->data->_hw[tn] > 0 &&
						this->data->_timetracker[tn] != this->time &&
						this->isInQ[tn] == false) {
					dynastack.emplace(
						WaCell<i_t, f_t>(tn, this->data->_surface[tn], 0., 0.));
				}
			}

			// nook.tok("disco took ");
			// nook.tik();

			// ###################################
			// Computing Water fluxes
			// ###################################

			// Flow velocity u using manning Equation
			// Note I am recasting the slope to a minimum value to avoid numerical
			// instabilities
			f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
								std::sqrt(std::max(1e-6, SS));

			// Converting to output volumetric discharge
			f_t tQwout = thw * u_w * SSdy;

			// Registering local output Volumetric discharge
			this->data->_Qwout[ctx.node] = tQwout;

			// Calculating the input Volumetric Discharge
			// Some water is transmitted via the dynamic cells (for example if they
			// have been in a local minima) rest of the water is in a temp vector
			// storing them
			this->data->_Qwin[ctx.node] = next.Qw + this->xtraQwin[ctx.node];
			this->xtraQwin[ctx.node] = 0.;

			// nook.tok("calc took ");
			// nook.tik();

			// ###################################
			// Optional Computations
			// ###################################

			if (this->computor == RUN_GF2::COMPUTEqr)
				qr_out[ctx.node] = thw * u_w;

			// ###################################
			// Transferring water to the neighbours
			// ###################################

			// For clarity
			f_t baseQw = this->data->_Qwin[ctx.node];

			// Iterating through the receivers
			for (int j = 0; j < nr; ++j) {

				// local id in the array of neighbours in the context
				int i = receivers[j];
				// actual index of the receivers
				int rec = ctx.neighbours[i];

				// Skipping if 0 to avoid  adding useless nodes to the Q
				if (receiversWeights[j] == 0 && receiversWeightsQs[j] == 0)
					continue;

				// Debug checker catching cases where I have a receivers that shoud not
				// be one. Keeping it cause I am still playing with receivers selection
				if (nodata(this->data->_boundaries[rec]))
					throw std::runtime_error(
						"Fraphflood2::run::run_hydro_mfd::receiver is no data");

				// Checking if my receiver has already been processed or not (Am I in
				// the process of solving a LM)
				bool tizdone = this->data->_timetracker[rec] == this->time;
				// Checking if my receiver is already in the Q
				bool inDaQ = this->isInQ[rec];

				if (inDaQ) {
					// Then I add the extra water to the bucket vector as it will come
					// later
					this->xtraQwin[rec] += baseQw * receiversWeights[j];
				}
				// If this is the sace
				else if (tizdone) {
					// Emplcing the cell and storing the watter within to be transmitted
					// out of the area
					dynastack.emplace(WaCell<i_t, f_t>(
						rec, this->data->_surface[rec], baseQw * receiversWeights[j]));

					// node is in the Q
					this->isInQ[rec] = true;

				} else {
					// Otherwise I just need to check if I don't fall out of the node
					// masks (Is it a part of the landscape I want to process)
					if (this->param->TSG_dist == false || this->tsg.xtraMask[rec]) {
						dynastack.emplace(WaCell<i_t, f_t>(
							rec, this->data->_surface[rec], baseQw * receiversWeights[j]));
						this->isInQ[rec] = true;
					}
				}
			}
			// nook.tok("transf took ");
			// nook.tik();
		}

		// #######################################
		//  End of the traversal loop
		//  Beginning of the post-process stage
		// #######################################

		// Registering the discharge per unit width if required to compute
		if (this->computor == RUN_GF2::COMPUTEqr)
			this->data->fbag["qr_out"] = std::move(qr_out);

		// Stopping the process here if computing a given metric
		if (this->computor != RUN_GF2::NORMAL)
			return;

		// nook.tok("PQ took ");
		// nook.tik();

		// #######################################
		//  Here there used to be a section checking and processing disconnected
		//  nodes These were nodes not belonging to the local stack yet with water
		//  For example nodes on a less frequent path.
		//  These are now caught when looping local neighbours and adding
		//  unprocessed upper ones to the queue a tiny number of nodes might still
		//  be sometiems uncaught but the speed gain is worth it
		// #######################################

		// #######################################
		//  Last step: incrementing flow depth
		// #######################################

		// Tacking N nodes have been incremented
		int NN = 0;

		// Looping through local stack
		for (int i = 0; i < tempnodes.size(); ++i) {

			// Local node
			int node = tempnodes[i];

			// double checking if I need to process it
			if (this->is_node_active(node) == false)
				continue;

			// Ref to water height and hydraulic surface
			f_t& thw = this->data->_hw[node];
			f_t& tsurf = this->data->_surface[node];

			// calculating increment
			f_t dhw = 0.;
			// Increment is equal to the divergence of the discharge...
			dhw = this->data->_Qwin[node] -		// input discharge adds water
						(this->data->_Qwout[node] * // output discharge removes water
						 this->param
							 ->capacityFacQw); // playing with the output discahrge can be
																 // used to speed up the iterative process
			dhw *= this->dt;					 // rate to actual increment
			dhw /= this->con->area(i); // Volume to height

			// rate of incrementation
			this->meandhstar += dhw / this->dt;
			++NN;

			// Applying it
			// # To the flow depth
			thw += dhw;
			// and the hydraulic surface
			tsurf += dhw;

			// Correcting if hw < 0 (rare, but can happen in small spots where dt is
			// slightly too high)
			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}
		}

		// Computing the average increment
		if (NN > 0)
			this->meandhstar /= NN;

		// nook.tok("post took ");

		// ##########################
		//  Done
		// ##########################
	}



	void run_hydro_sfd()
	{

		// ###################################
		// Initialising the data structures if needed
		// ###################################

		// In case I am computing the discharge per unit width
		std::vector<f_t> qr_out;
		if (this->computor == RUN_GF2::COMPUTEqr)
			qr_out = std::vector<f_t>(this->con->nxy(), 0.);

		// In case it is the first run
		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
		}

		// The priority queue is the dynamic stack popping/storing nodes in the
		// right way
		std::priority_queue<WaCell<i_t, f_t>,
												std::vector<WaCell<i_t, f_t>>,
												std::less<WaCell<i_t, f_t>>>
			dynastack;

		// ctx is the context neighbourer (helps nabigating through neighbours as
		// the stack is not computed)
		CT_neighbours<i_t, f_t> ctx;

		// Creating local stacks to speed up some operations
		std::vector<i_t> tempnodes, temprec;
		tempnodes.reserve(this->con->_nx);
		temprec.reserve(this->con->_nx);

		// placeholder for receiver data
		std::array<i_t, 8> receivers; // id of rec in the contextual neighbours

		// ###################################
		// Reinitialising the vectors to 0 if not furst run
		// ###################################

		// If active_nodes is activated, I only work on a subset of nodes
		if (this->active_nodes.size() == 0) {

			fillvec(this->data->_Qwin, 0.);
			fillvec(this->xtraQwin, 0.);
			fillvec(this->isInQ, false);
			// else I work on all the nodes
		} else {
			fillvec(this->data->_Qwin, 0., this->active_nodes);
			fillvec(this->xtraQwin, 0., this->active_nodes);
			fillvec(this->isInQ, false, this->active_nodes);
		}

		// ###################################
		// State variables
		// ###################################

		// Mass balance checker
		this->MB_Qwin_out = 0.;
		// increment water checker
		this->meandhstar = 0.;
		f_t sumout = 0.;

		// ###################################
		// Step 1: Preparing the inputs
		// ###################################

		// initialising the dynamic stack and the input points of water
		this->init_dstack<WaCell<i_t, f_t>, decltype(dynastack)>(dynastack);

		// Incrementing the timer
		this->time += this->dt;

		// ###################################
		// Step 2: Traversing the landscape
		// ###################################

		// Timer stuff
		// nook.tok("init took ");
		// nook.tik();

		// the main loop is running as long as there are still nodes in the priority
		// queue
		while (dynastack.empty() == false) {

			// Getting the next node
			auto next = this->_dstack_next<WaCell<i_t, f_t>>(dynastack);

			// deregistering it as being in the PQ stack
			this->isInQ[next.node] = false;

			// ispast is true if the node has not been processed yet
			bool ispast = this->data->_timetracker[next.node] != this->time;

			// Updating the timer
			this->data->_timetracker[next.node] = this->time;

			// Updating the context (fetching local neighbours index and otehr info)
			ctx.update(next.node, *this->con);

			// If the node is outletting the model I skip it
			if (can_out(ctx.boundary)) {
				continue;
			}

			// If the node is no_data I skip it
			if (nodata(ctx.boundary)) {
				continue;
			}

			// ###################################
			// Local state variables
			// ###################################

			// Number of receivers
			int nr = 0;
			// Index of the steepest receiver
			int SSi = 0;
			// Steepest slope
			f_t SS = 0;
			// dx in the steepest direction
			f_t SSdx = 1.;
			// dy in the steepest direction
			f_t SSdy = 1.;

			// Fectching local hw and hydraulic surface
			f_t& thw = this->data->_hw[next.node];
			f_t& tsurf = this->data->_surface[next.node];

			// ###################################
			// Computing the local graph
			// ###################################

			// This function takes care of filtering, selecting and calculating the
			// receivers' characteristics
			this->update_receivers_SFD(ctx,
														 receivers,
														 nr,
														 SS,
														 SSdx,
														 SSdy,
														 SSi);

			// Local stack
			if (ispast) {
				temprec.emplace_back(nr >= 0 ? SSi : next.node);
				tempnodes.emplace_back(next.node);
			}

			// Adding unprocessed upstream neighbours to the Queue to connect
			// disconnected node
			for (int j = 0; j < ctx.nn; ++j) {
				int tn = ctx.neighbours[j];
				if (this->data->_surface[tn] > tsurf && this->data->_hw[tn] > 0 &&
						this->data->_timetracker[tn] != this->time &&
						this->isInQ[tn] == false) {
					dynastack.emplace(
						WaCell<i_t, f_t>(tn, this->data->_surface[tn], 0., 0.));
				}
			}

			// ###################################
			// Computing Water fluxes
			// ###################################

			// Flow velocity u using manning Equation
			// Note I am recasting the slope to a minimum value to avoid numerical
			// instabilities
			f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
								std::sqrt(std::max(1e-6, SS));

			// Converting to output volumetric discharge
			f_t tQwout = thw * u_w * SSdy;

			// Registering local output Volumetric discharge
			this->data->_Qwout[ctx.node] = tQwout;

			// Calculating the input Volumetric Discharge
			// Some water is transmitted via the dynamic cells (for example if they
			// have been in a local minima) rest of the water is in a temp vector
			// storing them
			this->data->_Qwin[ctx.node] = next.Qw + this->xtraQwin[ctx.node];
			this->xtraQwin[ctx.node] = 0.;

			// ###################################
			// Optional Computations
			// ###################################

			if (this->computor == RUN_GF2::COMPUTEqr)
				qr_out[ctx.node] = thw * u_w;

			// ###################################
			// Transferring water to the neighbours
			// ###################################

			// For clarity
			f_t baseQw = this->data->_Qwin[ctx.node];

			// actual index of the receivers
			int rec = SSi;


			// Debug checker catching cases where I have a receivers that shoud not
			// be one. Keeping it cause I am still playing with receivers selection
			if (nodata(this->data->_boundaries[rec]))
				throw std::runtime_error(
					"Fraphflood2::run::run_hydro_sfd::receiver is no data");

			// Checking if my receiver has already been processed or not (Am I in
			// the process of solving a LM)
			bool tizdone = this->data->_timetracker[rec] == this->time;
			// Checking if my receiver is already in the Q
			bool inDaQ = this->isInQ[rec];

			if (inDaQ) {
				// Then I add the extra water to the bucket vector as it will come
				// later
				this->xtraQwin[rec] += baseQw;
			}
			// If this is the sace
			else if (tizdone) {
				// Emplcing the cell and storing the watter within to be transmitted
				// out of the area
				dynastack.emplace(WaCell<i_t, f_t>(
					rec, this->data->_surface[rec], baseQw));

				// node is in the Q
				this->isInQ[rec] = true;

			} else {
				// Otherwise I just need to check if I don't fall out of the node
				// masks (Is it a part of the landscape I want to process)
				if (this->param->TSG_dist == false || this->tsg.xtraMask[rec]) {
					dynastack.emplace(WaCell<i_t, f_t>(
						rec, this->data->_surface[rec], baseQw));
					this->isInQ[rec] = true;
				}
			}
			
		}

		// #######################################
		//  End of the traversal loop
		//  Beginning of the post-process stage
		// #######################################

		// Registering the discharge per unit width if required to compute
		if (this->computor == RUN_GF2::COMPUTEqr)
			this->data->fbag["qr_out"] = std::move(qr_out);

		// Stopping the process here if computing a given metric
		if (this->computor != RUN_GF2::NORMAL)
			return;

		// nook.tok("PQ took ");
		// nook.tik();

		// #######################################
		//  Here there used to be a section checking and processing disconnected
		//  nodes These were nodes not belonging to the local stack yet with water
		//  For example nodes on a less frequent path.
		//  These are now caught when looping local neighbours and adding
		//  unprocessed upper ones to the queue a tiny number of nodes might still
		//  be sometiems uncaught but the speed gain is worth it
		// #######################################

		// #######################################
		//  Last step: incrementing flow depth
		// #######################################

		// Tacking N nodes have been incremented
		int NN = 0;

		// Looping through local stack
		for (int i = 0; i < tempnodes.size(); ++i) {

			// Local node
			int node = tempnodes[i];

			// double checking if I need to process it
			if (this->is_node_active(node) == false)
				continue;

			// Ref to water height and hydraulic surface
			f_t& thw = this->data->_hw[node];
			f_t& tsurf = this->data->_surface[node];

			// calculating increment
			f_t dhw = 0.;
			// Increment is equal to the divergence of the discharge...
			dhw = this->data->_Qwin[node] -		// input discharge adds water
						(this->data->_Qwout[node] * // output discharge removes water
						 this->param
							 ->capacityFacQw); // playing with the output discahrge can be
																 // used to speed up the iterative process
			dhw *= this->dt;					 // rate to actual increment
			dhw /= this->con->area(i); // Volume to height

			// rate of incrementation
			this->meandhstar += dhw / this->dt;
			++NN;

			// Applying it
			// # To the flow depth
			thw += dhw;
			// and the hydraulic surface
			tsurf += dhw;

			// Correcting if hw < 0 (rare, but can happen in small spots where dt is
			// slightly too high)
			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}
		}

		// Computing the average increment
		if (NN > 0)
			this->meandhstar /= NN;

		// nook.tok("post took ");

		// ##########################
		//  Done
		// ##########################
	}

	void run_morpho_sfd()
	{

		// ###################################
		// Initialising the data structures if needed
		// ###################################

		// In case I am computing the discharge per unit width
		std::vector<f_t> qr_out;
		if (this->computor == RUN_GF2::COMPUTEqr)
			qr_out = std::vector<f_t>(this->con->nxy(), 0.);

		// In case it is the first run
		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
		}

		// In case it is the first run
		if (this->data->_Qsin.size() == 0) {
			this->data->_Qsin = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_Qsout = std::vector<f_t>(this->con->nxy(), 0.);
		}

		// The priority queue is the dynamic stack popping/storing nodes in the
		// right way
		std::priority_queue<WaCell<i_t, f_t>,
												std::vector<WaCell<i_t, f_t>>,
												std::less<WaCell<i_t, f_t>>>
			dynastack;

		// ctx is the context neighbourer (helps nabigating through neighbours as
		// the stack is not computed)
		CT_neighbours<i_t, f_t> ctx;

		// Creating local stacks to speed up some operations
		std::vector<i_t> tempnodes, temprec;
		tempnodes.reserve(this->con->_nx);
		temprec.reserve(this->con->_nx);

		// placeholder for receiver data
		std::array<i_t, 8> receivers; // id of rec in the contextual neighbours

		// ###################################
		// Reinitialising the vectors to 0 if not furst run
		// ###################################

		// If active_nodes is activated, I only work on a subset of nodes
		if (this->active_nodes.size() == 0) {

			fillvec(this->data->_Qsin, 0.);
			fillvec(this->data->_Qsout, 0.);
			fillvec(this->isInQ, false);
			fillvec(this->xtraQsout, 0.);
			// else I work on all the nodes
		} else {
			fillvec(this->data->_Qsin, 0., this->active_nodes);
			fillvec(this->data->_Qsout, 0., this->active_nodes);
			fillvec(this->isInQ, false, this->active_nodes);
			fillvec(this->xtraQsout, 0., this->active_nodes);
		}

		// ###################################
		// State variables
		// ###################################


		// ###################################
		// Step 1: Preparing the inputs
		// ###################################

		// initialising the dynamic stack and the input points of water
		this->init_dstack<WaCell<i_t, f_t>, decltype(dynastack)>(dynastack);

		// Incrementing the timer
		this->time += this->dt;

		// ###################################
		// Step 2: Traversing the landscape
		// ###################################

		// Timer stuff
		// nook.tok("init took ");
		// nook.tik();

		// the main loop is running as long as there are still nodes in the priority
		// queue
		while (dynastack.empty() == false) {

			// Getting the next node
			auto next = this->_dstack_next<WaCell<i_t, f_t>>(dynastack);

			// deregistering it as being in the PQ stack
			this->isInQ[next.node] = false;

			// ispast is true if the node has not been processed yet
			bool ispast = this->data->_timetracker[next.node] != this->time;

			// Updating the timer
			this->data->_timetracker[next.node] = this->time;

			// Updating the context (fetching local neighbours index and otehr info)
			ctx.update(next.node, *this->con);

			// If the node is outletting the model I skip it
			if (can_out(ctx.boundary)) {
				continue;
			}

			// If the node is no_data I skip it
			if (nodata(ctx.boundary)) {
				continue;
			}

			// ###################################
			// Local state variables
			// ###################################

			// Number of receivers
			int nr = 0;
			// Index of the steepest receiver
			int SSi = 0;
			// Steepest slope
			f_t SS = 0;
			// dx in the steepest direction
			f_t SSdx = 1.;
			// dy in the steepest direction
			f_t SSdy = 1.;

			// Fectching local hw and hydraulic surface
			f_t& thw = this->data->_hw[next.node];
			f_t& tsurf = this->data->_surface[next.node];

			// ###################################
			// Computing the local graph
			// ###################################

			// This function takes care of filtering, selecting and calculating the
			// receivers' characteristics
			this->update_receivers_SFD(ctx,
														 receivers,
														 nr,
														 SS,
														 SSdx,
														 SSdy,
														 SSi);

			// Local stack
			if (ispast) {
				temprec.emplace_back(nr >= 0 ? SSi : next.node);
				tempnodes.emplace_back(next.node);
			}

			// Adding unprocessed upstream neighbours to the Queue to connect
			// disconnected node
			for (int j = 0; j < ctx.nn; ++j) {
				int tn = ctx.neighbours[j];
				if (this->data->_surface[tn] > tsurf && this->data->_hw[tn] > 0 &&
						this->data->_timetracker[tn] != this->time &&
						this->isInQ[tn] == false) {
					dynastack.emplace(
						WaCell<i_t, f_t>(tn, this->data->_surface[tn], 0., 0.));
				}
			}

			// ###################################
			// Computing Sed fluxes
			// ###################################

			if(this->param->morphomode == MORPHOMODE::MPM){

				f_t tau = thw * SS * this->param->GRAVITY * this->param->rho_water;
				f_t capacity = 0.;
				if(tau > this->param->tau_c){
					capacity = this->param->E * std::pow(tau - this->param->tau_c,this->param->alpha);
				}
				this->data->_Qsout[ctx.node] = capacity;
				this->data->_Qsin[SSi] += capacity;


			}

			// ###################################
			// Transferring water to the neighbours
			// ###################################

			// actual index of the receivers
			int rec = SSi;


			// Debug checker catching cases where I have a receivers that shoud not
			// be one. Keeping it cause I am still playing with receivers selection
			if (nodata(this->data->_boundaries[rec]))
				throw std::runtime_error(
					"Fraphflood2::run::run_morpho_sfd::receiver is no data");

			// Checking if my receiver has already been processed or not (Am I in
			// the process of solving a LM)
			bool tizdone = this->data->_timetracker[rec] == this->time;
			// Checking if my receiver is already in the Q
			bool inDaQ = this->isInQ[rec];

			if (inDaQ) {
				continue;
			}
			// If this is the sace
			else if (tizdone) {
				// Emplcing the cell and storing the watter within to be transmitted
				// out of the area
				dynastack.emplace(WaCell<i_t, f_t>(
					rec, this->data->_surface[rec], 0.));

				// node is in the Q
				this->isInQ[rec] = true;

			} else {
				// Otherwise I just need to check if I don't fall out of the node
				// masks (Is it a part of the landscape I want to process)
				if (this->param->TSG_dist == false || this->tsg.xtraMask[rec]) {
					dynastack.emplace(WaCell<i_t, f_t>(
						rec, this->data->_surface[rec], 0.));
					this->isInQ[rec] = true;
				}
			}
			
		}

		// #######################################
		//  End of the traversal loop
		//  Beginning of the post-process stage
		// #######################################


		// nook.tok("PQ took ");
		// nook.tik();


		// #######################################
		//  Last step: incrementing flow depth
		// #######################################

		// Tacking N nodes have been incremented
		int NN = 0;

		// Looping through local stack
		for (int i = 0; i < tempnodes.size(); ++i) {

			// Local node
			int node = tempnodes[i];

			// double checking if I need to process it
			if (this->is_node_active(node) == false)
				continue;

			// Ref to water height and hydraulic surface
			f_t& thw = this->data->_hw[node];
			f_t& tsurf = this->data->_surface[node];

			// calculating increment
			f_t dhs = 0.;
			// Increment is equal to the divergence of the discharge...
			dhs = this->data->_Qsin[node] -		// input discharge adds water
						((this->data->_Qsout[node] + this->xtraQsout[node]) * // output discharge removes water
						 this->param
							 ->capacityFacQs); // playing with the output discahrge can be
																 // used to speed up the iterative process
			dhs *= this->dt;					 // rate to actual increment
			dhs /= this->con->area(i); // Volume to height


			// Applying it
			// # To the opposite flow depth as adding sediments decreases flow depth yo
			thw -= dhs;


			// Correcting if hw < 0 (rare, but can happen in small spots where dt is
			// slightly too high)
			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}
		}


		// nook.tok("post took ");

		// ##########################
		//  Done
		// ##########################
	}

	void morphoton(int N, f_t edt)
	{

		CT_neighbours<i_t, f_t> ctx;
		std::array<i_t, 8> receivers;
		std::array<f_t, 8> receiversWeights;

		// std::cout << "START" << std::endl;

		for (int __ = 0; __ < N; ++__) {

			this->time += edt;

			int node = -1;
			f_t tqs = 0;
			while (true) {
				i_t tid =
					std::floor(this->data->randu->get() * this->input_node_Qw.size());

				if (tid >= this->input_node_Qw.size())
					throw std::runtime_error("morphoton but no input node Qw");

				node = this->input_node_Qw[tid];
				if (this->is_node_active(node)) {

					if (tid >= this->input_Qs.size())
						throw std::runtime_error("morphoton but no input qs");

					tqs = this->input_Qs[tid] / this->con->_dy;
					break;
				}
			}

			while (true) {

				if (this->is_node_active(node) == false)
					break;

				ctx.update(node, *this->con);
				int nr = 0;
				int SSi = 0;
				f_t SS = 0;
				f_t SSdx = 1.;
				f_t SSdy = 1.;
				this->update_receivers(
					ctx, receivers, receiversWeights, nr, SS, SSdx, SSdy, SSi, true);

				f_t tdt = this->time - this->data->_timetracker[node];
				this->data->_timetracker[node] = this->time;

				f_t thw = this->data->_hw[node];

				// Calculating shear stress
				f_t tau = this->param->rho_water * this->param->GRAVITY * SS * thw;

				// rates
				// # basal erosion
				f_t edot = 0.;
				// # lateral erosion
				f_t eldot = 0;
				// # basal dep
				f_t ddot = 0.;
				// # lateral dep
				f_t dldot = 0.;

				// erosion only happens if critical shear stress increases
				if (tau > this->param->tau_c) {

					// calculating basal erosion: ke x (shear - critshear) ^ alpha
					edot = this->param->ke(node) *
								 std::pow(tau - this->param->tau_c, this->param->alpha);

					// lateral erosion
					if (this->param->bank_erosion) {
						// Getting the highest neighbour to erode
						int hni = ctx.idxHighestNeighbour(
							this->data->_surface, this->data->_hw, SSi);

						// hni is -1 if there is no higher neighbours
						if (hni >= 0) {

							// neihbour ID
							int hn = ctx.neighbours[hni];

							// difference in substrate elevation defining the max lateral
							// erosion rate
							f_t delta_Z =
								(this->data->_surface[hn] - this->data->_hw[hn] -
								 this->data->_surface[ctx.node] + this->data->_hw[ctx.node]);

							// Lateral gradient
							f_t latslope = delta_Z / ctx.neighboursDx[hni];

							// Actual rate
							eldot = this->param->kel * latslope * edot;

							// checking we do not erode too much
							if (eldot * this->dt > delta_Z)
								eldot = delta_Z / tdt;

							this->data->_hw[hn] += eldot * tdt;

							if (this->data->_hw[hn] < 0)
								this->data->_hw[hn] = 0.;
						}
					}
				}

				// TEMP DEACTIVATION
				// ddot = tqs / (this->param->kd * SSdy); //;
				f_t latslope = 0.;
				if (this->param->kdl > 0) {
					int hni =
						ctx.idxLowestNeighbour(this->data->_surface, this->data->_hw, SSi);

					// hni is -1 if there is no higher neighbours
					if (hni >= 0) {
						// neihbour ID
						int hn = ctx.neighbours[hni];

						// difference in substrate elevation defining the max lateral
						// erosion rate
						f_t delta_Z =
							(this->data->_surface[ctx.node] - this->data->_hw[ctx.node] -
							 this->data->_surface[hn] + this->data->_hw[hn]);

						latslope = delta_Z / ctx.neighboursDx[hni];

						dldot = latslope * this->param->kdl * tqs;

						this->data->_hw[hn] -= dldot * tdt;
						if (this->data->_hw[node] < 0)
							this->data->_hw[node] = 0.;
					}
				}

				// this->data->_debug[ctx.node] = ddot;
				// ddot = std::pow(this->data->_Qsin[ctx.node] / (this->param->kd *
				// SSdy), 0.5);

				// this->data->_Qsout[ctx.node] =
				// 	std::max(static_cast<f_t>(0.),
				// 					 this->data->_Qsin[ctx.node] +
				// 						 (edot + eldot - ddot - dldot) *
				// this->con->area(ctx.node));

				double K = (1. / this->param->kd + this->param->kdl * latslope);

				double edotpsy = (edot + eldot) / K;
				double C1 = tqs - edotpsy;

				f_t qs_out = (edotpsy + C1 * std::exp(-SSdx * K));

				if (qs_out < 0)
					throw std::runtime_error("qs < 0");

				this->data->_hw[node] -= (tqs - qs_out) / SSdx * tdt;
				if (this->data->_hw[node] < 0)
					this->data->_hw[node] = 0.;

				node = SSi;
				tqs = qs_out;
			}
		}
	}

	bool is_node_active(int i)
	{
		if (nodata(this->data->_boundaries[i]) ||
				can_out(this->data->_boundaries[i])) {
			return false;
		}

		if (this->param->TSG_dist && this->tsg.xtraMask[i] == false) {
			return false;
		}

		return true;
	}

	void solve_analytically_if(f_t threshold)
	{

		CT_neighbours<i_t, f_t> ctx;

		std::array<i_t, 8> receivers;
		std::array<f_t, 8> receiversWeights;

		for (int i = 0; i < this->con->nxy(); ++i) {

			if (this->is_node_active(i) == false)
				continue;

			if (this->data->_Qwin[i] == 0 ||
					abs(1 - this->data->_Qwout[i] / this->data->_Qwin[i]) < threshold)
				continue;

			// std::cout << next.node << std::endl;
			ctx.update(i, *this->con);
			;

			int nr = 0;
			int SSi = 0;
			f_t SS = 0;
			f_t SSdx = 1.;
			f_t SSdy = 1.;

			// std::cout << "A1" << std::endl;
			this->update_receivers(
				ctx, receivers, receiversWeights, nr, SS, SSdx, SSdy, SSi);

			f_t hw = this->data->_hw[i];

			f_t nhw = std::pow((this->mannings * this->data->_Qwin[i]) /
													 (SSdy * std::sqrt(SS)),
												 3. / 5.);
			this->data->_hw[i] = nhw;
			this->data->_surface[i] += nhw - hw;
		}
	}

	void run_dynamic()
	{

		// if (this->data->_timetracker.size() == 0) {
		// 	throw std::runtime_error("Run 1 by 1 cannot work is ");
		// }

		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_debug = std::vector<f_t>(this->con->nxy(), 0.);
		}

		if (this->param->gf2_morpho && this->data->_Qsin.size() == 0) {
			this->data->_Qsin = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_Qsout = std::vector<f_t>(this->con->nxy(), 0.);
		}

		std::vector<f_t> Qwinadd(this->data->_Qwin.size(), 0.);
		std::vector<f_t> Qsinadd(this->data->_Qwin.size(), 0.);

		this->add_init_flux();
		this->add_init_flux(Qwinadd, Qsinadd);

		CT_neighbours<i_t, f_t> ctx;

		this->time += this->dt;

		std::array<i_t, 8> receivers;
		std::array<f_t, 8> receiversWeights;

		f_t sumout = 0.;

		for (int i = 0; i < this->con->nxy(); ++i) {

			// Updating the timer
			this->data->_timetracker[i] = this->time;

			// std::cout << i << std::endl;
			ctx.update(i, *this->con);

			if (can_out(ctx.boundary))
				continue;

			if (nodata(ctx.boundary))
				continue;

			int nr = 0;
			int SSi = 0;
			f_t SS = 0;
			f_t SSdx = 1.;
			f_t SSdy = 1.;

			this->update_receivers(
				ctx, receivers, receiversWeights, nr, SS, SSdx, SSdy, SSi);

			f_t& thw = this->data->_hw[i];
			f_t& tsurf = this->data->_surface[i];
			// Actual flux calculation

			f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
								std::sqrt(std::max(1e-6, SS));
			f_t tQwout = thw * u_w * SSdy;

			this->data->_Qwout[ctx.node] = tQwout;

			if (this->param->gf2_morpho) {
				f_t tau = this->param->rho_water * this->param->GRAVITY * SS * thw;
				f_t edot = 0.;
				f_t ddot = 0.;
				if (tau > this->param->tau_c) {
					edot = this->param->ke(ctx.node) *
								 std::pow(tau - this->param->tau_c, this->param->alpha);
					// std::pow(tau - this->param->tau_c, this->param->alpha) *
					// (this->data->randu->get()*0.5 + 1.);
				}
				// ddot = this->data->_Qsin[ctx.node] / (this->param->kd * SSdy *
				// (this->data->randu->get()*0.5 + 1.));
				ddot = this->data->_Qsin[ctx.node] / (this->param->kd * SSdy); //;
				// ddot = std::pow(this->data->_Qsin[ctx.node] / (this->param->kd *
				// SSdy), 0.5);

				this->data->_Qsout[ctx.node] =
					std::max(static_cast<f_t>(0.),
									 this->data->_Qsin[ctx.node] + edot * SSdx * SSdy -
										 ddot * SSdx * SSdy);
			}

			// if (!ispast) {
			// 	thw += this->hw_increment_LM;
			// 	tsurf += this->hw_increment_LM;
			// }

			// NEED TO CARY ON ADDING MORPHO HERE

			f_t baseQw = this->data->_Qwin[i];
			f_t baseQs;
			if (this->param->gf2_morpho)
				baseQs = this->data->_Qsout[i];

			// if(nr==0){
			// 	Qwinadd[i] += this->data->_Qwin[i];
			// 	if (this->param->gf2_morpho) Qsinadd[i] += this->data->_Qsin[i];
			// }

			for (int j = 0; j < nr; ++j) {
				int tn = receivers[j];
				int rec = ctx.neighbours[tn];
				Qwinadd[rec] += (this->param->transient_flow)
													? tQwout * receiversWeights[j]
													: baseQw * receiversWeights[j];
				if (this->param->gf2_morpho)
					Qsinadd[rec] += baseQs * receiversWeights[j];
			}
		}

		this->data->_Qwin = std::move(Qwinadd);
		if (this->param->gf2_morpho)
			this->data->_Qsin = std::move(Qsinadd);

		for (int i = 0; i < this->con->nxy(); ++i) {

			if (nodata(this->data->_boundaries[i]) ||
					can_out(this->data->_boundaries[i]))
				continue;

			f_t& thw = this->data->_hw[i];
			f_t& tsurf = this->data->_surface[i];

			f_t dhw = 0.;
			if (this->sbg_method == SUBGRAPHMETHOD::V1) {
				dhw = this->data->_Qwin[i] - this->data->_Qwout[i];
			} else if (this->sbg_method == SUBGRAPHMETHOD::FILLONLY) {
				dhw = this->data->_Qwin[i];
			}

			dhw *= this->dt;
			dhw /= this->con->area(i);

			thw += dhw;
			tsurf += dhw;

			if (this->param->gf2_morpho) {
				f_t dz = this->data->_Qsin[i] - this->data->_Qsout[i];
				dz *= this->dt * this->param->time_dilatation_morpho;
				dz /= this->con->area(i);
				tsurf += dz;
			}

			if (std::isfinite(tsurf) == false) {
				std::cout << this->data->_Qwout[i] << "|" << this->data->_Qwin[i]
									<< std::endl;
				;
				throw std::runtime_error("blug");
			}

			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}
		}

		// std::cout << "Sumout: " << sumout << std::endl; ;
	}

	f_t evaluate_convergence(f_t tol = 0.05, bool onMask = true)
	{

		int Neq = 0;
		int N = 0;

		for (int i = 0; i < this->con->nxy(); ++i) {

			if (this->data->_Qwin[i] == 0)
				continue;

			if (onMask) {
				if (this->tsg.xtraMask[i]) {
					if (abs(1 - this->data->_Qwout[i] / this->data->_Qwin[i]) < tol)
						++Neq;
					++N;
				}
			} else {
				if (abs(1 - this->data->_Qwout[i] / this->data->_Qwin[i]) < tol)
					++Neq;
				++N;
			}
		}

		return (N > 0) ? static_cast<f_t>(Neq) / N : 0.;
	}

	/// Periodic refeeding of sediments
	/// The sediments outputting the system are refed to every entry points
	void feed_inputQs_with_out()
	{

		int n_inputs = static_cast<int>(this->input_Qs.size());

		if (n_inputs == 0 || this->param->gf2_morpho == false)
			return;

		f_t tot_QS = 0.;
		for (int i = 0; i < this->con->nxy(); ++i) {

			if (can_out(this->data->_boundaries[i]))
				tot_QS += this->data->_Qsin[i];
		}

		f_t sumsum = 0;
		for (int j = 0; j < n_inputs; ++j) {
			this->input_Qs[j] = tot_QS / n_inputs;
			sumsum += this->input_Qs[j];
		}

		// std::cout << "tot_QS:" << tot_QS << " vs " << sumsum << std::endl;
	}

	// experimental

	void anarun_subgraphflood()
	{

		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_debug = std::vector<f_t>(this->con->nxy(), 0.);
		}

		if (this->param->gf2_morpho && this->data->_Qsin.size() == 0) {
			this->data->_Qsin = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_Qsout = std::vector<f_t>(this->con->nxy(), 0.);
		}

		std::priority_queue<WaCell<i_t, f_t>,
												std::vector<WaCell<i_t, f_t>>,
												std::less<WaCell<i_t, f_t>>>
			dynastack;

		fillvec(this->data->_Qwin, 0.);
		fillvec(this->isInQ, false);
		if (this->param->gf2_morpho) {
			fillvec(this->data->_Qsin, 0.);
			fillvec(this->data->_Qsout, 0.);
		}

		this->init_dstack<WaCell<i_t, f_t>, decltype(dynastack)>(dynastack);

		// CT_neighbourer_WaCell<i_t, f_t> ctx;
		CT_neighbours<i_t, f_t> ctx;

		// int nndt = 0;
		// for(auto v:this->data->_boundaries){
		// 	if(nodata(v)) ++nndt;
		// }
		// std::cout << "I have " << nndt << "no data " << std::endl;

		// fillvec(this->data->_vmot_hw,0.);
		this->time += this->dt;
		// std::vector<bool> isdone(this->con->nxy(), false);
		// int ndone = 0;
		// int nredone = 0;

		std::array<i_t, 8> receivers;
		std::array<f_t, 8> receiversWeights;

		f_t sumout = 0.;

		// std::cout << "Starting the process" << std::endl;

		while (dynastack.empty() == false) {

			// Getting the next node
			auto next = this->_dstack_next<WaCell<i_t, f_t>>(dynastack);
			// auto next = dynastack.top();
			// dynastack.pop();
			// std::cout << next.node << "|" << this->data->_surface[next.node] <<
			// "||";

			this->isInQ[next.node] = false;

			bool ispast = this->data->_timetracker[next.node] != this->time;

			// if(isdone[next.node] == false){
			// 	ndone++;
			// 	isdone[next.node] = true;
			// 	if(ndone % 100 == 0)
			// 		std::cout << ndone << " vs " << nredone << " PQsizzla: " <<
			// dynastack.size() << " this->debugyolo " << this->debugyolo <<
			// std::endl; }else{ 	nredone++;
			// }

			// Updating the timer
			this->data->_timetracker[next.node] = this->time;

			// std::cout << next.node << std::endl;
			ctx.update(next.node, *this->con);
			;

			// std::cout << BC2str(ctx.boundary) << std::endl;

			if (ispast) {
				// this->data->_Qwin[next.node] = std::max(this->data->_Qwin[next.node],
				// next.Qw);
				this->data->_Qwin[next.node] += next.Qw;
			}

			if (can_out(ctx.boundary)) {
				if (ispast)
					sumout += this->data->_Qwin[ctx.node];
				continue;
			}

			if (nodata(ctx.boundary)) {
				// std::cout << "nodata reached" << std::endl;
				continue;
			}

			int nr = 0;
			int SSi = 0;
			f_t SS = 0;
			f_t SSdx = 1.;
			f_t SSdy = 1.;

			// std::cout << "A1" << std::endl;
			this->update_receivers(
				ctx, receivers, receiversWeights, nr, SS, SSdx, SSdy, SSi);
			// std::cout << "A2" << std::endl;

			f_t& thw = this->data->_hw[next.node];
			f_t& tsurf = this->data->_surface[next.node];
			// Actual flux calculation

			f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
								std::sqrt(std::max(1e-6, SS));
			f_t tQwout = thw * u_w * SSdy;

			this->data->_Qwout[ctx.node] = tQwout;

			if (this->param->gf2_morpho) {
				f_t tau = this->param->rho_water * this->param->GRAVITY * SS * thw;
				f_t edot = 0.;
				f_t ddot = 0.;
				if (tau > this->param->tau_c) {
					edot = this->param->ke(ctx.node) *
								 std::pow(tau - this->param->tau_c, this->param->alpha);
				}
				ddot = this->data->_Qsin[ctx.node] / (this->param->kd * SSdy);

				this->data->_Qsout[ctx.node] =
					std::max(static_cast<f_t>(0.),
									 this->data->_Qsin[ctx.node] + edot * SSdx * SSdy -
										 ddot * SSdx * SSdy);
			}

			if (!ispast) {
				thw += this->hw_increment_LM;
				tsurf += this->hw_increment_LM;
			}

			// NEED TO CARY ON ADDING MORPHO HERE

			f_t baseQw = (ispast) ? this->data->_Qwin[next.node] : next.Qw;

			f_t nhw = std::pow(this->mannings * baseQw /
													 (std::sqrt(std::max(1e-6, SS)) * SSdy),
												 3. / 5.);
			this->data->_Qwout[ctx.node] = nhw;

			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				int rec = ctx.neighbours[i];
				bool tizdone = this->data->_timetracker[rec] == this->time;
				if (tizdone) {
					if (this->isInQ[rec])
						throw std::runtime_error("happens45");
					dynastack.emplace(WaCell<i_t, f_t>(
						rec, this->data->_surface[rec], baseQw * receiversWeights[j]));
				} else {
					if (this->isInQ[rec] == false) {
						dynastack.emplace(WaCell<i_t, f_t>(rec, this->data->_surface[rec]));
						this->isInQ[rec] = true;
					}
					this->data->_Qwin[rec] += baseQw * receiversWeights[j];
				}
			}
		}

		for (int i = 0; i < this->con->nxy(); ++i) {
			f_t nhw = this->data->_Qwout[i] + this->data->_hw[i];
			nhw /= 2;
			f_t dhw = nhw - this->data->_hw[i];
			this->data->_hw[i] += dhw;
			this->data->_surface[i] += dhw;
			this->data->_Qwout[i] = 0;
		}

		// std::cout << "Sumout: " << sumout << std::endl; ;
	}

	template<class CTX>
	void update_receivers(CTX& ctx,
												std::array<i_t, 8>& receivers,
												std::array<f_t, 8>& receiversWeights,
												int& nr,
												f_t& SS,
												f_t& SSdx,
												f_t& SSdy,
												i_t& SSn,
												bool stochaslope = false)
	{
		nr = 0;
		bool recout;
		bool recdone;

		bool allout;
		bool alldone;

		f_t sumSdw = 0.;
		SS = 0.;
		f_t SSstoch = 0.;

		// pass 1: check the receivers, if no receivers fill up
		while (nr == 0) {
			recout = false;
			recdone = false;
			allout = true;
			alldone = true;

			for (int i = 0; i < ctx.nn; ++i) {

				if (nodata(ctx.neighboursCode[i]))
					continue;

				if (this->data->_surface[ctx.node] >
						this->data->_surface[ctx.neighbours[i]]) {
					if (can_receive(ctx.neighboursCode[i])) {
						receivers[nr] = i;
						if (can_out(ctx.neighboursCode[i])) {
							recout = true;
						} else {
							allout = false;
						}
						if (this->time == this->data->_timetracker[ctx.neighbours[i]]) {
							recdone = true;
						} else {
							alldone = false;
						}
						++nr;
					}
				}
			}
			if (nr == 0) {
				f_t tinc = this->hw_increment_LM +
									 (this->data->randu->get() * this->hw_increment_LM) / 2;
				this->data->_surface[ctx.node] += tinc;
				this->data->_hw[ctx.node] += tinc;
			}
			// std::cout << this->data->_surface[ctx.node] << "|";
		}

		// pass2: recast to only the that can out
		if (recout) {
			int nnr = 0;
			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				if (can_out(ctx.neighboursCode[i])) {
					receivers[nnr] = receivers[j];
					++nnr;
				}
			}
			nr = nnr;
		}
		// Pass 3: if no out, some done but not all done
		else if (recdone && !alldone) {
			int nnr = 0;
			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				if (this->time != this->data->_timetracker[ctx.neighbours[i]]) {
					receivers[nnr] = receivers[j];
					++nnr;
				}
			}
			nr = nnr;
		}
		// pass 4: all recs are done, I then choose a single one randomly
		else if (alldone) {
			int ti = static_cast<int>(this->data->randu->get() * nr);
			nr = 1;
			receivers[0] = receivers[ti];
		}

		// pass 5 precalculate SS and weights
		for (int j = 0; j < nr; ++j) {
			int i = receivers[j];
			f_t tS = this->get_Sw(ctx.node, ctx.neighbours[i], ctx.neighboursDx[i]);
			f_t tSdwinc = tS * ctx.neighboursDy[i];
			if (stochaslope)
				tSdwinc *= this->data->randu->get();
			// if(this->param->gf2_morpho)
			// 	tSdwinc *= this->data->randu->get();
			receiversWeights[j] = tSdwinc;

			sumSdw += tSdwinc;
			if (tS > SS && stochaslope == false) {
				SS = tS;
				SSdx = ctx.neighboursDx[i];
				SSdy = ctx.neighboursDy[i];
				SSn = ctx.neighbours[i];
			} else if (stochaslope) {
				if (tS > SS) {
					SS = tS;
				}

				if (SSstoch < tSdwinc) {
					SSstoch = tSdwinc;
					SSdx = ctx.neighboursDx[i];
					SSdy = ctx.neighboursDy[i];
					SSn = ctx.neighbours[i];
				}
			}
		}

		// pass 6: calculate weights
		if (sumSdw > 0) {
			for (int j = 0; j < nr; ++j) {
				receiversWeights[j] /= sumSdw;
			}
		} else {
			for (int j = 0; j < nr; ++j) {
				receiversWeights[j] /= nr;
			}
		}

		if (SSn >= 0)
			if (can_out(this->data->_boundaries[SSn]) &&
					this->param->gf2Bmode == BOUNDARY_HW::FIXED_SLOPE) {
				SS = std::min(this->param->gf2Bbval, SS);
			}

		// done
	}

	template<class CTX>
	void update_receivers(CTX& ctx,
												std::array<i_t, 8>& receivers,
												std::array<f_t, 8>& receiversWeights,
												int& nr,
												f_t& SS,
												f_t& SSdx,
												f_t& SSdy,
												i_t& SSn,
												std::array<f_t, 8>& receiversSlopes,
												bool stochaslope = false)
	{
		nr = 0;
		bool recout;
		bool recdone;

		bool allout;
		bool alldone;

		f_t sumSdw = 0.;
		SS = 0.;
		f_t SSstoch = 0.;

		// pass 1: check the receivers, if no receivers fill up
		while (nr == 0) {
			recout = false;
			recdone = false;
			allout = true;
			alldone = true;

			for (int i = 0; i < ctx.nn; ++i) {

				if (nodata(ctx.neighboursCode[i]))
					continue;

				if (this->data->_surface[ctx.node] >
						this->data->_surface[ctx.neighbours[i]]) {
					if (can_receive(ctx.neighboursCode[i])) {
						receivers[nr] = i;
						if (can_out(ctx.neighboursCode[i])) {
							recout = true;
						} else {
							allout = false;
						}
						if (this->time == this->data->_timetracker[ctx.neighbours[i]]) {
							recdone = true;
						} else {
							alldone = false;
						}
						++nr;
					}
				}
			}
			if (nr == 0) {
				f_t tinc = this->hw_increment_LM +
									 (this->data->randu->get() * this->hw_increment_LM) / 2;
				this->data->_surface[ctx.node] += tinc;
				this->data->_hw[ctx.node] += tinc;
			}
			// std::cout << this->data->_surface[ctx.node] << "|";
		}

		// pass2: recast to only the that can out
		if (recout) {
			int nnr = 0;
			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				if (can_out(ctx.neighboursCode[i])) {
					receivers[nnr] = receivers[j];
					++nnr;
				}
			}
			nr = nnr;
		}
		// Pass 3: if no out, some done but not all done
		else if (recdone && !alldone) {
			int nnr = 0;
			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				if (this->time != this->data->_timetracker[ctx.neighbours[i]]) {
					receivers[nnr] = receivers[j];
					++nnr;
				}
			}
			nr = nnr;
		}
		// pass 4: all recs are done, I then choose a single one randomly
		else if (alldone) {
			int ti = static_cast<int>(this->data->randu->get() * nr);
			nr = 1;
			receivers[0] = receivers[ti];
		}

		// pass 5 precalculate SS and weights
		for (int j = 0; j < nr; ++j) {
			int i = receivers[j];
			f_t tS = this->get_Sw(ctx.node, ctx.neighbours[i], ctx.neighboursDx[i]);
			receiversSlopes[j] = tS;
			f_t tSdwinc = tS * ctx.neighboursDy[i];
			if (stochaslope)
				tSdwinc *= this->data->randu->get();
			// if(this->param->gf2_morpho)
			// 	tSdwinc *= this->data->randu->get();
			receiversWeights[j] = tSdwinc;

			sumSdw += tSdwinc;
			if (tS > SS && stochaslope == false) {
				SS = tS;
				SSdx = ctx.neighboursDx[i];
				SSdy = ctx.neighboursDy[i];
				SSn = ctx.neighbours[i];
			} else if (stochaslope) {
				if (tS > SS) {
					SS = tS;
				}

				if (SSstoch < tSdwinc) {
					SSstoch = tSdwinc;
					SSdx = ctx.neighboursDx[i];
					SSdy = ctx.neighboursDy[i];
					SSn = ctx.neighbours[i];
				}
			}
		}

		// pass 6: calculate weights

		bool needCorrection = false;
		if (sumSdw > 0) {
			for (int j = 0; j < nr; ++j) {
				receiversWeights[j] /= sumSdw;
				if (receiversWeights[j] < this->param->minWeightQw) {
					needCorrection = true;
					receiversWeights[j] = 0.;
				}
			}
			if (needCorrection) {
				sumSdw = 0;
				for (int j = 0; j < nr; ++j)
					sumSdw += receiversWeights[j];
				for (int j = 0; j < nr; ++j)
					receiversWeights[j] /= sumSdw;
			}

		} else {
			for (int j = 0; j < nr; ++j) {
				receiversWeights[j] /= nr;
			}
		}

		if (SSn >= 0)
			if (can_out(this->data->_boundaries[SSn]) &&
					this->param->gf2Bmode == BOUNDARY_HW::FIXED_SLOPE) {
				SS = std::min(this->param->gf2Bbval, SS);
			}

		// done
	}


	template<class CTX>
	void update_receivers_SFD(CTX& ctx,
												std::array<i_t, 8>& receivers,
												int& nr,
												f_t& SS,
												f_t& SSdx,
												f_t& SSdy,
												i_t& SSn,
												bool stochaslope = false)
	{
		nr = 0;
		bool recout;
		bool recdone;

		bool allout;
		bool alldone;

		f_t sumSdw = 0.;
		SS = 0.;
		f_t SSstoch = 0.;

		// pass 1: check the receivers, if no receivers fill up
		while (nr == 0) {
			recout = false;
			recdone = false;
			allout = true;
			alldone = true;

			for (int i = 0; i < ctx.nn; ++i) {

				if (nodata(ctx.neighboursCode[i]))
					continue;

				if (this->data->_surface[ctx.node] >
						this->data->_surface[ctx.neighbours[i]]) {
					if (can_receive(ctx.neighboursCode[i])) {
						receivers[nr] = i;
						if (can_out(ctx.neighboursCode[i])) {
							recout = true;
						} else {
							allout = false;
						}
						if (this->time == this->data->_timetracker[ctx.neighbours[i]]) {
							recdone = true;
						} else {
							alldone = false;
						}
						++nr;
					}
				}
			}
			if (nr == 0) {
				f_t tinc = this->hw_increment_LM +
									 (this->data->randu->get() * this->hw_increment_LM) / 2;
				this->data->_surface[ctx.node] += tinc;
				this->data->_hw[ctx.node] += tinc;
			}
			// std::cout << this->data->_surface[ctx.node] << "|";
		}

		// pass2: recast to only the that can out
		if (recout) {
			int nnr = 0;
			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				if (can_out(ctx.neighboursCode[i])) {
					receivers[nnr] = receivers[j];
					++nnr;
				}
			}
			nr = nnr;
		}
		// Pass 3: if no out, some done but not all done
		else if (recdone && !alldone) {
			int nnr = 0;
			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				if (this->time != this->data->_timetracker[ctx.neighbours[i]]) {
					receivers[nnr] = receivers[j];
					++nnr;
				}
			}
			nr = nnr;
		}
		// pass 4: all recs are done, I then choose a single one randomly
		else if (alldone) {
			int ti = static_cast<int>(this->data->randu->get() * nr);
			nr = 1;
			receivers[0] = receivers[ti];
		}

		// pass 5 precalculate SS and weights
		for (int j = 0; j < nr; ++j) {
			int i = receivers[j];
			f_t tS = this->get_Sw(ctx.node, ctx.neighbours[i], ctx.neighboursDx[i]);
			f_t tSdwinc = tS * ctx.neighboursDy[i];
			if (stochaslope)
				tSdwinc *= this->data->randu->get();

			sumSdw += tSdwinc;
			if (tS > SS && stochaslope == false) {
				SS = tS;
				SSdx = ctx.neighboursDx[i];
				SSdy = ctx.neighboursDy[i];
				SSn = ctx.neighbours[i];
			} else if (stochaslope) {
				if (tS > SS) {
					SS = tS;
				}

				if (SSstoch < tSdwinc) {
					SSstoch = tSdwinc;
					SSdx = ctx.neighboursDx[i];
					SSdy = ctx.neighboursDy[i];
					SSn = ctx.neighbours[i];
				}
			}
		}

		if (SSn >= 0)
			if (can_out(this->data->_boundaries[SSn]) &&
					this->param->gf2Bmode == BOUNDARY_HW::FIXED_SLOPE) {
				SS = std::min(this->param->gf2Bbval, SS);
			}

		// done
	}

	f_t get_Sw(int d, int r, f_t dx)
	{
		return (this->data->_surface[d] - this->data->_surface[r]) / dx;
	}

	template<class CELL>
	CELL _dstack_next(
		std::priority_queue<CELL, std::vector<CELL>, std::less<CELL>>& dynastack)
	{
		auto next = dynastack.top();
		// std::cout << next.Qw;

		dynastack.pop();

		if (dynastack.empty() == false) {
			while (dynastack.top().node == next.node) {
				// std::cout << next.Qw;
				next.ingest(dynastack.top());
				// std::cout << " | " << next.Qw << std::endl;
				dynastack.pop();
				if (dynastack.empty())
					break;
			}
		}
		// std::cout << " | " << next.Qw << std::endl;;

		return next;
	}

	template<class CELL, class CTX, class Q_t>
	void _subGF_process_node(CELL& next, CTX& ctx, Q_t& dynastack)
	{

		f_t& thw = this->data->_hw[next.node];
		f_t& tsurf = this->data->_surface[next.node];

		// Regulating weights to avoid useless spreading
		this->regulate_weights(ctx);

		// Actual flux calculation

		f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
							std::sqrt(ctx.receiversSlopes[ctx.SSj]);
		f_t tQwout = thw * u_w * ctx.receiversDy[ctx.SSj];

		// std::cout << ctx.receiversSlopes[ctx.SSj] << "|";

		// f_t dhw = next.Qw - tQwout;
		// dhw *= this->dt;
		// dhw /= this->con->area(next.node);

		// thw += dhw;
		// tsurf += dhw;

		// if(dhw>1){
		// 	std::cout << dhw << std::endl;

		// // if(std::isfinite(thw) == false){
		// 	std::cout <<  next.Qw <<  "|" << tQwout << std::endl;
		// 	throw std::runtime_error("bite");
		// }

		// this->data->_Qwin[ctx.node] =
		// 	next.Qw; // std::max( this->data->_Qwin[ctx.node], next.Qw);
		this->data->_Qwout[ctx.node] = tQwout;

		// if (thw < 0) {
		// 	tsurf -= thw;
		// 	thw = 0;
		// }

		this->emplace_transfer(next, ctx, dynastack);
	}

	template<class CTX>
	void _calculate_Qwout_for_disconnected_nodes(CTX& ctx)
	{

		f_t& thw = this->data->_hw[ctx.node];
		f_t& tsurf = this->data->_surface[ctx.node];

		// Actual flux calculation

		if (ctx.nn <= 0)
			return;

		int SSj = ctx.idxLowestNeighbour(this->data->_surface);
		if (SSj == -1)
			return;
		int SSi = ctx.neighbours[SSj];
		f_t SSdx = ctx.neighboursDx[SSj];
		f_t SSdy = ctx.neighboursDy[SSj];

		f_t SS = this->data->_surface[ctx.node] - this->data->_surface[SSi];
		SS /= SSdx;

		f_t u_w =
			std::pow(thw, (2. / 3.)) / this->mannings * std::sqrt(std::max(SS, 1e-6));
		// if(std::isfinite(u_w) == false) throw std::runtime_error("here" +
		// std::to_string(SSdx));
		f_t tQwout = thw * u_w * SSdy;
		this->data->_Qwout[ctx.node] = tQwout;
	}

	template<class CELL, class CTX, class Q_t>
	void _subGF_reprocess_node(CELL& next, CTX& ctx, Q_t& dynastack)
	{

		f_t& thw = this->data->_hw[next.node];
		f_t& tsurf = this->data->_surface[next.node];

		// Regulating weights to avoid useless spreading
		this->regulate_weights(ctx);
		f_t tinc = this->hw_increment_LM * (this->data->randu->get() + 0.5);
		thw += tinc;
		tsurf += tinc;

		// this->data->_Qwin[ctx.node] =
		// 	next.Qw; // std::max( this->data->_Qwin[ctx.node], next.Qw);
		// this->data->_Qwin[ctx.node] = 0.;
		// this->data->_Qwout[ctx.node] = 0.;

		this->emplace_transfer(next, ctx, dynastack);
	}

	template<class CTX>
	void regulate_weights(CTX& ctx)
	{
		if (ctx.nr <= 1)
			return;

		if (this->min_part_Qw > 0) {
			f_t nSumSlopesDw = 0.;
			for (int i = 0; i < ctx.nr; ++i) {
				if (ctx.receiversWeights[i] < this->min_part_Qw)
					continue;
				nSumSlopesDw += ctx.receiversSlopes[i] * ctx.receiversDy[i];
			}

			if (nSumSlopesDw != 0. && nSumSlopesDw != ctx.sumslopesdw) {
				for (int i = 0; i < ctx.nr; ++i) {
					if (ctx.receiversWeights[i] < this->min_part_Qw) {
						ctx.receiversWeights[i] = 0.;
						ctx.receiversSlopes[i] = 0.;
					} else
						ctx.receiversWeights[i] =
							ctx.receiversSlopes[i] * ctx.receiversDy[i] / nSumSlopesDw;
				}
			}
			ctx.sumslopesdw = nSumSlopesDw;
		}
	}

	template<class CELL, class CTX, class Q_t>
	void emplace_transfer(CELL& next, CTX& ctx, Q_t& dynastack)
	{
		if (ctx.canout)
			return;

		if (ctx.nr == 0) {
			next.topo = this->data->_surface[next.node];
			dynastack.emplace(next);
			++this->debugyolo;
		} else {
			// std::cout << ctx.nr << "|";
			f_t sumW = 0;
			for (int i = 0; i < ctx.nr; ++i) {
				if (ctx.receiversWeights[i] > 0) {
					sumW += ctx.receiversWeights[i];
					WaCell<i_t, f_t> tnext;
					if (this->time != this->data->_timetracker[ctx.receivers[i]]) {
						tnext = WaCell<i_t, f_t>(ctx.receivers[i],
																		 this->data->_surface[ctx.receivers[i]]);

						this->data->_Qwin[ctx.receivers[i]] +=
							this->data->_Qwin[next.node] * ctx.receiversWeights[i];
					} else {
						tnext = WaCell<i_t, f_t>(ctx.receivers[i],
																		 this->data->_surface[ctx.receivers[i]],
																		 this->data->_Qwin[next.node] *
																			 ctx.receiversWeights[i]);
					}
					// if(this->data->_hw[tnext.node] == 0) std::cout << tnext.Qw <<
					// std::endl;
					// this->data->_debug[tnext.node] = 1;
					dynastack.emplace(tnext);
				}
			}

			// if (abs(sumW - 1) > 1e-5)
			// 	std::cout << sumW << std::endl;
		}
	}

	template<class CTX, class Q_t>
	void emplace_transfer(ExpCell<i_t, f_t>& next, CTX& ctx, Q_t& dynastack)
	{
		if (ctx.canout)
			return;

		if (ctx.nr == 0) {
			next.topo = this->data->_surface[next.node];
			dynastack.emplace(next);
		} else {
			// std::cout << ctx.nr << "|";
			f_t sumW = 0;
			for (int i = 0; i < ctx.nr; ++i) {
				if (ctx.receiversWeights[i] > 0) {
					sumW += ctx.receiversWeights[i];
					auto tnext = ExpCell<i_t, f_t>(ctx.receivers[i],
																				 this->data->_surface[ctx.receivers[i]],
																				 next.Qw * ctx.receiversWeights[i],
																				 next.Qs * ctx.receiversWeights[i]);
					dynastack.emplace(tnext);
				}
			}

			// if (abs(sumW - 1) > 1e-5)
			// 	std::cout << sumW << std::endl;
		}
	}

	void run_subgraphflood_expA() { std::cout << "DEPRECATED" << std::endl; }

	f_t bedrockZatI(int i)
	{
		return this->data->_surface[i] - this->data->_hw[i];
	}

	// Access to steepest rec of node i but reprocessed on the spot (as opposed to
	// preprocessed from a compute() operation)
	int Sreceivers_bedrock_raw(i_t i,
														 std::array<i_t, 8>& arr,
														 std::array<f_t, 8>& arrdx,
														 i_t& trec,
														 f_t& tdx)
	{
		i_t nn = this->con->Neighbours(i, arr);
		nn = this->con->NeighboursDx(i, arrdx);
		trec = -1;
		tdx = this->con->_dx;
		f_t tSS = 0.;
		for (int j = 0; j < nn; ++j) {
			// arr[j] += i;

			f_t ttSS = (this->bedrockZatI(i) - this->bedrockZatI(arr[j])) / arrdx[j];

			if (ttSS > tSS) {
				tSS = ttSS;
				trec = arr[j];
				tdx = arrdx[j];
			}
		}

		return nn;
	}

	f_t _quick_slipos_from_point(i_t node,
															 f_t friction_slope,
															 f_t S_c,
															 f_t Kdep,
															 f_t Hpart,
															 f_t minH)
	{

		std::unordered_set<int> alreadyIn;

		std::queue<int> tQ;
		tQ.emplace(node);
		alreadyIn.insert(node);

		std::array<int, 8> tneigh;
		std::array<f_t, 8> tneighdx;

		f_t HLS = 0;

		// determine rupture_slope
		f_t rupture_slope = 0.;

		i_t nn = this->con->Neighbours(node, tneigh);
		this->con->NeighboursDx(node, tneighdx);

		f_t tsurf = this->bedrockZatI(node);

		for (int j = 0; j < nn; ++j) {
			i_t tn = tneigh[j];
			f_t tdx = tneighdx[j];
			std::cout << this->bedrockZatI(tn) << "vs" << tsurf << std::endl;
			if ((this->bedrockZatI(tn) - tsurf) / tdx > rupture_slope) {
				rupture_slope = (this->bedrockZatI(tn) - tsurf) / tdx;
			}
		}

		if (rupture_slope < friction_slope) {
			std::cout << "No ldsl" << std::endl;
			return 0.;
		}

		rupture_slope += friction_slope;
		rupture_slope /= 2;

		// Gather operation
		while (tQ.empty() == false) {
			i_t next = tQ.front();
			tQ.pop();
			// get donors
			i_t nn = this->con->Neighbours(next, tneigh);
			this->con->NeighboursDx(next, tneighdx);

			tsurf = this->data->_surface[next];

			for (int j = 0; j < nn; ++j) {
				i_t tn = tneigh[j];
				if (alreadyIn.find(tn) != alreadyIn.end())
					continue;
				f_t tdx = tneighdx[j];

				f_t tS = (this->bedrockZatI(tn) - tsurf) / tdx;
				if (tS > rupture_slope) {
					alreadyIn.insert(tn);
					tQ.emplace(tn);
					f_t tgz = tsurf + rupture_slope * tdx;
					HLS += (this->bedrockZatI(tn) - tgz);
					this->data->_surface[tn] = tgz + this->data->_hw[tn];
				}
			}
		}
		// return 0;

		i_t trec;
		f_t tdx;

		i_t Npart = std::floor(HLS / Hpart);
		std::cout << "n parts: " << Npart << std::endl;
		for (int i = 0; i < Npart; ++i) {
			f_t tH = Hpart;
			int nexnode = node;
			// this->Sreceivers_bedrock_raw(nexnode, tneigh, tneighdx, trec, tdx);
			// nexnode = trec;

			while (true) {
				this->Sreceivers_bedrock_raw(nexnode, tneigh, tneighdx, trec, tdx);

				if (can_out(this->data->_boundaries[nexnode]))
					break;

				if (trec == -1) {
					f_t delta_hk = tH * Kdep;
					this->data->_surface[nexnode] += delta_hk;
					tH -= delta_hk;

					if (tH < minH) {
						this->data->_surface[nexnode] += tH;
						break;
					}
					continue;
				}
				if (node == nexnode) {
					nexnode = trec;

					continue;
				}

				// std::cout << "YOLO:" << trec << std::endl;

				f_t tS = (this->bedrockZatI(nexnode) - this->bedrockZatI(trec)) / tdx;

				if (tS < S_c) {
					f_t delta_h = std::min((S_c - tS) * tdx, tH);
					this->data->_surface[nexnode] += delta_h;
					tH -= delta_h;
					// std::cout << delta_h << std::endl;
				}

				f_t delta_hk = tH * Kdep;
				this->data->_surface[nexnode] += delta_hk;
				tH -= delta_hk;
				// std::cout << delta_hk << "::" << std::endl;

				if (tH < minH) {
					this->data->_surface[nexnode] += tH;
					break;
				}

				nexnode = trec;
			}
		}

		return HLS;

		// Sreceivers_raw(i_t i, std::array<i_t, 8>& arr, std::array<f_t, 8>& arrdx,
		// i_t& trec, f_t& tdx)
	}

	void standalone_Qwin() { this->data->_Qwin = this->_standalone_Qwin(); }

	template<class out_t>
	out_t standalone_Qwin_D8()
	{
		auto out = this->_standalone_Qwin_D8();
		return format_output<std::vector<f_t>, out_t>(out);
	}

	std::vector<f_t> _standalone_Qwin()
	{

		this->con->PFcompute_all(false);

		std::vector<f_t> tQw(this->con->nxy(), 0.);

		if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_QW) {
			for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
				tQw[this->input_node_Qw[i]] += this->input_Qw[i];
			}
		}

		std::array<i_t, 8> recs;
		std::array<f_t, 8> dxs;
		std::array<f_t, 8> dys;
		std::array<f_t, 8> slopes;
		for (int i = static_cast<int>(this->data->_stack.size()) - 1; i >= 0; --i) {
			int node = this->data->_stack[i];

			if (nodata(this->data->_boundaries[node]) ||
					can_out(this->data->_boundaries[node]))
				continue;

			if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT)
				tQw[node] += this->Prate * this->con->area(node);

			int nr = this->con->Receivers(node, recs);

			if (nr == 0)
				continue;

			this->con->ReceiversDx(node, dxs);
			this->con->ReceiversDy(node, dys);

			f_t sumslopes = 0.;

			for (int j = 0; j < nr; ++j) {
				slopes[j] =
					std::max(this->data->_surface[node] - this->data->_surface[recs[j]],
									 1e-6) /
					dxs[j] * dys[j];
				sumslopes += slopes[j];
			}

			if (sumslopes <= 0)
				sumslopes = nr;

			for (int j = 0; j < nr; ++j) {
				slopes[j] /= sumslopes;
				tQw[recs[j]] += slopes[j] * tQw[node];
			}
		}

		return tQw;
	}

	std::vector<f_t> _standalone_Qwin_D8()
	{

		this->con->PFcompute_all(false);

		std::vector<f_t> tQw(this->con->nxy(), 0.);

		if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_QW) {
			for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
				tQw[this->input_node_Qw[i]] += this->input_Qw[i];
			}
		}

		for (int i = static_cast<int>(this->data->_stack.size()) - 1; i >= 0; --i) {
			int node = this->data->_stack[i];

			if (nodata(this->data->_boundaries[node]) ||
					can_out(this->data->_boundaries[node]))
				continue;

			if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT)
				tQw[node] += this->Prate * this->con->area(node);

			int trec = con->Sreceivers(node);
			if (trec != node)
				tQw[trec] += tQw[node];
		}

		return tQw;
	}

	void diffuse_topo(f_t mag)
	{
		this->data->_surface = On_gaussian_blur(
			mag, this->data->_surface, this->con->_nx, this->con->_ny);
	}

	void chunk_by_distance_to_outlet(f_t mindist, f_t maxdist)
	{

		// First computing the bloody graph stuffies yolo
		this->con->PFcompute_all(false);

		this->active_nodes = std::vector<i_t>();
		this->active_nodes.reserve(this->con->nxy());

		this->tsg.xtraMask = std::vector<std::uint8_t>(this->con->nxy(), false);

		// Reset the entry points
		this->input_node_Qw_tsg.clear();
		this->input_Qw_tsg.clear();

		// Calculating the distance from outlet
		std::vector<f_t> dist2out =
			_compute_min_distance_from_outlets<i_t, f_t, CONNECTOR_T>(*this->con);

		normalise_vector(dist2out);

		// and masking
		for (int i = 0; i < this->con->nxy(); ++i) {
			// std::cout << dist2out[i] << " ";
			if (dist2out[i] > mindist && dist2out[i] < maxdist) {
				// std::cout << "???" << std::endl;
				this->tsg.xtraMask[i] = true;
				this->active_nodes.emplace_back(i);
			} else
				this->tsg.xtraMask[i] = false;
		}

		std::vector<int> temp(this->con->nxy());
		std::copy(
			this->tsg.xtraMask.begin(), this->tsg.xtraMask.end(), temp.begin());

		this->data->ibag["tsgMask"] = temp;
	}

	void prepare_tsg()
	{

		// process non-local fluxes entries
		// # Get the right Qwin
		// # Note it also computes the TopoSort
		std::vector<f_t> tQw = this->_standalone_Qwin();
		std::vector<std::uint8_t> isDone(this->con->nxy(), false);

		std::unordered_map<i_t, f_t> i2Qw;
		for (int i = 0; i < static_cast<int>(this->input_node_Qw.size()); ++i) {
			i2Qw[this->input_node_Qw[i]] = this->input_Qw[i];
		}

		// flush the old entry points
		this->input_node_Qw.clear();
		this->entry_node_PQ.clear();
		this->input_Qw.clear();

		// reserve new space (let's say 1% of the total number of nodes, to be
		// tested but should not be performance critical)
		// this->input_node_Qw.reserve(static_cast<int>(this->con->nxy() / 100));
		// this->input_Qw.reserve(static_cast<int>(this->con->nxy() / 100));

		// # iterates upstream to downstream to gather the entry points
		std::array<i_t, 8> recs;
		std::array<f_t, 8> dxs, dys, slopes;
		for (int i = static_cast<int>(this->data->_stack.size()) - 1; i >= 0; --i) {
			// ## Getting next upstreamest unprocessed node
			int node = this->data->_stack[i];

			// ## ignoring if already processed or already in the tsg
			if (nodata(this->data->_boundaries[node]) ||
					can_out(this->data->_boundaries[node]) ||
					(this->tsg.xtraMask[node]) || isDone[node])
				continue;

			// if()

			// ## Getting receivers
			int nr = this->con->Receivers(node, recs);
			if (nr == 0)
				continue;

			this->con->ReceiversDx(node, dxs);
			this->con->ReceiversDy(node, dys);

			f_t sumslopes = 0.;

			for (int j = 0; j < nr; ++j) {
				slopes[j] =
					std::max(this->data->_surface[node] - this->data->_surface[recs[j]],
									 1e-6) /
					dxs[j] * dys[j];
				sumslopes += slopes[j];
			}

			if (sumslopes <= 0)
				sumslopes = nr;

			for (int j = 0; j < nr; ++j) {
				slopes[j] /= sumslopes;
			}

			bool is2add = false;
			f_t sumQW = 0;
			// ## Iterating through receivers
			for (int j = 0; j < nr; ++j) {
				int trec = recs[j];
				if (this->tsg.xtraMask[trec]) {

					is2add = true;

					// ## and labelling the ones that are in the tsg while current node
					// isn't
					isDone[trec] = true;
					sumQW += tQw[node] * slopes[j];
				}
			}

			if (is2add && sumQW > 0) {
				this->input_node_Qw.emplace_back(node);
				this->entry_node_PQ.emplace_back(node);
				this->input_Qw.emplace_back(sumQW);
			}
		}

		for (auto& pair : i2Qw) {
			if (this->tsg.xtraMask[pair.first] && isDone[pair.first] == false) {
				// tQw[pair.first] = pair.second;
				this->input_node_Qw.emplace_back(pair.first);
				this->entry_node_PQ.emplace_back(pair.first);
				this->input_Qw.emplace_back(pair.second);
				// isDone[pair.first] = true;
			}
		}

		// // Saving the labelled entry points to the subgraph
		// for (int i = 0; i < this->con->nxy(); ++i) {
		// 	// # First the ones labelled non=locally above
		// 	if (isDone[i]) {
		// 		this->input_node_Qw.emplace_back(i);
		// 		this->entry_node_PQ.emplace_back(i);
		// 		this->input_Qw.emplace_back(tQw[i]);
		// 		// # but also the local ones in case a node is not labelled, in the
		// tsg,
		// 		// and if precipitations rates are global
		// 	} else if (this->tsg.xtraMask[i] &&
		// 						 this->water_input_mode ==
		// 							 WATER_INPUT::PRECIPITATIONS_CONSTANT) {
		// 		this->input_node_Qw.emplace_back(i);
		// 		this->entry_node_PQ.emplace_back(i);
		// 		this->input_Qw.emplace_back(this->Prate * this->con->area(i));
		// 	}
		// }

		// dealing now with the local inputs
		// if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_QW) {
		// 	for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
		// 		if (this->tsg.xtraMask[this->input_node_Qw[i]] &&
		// 				isDone[this->input_node_Qw[i]] == false) {
		// 			if (this->input_Qw[i] > 0) {
		// 				this->input_node_Qw.emplace_back(this->input_node_Qw[i]);
		// 				this->input_Qw.emplace_back(this->input_Qw[i]);
		// 			}
		// 		}
		// 	}
		// }

		// std::cout << input_node_Qw_tsg.size() << " in da Q" << std::endl;
	}

	void run_tinysubgraph()
	{

		fillvec(this->data->_Qwin, 0.);
		fillvec(this->data->_Qwout, 0.);

		this->tsg.build_from_donor_sources(this->input_node_Qw_tsg);

		// std::cout << "Node in tsg::" << this->tsg.nodes.size() << std::endl;

		for (int i = 0; i < static_cast<int>(this->input_node_Qw_tsg.size()); ++i) {
			this->data->_Qwin[this->input_node_Qw_tsg[i]] += this->input_Qw_tsg[i];
		}

		// throw std::runtime_error("StopA");
		std::array<i_t, 8> recs;
		std::array<f_t, 8> dxs;
		std::array<f_t, 8> dys;
		std::array<f_t, 8> slopes;
		for (int i = static_cast<int>(this->tsg.stack.size()) - 1; i >= 0; --i) {
			int node = this->tsg.stack[i];

			int nr = this->con->Receivers(node, recs);
			if (nr == 0) {
				f_t oSS = this->param->gf2Bbval;
				f_t u_w = std::pow(this->data->_hw[node], (2. / 3.)) / this->mannings *
									std::sqrt(std::max(1e-6, oSS));
				this->data->_Qwout[node] = this->data->_hw[node] * u_w * this->con->_dy;
				continue;
			}

			this->con->ReceiversDx(node, dxs);
			this->con->ReceiversDy(node, dys);

			f_t sumslopes = 0.;

			for (int j = 0; j < nr; ++j) {
				slopes[j] =
					std::max(this->data->_surface[node] - this->data->_surface[recs[j]],
									 1e-6) /
					dxs[j] * dys[j];
				sumslopes += slopes[j];
			}

			if (sumslopes <= 0)
				sumslopes = nr;

			f_t sumweight = 0;

			for (int j = 0; j < nr; ++j) {
				slopes[j] /= sumslopes;
				this->data->_Qwin[recs[j]] += slopes[j] * this->data->_Qwin[node];
				sumweight += slopes[j];
			}

			if (std::abs(sumweight - 1.) > 1e-3 && sumweight != 0.)
				throw std::runtime_error("WUf_t");

			f_t SS = (this->data->_surface[node] -
								this->data->_surface[this->con->Sreceivers(node)]) /
							 this->con->SreceiversDx(node);
			f_t SSdy = this->con->SreceiversDy(node);

			f_t u_w = std::pow(this->data->_hw[node], (2. / 3.)) / this->mannings *
								std::sqrt(std::max(1e-6, SS));

			this->data->_Qwout[node] = this->data->_hw[node] * u_w * SSdy;
		}

		for (auto i : this->tsg.nodes) {

			if (can_out(this->data->_boundaries[i]))
				continue;

			f_t& thw = this->data->_hw[i];
			f_t& tsurf = this->data->_surface[i];

			f_t dhw = 0.;
			if (this->sbg_method == SUBGRAPHMETHOD::V1) {
				dhw = this->data->_Qwin[i] - this->data->_Qwout[i];
			} else if (this->sbg_method == SUBGRAPHMETHOD::FILLONLY) {
				dhw = this->data->_Qwin[i];
			}

			dhw *= this->dt;
			dhw /= this->con->area(i);

			thw += dhw;
			tsurf += dhw;

			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}
		}
	}

	void run_tinysubgraph_dyn() { ; }
	// {

	// 	if (this->data->_timetracker.size() == 0) {
	// 		this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
	// 		this->data->_debug = std::vector<f_t>(this->con->nxy(), 0.);
	// 	}

	// 	if (this->param->gf2_morpho && this->data->_Qsin.size() == 0) {
	// 		this->data->_Qsin = std::vector<f_t>(this->con->nxy(), 0.);
	// 		this->data->_Qsout = std::vector<f_t>(this->con->nxy(), 0.);
	// 	}

	// 	std::priority_queue<WaCell<i_t, f_t>,
	// 											std::vector<WaCell<i_t, f_t>>,
	// 											std::less<WaCell<i_t, f_t>>>
	// 		dynastack;

	// 	fillvec(this->data->_Qwin, 0.);
	// 	fillvec(this->isInQ, false);
	// 	if (this->param->gf2_morpho) {
	// 		fillvec(this->data->_Qsin, 0.);
	// 		fillvec(this->data->_Qsout, 0.);
	// 	}

	// 	this->init_dstack_tsbdyn<WaCell<i_t, f_t>,
	// decltype(dynastack)>(dynastack);

	// 	// CT_neighbourer_WaCell<i_t, f_t> ctx;
	// 	CT_neighbours<i_t, f_t> ctx;

	// 	// int nndt = 0;
	// 	// for(auto v:this->data->_boundaries){
	// 	// 	if(nodata(v)) ++nndt;
	// 	// }
	// 	// std::cout << "I have " << nndt << "no data " << std::endl;

	// 	// fillvec(this->data->_vmot_hw,0.);
	// 	this->time += this->dt;
	// 	// std::vector<bool> isdone(this->con->nxy(), false);
	// 	// int ndone = 0;
	// 	// int nredone = 0;

	// 	std::array<i_t, 8> receivers;
	// 	std::array<f_t, 8> receiversWeights;

	// 	f_t sumout = 0.;

	// 	// std::cout << "Starting the process" << std::endl;

	// 	while (dynastack.empty() == false) {

	// 		// Getting the next node
	// 		auto next = this->_dstack_next<WaCell<i_t, f_t>>(dynastack);

	// 		if (this->tsg.xtraMask[next.node] == false) {
	// 			f_t oSS = this->param->gf2Bbval;
	// 			f_t u_w = std::pow(this->data->_hw[next.node], (2. / 3.)) /
	// 								this->mannings * std::sqrt(std::max(1e-6, oSS));
	// 			this->data->_Qwout[next.node] =
	// 				this->data->_hw[next.node] * u_w * this->con->_dy;
	// 			continue;
	// 		}

	// 		this->isInQ[next.node] = false;

	// 		bool ispast = this->data->_timetracker[next.node] != this->time;

	// 		// if(isdone[next.node] == false){
	// 		// 	ndone++;
	// 		// 	isdone[next.node] = true;
	// 		// 	if(ndone % 100 == 0)
	// 		// 		std::cout << ndone << " vs " << nredone << " PQsizzla: " <<
	// 		// dynastack.size() << " this->debugyolo " << this->debugyolo <<
	// 		// std::endl; }else{ 	nredone++;
	// 		// }

	// 		// Updating the timer
	// 		this->data->_timetracker[next.node] = this->time;

	// 		// std::cout << next.node << std::endl;
	// 		ctx.update(next.node, *this->con);
	// 		;

	// 		// std::cout << BC2str(ctx.boundary) << std::endl;

	// 		if (ispast) {
	// 			// this->data->_Qwin[next.node] =
	// std::max(this->data->_Qwin[next.node],
	// 			// next.Qw);
	// 			if (this->data->_Qwin[next.node] > next.Qw) {
	// 				this->data->_Qwin[next.node] += next.Qw;
	// 			} else {
	// 				this->data->_Qwin[next.node] = next.Qw;
	// 			}
	// 		}

	// 		if (can_out(ctx.boundary)) {
	// 			if (ispast)
	// 				sumout += this->data->_Qwin[ctx.node];
	// 			continue;
	// 		}

	// 		if (nodata(ctx.boundary)) {
	// 			// std::cout << "nodata reached" << std::endl;
	// 			continue;
	// 		}

	// 		int nr = 0;
	// 		int SSi = 0;
	// 		f_t SS = 0;
	// 		f_t SSdx = 1.;
	// 		f_t SSdy = 1.;

	// 		// std::cout << "A1" << std::endl;
	// 		this->update_receivers(
	// 			ctx, receivers, receiversWeights, nr, SS, SSdx, SSdy, SSi);
	// 		// std::cout << "A2" << std::endl;

	// 		f_t& thw = this->data->_hw[next.node];
	// 		f_t& tsurf = this->data->_surface[next.node];
	// 		// Actual flux calculation

	// 		f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
	// 							std::sqrt(std::max(1e-6, SS));
	// 		f_t tQwout = thw * u_w * SSdy;

	// 		this->data->_Qwout[ctx.node] = tQwout;

	// 		if (this->param->gf2_morpho) {
	// 			f_t tau = this->param->rho_water * this->param->GRAVITY * SS * thw;
	// 			f_t edot = 0.;
	// 			f_t ddot = 0.;
	// 			if (tau > this->param->tau_c) {
	// 				edot = this->param->ke(ctx.node) *
	// 							 std::pow(tau - this->param->tau_c, this->param->alpha);
	// 			}
	// 			ddot = this->data->_Qsin[ctx.node] / (this->param->kd * SSdy);

	// 			this->data->_Qsout[ctx.node] =
	// 				std::max(static_cast<f_t>(0.),
	// 								 this->data->_Qsin[ctx.node] + edot * SSdx * SSdy -
	// 									 ddot * SSdx * SSdy);
	// 		}

	// 		// if (!ispast) {
	// 		// 	thw += this->hw_increment_LM;
	// 		// 	tsurf += this->hw_increment_LM;
	// 		// }

	// 		// NEED TO CARY ON ADDING MORPHO HERE

	// 		f_t baseQw = (ispast) ? this->data->_Qwin[next.node] : next.Qw;
	// 		f_t baseQs;
	// 		if (this->param->gf2_morpho)
	// 			baseQs = (ispast) ? this->data->_Qsin[next.node] : next.Qs;

	// 		for (int j = 0; j < nr; ++j) {
	// 			int i = receivers[j];
	// 			int rec = ctx.neighbours[i];
	// 			bool tizdone = this->data->_timetracker[rec] == this->time;
	// 			if (tizdone) {
	// 				if (this->param->gf2_morpho)
	// 					dynastack.emplace(WaCell<i_t, f_t>(rec,
	// 																						 this->data->_surface[rec],
	// 																						 baseQw * receiversWeights[j],
	// 																						 baseQs * receiversWeights[j]));
	// 				else
	// 					dynastack.emplace(WaCell<i_t, f_t>(
	// 						rec, this->data->_surface[rec], baseQw * receiversWeights[j]));

	// 			} else {
	// 				if (this->isInQ[rec] == false) {
	// 					dynastack.emplace(WaCell<i_t, f_t>(rec,
	// this->data->_surface[rec])); 					this->isInQ[rec] = true;
	// 				}
	// 				this->data->_Qwin[rec] += baseQw * receiversWeights[j];
	// 				if (this->param->gf2_morpho)
	// 					this->data->_Qsin[rec] += baseQs * receiversWeights[j];
	// 			}
	// 		}
	// 	}
	// 	// std::cout << "done" << std::endl;

	// 	CT_neighbourer_WaCell<i_t, f_t> ctx2;

	// 	for (int i = 0; i < this->con->nxy(); ++i) {

	// 		if (nodata(this->data->_boundaries[i]) ||
	// 				can_out(this->data->_boundaries[i]) ||
	// 				this->sbg_method == SUBGRAPHMETHOD::FILLONLY ||
	// 				this->tsg.xtraMask[i] == false)
	// 			continue;
	// 		if (this->data->_timetracker[i]<this->time&& this->data->_hw[i]> 0) {
	// 			this->data->_timetracker[i] = this->time;
	// 			this->data->_Qwin[i] = 0.;
	// 			ctx.update(i, *this->con);
	// 			this->_calculate_Qwout_for_disconnected_nodes(ctx2);
	// 		}
	// 	}

	// 	if (this->param->gf2_diffuse_Qwin)
	// 		this->data->_Qwin =
	// 			On_gaussian_blur(1., this->data->_Qwin, this->con->_nx,
	// this->con->_ny);

	// 	for (int i = 0; i < this->con->nxy(); ++i) {

	// 		if (nodata(this->data->_boundaries[i]) ||
	// 				can_out(this->data->_boundaries[i]) || this->tsg.xtraMask[i] ==
	// false) 			continue;

	// 		f_t& thw = this->data->_hw[i];
	// 		f_t& tsurf = this->data->_surface[i];

	// 		f_t dhw = 0.;
	// 		if (this->sbg_method == SUBGRAPHMETHOD::V1) {
	// 			dhw = this->data->_Qwin[i] - this->data->_Qwout[i];
	// 		} else if (this->sbg_method == SUBGRAPHMETHOD::FILLONLY) {
	// 			dhw = this->data->_Qwin[i];
	// 		}

	// 		dhw *= this->dt;
	// 		dhw /= this->con->area(i);

	// 		thw += dhw;
	// 		tsurf += dhw;

	// 		if (this->param->gf2_morpho) {
	// 			f_t dz = this->data->_Qsin[i] - this->data->_Qsout[i];
	// 			dz *= this->dt;
	// 			dz /= this->con->area(i);
	// 			if (thw - dz > 0) {
	// 				thw -= dz;
	// 			} else {
	// 				tsurf -= thw - dz;
	// 				thw = 1e-4;
	// 			}
	// 		}

	// 		if (std::isfinite(tsurf) == false) {
	// 			std::cout << this->data->_Qwout[i] << "|" << this->data->_Qwin[i]
	// 								<< std::endl;
	// 			;
	// 			throw std::runtime_error("blug");
	// 		}

	// 		if (thw < 0) {
	// 			tsurf -= thw;
	// 			thw = 0;
	// 		}
	// 	}

	// 	// std::cout << "Sumout: " << sumout << std::endl; ;
	// }

	// Test tinygraph
	void _run_tinysubgraph_v1()
	{

		fillvec(this->data->_Qwin, 0.);
		fillvec(this->data->_Qwout, 0.);

		this->tsg.build_from_donor_sources(this->entry_node_PQ);
		this->data->ibag["debugTSG"] = this->tsg.nodes;

		for (int i = 0; i < static_cast<int>(this->input_node_Qw.size()); ++i) {
			this->data->_Qwin[this->input_node_Qw[i]] += this->input_Qw[i];
		}

		std::array<i_t, 8> recs;
		std::array<f_t, 8> dxs;
		std::array<f_t, 8> dys;
		std::array<f_t, 8> slopes;
		for (int i = static_cast<int>(this->tsg.stack.size()) - 1; i >= 0; --i) {
			int node = this->tsg.stack[i];

			int nr = this->con->Receivers(node, recs);
			if (nr == 0)
				continue;

			this->con->ReceiversDx(node, dxs);
			this->con->ReceiversDy(node, dys);

			f_t sumslopes = 0.;

			for (int j = 0; j < nr; ++j) {
				slopes[j] =
					(this->data->_surface[node] - this->data->_surface[recs[j]]) /
					dxs[j] * dys[j];
				sumslopes += slopes[j];
			}

			if (sumslopes <= 0)
				sumslopes = nr;

			f_t sumweight = 0;

			for (int j = 0; j < nr; ++j) {
				slopes[j] /= sumslopes;
				this->data->_Qwin[recs[j]] += slopes[j] * this->data->_Qwin[node];
				sumweight += slopes[j];
			}

			// if(std::abs(sumweight - 1.) > 1e-3) throw std::runtime_error("WUf_t");

			f_t SS = (this->data->_surface[node] -
								this->data->_surface[this->con->Sreceivers(node)]) /
							 this->con->SreceiversDx(node);
			f_t SSdy = this->con->SreceiversDy(node);

			f_t u_w = std::pow(this->data->_hw[node], (2. / 3.)) / this->mannings *
								std::sqrt(std::max(1e-6, SS));

			this->data->_Qwout[node] = this->data->_hw[node] * u_w * SSdy;
		}

		for (auto i : this->tsg.nodes) {

			if (can_out(this->data->_boundaries[i]))
				continue;

			f_t& thw = this->data->_hw[i];
			f_t& tsurf = this->data->_surface[i];

			f_t dhw = 0.;
			if (this->sbg_method == SUBGRAPHMETHOD::V1) {
				dhw = this->data->_Qwin[i] - this->data->_Qwout[i];
			} else if (this->sbg_method == SUBGRAPHMETHOD::FILLONLY) {
				dhw = this->data->_Qwin[i];
			}

			dhw *= this->dt;
			dhw /= this->con->area(i);

			thw += dhw;
			tsurf += dhw;

			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}
		}
	}

	void bluplift(f_t rate)
	{
		for (int i = 0; i < this->con->nxy(); ++i) {
			if (can_out(this->data->_boundaries[i]) == false) {
				this->data->_surface[i] += rate * this->dt;
			}
		}
	}


	void smooth_bedrock(f_t R){

		this->data->_surface = On_gaussian_blur(R, this->data->_surface, this->con->_nx, this->con->_ny);
	
	}

}; // end of class

} // end of namespace DAGGER
