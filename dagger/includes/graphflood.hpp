#pragma once

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
#include "D8connector.hpp"

#include "graphflood_enums.hpp"

// defines all the format_input depnding on the eventual wrapper
#include "fastflood_recorder.hpp"
#include "wrap_helper.hpp"

namespace DAGGER {

template<class T, class U>
class dynanode
{
public:
	// empty constructor
	dynanode() = default;
	// Constructor by default
	dynanode(T node, U score, U Qw)
	{
		this->node = node;
		this->topo = score;
		this->Qw = Qw;
	};

	dynanode(T node, U score, U Qw, U Qs)
	{
		this->node = node;
		this->topo = score;
		this->Qw = Qw;
		this->Qs = Qs;
	};

	// Node index
	T node;
	// Score data
	U topo;
	U Qw;
	U Qs = 0.;

	// void ingest(dynanode<T,U>& other){this->Qw += other.Qw;}
	void ingest(dynanode<T, U> other)
	{
		this->Qw += other.Qw;
		this->Qs += other.Qs;
	}
};
;

// Custom operator sorting the nodes by scores
template<class T, class U>
inline bool
operator>(const dynanode<T, U>& lhs, const dynanode<T, U>& rhs)
{
	if (lhs.topo != rhs.topo)
		return lhs.topo > rhs.topo;
	else
		return lhs.node > rhs.node;
}

// Custom operator sorting the nodes by topos
template<class T, class U>
inline bool
operator<(const dynanode<T, U>& lhs, const dynanode<T, U>& rhs)
{
	if (lhs.topo != rhs.topo)
		return lhs.topo < rhs.topo;
	else
		return lhs.node < rhs.node;
}

constexpr double GRAVITY = 9.81, FIVETHIRD = 5. / 3., TWOTHIRD = 2. / 3.;

template<class fT, class Graph_t, class Connector_t>
class graphflood
{
public:
	// Underlying connector and graph objects:
	Graph_t* graph;
	Connector_t* connector;

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~ Main data holders ~=~=~=~=~=~~=~=~=~=~=~~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # Hydraulic Surface elevation (bedrock + sed + water)
	std::vector<fT> _surface;

	// # Water depth
	std::vector<fT> _hw;

	// # Water discharge
	std::vector<fT> _Qw;

	// # Sediment discahrge
	std::vector<fT> _Qs;

	// # Sediment height
	std::vector<fT> _hs;

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~= Hydro params and related functions ~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// Model options
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # Is the graph SFD or MFD (or else)
	HYDRO hydromode = HYDRO::GRAPH_MFD;

	// # Is morpho activated (potentially deprecated)
	MORPHO morphomode = MORPHO::NONE;
	// # How to deal with local minima, from graphflood point of view (the
	// rerouting is managed by the graph)
	HYDROGRAPH_LM depression_management = HYDROGRAPH_LM::FILL;

	// # How to partition flow in MFD: proportional to the slope, sqrt, ...
	MFD_PARTITIONNING weight_management = MFD_PARTITIONNING::PROPOSLOPE;

	// # How to manage tthe time step for hydrological calculation
	PARAM_DT_HYDRO mode_dt_hydro = PARAM_DT_HYDRO::COURANT;

	// # How to manage Water inputs
	WATER_INPUT water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;

	// # How to manage boundary conditions
	BOUNDARY_HW boundhw = BOUNDARY_HW::FIXED_HW;

	// # method to check convergence (Experimental)
	CONVERGENCE convergence_mode = CONVERGENCE::NONE;

	// Param value holder
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # Hydro timestep
	std::vector<fT> _dt_hydro = { 1e-3 };

	// # courant stuff
	// ## actual courant number
	fT courant_number = 5e-4;
	// ## maximum dt under courant conditions
	fT max_courant_dt_hydro = 1e3;
	// ## minimum dt under courant conditions
	fT min_courant_dt_hydro = 1e-4;
	// ## actual time step (default -1 means not calculated yet)
	fT courant_dt_hydro = -1;

	// # Hydro Dt regulator
	fT maxdHw = 1e2; // set to very high values by default, as we do not want to
									 // impose artificial regulation
	void set_maxdHw(fT val) { this->maxdHw = val; }
	fT mindHw = -1e2; // set to very high values by default, as we do not want to
										// impose artificial regulation
	void set_mindHw(fT val)
	{
		if (val < 0) {
			this->mindHw = val;
		} else {
			throw std::runtime_error("mindHw needs to be negative (decrement)");
		}
	}

	// # minimum slope for manning's calculation
	fT minslope = 0.;

	// # Stochasticity management (experimental)
	// ## Enable stochasticity added to the slope
	bool stochaslope = false;
	// # Stochasticity magnitude
	fT stochaslope_coeff = 1.;
	// # Stochasticity magnitude
	fT bou_fixed_val = 0.;

	// # Convergence checekrs
	// ## nodes used to check convergence
	std::vector<int> conv_nodes;
	// ## Initial discharge at convergence nodes
	std::vector<fT> conv_ini_Qw;
	// ## tracking delta hw at convergence nodes
	std::vector<std::vector<fT>> conv_dhw;
	// ## tracking discharge ratio at convergence nodes
	std::vector<std::vector<fT>> conv_Qr;
	// ## tracking which nodes have converged
	std::vector<std::uint8_t> converged;

	// # Manning's constant
	// ## is manning spatially constant
	bool mode_mannings = false;
	// ## Data holder
	std::vector<fT> _mannings = { 0.033 };

	// # Input discharge water
	// ## As precipitations
	std::vector<fT> _precipitations = { 1e-4 };
	// ## As discrete entry points
	// ### Entry points
	std::vector<int> _water_entry_nodes;
	// ### Entry discahrge
	std::vector<fT> _water_entries;

	// # Topological number correcting for MFD (see paper)
	fT topological_number = 4. / 8.;

	// # Use flow depth sensu Caesar Lisflood (only taking the height where water
	// can flow)
	bool hflow = false;

	// # Is the model in hydrostationary mode
	bool hydrostationary = true;
	// # Experimental
	fT Qwin_crit = 0.;

	// # Debugging temp holder
	bool debugntopo = true;
	std::vector<fT> DEBUGNTOPO;
	fT debug_CFL = 0.;

	// # Create random device and Mersenne Twister engine
	std::random_device rd;
	std::mt19937 gen;
	std::uniform_int_distribution<> dis;

	// Randomiser helper
	DAGGER::easyRand randu;

	void gen_rdid() { this->graph->gen_rdid(); };
	fT get_rdid() { return this->graph->get_rdid(); };

	// # Precipitons stuff (will wvolve or get remove or remain experimental
	// anyway)
	std::vector<fT> last_Smax;
	std::vector<fT> last_dt_prec;
	std::vector<fT> last_dt_prec_e;
	std::vector<fT> last_sw_prec;
	std::vector<fT> last_dx_prec;
	fT current_dt_prec = 0.;
	fT Vp = 0.;

	// Setters and Getters
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # Export flow depth
	template<class out_t>
	out_t get_hw()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_hw);
	}

	// # Export Surface
	template<class out_t>
	out_t get_surface_topo()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_surface);
	}

	// # Export topographic surface
	template<class out_t>
	out_t get_bedrock_topo()
	{
		std::vector<fT> diff(this->_surface);

		for (int i = 0; i < this->graph->nnodes; ++i)
			diff[i] -= this->_hw[i];

		return DAGGER::format_output<std::vector<fT>, out_t>(diff);
	}

	// # Export Qwin
	template<class out_t>
	out_t get_Qwin()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_Qw);
	}

	// # Export Qwin
	template<class out_t>
	out_t get_Qs()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_Qs);
	}

	// # Export Debug info (changes from time to time)
	template<class out_t>
	out_t get_SSTACKDEBUG()
	{
		return DAGGER::format_output<std::vector<size_t>, out_t>(
			this->graph->Sstack);
	}

	// # Manning's setter
	void set_mannings(fT val) { this->_mannings = { val }; }

	// # Hydro constant timestep getter
	fT get_dt_hydro() const { return this->_dt_hydro[0]; }
	// # Setting the hydrologic time step to courant number
	void enable_courant_dt_hydro()
	{
		this->mode_dt_hydro = PARAM_DT_HYDRO::COURANT;
	}
	// # Setting the hydrologic time step to courant number
	void disable_courant_dt_hydro()
	{
		this->mode_dt_hydro = PARAM_DT_HYDRO::CONSTANT;
	}

	// # Setting the courant number
	void set_courant_number(fT val) { this->courant_number = val; };

	// # Setting the max dt courant can set
	void set_max_courant_dt_hydro(fT val) { this->max_courant_dt_hydro = val; };

	// # Setting the min dt courant can set
	void set_min_courant_dt_hydro(fT val) { this->min_courant_dt_hydro = val; };

	// # Getting the courant dt used for the last time step
	fT get_courant_dt_hydro() const { return this->courant_dt_hydro; }

	// # Experimental
	int n_nodes_convergence() const
	{
		return static_cast<int>(conv_nodes.size());
	}
	int n_stack_convergence() const
	{
		if (this->convergence_mode == CONVERGENCE::ALL ||
				this->convergence_mode == CONVERGENCE::QWR)
			return static_cast<int>(this->conv_Qr[0].size());
		if (this->convergence_mode == CONVERGENCE::ALL ||
				this->convergence_mode == CONVERGENCE::DHW)
			return static_cast<int>(this->conv_dhw[0].size());
		return 0;
	}

	// # Set topological number (deprecated, now dynamically calculated)
	void set_topological_number(fT val) { this->topological_number = val; };

	// # get topological number (deprecated, now dynamically calculated)
	fT get_topological_number() const { return this->topological_number; };

	// # Set partition method
	void set_partition_method(MFD_PARTITIONNING& tmffmeth)
	{
		this->weight_management = tmffmeth;
	}

	// # Set the stochasticity coeff (de facto enable stochaslope stuff)
	void set_stochaslope(fT val)
	{
		this->stochaslope = true;
		this->stochaslope_coeff = val;
		this->connector->set_stochaticiy_for_SFD(val);
	}

	// # Disable the stochasticity on slope calculations
	void disable_stochaslope() { this->stochaslope = false; }

	// # Set the model to hydrostationary mode
	void enable_hydrostationary() { this->hydrostationary = true; };

	// # Set the model to dynamic mode
	void disable_hydrostationary() { this->hydrostationary = false; };

	// # (Past experiment, does not work) Set the model to Qwincrit
	void set_Qwin_crit(fT val)
	{
		this->Qwin_crit = val;
		this->hydromode = HYDRO::GRAPH_HYBRID;
	}

	// # set the boundaries to a constant flow depth
	void set_fixed_hw_at_boundaries(fT val)
	{
		this->boundhw = BOUNDARY_HW::FIXED_HW;
		this->bou_fixed_val = val;
	}

	// # set the boundaries to a constant hydraulic slope
	void set_fixed_slope_at_boundaries(fT val)
	{
		this->boundhw = BOUNDARY_HW::FIXED_SLOPE;
		this->bou_fixed_val = val;
	}

	// # get the array of topological numbers
	std::vector<fT> get_nT() { return this->DEBUGNTOPO; }

	// # Setting the timestep for constant hydro (does not have any effect if
	// courant is on)
	void set_dt_hydro(fT tdt) { this->_dt_hydro = { tdt }; }

	// # Setting dicrete entry point for water discharge
	template<class out_t, class in_t>
	void set_water_input_by_entry_points(out_t& hw_entry, in_t& hw_indices)
	{
		// preformatting the inputs
		auto tin = DAGGER::format_input(hw_entry);
		auto tin_idx = DAGGER::format_input(hw_indices);

		// Setting the general mode
		this->water_input_mode = WATER_INPUT::ENTRY_POINTS_H;
		this->_water_entry_nodes = DAGGER::to_vec(tin_idx);
		this->_water_entries = DAGGER::to_vec(tin);

		// Correcting boundary conditions
		for (size_t i = 0; i < this->_water_entry_nodes.size(); ++i) {
			this->connector->boundaries.codes[this->_water_entry_nodes[i]] =
				BC::FORCE_IN;
		}

		this->_Qw_adder = std::vector<fT>(this->connector->nxy(), 0.);
		for (int i = 0; i < this->_water_entry_nodes.size(); ++i) {
			this->_Qw_adder[this->_water_entry_nodes[i]] =
				this->_water_entries[i] *
				this->connector->get_area_at_node(this->_water_entry_nodes[i]);
		}
		this->init_Qw();
	}

	// # Setting global constant precipitation rates
	void set_water_input_by_constant_precipitation_rate(fT precipitations)
	{
		this->water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;
		this->_precipitations = { precipitations };
	}

	// # Setting global constant precipitation rates
	template<class out_t>
	void set_water_input_by_variable_precipitation_rate(out_t& precipitations)
	{
		this->water_input_mode = WATER_INPUT::PRECIPITATIONS_VARIABLE;
		auto tin = format_input(precipitations);
		this->_precipitations = DAGGER::to_vec(tin);
	}

	// # Setting Flow partition to MFD
	void enable_MFD() { this->hydromode = HYDRO::GRAPH_MFD; }

	// # Setting Flow partition to SFD
	void enable_SFD() { this->hydromode = HYDRO::GRAPH_SFD; }

	// # Local Minima will be filled with water automatically
	void fill_minima() { this->depression_management = HYDROGRAPH_LM::FILL; }

	// # Flow will route through local minimas, negative slopes will be treated as
	// 1e-6
	void reroute_minima()
	{
		this->depression_management = HYDROGRAPH_LM::REROUTE;
	}

	// # Flow stops at local minimas
	void ignore_minima() { this->depression_management = HYDROGRAPH_LM::IGNORE; }

	// Optional hydraulic model monitoring
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # Tracking the total input of water added to the model
	fT tot_Qw_input = 0;
	fT get_tot_Qw_input() const { return this->tot_Qw_input; }

	// # Tracking the total Qwin exiting the model (for mass balance cheking for
	// example)
	fT tot_Qwin_output = 0;
	fT get_tot_Qwin_output() const { return this->tot_Qwin_output; }

	// # Tracking the total Qwout exiting the model (for mass balance cheking for
	// example)
	fT tot_Qw_output = 0;
	fT get_tot_Qw_output() const { return this->tot_Qw_output; }

	// # Tracking the total Qs exiting the model (for mass balance cheking for
	// example)
	fT tot_Qs_output = 0;
	fT get_tot_Qs_output() const { return this->tot_Qs_output; }

	// # Qw_out recorder for the whole landscape
	// ## switch activating the recording
	bool record_Qw_out = false;
	// ## Data holder
	std::vector<fT> _rec_Qwout;
	// ## Switching on the recording
	void enable_Qwout_recording() { this->record_Qw_out = true; };
	// ## Switching off the recording
	void disable_Qwout_recording()
	{
		this->record_Qw_out = false;
		this->_rec_Qwout.clear();
	};
	// ## Exporting the output
	template<class out_t>
	out_t get_Qwout_recording()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_rec_Qwout);
	}

	// # Hydraulic slope recorder for the whole landscape
	// ## switch activating the recording
	bool record_Sw = false;
	// ## Data holder
	std::vector<fT> _rec_Sw;
	// ## Switching on the recording
	void enable_Sw_recording() { this->record_Sw = true; };
	// ## Switching off the recording
	void disable_Sw_recording()
	{
		this->record_Sw = false;
		this->_rec_Sw.clear();
	};
	// ## Exporting the output
	template<class out_t>
	out_t get_Sw_recording()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_rec_Sw);
	}

	// # Flow depth increment recorder for the whole landscape
	// ## switch activating the recording
	bool record_dhw = false;
	// ## Data holder
	std::vector<fT> _rec_dhw;
	// ## Switching on the recording
	void enable_dhw_recording() { this->record_dhw = true; };
	// ## Switching off the recording
	void disable_dhw_recording()
	{
		this->record_dhw = false;
		this->_rec_dhw.clear();
	};
	// ## Exporting the output
	template<class out_t>
	out_t get_dhw_recording()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_rec_dhw);
	}

	// # Local minima initial filling recorder for the whole landscape
	// ## switch activating the recording
	bool record_filling = false;
	// ## Data holder
	std::vector<fT> _rec_filling;
	// ## Switching on the recording
	void enable_filling_recording() { this->record_filling = true; };
	// ## Switching off the recording
	void disable_filling_recording()
	{
		this->record_filling = false;
		this->_rec_filling.clear();
	};
	// ## Exporting the output
	template<class out_t>
	out_t get_filling_recording()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_rec_filling);
	}

	// # FLow vector recorder for the whole landscape
	// ## switch activating the recording
	bool record_flowvec = false;
	// ## Data holder
	std::vector<fT> _rec_flowvec;
	// ## Switching on the recording
	void enable_flowvec_recording() { this->record_flowvec = true; };
	// ## Switching off the recording
	void disable_flowvec_recording()
	{
		this->record_flowvec = false;
		this->_rec_flowvec.clear();
	};
	// ## Exporting the output
	template<class out_t>
	out_t get_flowvec_recording()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_rec_flowvec);
	}

	// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	// =~=~=~=~=~= Morpho params and related functions =~=~=~=~=~=~=~
	// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	//                   . - ~ ~ ~ - .
	//         _     .-~               ~- .
	//         \ `..~                       ` .
	//           }  }              /       \   \
// (\   \\ \~^..'                 |       }  \
// \`.-~  o      /       }       |        /  \
// (__          |       /        |       /    `.
	// `- - ~ ~ -._|      /_ - ~ ~ ^|   ____/- _    `.
	//         |     /          |     /     ~-.     ~- _
	//         |_____|          |_____|         ~ - . _ _~_-_
	// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
	// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

	// The whole morpho section is highly WIP and not really detailed yet as it
	// will really evolve
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # method to manage time steps for the morpho (WIP)
	PARAM_DT_MORPHO mode_dt_morpho = PARAM_DT_MORPHO::HYDRO;

	// # method to manage sediment input (WIP)
	SED_INPUT sed_input_mode = SED_INPUT::NONE;

	fT tau_max = 0.;
	fT current_dt_prec_e = 0.;
	fT Vps = 0.;

	std::vector<int> _sed_entry_nodes;
	std::vector<fT> _sed_entries;

	bool record_edot = false;
	std::vector<fT> _rec_edot;
	void enable_edot_recording() { this->record_edot = true; };
	void disable_edot_recording()
	{
		this->record_edot = false;
		this->_rec_edot.clear();
	};
	template<class out_t>
	out_t get_edot_recording()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_rec_edot);
	}

	bool record_ddot = false;
	std::vector<fT> _rec_ddot;
	void enable_ddot_recording() { this->record_ddot = true; };
	void disable_ddot_recording()
	{
		this->record_ddot = false;
		this->_rec_ddot.clear();
	};
	template<class out_t>
	out_t get_ddot_recording()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->_rec_ddot);
	}

	// # a exponent for erosion
	bool mode_aexp = false;
	std::vector<fT> _aexp = { 1.5 };

	// # ke (coefficient for erosion)
	PARAM_KE mode_ke = PARAM_KE::CONSTANT;
	std::vector<fT> _ke = { 1e-4 };

	// # ke (coefficient for erosion)
	bool mode_ke_lateral = false;
	std::vector<fT> _ke_lateral = { 0.1 };

	// # ke (coefficient for erosion)
	bool mode_kd = false;
	std::vector<fT> _kd = { 100 };
	;

	// # kd (coefficient for erosion)
	bool mode_kd_lateral = false;
	std::vector<fT> _kd_lateral = { 0.1 };
	;

	// # ke (coefficient for erosion)
	bool mode_rho = false;
	std::vector<fT> _rho = { 1000 };
	;

	// # ke (coefficient for erosion)
	bool mode_tau_c = false;
	std::vector<fT> _tau_c = { 6 };
	;

	// # ke (coefficient for erosion)
	fT dt_morpho_multiplier = 1.;
	void set_dt_morpho_multiplier(fT val) { this->dt_morpho_multiplier = val; }
	std::vector<fT> _dt_morpho = { 1e-3 };

	void set_dt_morpho(fT tdt)
	{
		this->mode_dt_morpho = PARAM_DT_MORPHO::CONSTANT;
		this->_dt_morpho = { tdt };
	}

	template<class out_t, class in_t>
	void set_sed_input_by_entry_points(out_t& sed_entry, in_t& sed_indices)
	{
		// preformatting the inputs
		auto tin = DAGGER::format_input(sed_entry);
		auto tin_idx = DAGGER::format_input(sed_indices);

		// Setting the general mode
		this->sed_input_mode = SED_INPUT::ENTRY_POINTS_Q;

		this->_sed_entry_nodes = DAGGER::to_vec(tin_idx);
		this->_sed_entries = DAGGER::to_vec(tin);

		this->_Qs_adder = std::vector<fT>(this->connector->nxy(), 0.);
		for (int i = 0; i < this->_sed_entries.size(); ++i) {
			this->_Qs_adder[this->_sed_entry_nodes[i]] =
				this->_sed_entries[i] * this->connector->dy;
		}
	}

	void enable_morpho() { this->morphomode = MORPHO::TL; }
	void disable_morpho() { this->morphomode = MORPHO::NONE; }

	void set_single_aexp(fT ta)
	{
		this->mode_aexp = false;
		this->_aexp = { ta };
	}
	void set_single_ke(fT ta)
	{
		this->mode_ke = PARAM_KE::CONSTANT;
		this->_ke = { ta };
	}
	void set_single_ke_lateral(fT ta)
	{
		this->mode_ke_lateral = false;
		this->_ke_lateral = { ta };
	}
	void set_single_kd(fT ta)
	{
		this->mode_kd = false;
		this->_kd = { ta };
	}
	void set_single_kd_lateral(fT ta)
	{
		this->mode_kd_lateral = false;
		this->_kd_lateral = { ta };
	}
	void set_single_tau_c(fT ta)
	{
		this->mode_tau_c = false;
		this->_tau_c = { ta };
	}

	template<class topo_t>
	void set_variable_ke(topo_t& variable_ke)
	{
		auto tke = format_input(variable_ke);
		this->_ke = to_vec(tke);
		this->mode_ke = PARAM_KE::VARIABLE;
	}

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~= dynamic graph experiment ~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	std::priority_queue<dynanode<int, fT>,
											std::vector<dynanode<int, fT>>,
											std::less<dynanode<int, fT>>>
		dynastack;

	std::queue<int> dynaqueue;

	// incrementing topography in local minimas
	fT dynincr_LM = 1e-4;

	fT glob_dynatime = 0.;
	void set_glob_dynatime(fT val) { this->glob_dynatime = val; };
	fT get_glob_dynatime(fT val) { return this->glob_dynatime; };
	std::vector<fT> dttracker;
	std::vector<fT> _Qw_adder;
	std::vector<fT> _Qs_adder;

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~= Functions setting up the model ~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// Constructors
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # Model constructor without any params (need to seed the randomiser though)
	graphflood() { this->gen = std::mt19937(this->rd()); };

	// # Model classic constructor from existing graph and connector objects
	graphflood(Graph_t& graph, Connector_t& connector)
	{
		// Ingesting graph and connectors
		this->graph = &graph;
		this->connector = &connector;

		// Seeding the random engine
		this->gen = std::mt19937(this->rd());
	}

	// # Setting initial topo (or updating topo from external sources)
	template<class topo_t>
	void set_topo(topo_t& topo)
	{
		// formatting input from higher level languages to c++ veclike stuff
		auto tin = DAGGER::format_input(topo);
		std::vector<fT> temp = DAGGER::to_vec(tin);

		// Initialising flow depth if not done yet
		bool inithw = false;
		if (this->_hw.size() == 0) {
			inithw = true;
			this->_hw = std::vector<fT>(this->graph->nnodes, 0);
		}

		// Fedding the surface (and correcting to eventual existing flow depth)
		this->_surface = std::vector<fT>(temp);
		if (inithw == false) {
			for (int i = 0; i < this->graph->nnodes; ++i) {
				this->_surface[i] += this->_hw[i];
			}
		}

		// Done
	}

	// # Setting the flow depth to an external value
	template<class topo_t>
	void set_hw(topo_t& thw)
	{
		// formatting input from higher level languages to c++ veclike stuff
		auto tin = DAGGER::format_input(thw);
		std::vector<fT> temp = DAGGER::to_vec(tin);

		// case one: no preexisting flow depth, adding hw to hydraulic surface
		if (this->_hw.size() == 0) {
			this->_hw = std::move(temp);
			for (int i = 0; i < this->graph->nnodes; ++i)
				this->_surface[i] += this->_hw[i];
		}
		// case two: preexisting flow depth, correcting to hw
		else {
			for (int i = 0; i < this->graph->nnodes; ++i)
				this->_surface[i] += temp[i] - this->_hw[i];
			this->_hw = std::move(temp);
		}
	}

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~= Generic accessors for internal use ~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// Surfaces and Hydraulic components
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # get local flow depth at node
	fT hw(int i) { return this->_hw[i]; }

	// # get local input flow discharge at node
	fT Qw(int i) { return this->_Qw[i]; }

	// # get local sediment input discharge at node
	fT Qs(int i) { return this->_Qs[i]; }

	// # get local topographic + water surface at node
	fT surface(int i) { return this->_surface[i]; }

	// # Manning's coefficient at node
	fT mannings(int i)
	{
		if (this->mode_mannings)
			return this->_mannings[i];
		else
			return this->_mannings[0];
	}

	// # Precipitation rate at node
	fT precipitations(int i)
	{
		if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE)
			return this->_precipitations[i];
		else
			return this->_precipitations[0];
	}

	// small helper function returning node index from the right stack
	int get_istack_node(int i)
	{
		if (this->hydromode != HYDRO::GRAPH_SFD)
			return this->graph->stack[i];
		else
			return this->graph->Sstack[i];
	}

	fT get_Sw(int lix, fT minslope)
	{
		int from, to;
		this->connector->from_to_from_link_index(lix, from, to);
		return std::max((this->_surface[from] - this->_surface[to]) /
											this->connector->get_dx_from_links_idx(lix),
										minslope);
	}

	fT get_Sw(int node, int rec, fT dx, fT minslope)
	{
		return std::max((this->_surface[node] - this->_surface[rec]) / dx,
										minslope);
	}

	fT get_Sw(int node, int rec, fT dx)
	{
		return (this->_surface[node] - this->_surface[rec]) / dx;
	}

	fT get_Stopo(int node, int rec, fT dx)
	{
		return (this->_surface[node] - this->_hw[node] -
						(this->_surface[rec] - this->_hw[rec])) /
					 dx;
	}

	// Morpho params (WIP, no working version yet)
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # get topographic + water surface
	fT aexp(int i)
	{
		if (this->mode_aexp)
			return this->_aexp[i];
		else
			return this->_aexp[0];
	}

	fT ke(int i)
	{
		if (PARAM_KE::VARIABLE == this->mode_ke)
			return this->_ke[i];
		else
			return this->_ke[0];
	}

	fT ke_lateral(int i)
	{
		if (this->mode_ke_lateral)
			return this->_ke_lateral[i];
		else
			return this->_ke_lateral[0];
	}

	fT kd(int i)
	{
		if (this->mode_kd)
			return this->_kd[i];
		else
			return this->_kd[0];
	}

	fT kd_lateral(int i)
	{
		if (this->mode_kd_lateral)
			return this->_kd_lateral[i];
		else
			return this->_kd_lateral[0];
	}

	fT dt_hydro(int i)
	{
		if (this->mode_dt_hydro == PARAM_DT_HYDRO::VARIABLE)
			return this->_dt_hydro[i];
		else if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT) {
			return this->courant_dt_hydro;
		} else
			return this->_dt_hydro[0];
	}

	fT dt_morpho(int i)
	{
		return this->dt_hydro(i) * this->dt_morpho_multiplier;
		// if (this->mode_dt_morpho == PARAM_DT_MORPHO::VARIABLE)
		// 	return this->_dt_morpho[i];
		// else if (this->mode_dt_morpho == PARAM_DT_MORPHO::HYDRO)
		// 	return this->dt_hydro(i);
		// else
		// 	return this->_dt_morpho[0];
	}

	fT tau_c(int i)
	{
		if (this->mode_tau_c)
			return this->_tau_c[i];
		else
			return this->_tau_c[0];
	}

	fT rho(int i)
	{
		if (this->mode_rho)
			return this->_rho[i];
		else
			return this->_rho[0];
	}

	void block_uplift(fT val)
	{
		for (int i = 0; i < this->connector->nnodes; ++i) {
			if (this->connector->boundaries.can_out(i) == false)
				this->_surface[i] += this->dt_morpho(i) * val;
		}
	}

	template<class topo_t>
	void variable_uplift(topo_t& ival)
	{
		auto val = format_input(ival);
		for (int i = 0; i < this->connector->nnodes; ++i) {
			// if(this->connector->boundaries.can_out(i) == false)
			this->_surface[i] += this->dt_morpho(i) * val[i];
		}
	}

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~= Running and helper functions ~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// # MAin running function
	void run()
	{
		// Saving the topological number if needed
		if (this->debugntopo)
			this->DEBUGNTOPO = std::vector<fT>(this->connector->nnodes, 0);

		this->debug_CFL = 0.;

		this->tau_max = 0.;

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition
		// and other monitoring features:
		this->init_Qw();

		// reinitialising the sediments if needed
		if (this->morphomode != MORPHO::NONE)
			this->init_Qs();

		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// To be used if courant dt hydro is selected
		fT tcourant_dt_hydro = std::numeric_limits<fT>::max();

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> vmot, vmot_hw(this->graph->nnodes, 0.);

		// -> Only initialising vertical motions for the bedrock if morpho is on
		if (this->morphomode != MORPHO::NONE)
			vmot = std::vector<fT>(this->graph->nnodes, 0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT, 8> weights, slopes;

		// main loop
		for (int i = this->graph->nnodes - 1; i >= 0; --i) {

			// Getting next node in line
			int node = this->get_istack_node(i);

			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need
			// to be processed THis is where all the boundary treatment happens, if
			// you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers, vmot_hw))
				continue;

			// if(this->hydrostationary == false)
			// {
			// 	if(this->dt_hydro(node) > 0)
			// 		this->_hw[node] += this->dt_hydro(node) *
			// this->precipitations(node); 	else 		this->_hw[node] += 1e-3
			// * this->precipitations(node);

			// }

			// CFL calculator
			// fT sum_ui_over_dxi = 0.;

			// Deprecated test to switch dynamically between SFD and MFD (does not
			// really add anything and is buggy in rivers) if(this->hydromode ==
			// HYDRO::GRAPH_HYBRID) 	SF = this->_Qw[node] < this->Qwin_crit;

			// Getting the receivers
			int nrecs;
			if (SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node, receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			// NOTE:
			//  No need to calculate the topological number anymore, but keeping it
			//  for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al(node,
																					SF,
																					Smax,
																					slopes,
																					weights,
																					nrecs,
																					receivers,
																					recmax,
																					dx,
																					dw0max,
																					topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], TWOTHIRD);
			// std:: cout << "pohw:" << pohw << "|" << this->_hw[node] << std::endl;

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax / this->mannings(node);
			// Volumetric discahrge
			fT Qwout = dw0max * this->_hw[node] * u_flow;
			// std::cout << Qwout << "|";

			// Eventually recording Smax
			if (this->record_Sw)
				this->_rec_Sw[node] = Smax;

			// temp calc for courant
			// if(this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && this->_hw[node] >
			// 0 && dx >0) 	sum_ui_over_dxi = u_flow/dx;

			// Automates the computation of morpho LEM if needed
			this->_compute_morpho(node, recmax, dx, Smax, vmot);

			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(
				nrecs, recmax, node, SF, receivers, weights, Qwin, Qwout);

			// Computing courant based dt
			if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && u_flow > 0) {
				fT provisional_dt = (this->courant_number * dx) / u_flow;
				provisional_dt = std::min(provisional_dt, this->max_courant_dt_hydro);
				provisional_dt = std::max(provisional_dt, this->min_courant_dt_hydro);
				tcourant_dt_hydro = std::min(provisional_dt, tcourant_dt_hydro);
			}

			// computing hydro vertical motion changes for next time step
			if (this->hydrostationary)
				vmot_hw[node] +=
					(this->_Qw[node] - Qwout) / this->connector->get_area_at_node(node);
			else
				vmot_hw[node] +=
					(this->_Qw[node] - Qwout) / this->connector->get_area_at_node(node);

			if (this->record_Qw_out)
				this->_rec_Qwout[node] += Qwout;
		}

		// if(this->tau_max > 500)
		// 	std::cout << "WARNING::tau_max is " << this->tau_max << std::endl;

		// END OF MAIN LOOP

		// Computing final courant based dt

		if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT) {
			if (this->courant_dt_hydro == -1)
				this->courant_dt_hydro = 1e-3;
			else if (tcourant_dt_hydro > 0 &&
							 tcourant_dt_hydro != std::numeric_limits<fT>::max())
				this->courant_dt_hydro = tcourant_dt_hydro;
		}

		if (this->convergence_mode == CONVERGENCE::ALL ||
				this->convergence_mode == CONVERGENCE::DHW) {
			for (int i = 0; i < this->n_nodes_convergence(); ++i)
				this->conv_dhw[i].emplace_back(
					vmot_hw[this
										->conv_nodes[i]]); // /this->dt_hydro(this->conv_nodes[i]));
		}

		if (this->convergence_mode == CONVERGENCE::ALL ||
				this->convergence_mode == CONVERGENCE::QWR) {
			for (int i = 0; i < this->n_nodes_convergence(); ++i)
				this->conv_Qr[i].emplace_back(this->_rec_Qwout[this->conv_nodes[i]] /
																			this->_Qw[this->conv_nodes[i]]);
		}

		// Applying vmots with the right dt
		this->_compute_vertical_motions(vmot_hw, vmot);
	}

	void init_Qw()
	{
		// resetting Qwin (needed to be stored at all time)
		this->_Qw = std::vector<fT>(this->graph->nnodes, 0.);

		// resetting global monitors (low cost)
		this->tot_Qw_input = 0;
		this->tot_Qw_output = 0;
		this->tot_Qwin_output = 0;

		if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT ||
				this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE) {
			for (int i = 0; i < this->graph->nnodes; ++i) {
				if (this->connector->boundaries.can_give(i) &&
						this->connector->flow_out_or_pit(i) == false) {
					this->_Qw[i] +=
						this->precipitations(i) * this->connector->get_area_at_node(i);
					this->tot_Qw_input += this->_Qw[i];
				}
			}
		} else if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_H) {
			for (size_t i = 0; i < this->_water_entries.size(); ++i) {
				int node = this->_water_entry_nodes[i];
				this->_Qw[node] +=
					this->_water_entries[i] * this->connector->get_area_at_node(node);
				this->tot_Qw_input += this->_Qw[node];
				// std::cout << node << " is given " << this->_water_entries[i] *
				// this->connector->get_area_at_node(node);
			}
		}

		if (this->record_Qw_out)
			this->_rec_Qwout = std::vector<fT>(this->graph->nnodes, 0.);

		if (this->record_Sw)
			this->_rec_Sw = std::vector<fT>(this->graph->nnodes, 0.);
		if (this->record_dhw)
			this->_rec_dhw = std::vector<fT>(this->graph->nnodes, 0.);

		if (record_flowvec)
			this->_rec_flowvec = std::vector<fT>(this->graph->nnodes * 2, 0.);
	}

	void init_Qs()
	{
		this->_Qs = std::vector<fT>(this->graph->nnodes, 0.);
		if (this->sed_input_mode == SED_INPUT::ENTRY_POINTS_Q) {
			for (size_t i = 0; i < this->_sed_entries.size(); ++i) {
				int node = this->_sed_entry_nodes[i];
				this->_Qs[node] += this->_sed_entries[i] * this->connector->dy;
				// this->tot_Qw_input += this->_Qs[node];
			}
		}

		// Initialising the Qs recorders

		if (this->record_edot)
			this->_rec_edot = std::vector<fT>(this->connector->nnodes, 0.);

		this->tot_Qs_output = 0.;
	}

	/// Automates the processing of the graph
	/// Function of global parameters for flow topology, local minima management,
	/// ...
	void graph_automator()
	{

		// is SS?
		bool only_SD = (this->hydromode == HYDRO::GRAPH_SFD);

		// making sure it has the right depression solver (SHOULD BE MOVED TO THE
		// GRAPH MANAGEMENT LATER)
		if (this->depression_management == HYDROGRAPH_LM::IGNORE)
			this->graph->set_LMR_method(DEPRES::none);

		if (this->record_filling)
			this->_rec_filling = std::vector<fT>(this->graph->nnodes, 0.);

		// preformatting post-topo
		std::vector<fT> post_topo(this->_surface.size(), 0);
		for (int i = 0; i < this->graph->nnodes; ++i) {
			post_topo[i] = this->_surface[i];
		}

		this->graph->_compute_graph(post_topo, only_SD, false);

		// fill water where depressions have been solved
		if (this->depression_management == HYDROGRAPH_LM::FILL) {
			for (int i = 0; i < this->graph->nnodes; ++i) {
				if (this->connector->boundaries.no_data(i) ||
						this->connector->boundaries.forcing_io(i))
					continue;

				if (this->_surface[i] < post_topo[i]) {

					fT ddhhee = post_topo[i] - this->_surface[i];

					if (this->record_filling)
						this->_rec_filling[i] = ddhhee;

					this->_hw[i] += ddhhee;

					this->_surface[i] = post_topo[i];
				}
			}
		}
	}

	// initial check for boundary conditions and eventually applying relevant
	// changes
	bool _initial_check_boundary_pit(int& node,
																	 std::array<int, 8>& receivers,
																	 std::vector<fT>& vmot_hw)
	{
		if (this->connector->boundaries.no_data(node) ||
				this->connector->flow_out_or_pit(node) ||
				(this->hydrostationary == false && this->_hw[node] == 0 &&
				 this->_Qw[node] == 0)) {
			// Checking mass conservations
			this->tot_Qwin_output += this->_Qw[node];

			// I encountered a pit, that happened when LM are not preprocessed
			// Then I fill it slightly, but it does not really work
			if (this->connector->flow_out_model(node) == false &&
					this->connector->boundaries.no_data(node) == false) {
				// int nn = this->connector->get_neighbour_idx(node,receivers);
				// fT hzh = this->_surface[node];
				// for(int j=0; j < nn; ++j)
				// {
				// 	if(this->_surface[receivers[j]] > hzh)
				// 		hzh = this->_surface[receivers[j]];
				// }
				vmot_hw[node] =
					this->_Qw[node] / this->connector->get_area_at_node(node);
				// if(vmot_hw[node]>0)
				// {
				// 	std::cout << this->_Qw[node] << "??::" << vmot_hw[node] *
				// this->dt_hydro(node) << " ||| " << this->dt_hydro(node) <<
				// std::endl;; 	if(vmot_hw[node]> 1e6) 		throw
				// std::runtime_error("flkflsdfkjdsf");
				// }
				// vmot_hw[node] += hzh - this->_surface[node];
			}
			return true;
		}
		return false;
	}

	// offsetting some of the calculations from the run to make it more
	// undertandable adn reusable
	void _compute_slopes_weights_et_al(int& node,
																		 bool& SF,
																		 fT& Smax,
																		 std::array<fT, 8>& slopes,
																		 std::array<fT, 8>& weights,
																		 int& nrecs,
																		 std::array<int, 8>& receivers,
																		 int& recmax,
																		 fT& dx,
																		 fT& dw0max,
																		 fT& topological_number_v2)
	{
		// Calculating the slopes/weights and topological numbers
		if (SF == false) {

			// ROUTINES FOR Multiple FLow Directions
			// fill the arrays of slopes/weights/.. and calculate Smax, dmax, recmas
			// and all
			Smax = this->weights_automator_v2(receivers,
																				weights,
																				slopes,
																				node,
																				nrecs,
																				topological_number_v2,
																				dw0max,
																				recmax,
																				dx);

			// Debugging, I guess deprecated TODO::CHECK
			if (topological_number_v2 == 0) {
				topological_number_v2 = 1.;
			}

			// Check the recording of topological number
			if (this->debugntopo) {
				this->DEBUGNTOPO[node] = topological_number_v2;
			}

		} else {

			// ROUTINES FOR Single FLow Directions
			// -> Slope
			Smax = this->get_Sw(node,
													this->connector->_Sreceivers[node],
													this->connector->Sdistance2receivers[node],
													this->minslope);
			// -> rec
			recmax = this->connector->_Sreceivers[node];
			// -> dy (integrated width)
			dw0max = this->connector->get_travers_dy_from_dx(
				this->connector->Sdistance2receivers[node]);
			// -> dx (flow distance)
			dx = this->connector->Sdistance2receivers[node];

			// boundary case
			if (this->connector->flow_out_model(recmax) &&
					this->boundhw == BOUNDARY_HW::FIXED_SLOPE) {
				dw0max = this->connector->dy;
				Smax = this->bou_fixed_val;
				dx = this->connector->dx;
			}
		}

		// I can only have an hydraulic slope if I have water
		if (this->_hw[node] == 0)
			Smax = 0;

		// Done
	}

	void _compute_morpho(int& node,
											 int& recmax,
											 fT& dx,
											 fT& Smax,
											 std::vector<fT>& vmot)
	{
		if (this->morphomode != MORPHO::NONE &&
				this->connector->boundaries.forcing_io(node) == false) {

			// precaching rates
			fT edot = 0., ddot = 0., eldot_A = 0., eldot_B = 0., dldot_A = 0.,
				 dldot_B = 0.;

			// Lateral dx for lat e/d
			fT dy = this->connector->get_travers_dy_from_dx(dx);

			// And gathering the orthogonal nodes
			std::pair<int, int> orthonodes =
				this->connector->get_orthogonal_nodes(node, recmax);

			// Calculating sheer stress
			fT tau = this->rho(node) * this->_hw[node] * GRAVITY * Smax;

			if (tau > tau_max)
				this->tau_max = tau;

			// if(tau>150)
			// 	tau = 150;

			// Double checking the orthogonal nodes and if needs be to process them
			// (basically if flow outs model or has boundary exceptions)
			int oA = orthonodes.first;
			if (this->connector->boundaries.forcing_io(oA) ||
					this->connector->is_in_bound(oA) == false ||
					this->connector->boundaries.no_data(oA) ||
					this->connector->flow_out_model(oA))
				oA = -1;
			int oB = orthonodes.second;
			if (this->connector->boundaries.forcing_io(oB) ||
					this->connector->is_in_bound(oB) == false ||
					this->connector->boundaries.no_data(oB) ||
					this->connector->flow_out_model(oB))
				oB = -1;

			// Calculating erosion if sheer stress is above critical
			if (tau > this->tau_c(node)) {
				// Basal erosion rates
				edot =
					this->ke(node) * std::pow(tau - this->tau_c(node), this->aexp(node));
				// if recording stuff
				if (this->record_edot)
					this->_rec_edot[node] += edot;
			}

			// Claculating the deposition rates
			ddot = this->_Qs[node] / this->kd(node);

			// Dealing with lateral deposition if lS > 0 and erosion if lS <0
			if (oA >= 0) {
				fT tSwl = this->get_Stopo(node, oA, dy);
				if (tSwl > 0) {
					dldot_A = tSwl * this->kd_lateral(node) * ddot;
				} else {
					eldot_A = std::abs(tSwl) * this->ke_lateral(node) * edot;
				}
			}

			if (oB >= 0) {
				fT tSwl = this->get_Stopo(node, oB, dy);
				if (tSwl > 0) {
					dldot_B = tSwl * this->kd_lateral(node) * ddot;
				} else {
					eldot_B = std::abs(tSwl) * this->ke_lateral(node) * edot;
				}
			}

			// Am I depositing more than I can chew?
			fT totdqs = dx * (dldot_B + dldot_A + ddot);
			if (totdqs > this->_Qs[node]) {
				// std::cout << " happens??? " << totdqs;
				fT corrqs = this->_Qs[node] / totdqs;
				dldot_B *= corrqs;
				dldot_A *= corrqs;
				ddot *= corrqs;
			}
			fT sqs = this->_Qs[node];
			fT fbatch = (ddot + dldot_B + dldot_A - edot - eldot_A - eldot_B) * dx;
			this->_Qs[node] -= fbatch;
			if (std::isfinite(this->_Qs[node]) == false) {
				std::cout << "QS NAN:" << this->_Qs[node] << " vs " << sqs << std::endl;
				throw std::runtime_error("BITE");
			}

			if (this->_Qs[node] < 0) {
				this->_Qs[node] = 0;
			}

			vmot[node] += ddot - edot;
			if (oA >= 0) {
				vmot[oA] += dldot_A;
				vmot[oA] -= eldot_A;
			}
			if (oB >= 0) {
				vmot[oB] += dldot_B;
				vmot[oB] -= eldot_B;
			}

			if (std::isfinite(vmot[node]) == false) {
				std::cout << "edot:" << edot << " ddot" << ddot << std::endl;
				std::cout << "qs:" << sqs << " tau" << tau << std::endl;
				throw std::runtime_error("Non finite vmot gaft");
			}
		}
	}

	void _compute_transfers(int& nrecs,
													int& recmax,
													int& node,
													bool& SF,
													std::array<int, 8>& receivers,
													std::array<fT, 8>& weights,
													fT& Qwin,
													fT& Qwout)
	{
		// going through the receiver(s)
		for (int j = 0; j < nrecs; ++j) {

			// Hydro for the link
			// universal params
			int rec;
			if (SF)
				rec = recmax;
			else
				rec = this->connector->get_to_links(receivers[j]);

			if (rec < 0)
				continue;

			if (this->connector->flow_out_model(rec)) {
				if (SF) {
					this->tot_Qw_output += Qwout;
				} else if (weights[j] > 0 && Qwin > 0) {
					this->tot_Qw_output += weights[j] * Qwout;
				}
			}

			if (this->hydrostationary) {
				if (SF) {
					this->_Qw[rec] += Qwin;
				} else if (weights[j] > 0 && Qwin > 0) {
					this->_Qw[rec] += weights[j] * Qwin;
				}

			} else {
				// this->_Qw[rec] += std::min(Qwout * weights[j], Qwin * weights[j]);
				if (SF)
					this->_Qw[rec] += Qwout;
				else
					this->_Qw[rec] += Qwout * weights[j];
				// std::cout << "transferring:" << Qwout << std::endl;
			}

			if (this->morphomode != MORPHO::NONE) {
				if (std::isfinite(this->_Qs[node]) == false)
					throw std::runtime_error("QS NAN");
				if (std::isfinite(this->_Qs[rec]) == false)
					throw std::runtime_error("QSREC NAN");
				this->_Qs[rec] +=
					(SF == false) ? weights[j] * this->_Qs[node] : this->_Qs[node];
				if (std::isfinite(this->_Qs[rec]) == false) {
					std::cout << weights[j] << std::endl;
					;
					throw std::runtime_error("QSREC NAN AFTER");
				}
			}
		}
	}

	void _compute_vertical_motions(std::vector<fT>& vmot_hw,
																 std::vector<fT>& vmot,
																 bool use_dt = true)
	{

		for (int i = 0; i < this->graph->nnodes; ++i) {

			if (this->connector->flow_out_model(i) &&
					this->boundhw == BOUNDARY_HW::FIXED_HW)
				this->_hw[i] = this->bou_fixed_val;

			if (this->connector->boundaries.forcing_io(i))
				continue;

			fT tvh;
			if (tvh > 0)
				tvh = std::min(vmot_hw[i], this->maxdHw);
			else
				tvh = std::max(vmot_hw[i], this->mindHw);

			if (use_dt) {
				tvh *= this->dt_hydro(i);
				// std::cout << this->dt_hydro(i) << std::endl;
			}
			if (tvh < -this->_hw[i]) {
				tvh = -this->_hw[i];
			}

			this->_hw[i] += tvh;
			if (this->record_dhw)
				this->_rec_dhw[i] = tvh;

			this->_surface[i] += tvh;

			if (this->morphomode != MORPHO::NONE) {

				if (std::isfinite(vmot[i]) == false) {
					throw std::runtime_error("non finite vmot prabul");
				}

				this->_surface[i] += vmot[i] * this->dt_morpho(i);
			}

			if (this->_hw[i] < 0)
				throw std::runtime_error("hw < 0???");
		}
	}

	void _compute_vertical_motions_hw(std::vector<fT>& vmot_hw,
																		bool use_dt = true)
	{

		for (int i = 0; i < this->graph->nnodes; ++i) {

			if (this->connector->flow_out_model(i) &&
					this->boundhw == BOUNDARY_HW::FIXED_HW)
				this->_hw[i] = this->bou_fixed_val;

			if (this->connector->boundaries.forcing_io(i) && this->hydrostationary)
				continue;

			fT tvh;
			if (tvh > 0)
				tvh = std::min(vmot_hw[i], this->maxdHw);
			else
				tvh = std::max(vmot_hw[i], this->mindHw);

			if (use_dt) {
				tvh *= this->dt_hydro(i);
				// std::cout << this->dt_hydro(i) << std::endl;
			}
			if (tvh < -this->_hw[i]) {
				tvh = -this->_hw[i];
			}

			this->_hw[i] += tvh;
			if (this->record_dhw)
				this->_rec_dhw[i] = tvh;

			this->_surface[i] += tvh;
		}
	}

	void _compute_vertical_motions_averaged_test(std::vector<fT>& vmot_hw,
																							 std::vector<fT>& vmot,
																							 bool use_dt = true)
	{

		for (int i = 0; i < this->graph->nnodes; ++i) {

			if (this->connector->flow_out_model(i) &&
					this->boundhw == BOUNDARY_HW::FIXED_HW)
				this->_hw[i] = this->bou_fixed_val;

			if (this->connector->boundaries.forcing_io(i))
				continue;

			fT tvh = vmot_hw[i];
			if (use_dt) {
				tvh *= this->dt_hydro(i);
				// std::cout << this->dt_hydro(i) << std::endl;
			}
			if (tvh < -this->_hw[i]) {
				tvh = -this->_hw[i];
			}

			fT newhw = this->_hw[i] + tvh;
			fT prevhw = this->_hw[i];
			this->_hw[i] = (this->_hw[i] + newhw) / 2;

			if (this->record_dhw)
				this->_rec_dhw[i] = tvh;

			this->_surface[i] += newhw - prevhw;

			if (this->morphomode != MORPHO::NONE)
				this->_surface[i] += vmot[i] * this->dt_morpho(i);

			if (this->_hw[i] < 0)
				throw std::runtime_error("hw < 0???");
		}
	}

	fT weights_automator_v2(std::array<int, 8>& receivers,
													std::array<fT, 8>& weights,
													std::array<fT, 8>& slopes,
													int& node,
													int& nrecs,
													fT& topological_number_v2,
													fT& dw0max,
													int& recmax,
													fT& dx)
	{

		// Placeholders
		fT sumw = 0., Smax = this->minslope, sumSdw = 0.;
		;
		std::pair<fT, fT> fvech;
		int nval = 0;

		for (int i = 0; i < nrecs; ++i) {
			int lix = receivers[i];

			if (this->connector->is_link_valid(lix) == false) {
				continue;
			}
			++nval;

			// int_dw += this->connector->get_dx_from_links_idx(lix);
			fT tdw = this->connector->get_traverse_dx_from_links_idx(lix);
			int rec = this->connector->get_to_links(lix);
			bool isboud = false;
			if (this->connector->flow_out_or_pit(rec) &&
					this->boundhw == BOUNDARY_HW::FIXED_SLOPE) {
				slopes[i] = this->bou_fixed_val;
				tdw = this->connector->dy;
				isboud = true;
			} else
				slopes[i] = this->get_Sw(lix, static_cast<fT>(1e-6));

			sumSdw += slopes[i] * tdw;

			if (this->record_flowvec) {
				this->connector->get_dxdy_from_links_idx(lix, node, fvech, false);
				this->_rec_flowvec[node * 2] += slopes[i] * tdw * fvech.first;
				this->_rec_flowvec[node * 2 + 1] += slopes[i] * tdw * fvech.second;
				// if(this->_rec_flowvec[node * 2] > 0)
				// std::cout << this->_rec_flowvec[node * 2] << std::endl;
			}

			if (slopes[i] > Smax) {
				Smax = slopes[i];
				dw0max = tdw;
				dx = (isboud) ? this->connector->dx
											: this->connector->get_dx_from_links_idx(lix);
				recmax = rec;
			}

			// // slopes[i] = slope;
			// if(this->stochaslope)
			// 	weights[i] = this->randu.get();

			if (this->weight_management == MFD_PARTITIONNING::PROPOSLOPE)
				weights[i] = slopes[i] * tdw;

			else if (this->weight_management == MFD_PARTITIONNING::SQRTSLOPE)
				weights[i] = std::sqrt(slopes[i] * tdw);

			else if (this->weight_management == MFD_PARTITIONNING::PROPOREC)
				weights[i] = 1.;

			else if (this->weight_management == MFD_PARTITIONNING::PROPOSLOPE_NODIAG)
				weights[i] = slopes[i] * tdw;

			if (this->stochaslope)
				weights[i] += this->stochaslope_coeff * this->randu.get();

			sumw += weights[i];
		}

		fT sumf = 0.;
		for (int i = 0; i < nrecs; ++i) {
			int lix = receivers[i];

			if (this->connector->is_link_valid(lix) == false)
				continue;

			if (sumw > 0)
				weights[i] = weights[i] / sumw;
			else
				weights[i] = 1. / nval;

			sumf += weights[i];
		}

		if (std::abs(1 - sumf) > 1e-3)
			std::cout << "|||"
								<< "WARNING::MASSLOSS at node " << node << " -> " << sumf << "/"
								<< std::to_string(sumw == 0) << " ? " << nval << " vs " << nrecs
								<< " |||";

		if (this->record_flowvec) {

			fT length = std::sqrt(std::pow(this->_rec_flowvec[node * 2], 2) +
														std::pow(this->_rec_flowvec[node * 2 + 1], 2));
			if (length != 0) {
				this->_rec_flowvec[node * 2] /= length;
				this->_rec_flowvec[node * 2 + 1] /= length;
			}
		}

		if (this->weight_management == MFD_PARTITIONNING::PROPOSLOPE_NODIAG) {
			dw0max = this->connector->dx;
		}

		topological_number_v2 = (Smax * dw0max) / sumSdw;

		return Smax;
	}

	void init_convergence_checker(int tN, CONVERGENCE conv)
	{

		this->conv_nodes.clear();
		this->conv_ini_Qw.clear();
		this->conv_Qr.clear();
		this->conv_dhw.clear();
		this->converged.clear();

		this->convergence_mode = conv;

		// init convergence vectors to the right size
		this->conv_nodes.reserve(tN);
		this->conv_ini_Qw.reserve(tN);
		this->converged.reserve(tN);

		if (this->convergence_mode == CONVERGENCE::QWR ||
				this->convergence_mode == CONVERGENCE::ALL)
			this->record_Qw_out = true;

		// Calculating relevant metrics
		auto DA = this->calculate_drainage_area();
		auto FD = this->calculate_flow_distance();

		// finding the longest flow line
		auto maxElement = std::max_element(FD.begin(), FD.end());
		int ti = std::distance(FD.begin(), maxElement);

		std::vector<int> flow_line;
		flow_line.reserve(this->connector->nnodes);

		while (this->connector->flow_out_or_pit(ti) == false) {
			flow_line.emplace_back(ti);
			ti = this->connector->_Sreceivers[ti];
		}
		// std::cout << "Last node is " << ti << std::endl;

		if (this->connector->flow_out_model(ti) == false)
			std::cout << "WARNING::convergence checker flow line ends in a pit, this "
									 "can impact its liability"
								<< std::endl;

		int nflodes = static_cast<int>(flow_line.size());
		fT step = static_cast<fT>(nflodes) / tN;
		if (step == 0)
			throw std::runtime_error("cannot init convergence checker: Nnodes bigger "
															 "than the size of the biggest flow line");

		// std::cout << "nflodes is " << nflodes << " step is " << step <<
		// std::endl;

		fT tstep = 0;
		while (conv_nodes.size() < tN &&
					 tstep < static_cast<fT>(flow_line.size())) {
			// std::cout << tstep << " vs " << flow_line.size();
			int node = flow_line[std::floor(tstep)];
			// std::cout << " and node is " << node << std::endl;
			this->conv_nodes.emplace_back(node);
			this->conv_ini_Qw.emplace_back(DA[node]);
			tstep += step;
		}

		std::cout << "Lest tstep is " << tstep << std::endl;

		if (this->n_nodes_convergence() == tN - 1) {
			this->conv_nodes.emplace_back(flow_line.back());
			this->conv_ini_Qw.emplace_back(DA[flow_line.back()]);
		}

		if (this->convergence_mode == CONVERGENCE::ALL ||
				this->convergence_mode == CONVERGENCE::QWR)
			this->conv_Qr = std::vector<std::vector<fT>>(tN);
		if (this->convergence_mode == CONVERGENCE::ALL ||
				this->convergence_mode == CONVERGENCE::DHW)
			this->conv_dhw = std::vector<std::vector<fT>>(tN);
	}

	template<class out_t>
	out_t get_conv_nodes()
	{
		return DAGGER::format_output<std::vector<int>, out_t>(this->conv_nodes);
	}

	template<class out_t>
	out_t get_conv_ini_Qw()
	{
		return DAGGER::format_output<std::vector<fT>, out_t>(this->conv_ini_Qw);
	}

	std::vector<fT> calculate_drainage_area()
	{

		// Graph Processing (SFD will be calculated whatsoever)
		this->graph_automator();
		return this->graph->_accumulate_constant_downstream_SFD(
			this->connector->get_area_at_node(0));
	}

	std::vector<fT> calculate_flow_distance()
	{
		std::vector<fT> FD(this->graph->nnodes, 0.);
		this->graph->_get_SFD_distance_from_outlets(FD);
		return FD;
	}

	template<class out_t>
	out_t get_conv_mean_Qr(int n_n)
	{
		auto temp = this->_get_conv_mean_Qr(n_n);
		return DAGGER::format_output<std::vector<fT>, out_t>(temp);
	}

	std::vector<fT> _get_conv_mean_Qr(int n_n)
	{
		std::vector<fT> out(this->n_nodes_convergence(), 0.);
		n_n = std::min(static_cast<int>(this->n_stack_convergence()), n_n);
		if (n_n == 0)
			return out;

		for (int i = 0; i < this->n_nodes_convergence(); ++i) {
			out[i] = this->__get_conv_mean_Qr(i, n_n);
		}
		return out;
	}

	fT __get_conv_mean_Qr(int i, int n_n)
	{
		fT mean = 0;
		for (int j = static_cast<int>(this->conv_Qr[0].size()) - n_n;
				 j < static_cast<int>(this->conv_Qr[0].size());
				 ++j)
			mean += this->conv_Qr[i][j];
		mean /= n_n;
		return mean;
	}

	template<class out_t>
	out_t get_conv_mean_dhw(int n_n)
	{
		auto temp = this->_get_conv_mean_dhw(n_n);
		return DAGGER::format_output<std::vector<fT>, out_t>(temp);
	}

	std::vector<fT> _get_conv_mean_dhw(int n_n)
	{
		std::vector<fT> out(this->n_nodes_convergence(), 0.);
		n_n = std::min(static_cast<int>(this->n_stack_convergence()), n_n);
		if (n_n == 0)
			return out;
		for (int i = 0; i < this->n_nodes_convergence(); ++i) {
			out[i] = this->__get_conv_mean_dhw(i, n_n);
		}
		return out;
	}

	fT __get_conv_mean_dhw(int i, int n_n)
	{
		fT mean = 0;
		for (int j = static_cast<int>(this->conv_dhw[0].size()) - n_n;
				 j < static_cast<int>(this->conv_dhw[0].size());
				 ++j)
			mean += this->conv_dhw[i][j];
		mean /= n_n;
		return mean;
	}

	void catch_nan(fT testval, std::string error_message)
	{
		if (std::isfinite(testval) == false)
			throw std::runtime_error(error_message);
	}

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~= Topographic Analysis ~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// computes u, q or Q
	// dim: 1 for velocity u, 2 for discarge per unit width for q or 3 for
	// volumetric discharge Q
	template<class out_t>
	out_t compute_tuqQ(int dimensions)
	{
		if (dimensions <= 0 || dimensions > 3)
			throw std::runtime_error(
				"Invalid number of dimensions. Needs to be 0 for the topological "
				"number (dfimensionless coefficient of paritionning), 1 (velocity in "
				"m/s), 2 (discharge per unit width in m^2/s) or 3 (volumetric "
				"discarghe in m^3/s)");

		std::vector<fT> out = this->_compute_tuqQ(dimensions);

		return format_output<decltype(out), out_t>(out);
	}

	std::vector<fT> _compute_tuqQ(int dimensions)
	{

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition
		// and other monitoring features:
		// this->init_Qw();

		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> out(this->graph->nnodes, 0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT, 8> weights, slopes;

		// main loop
		for (int i = this->graph->nnodes - 1; i >= 0; --i) {

			// Getting next node in line
			int node = this->get_istack_node(i);

			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need
			// to be processed THis is where all the boundary treatment happens, if
			// you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers))
				continue;

			// Getting the receivers
			int nrecs;
			if (SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node, receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			// NOTE:
			//  No need to calculate the topological number anymore, but keeping it
			//  for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al(node,
																					SF,
																					Smax,
																					slopes,
																					weights,
																					nrecs,
																					receivers,
																					recmax,
																					dx,
																					dw0max,
																					topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], TWOTHIRD);

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Metric to calculate:
			fT metric = topological_number_v2;
			if (dimensions > 0)
				metric = pohw * sqrtSmax / this->mannings(node);
			if (dimensions > 1)
				metric *= this->_hw[node];
			if (dimensions > 2)
				metric *= dw0max;

			// transfer fluxes
			// Computing the transfer of water and sed
			// this->_compute_transfers(
			// 	nrecs, recmax, node, SF, receivers, weights, Qwin, metric);

			out[node] = metric;

			// Volumetric discahrge
			// fT Qwout = dw0max * this->_hw[node] * u_flow;
		}

		return out;
	}

	template<class out_t>
	out_t compute_AD8(fT exp_slope)
	{
		auto hydrocache = this->hydromode;

		this->hydromode = HYDRO::GRAPH_SFD;

		auto Aout = this->graph->_get_drainage_area_SFD();
		return format_output<decltype(Aout), out_t>(Aout);
	}

	template<class out_t>
	out_t compute_AD8_maxQw()
	{
		auto hydrocache = this->hydromode;

		this->hydromode = HYDRO::GRAPH_SFD;

		// getting the volumetric discharge
		std::vector<fT> tQwout = this->_compute_tuqQ(3);
		std::vector<fT> Aout(this->connector->nxy(), 0.);

		auto receivers = this->connector->get_empty_neighbour();
		for (int i = this->connector->nxy() - 1; i >= 0; --i) {

			int node = this->graph->stack[i];
			Aout[node] += this->connector->get_area_at_node(node);

			if (this->_initial_check_boundary_pit(node, receivers))
				continue;

			int nrecs = this->connector->get_receivers_idx_links(node, receivers);
			fT maxQw = -1;
			int tSrec = this->connector->Sreceivers(node);
			fT tSdx = this->connector->Sdistance2receivers[node];
			for (int j = 0; j < nrecs; ++j) {
				int tl = receivers[j];
				fT dx = this->connector->get_dx_from_links_idx(tl);
				int tn = this->connector->get_to_links(tl);
				if (tQwout[tn] > maxQw) {
					maxQw = tQwout[tn];
					tSrec = tn;
					tSdx = dx;
				}
			}

			this->connector->_Sreceivers[node] = tSrec;
			this->connector->Sdistance2receivers[node] = tSdx;

			Aout[tSrec] += Aout[node];
		}

		this->connector->recompute_SF_donors_from_receivers();
		this->connector->recompute_SF_donors_from_receivers();
		this->graph->topological_sorting_SF();

		this->hydromode = hydrocache;
		return format_output<decltype(Aout), out_t>(Aout);
	}

	template<class out_t>
	out_t compute_AD8_stochastic_Qw(fT exp)
	{

		auto hydrocache = this->hydromode;

		this->hydromode = HYDRO::GRAPH_SFD;

		// getting the volumetric discharge
		std::vector<fT> tQwout = this->_compute_tuqQ(3);
		std::vector<fT> Aout(this->connector->nxy(), 0.);

		auto receivers = this->connector->get_empty_neighbour();
		for (int i = this->connector->nxy() - 1; i >= 0; --i) {

			int node = this->graph->stack[i];
			Aout[node] += this->connector->get_area_at_node(node);

			if (this->_initial_check_boundary_pit(node, receivers))
				continue;

			int nrecs = this->connector->get_receivers_idx_links(node, receivers);
			fT maxQw = -1;
			int tSrec = this->connector->Sreceivers(node);
			fT tSdx = this->connector->Sdistance2receivers[node];
			for (int j = 0; j < nrecs; ++j) {
				int tl = receivers[j];
				fT dx = this->connector->get_dx_from_links_idx(tl);
				int tn = this->connector->get_to_links(tl);

				fT tmaxQw = std::pow(tQwout[tn], exp) * this->connector->randu->get();

				if (tmaxQw > maxQw) {
					maxQw = tmaxQw;
					tSrec = tn;
					tSdx = dx;
				}
			}

			this->connector->_Sreceivers[node] = tSrec;
			this->connector->Sdistance2receivers[node] = tSdx;

			Aout[tSrec] += Aout[node];
		}

		this->connector->recompute_SF_donors_from_receivers();
		this->connector->recompute_SF_donors_from_receivers();
		this->graph->topological_sorting_SF();

		this->hydromode = hydrocache;
		return format_output<decltype(Aout), out_t>(Aout);
	}

	template<class out_t>
	out_t compute_QW8_stochastic_Qw(fT exp)
	{

		auto hydrocache = this->hydromode;

		this->hydromode = HYDRO::GRAPH_SFD;

		// getting the volumetric discharge out
		std::vector<fT> tQwout = this->_compute_tuqQ(3);

		// initialising the ouwput
		std::vector<fT> Aout(this->connector->nxy(), 0.);
		if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT ||
				this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE) {
			for (int i = 0; i < this->graph->nnodes; ++i) {
				if (this->connector->boundaries.can_give(i) &&
						this->connector->flow_out_or_pit(i) == false) {
					Aout[i] +=
						this->precipitations(i) * this->connector->get_area_at_node(i);
				}
			}
		} else if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_H) {
			for (size_t i = 0; i < this->_water_entries.size(); ++i) {
				int node = this->_water_entry_nodes[i];
				Aout[node] +=
					this->_water_entries[i] * this->connector->get_area_at_node(node);
			}
		}

		auto receivers = this->connector->get_empty_neighbour();
		for (int i = this->connector->nxy() - 1; i >= 0; --i) {

			int node = this->graph->stack[i];

			if (this->_initial_check_boundary_pit(node, receivers))
				continue;

			int nrecs = this->connector->get_receivers_idx_links(node, receivers);
			fT maxQw = -1;
			int tSrec = this->connector->Sreceivers(node);
			fT tSdx = this->connector->Sdistance2receivers[node];
			for (int j = 0; j < nrecs; ++j) {
				int tl = receivers[j];
				fT dx = this->connector->get_dx_from_links_idx(tl);
				int tn = this->connector->get_to_links(tl);

				fT tmaxQw = std::pow(tQwout[tn], exp) * this->connector->randu->get();

				if (tmaxQw > maxQw) {
					maxQw = tmaxQw;
					tSrec = tn;
					tSdx = dx;
				}
			}

			this->connector->_Sreceivers[node] = tSrec;
			this->connector->Sdistance2receivers[node] = tSdx;

			Aout[tSrec] += Aout[node];
		}

		this->connector->recompute_SF_donors_from_receivers();
		this->graph->topological_sorting_SF();

		this->hydromode = hydrocache;
		return format_output<decltype(Aout), out_t>(Aout);
	}

	template<class out_t>
	out_t compute_AD8_stochastic_Sw(fT exp)
	{

		auto hydrocache = this->hydromode;

		this->hydromode = HYDRO::GRAPH_SFD;

		// getting the volumetric discharge
		std::vector<fT> Aout(this->connector->nxy(), 0.);

		auto receivers = this->connector->get_empty_neighbour();
		for (int i = this->connector->nxy() - 1; i >= 0; --i) {

			int node = this->graph->stack[i];
			Aout[node] += this->connector->get_area_at_node(node);

			if (this->_initial_check_boundary_pit(node, receivers))
				continue;

			int nrecs = this->connector->get_receivers_idx_links(node, receivers);
			fT maxQw = -1;
			int tSrec = this->connector->Sreceivers(node);
			fT tSdx = this->connector->Sdistance2receivers[node];
			for (int j = 0; j < nrecs; ++j) {
				int tl = receivers[j];
				fT dx = this->connector->get_dx_from_links_idx(tl);
				int tn = this->connector->get_to_links(tl);

				fT testSw =
					std::pow(this->get_Sw(tl, 1e-6), exp) * this->connector->randu->get();

				if (testSw > maxQw) {
					maxQw = testSw;
					tSrec = tn;
					tSdx = dx;
				}
			}

			this->connector->_Sreceivers[node] = tSrec;
			this->connector->Sdistance2receivers[node] = tSdx;

			Aout[tSrec] += Aout[node];
		}

		this->connector->recompute_SF_donors_from_receivers();
		this->graph->topological_sorting_SF();

		this->hydromode = hydrocache;
		return format_output<decltype(Aout), out_t>(Aout);
	}

	template<class in_t, class out_t>
	out_t compute_elemental_transfer(in_t& _production,
																	 fT total_time,
																	 fT concentration_max)
	{

		auto production = format_input(_production);

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition
		// and other monitoring features:
		this->init_Qw();

		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// Vertical motions are applied at the end of the timestep
		std::vector<fT> out(this->graph->nnodes, 0.),
			remaining_t(this->graph->nnodes, 0.);

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT, 8> weights, slopes;

		// main loop
		for (int i = this->graph->nnodes - 1; i >= 0; --i) {

			// Getting next node in line
			int node = this->get_istack_node(i);

			production[node] =
				std::min(production[node], this->_Qw[node] * concentration_max);

			if (production[node] > 0) {
				remaining_t[node] =
					(remaining_t[node] * out[node] + production[node] * total_time) /
					(out[node] + production[node]);
				out[node] += production[node] / this->_hw[node];
			}

			out[node] = std::min(out[node], this->_Qw[node] * concentration_max);

			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need
			// to be processed THis is where all the boundary treatment happens, if
			// you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers))
				continue;

			// Getting the receivers
			int nrecs;
			if (SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node, receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			// NOTE:
			//  No need to calculate the topological number anymore, but keeping it
			//  for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al(node,
																					SF,
																					Smax,
																					slopes,
																					weights,
																					nrecs,
																					receivers,
																					recmax,
																					dx,
																					dw0max,
																					topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], TWOTHIRD);

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Metric to calculate:
			fT transfer_time = dx / (pohw * sqrtSmax / this->mannings(node));
			remaining_t[node] -= transfer_time;

			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(
				nrecs, recmax, node, SF, receivers, weights, Qwin, transfer_time);

			if (out[node] > 0 && remaining_t[node] > 0) {
				for (int j = 0; j < nrecs; ++j) {
					int trec = this->connector->get_to_links(receivers[j]);
					fT tadd = out[node] * weights[j];

					remaining_t[trec] =
						(remaining_t[trec] * out[trec] + remaining_t[node] * tadd) /
						(tadd + out[trec]);

					out[trec] = (tadd + out[trec]);

					out[trec] = std::min(out[trec], this->_Qw[trec] * concentration_max);
				}
			}

			// out[node] = metric;

			// Volumetric discahrge
			// fT Qwout = dw0max * this->_hw[node] * u_flow;
		}

		for (int i = 0; i < this->connector->nnodes; ++i)
			out[i] *= remaining_t[i];

		return format_output<decltype(out), out_t>(out);
	}

	// Experimental drainage area stuff
	// template<class out_t

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~= Precipitests ~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	int spawn_precipition()
	{

		int node;

		do {
			bool OK = false;

			if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT ||
					this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE) {
				node = this->dis(this->gen);
				if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT)
					OK = true;
				else if (this->precipitations(node) > 0)
					OK = true;
			} else if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_H) {
				node = this->_water_entry_nodes[this->dis(this->gen)];
				OK = true;
			} else
				throw std::runtime_error("NOT DEFINED");

			if (this->connector->boundaries.no_data(node) == false && OK)
				break;

		} while (true);
		// std::cout << node << std::endl;

		return node;
	}

	// Run a basic version of floodos, for testing purposes
	// - N is the number of precipitons to launch
	// - precdt is the time step between 2 successive precipitons
	void run_precipitions(int N, fT precdt)
	{
		// check if this is the first lauch and allocate the vectors if needed
		if (this->last_dt_prec.size() == 0) {
			this->last_dt_prec = std::vector<fT>(this->connector->nnodes, 0.);
			this->last_dt_prec_e = std::vector<fT>(this->connector->nnodes, 0.);
			this->last_sw_prec = std::vector<fT>(this->connector->nnodes, 0.);
			this->last_dx_prec = std::vector<fT>(this->connector->nnodes, 0.);
		}

		// Computes the total volume of Qin.t for all entry points
		// -> if precipitation: sums all the precipitations * cellarea * dt
		// -> if entry Qwin: sum all the entryQwin*dt
		// this also sets the stochastic spawner
		this->compute_Vp(precdt);

		// current_dt_prec
		// caching neighbours and slopes
		auto neighbours = this->connector->get_empty_neighbour();

		// launching n precipitons
		for (int _ = 0; _ < N; ++_) {

			// Spawn location
			int node = this->spawn_precipition();
			// caching the starting point for debugging
			// int onode = node;

			// increment the global chronometer
			this->current_dt_prec += precdt;

			// and for erosion if activated
			if (this->morphomode != MORPHO::NONE)
				this->current_dt_prec_e += precdt;

			// sediment flux
			fT tQs = this->Vps * this->connector->dy;

			// adding a safeguard stopping the precipiton if it runs for too long
			// (mostly used to avoid infinite loop when divergence occurs)
			int NMAX = this->connector->nnodes * 2;
			while (true) {
				--NMAX;

				// stops the precipitons if it has reached model edge
				if (this->connector->boundaries.has_to_out(node) || NMAX == 0)
					break;

				// local time step
				fT tdt = this->current_dt_prec - this->last_dt_prec[node];
				// this->last_dt_prec[node] = this->current_dt_prec;

				// caching receivers/slope/...
				int rec = node; // will be rec
				fT Sw = 0.;			// will be hydraulic slope
				fT Sw_sel = 0.; // stochastic slope to select next receiver
				// spatial steps
				fT dx = 0.;
				fT dy = 0.;

				// getting the current neighbours
				int nn = this->connector->get_neighbour_idx_links(node, neighbours);
				for (int j = 0; j < nn; ++j) {
					// local neighbour
					int trec =
						this->connector->get_other_node_from_links(neighbours[j], node);
					// ignore no data
					if (this->connector->boundaries.no_data(trec))
						continue;

					// if(this->last_dt_prec[trec] != this->current_dt_prec)
					// decrease the water height of neighbour based on previous sw
					this->update_hw_prec_and_dt(trec);
				}

				// decrease node's hw based on previous sw
				this->update_hw_prec_and_dt(node);

				// Takes care of adding the precipiton's hw (making sure water is not
				// negative after the decrease too)
				fT deltah = 0.;
				if (this->connector->boundaries.can_out(node) == false) {
					deltah = this->Vp;
					this->_hw[node] += deltah;
					if (this->_hw[node] < 0.) {
						deltah += this->_hw[node];
						this->_hw[node] = 0.;
					}
				}

				// propagating to the surface
				this->_surface[node] += deltah;

				// Determining the next path and updating the slopw/time/...
				for (int j = 0; j < nn; ++j) {

					// testing next receiver
					int trec =
						this->connector->get_other_node_from_links(neighbours[j], node);

					if (this->connector->boundaries.no_data(trec)) {
						continue;
					}

					// if potential receiver
					if (this->_surface[trec] <= this->_surface[node]) {

						// calculating slope
						fT tsw = (this->_surface[node] - this->_surface[trec]) /
										 this->connector->get_dx_from_links_idx(neighbours[j]);
						// stochasticity f(sqrt(tsw))
						fT tsw_sel =
							std::sqrt(tsw) * this->randu.get() *
							this->connector->get_traverse_dx_from_links_idx(neighbours[j]);

						// Selecting these info if selected
						if (tsw_sel > Sw_sel) {
							Sw_sel = tsw_sel;
							rec = trec;
							dx = this->connector->get_dx_from_links_idx(neighbours[j]);
							dy = this->connector->get_travers_dy_from_dx(dx);
						}

						// finding Smax
						if (tsw > Sw)
							Sw = tsw;

						// managing preboundary conditions
						if (this->connector->boundaries.can_out(trec) &&
								this->boundhw == BOUNDARY_HW::FIXED_SLOPE) {
							Sw_sel = this->bou_fixed_val;
							Sw = this->bou_fixed_val;
							rec = trec;
							dx = this->connector->dx;
							dy = this->connector->dy;
							break;
						}
					}
				}

				// I have the Smax and other thingies, saving them
				this->last_sw_prec[node] = Sw;
				this->last_dx_prec[node] = dx;

				// Dealing with morpho if enabled
				if (this->morphomode != MORPHO::NONE &&
						this->connector->boundaries.forcing_io(node) == false) {
					if (Sw > 0 && rec != node && dx > 0) {
						tdt = (this->current_dt_prec_e - this->last_dt_prec_e[node]) *
									this->dt_morpho_multiplier;
						this->last_dt_prec_e[node] = this->current_dt_prec_e;
						fT tau = this->rho(node) * this->_hw[node] * GRAVITY * Sw;
						fT edot = 0., ddot = 0., edot_A = 0., edot_B = 0, ddot_A = 0,
							 ddot_B = 0.;
						if (tau > this->tau_c(node))
							edot += this->ke(node) *
											std::pow(tau - this->tau_c(node), this->aexp(node));

						ddot = tQs / this->kd(node);

						if (tQs < 0) {
							// std::cout << "happens";
							tQs = 0;
						}

						std::pair<int, int> orthonodes =
							this->connector->get_orthogonal_nodes(node, rec);
						int oA = orthonodes.first;
						if (this->connector->boundaries.forcing_io(oA) ||
								this->connector->is_in_bound(oA) == false ||
								this->connector->boundaries.no_data(oA) ||
								this->connector->flow_out_model(oA))
							oA = -1;
						int oB = orthonodes.second;
						if (this->connector->boundaries.forcing_io(oB) ||
								this->connector->is_in_bound(oB) == false ||
								this->connector->boundaries.no_data(oB) ||
								this->connector->flow_out_model(oB))
							oB = -1;

						// Dealing with lateral deposition if lS > 0 and erosion if lS <0
						if (oA >= 0) {
							fT tSwl = this->get_Stopo(node, oA, dy);
							if (tSwl > 0) {
								ddot_A = tSwl * this->kd_lateral(node) * ddot;
							} else {
								edot_A = std::abs(tSwl) * this->ke_lateral(node) * edot;
							}

							// std::cout << edot_A << "|" << ddot_A  << std::endl;
						}

						if (oB >= 0) {
							fT tSwl = this->get_Stopo(node, oB, dy);
							if (tSwl > 0) {
								ddot_B = tSwl * this->kd_lateral(node) * ddot;
							} else {
								edot_B = std::abs(tSwl) * this->ke_lateral(node) * edot;
							}
						}

						fT totd = ddot + ddot_B + ddot_A;
						totd *= dx;
						if (totd > tQs) {
							auto cor = tQs / totd;
							ddot *= cor;
							ddot_B *= cor;
							ddot_A *= cor;
						}

						tQs += (edot - ddot) * dx;
						tQs += (edot_B - ddot_B) * dy;
						tQs += (edot_A - ddot_A) * dy;
						if (oA >= 0)
							this->_surface[oA] += (ddot_A - edot_A) * tdt;
						if (oB >= 0)
							this->_surface[oB] += (ddot_B - edot_B) * tdt;

						this->_surface[node] += (ddot - edot) * tdt;
					}
					// else
					// {
					// 	this->_surface[node] += this->Vp;
					// }
				}

				// stopping the loop if node is leaving the model
				if (node == rec && this->connector->boundaries.can_out(node))
					break;

				// updating positions
				node = rec;

			} // end of while loop for a single precipiton

		} // end of the for loop for all the launches

	} // end of the precipiton function

	std::vector<fT> get_precipitations_vector()
	{
		std::vector<fT> out(this->connector->nnodes, 0.);
		for (int i = 0; i < this->connector->nnodes; ++i) {
			out[i] = this->precipitations(i);
		}
		return out;
	}

	void define_precipitations_Ath(fT Ath)
	{
		std::vector<BC> nBCs(this->connector->boundaries.codes);
		std::vector<fT> fake(this->_surface);
		this->graph->_compute_graph(fake, false, false);
		auto DA = this->graph->_accumulate_constant_downstream_SFD(
			this->connector->get_area_at_node(0));
		std::vector<std::uint8_t> vis(this->connector->nnodes, false);
		std::vector<fT> tprecs = this->get_precipitations_vector();

		for (int i = this->connector->nnodes - 1; i >= 0; --i) {

			int node = this->graph->Sstack[i];

			if (this->connector->boundaries.no_data(node))
				continue;

			// if(vis[node]) continue;

			int rec = this->connector->_Sreceivers[node];
			if (Ath < DA[node]) {
				vis[node] = true;
			} else {
				tprecs[rec] += tprecs[node];
			}

			vis[rec] = vis[node];
		}

		std::cout << "LKDFJDLKF" << std::endl;

		auto receivers = this->connector->get_empty_neighbour();
		for (int i = this->connector->nnodes - 1; i >= 0; --i) {
			int node = this->graph->stack[i];

			if (vis[node]) {
				int nn = this->connector->get_receivers_idx(node, receivers);
				for (int j = 0; j < nn; ++j)
					vis[receivers[j]] = true;
			}
		}

		for (int i = this->connector->nnodes - 1; i >= 0; --i) {
			if (vis[i] == false)
				nBCs[i] = BC::NO_FLOW;
		}

		std::cout << "??" << std::endl;
		this->connector->set_custom_boundaries(nBCs);
		this->set_water_input_by_variable_precipitation_rate(tprecs);
		std::cout << "!!" << std::endl;
	}

	void run_graphipiton(int N, fT precdt, fT Ath)
	{
		if (this->last_dt_prec.size() == 0) {
			this->last_dt_prec = std::vector<fT>(this->connector->nnodes, 0.);
			this->last_dt_prec_e = std::vector<fT>(this->connector->nnodes, 0.);
			this->last_sw_prec = std::vector<fT>(this->connector->nnodes, 0.);
			this->last_dx_prec = std::vector<fT>(this->connector->nnodes, 0.);
		}

		std::vector<fT> fake(this->_surface);
		this->graph->_compute_graph(fake, true, false);
		auto DA = this->graph->_accumulate_constant_downstream_SFD(
			this->connector->get_area_at_node(0));
		std::vector<std::uint8_t> vis(this->connector->nnodes, false);
		std::vector<int> starters;
		starters.reserve(1000);
		for (int i = this->connector->nnodes - 1; i >= 0; --i) {
			int node = this->graph->Sstack[i];
			if (vis[node])
				continue;
			int rec = this->connector->_Sreceivers[node];
			vis[node] = true;
			if (Ath < DA[node]) {
				vis[rec] = true;
				starters.emplace_back(node);
			}
		}

		this->compute_Vp(precdt);

		this->dis = std::uniform_int_distribution<>(0, starters.size() - 1);

		// current_dt_prec
		auto neighbours = this->connector->get_empty_neighbour();

		for (int _ = 0; _ < N; ++_) {

			// int node = this->spawn_precipition();
			int node = starters[this->dis(this->gen)];
			// int onode = node;
			// std::cout << "init at " << node << std::endl;

			this->current_dt_prec += precdt;

			if (this->morphomode != MORPHO::NONE)
				this->current_dt_prec_e += precdt;

			// fT tQs = this->Vps * this->connector->dy;

			int NMAX = this->connector->nnodes * 2;
			while (true) {
				--NMAX;
				if (this->connector->boundaries.has_to_out(node) || NMAX == 0)
					break;

				// // update step from Qwout
				// fT tdt = this->current_dt_prec - this->last_dt_prec[node];
				// // this->last_dt_prec[node] = this->current_dt_prec;

				int rec = node;
				fT Sw = 0.;
				fT Sw_sel = 0.;
				fT dx = 0.;
				// fT dy = 0.;
				int nn = this->connector->get_neighbour_idx_links(node, neighbours);

				for (int j = 0; j < nn; ++j) {
					int trec =
						this->connector->get_other_node_from_links(neighbours[j], node);
					if (this->connector->boundaries.no_data(trec))
						continue;

					// if(this->last_dt_prec[trec] != this->current_dt_prec)
					this->update_hw_prec_and_dt(trec);
				}
				this->update_hw_prec_and_dt(node);

				fT deltah = 0.;
				if (true) {
					deltah = this->Vp;
					this->_hw[node] += deltah;
					if (this->_hw[node] < 0.) {
						deltah += this->_hw[node];
						this->_hw[node] = 0.;
					}
				}
				this->_surface[node] += deltah;

				for (int j = 0; j < nn; ++j) {
					int trec =
						this->connector->get_other_node_from_links(neighbours[j], node);

					if (this->connector->boundaries.no_data(trec)) {
						continue;
					}

					if (this->_surface[trec] <= this->_surface[node]) {

						fT tsw = (this->_surface[node] - this->_surface[trec]) /
										 this->connector->get_dx_from_links_idx(neighbours[j]);
						fT tsw_sel = tsw * this->randu.get();
						// std::cout << tsw_sel << std::endl;
						if (tsw_sel > Sw_sel) {
							Sw_sel = tsw_sel;
							rec = trec;
							dx = this->connector->get_dx_from_links_idx(neighbours[j]);
							// dy = this->connector->get_travers_dy_from_dx(dx);
						}

						if (tsw > Sw)
							Sw = tsw;

						if (this->connector->boundaries.can_out(trec) &&
								this->boundhw == BOUNDARY_HW::FIXED_SLOPE) {
							Sw_sel = this->bou_fixed_val;
							Sw = this->bou_fixed_val;
							rec = trec;
							dx = this->connector->dx;
							// dy = this->connector->dy;
						}
					}
				}

				this->last_sw_prec[node] = Sw;
				this->last_dx_prec[node] = dx;

				// std::cout << onode << "||" << node << " vs " << rec << "||";

				if (node == rec && this->connector->boundaries.can_out(node))
					break;

				// if(node != rec)
				// 	throw std::runtime_error("DONE!");
				node = rec;
			}
			// throw std::runtime_error("BITE2");
		}
	}

	void update_hw_prec_and_dt(int node)
	{

		fT tdt = this->current_dt_prec - this->last_dt_prec[node];
		this->last_dt_prec[node] = this->current_dt_prec;

		if (tdt == 0 || this->last_sw_prec[node] == 0 ||
				this->last_dx_prec[node] == 0) {
			return;
		}

		auto delta = tdt * std::pow(this->_hw[node], 5. / 3.) *
								 std::sqrt(this->last_sw_prec[node]) / this->mannings(node) /
								 this->last_dx_prec[node];
		this->_hw[node] -= delta;
		this->_surface[node] -= delta;
	}

	void update_hw_prec_and_dt_exp_2(int node,
																	 std::vector<fT>& gHwin,
																	 fT tHwin,
																	 fT precdt)
	{

		fT tdt = this->current_dt_prec - this->last_dt_prec[node];
		this->last_dt_prec[node] = this->current_dt_prec;

		// SHAKING THINGS HERE
		tdt = precdt;

		if (tdt == 0 || this->last_sw_prec[node] == 0 ||
				this->last_dx_prec[node] == 0) {
			// fT dHwin =  std::max(tHwin, gHwin[node]);

			// this->_hw[node] += dtfill * dHwin ;
			// this->_surface[node] += dtfill * dHwin ;
			return;
		}

		fT dHwin = std::max(tHwin, gHwin[node]);
		// std::cout << "+:" << dHwin * tdt;

		auto delta = tdt * std::pow(this->_hw[node], 5. / 3.) *
								 std::sqrt(this->last_sw_prec[node]) / this->mannings(node) /
								 this->last_dx_prec[node];
		// std::cout << " vs - :" << delta << std::endl;;
		delta -= dHwin * tdt;
		if (delta < -this->_hw[node])
			delta = -this->_hw[node];

		this->_hw[node] -= delta;
		this->_surface[node] -= delta;
		// if(delta > 0)
		// 	std::cout <<"delta::" << delta << std::endl;
	}

	void update_hw_prec_and_dt_exp(int node,
																 std::vector<int>& baslab,
																 std::vector<fT>& basdt)
	{

		fT tdt = basdt[baslab[node]] - this->last_dt_prec[node];
		this->last_dt_prec[node] = basdt[baslab[node]];

		if (tdt == 0 || this->last_sw_prec[node] == 0 ||
				this->last_dx_prec[node] == 0) {
			return;
		}

		auto delta = tdt * std::pow(this->_hw[node], 5. / 3.) *
								 std::sqrt(this->last_sw_prec[node]) / this->mannings(node) /
								 this->last_dx_prec[node];
		this->_hw[node] -= delta;
		this->_surface[node] -= delta;
	}

	void compute_Vp(fT precdt)
	{
		this->Vp = 0.;
		if (this->water_input_mode == WATER_INPUT::PRECIPITATIONS_CONSTANT ||
				this->water_input_mode == WATER_INPUT::PRECIPITATIONS_VARIABLE) {

			for (int i = 0; i < this->connector->nnodes; ++i) {
				this->Vp += this->precipitations(i) * precdt;
			}
			this->dis =
				std::uniform_int_distribution<>(0, this->connector->nnodes - 1);
		} else if (this->water_input_mode == WATER_INPUT::ENTRY_POINTS_H) {
			for (size_t i = 0; i < this->_water_entries.size(); ++i) {
				this->Vp += this->_water_entries[i] * precdt;
			}
			this->dis =
				std::uniform_int_distribution<>(0, this->_water_entries.size() - 1);
		}

		this->Vps = 0;
		if (this->sed_input_mode == SED_INPUT::ENTRY_POINTS_Q) {
			for (size_t i = 0; i < this->_sed_entries.size(); ++i) {
				this->Vps += this->_sed_entries[i] * this->connector->dy * precdt *
										 this->dt_morpho_multiplier;
			}
			// this->dis = std::uniform_int_distribution<>
			// (0,this->_water_entries.size()-1);
		}
	}

	void run_precipitions_exp(int N, fT precdt, fT dtfill)
	{
		// place holder for debugging
		return;
	}

	void run_precipitions_exp_1(int N, fT precdt)
	{
		// place holder for debugging
		return;
	}

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~= Legacy archive ~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	// Note that A lot of experimental functions have been erased (see former
	// commit)
	void run_exp()
	{
		std::cout << "Run_exp is not available at the moment" << std::endl;
	}
	// initial check for boundary conditions WITHOUT applying relevant changes
	bool _initial_check_boundary_pit(int& node, std::array<int, 8>& receivers)
	{
		if (this->connector->boundaries.no_data(node) ||
				this->connector->flow_out_or_pit(node)) {
			// Checking mass conservations
			this->tot_Qwin_output += this->_Qw[node];
			return true;
		}
		return false;
	}

	fT weights_automator(std::array<int, 8>& receivers,
											 std::vector<fT>& weights,
											 std::vector<fT>& slopes,
											 int& node,
											 int& nrecs,
											 fT& topological_number_v2)
	{

		fT sumw = 0., Smax = this->minslope, dw0max = 0, sumSdw = 0.;
		;

		// int_dw = 0;
		// int yolo1 = node;
		// bool isit = yolo1 == 300;
		for (int i = 0; i < nrecs; ++i) {
			int lix = receivers[i];

			if (this->connector->is_link_valid(lix) == false) {
				continue;
			}

			// int_dw += this->connector->get_dx_from_links_idx(lix);
			fT tdw = this->connector->get_traverse_dx_from_links_idx(lix);
			int rec = this->connector->get_to_links(lix);
			if (this->connector->flow_out_or_pit(rec) &&
					this->boundhw == BOUNDARY_HW::FIXED_SLOPE) {
				slopes[i] = this->bou_fixed_val;
				tdw = this->connector->dy;
			} else
				slopes[i] = this->get_Sw(lix, this->minslope);

			sumSdw += slopes[i] * tdw;

			if (slopes[i] > Smax) {
				Smax = slopes[i];
				dw0max = tdw;
			}

			// slopes[i] = slope;
			if (this->weight_management == MFD_PARTITIONNING::PROPOSLOPE)
				weights[i] = slopes[i] * tdw;

			else if (this->weight_management == MFD_PARTITIONNING::SQRTSLOPE)
				weights[i] = std::sqrt(slopes[i] * tdw);

			else if (this->weight_management == MFD_PARTITIONNING::PROPOREC)
				weights[i] = 1.;

			if (this->stochaslope)
				weights[i] *= this->stochaslope_coeff * this->randu.get();

			sumw += weights[i];
		}

		fT sumf = 0.;
		for (int i = 0; i < nrecs; ++i) {
			int lix = receivers[i];

			if (this->connector->is_link_valid(lix) == false)
				continue;

			weights[i] = weights[i] / sumw;
			sumf += weights[i];
		}

		topological_number_v2 = (Smax * dw0max) / sumSdw;

		return Smax;
	}

	void _compute_transfers_exp(int& nrecs,
															int& recmax,
															int& node,
															bool& SF,
															std::array<int, 8>& receivers,
															std::array<fT, 8>& weights,
															fT& Qwin)
	{
		// going through the receiver(s)
		for (int j = 0; j < nrecs; ++j) {

			// Hydro for the link
			// universal params
			int rec;
			if (SF)
				rec = recmax;
			else
				rec = this->connector->get_to_links(receivers[j]);

			if (rec < 0)
				continue;

			// if(this->hydrostationary)
			// {
			if (SF) {
				this->_Qw[rec] += Qwin;
			} else if (weights[j] > 0 && Qwin > 0) {
				this->_Qw[rec] += weights[j] * Qwin;
			}

			// if(this->morphomode != MORPHO::NONE)
			// {
			// 	if(std::isfinite(this->_Qs[node]) == false)
			// 		throw std::runtime_error("QS NAN");
			// 	if(std::isfinite(this->_Qs[rec]) == false)
			// 		throw std::runtime_error("QSREC NAN");
			// 	this->_Qs[rec] += (SF == false)?weights[j] *
			// this->_Qs[node]:this->_Qs[node] ; if(std::isfinite(this->_Qs[rec]) ==
			// false)
			// 	{
			// 		std::cout << weights[j] << std::endl;;
			// 		throw std::runtime_error("QSREC NAN AFTER");
			// 	}
			// }
		}
	}

	// EXPERIMENTAL STUFF
	// # MAin running function
	void run_hydro_only()
	{
		// Saving the topological number if needed
		if (this->debugntopo)
			this->DEBUGNTOPO = std::vector<fT>(this->connector->nnodes, 0);

		this->debug_CFL = 0.;

		this->tau_max = 0.;

		if (this->courant_dt_hydro <= 0)
			this->courant_dt_hydro = 1e-3;

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition
		// and other monitoring features:
		this->init_Qw();

		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// To be used if courant dt hydro is selected
		fT tcourant_dt_hydro = std::numeric_limits<fT>::max();

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT, 8> weights, slopes;

		// main loop
		for (int i = this->graph->nnodes - 1; i >= 0; --i) {

			// Getting next node in line
			int node = this->get_istack_node(i);

			// Processing case where the node is a model edge, or no data
			// this function  returns true if the node was boundary and does not need
			// to be processed THis is where all the boundary treatment happens, if
			// you need to add something happening at the boundaries
			if (this->_initial_check_boundary_pit(node, receivers))
				continue;

			// Getting the receivers
			int nrecs;
			if (SF)
				nrecs = 1;
			else
				nrecs = this->connector->get_receivers_idx_links(node, receivers);

			// Caching Slope max
			fT Smax;
			fT dw0max;
			fT dx;
			int recmax = node;

			// NOTE:
			//  No need to calculate the topological number anymore, but keeping it
			//  for recording its value
			fT topological_number_v2 = 0.;

			this->_compute_slopes_weights_et_al(node,
																					SF,
																					Smax,
																					slopes,
																					weights,
																					nrecs,
																					receivers,
																					recmax,
																					dx,
																					dw0max,
																					topological_number_v2);

			// Initialising the total Qout
			fT Qwin = this->_Qw[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], TWOTHIRD);
			// std:: cout << "pohw:" << pohw << "|" << this->_hw[node] << std::endl;

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax / this->mannings(node);
			// Volumetric discahrge
			fT Qwout = dw0max * this->_hw[node] * u_flow;
			// std::cout << Qwout << "|";

			// Eventually recording Smax
			if (this->record_Sw)
				this->_rec_Sw[node] = Smax;

			// transfer fluxes
			// Computing the transfer of water and sed
			this->_compute_transfers(
				nrecs, recmax, node, SF, receivers, weights, Qwin, Qwout);

			// Computing courant based dt
			if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && u_flow > 0) {
				fT provisional_dt = (this->courant_number * dx) / u_flow;
				provisional_dt = std::min(provisional_dt, this->max_courant_dt_hydro);
				provisional_dt = std::max(provisional_dt, this->min_courant_dt_hydro);
				tcourant_dt_hydro = std::min(provisional_dt, tcourant_dt_hydro);
			}

			// computing hydro vertical motion changes for next time step
			fT dH = this->dt_hydro(node) * (this->_Qw[node] - Qwout) /
							this->connector->get_area_at_node(node);

			this->_hw[node] += dH;
			if (this->_hw[node] < 0) {
				dH -= this->_hw[node];
				this->_hw[node] = 0;
			}
			this->_surface[node] += dH;

			if (this->record_Qw_out)
				this->_rec_Qwout[node] += Qwout;
		}

		// if(this->tau_max > 500)
		//  std::cout << "WARNING::tau_max is " << this->tau_max << std::endl;

		// END OF MAIN LOOP

		// Computing final courant based dt

		if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT) {
			if (this->courant_dt_hydro == -1)
				this->courant_dt_hydro = 1e-3;
			else if (tcourant_dt_hydro > 0 &&
							 tcourant_dt_hydro != std::numeric_limits<fT>::max())
				this->courant_dt_hydro = tcourant_dt_hydro;
		}

		// if (this->convergence_mode == CONVERGENCE::ALL ||
		//     this->convergence_mode == CONVERGENCE::DHW) {
		//   for (int i = 0; i < this->n_nodes_convergence(); ++i)
		//     this->conv_dhw[i].emplace_back(
		//         vmot_hw[this->conv_nodes
		//                     [i]]); // /this->dt_hydro(this->conv_nodes[i]));
		// }

		if (this->convergence_mode == CONVERGENCE::ALL ||
				this->convergence_mode == CONVERGENCE::QWR) {
			for (int i = 0; i < this->n_nodes_convergence(); ++i)
				this->conv_Qr[i].emplace_back(this->_rec_Qwout[this->conv_nodes[i]] /
																			this->_Qw[this->conv_nodes[i]]);
		}
	}

	// EXPERIMENTAL STUFF
	// # MAin running function
	void run_experimental()
	{
		this->hydromode = HYDRO::GRAPH_SFD;
		// Saving the topological number if needed
		if (this->debugntopo)
			this->DEBUGNTOPO = std::vector<fT>(this->connector->nnodes, 0);

		this->debug_CFL = 0.;

		this->tau_max = 0.;

		if (this->courant_dt_hydro <= 0)
			this->courant_dt_hydro = 1e-3;

		// Graph Processing
		this->graph_automator();

		// Initialise the water discharge fields according to water input condition
		// and other monitoring features:
		this->init_Qw();

		std::vector<fT> tQwin =
			this->graph->_accumulate_constant_downstream_SFD(this->precipitations(0));
		std::vector<fT> tQwout(this->connector->nnodes, 0.);
		;

// main loop
#ifdef _OPENMP
#pragma omp parallel for num_threads(4)
#endif
		for (int node = 0; node < this->connector->nnodes; ++node) {
			// Caching Slope max
			fT Smax = this->connector->SS[node];
			fT dx = this->connector->Sdistance2receivers[node];
			fT dw0max = this->connector->get_travers_dy_from_dx(dx);
			int recmax = this->connector->_Sreceivers[node];

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], TWOTHIRD);
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax / this->mannings(node);
			// Volumetric discahrge
			tQwout[node] = dw0max * this->_hw[node] * u_flow;

			// // Computing courant based dt
			// if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && u_flow > 0) {
			//   fT provisional_dt = (this->courant_number * dx) / u_flow;
			//   provisional_dt = std::min(provisional_dt,
			//   this->max_courant_dt_hydro); provisional_dt =
			//   std::max(provisional_dt, this->min_courant_dt_hydro);
			//   tcourant_dt_hydro = std::min(provisional_dt, tcourant_dt_hydro);
			// }
		}

#ifdef _OPENMP
#pragma omp parallel for num_threads(4)
#endif
		for (int node = 0; node < this->connector->nnodes; ++node) {
			// computing hydro vertical motion changes for next time step
			fT dH = this->dt_hydro(node) * (this->_Qw[node] - tQwout[node]) /
							this->connector->get_area_at_node(node);
			this->_hw[node] += dH;
			if (this->_hw[node] < 0) {
				dH -= this->_hw[node];
				this->_hw[node] = 0;
			}
			this->_surface[node] += dH;
		}

		if (this->record_Qw_out)
			this->_rec_Qwout = tQwout;

		// if(this->tau_max > 500)
		//  std::cout << "WARNING::tau_max is " << this->tau_max << std::endl;

		// END OF MAIN LOOP

		// Computing final courant based dt

		// if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT) {
		//   if (this->courant_dt_hydro == -1)
		//     this->courant_dt_hydro = 1e-3;
		//   else if (tcourant_dt_hydro > 0 &&
		//            tcourant_dt_hydro != std::numeric_limits<fT>::max())
		//     this->courant_dt_hydro = tcourant_dt_hydro;
		// }

		// if (this->convergence_mode == CONVERGENCE::ALL ||
		//     this->convergence_mode == CONVERGENCE::DHW) {
		//   for (int i = 0; i < this->n_nodes_convergence(); ++i)
		//     this->conv_dhw[i].emplace_back(
		//         vmot_hw[this->conv_nodes
		//                     [i]]); // /this->dt_hydro(this->conv_nodes[i]));
		// }

		// if (this->convergence_mode == CONVERGENCE::ALL ||
		//     this->convergence_mode == CONVERGENCE::QWR) {
		//   for (int i = 0; i < this->n_nodes_convergence(); ++i)
		//     this->conv_Qr[i].emplace_back(this->_rec_Qwout[this->conv_nodes[i]] /
		//                                   this->_Qw[this->conv_nodes[i]]);
		// }
	}

	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~= Hydro params and related functions ~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	// _      _      _      _      _      _      _      _
	// )`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_)`'-.,_
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
	//~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=

	void setup_entry_points_dynagraph_farea(fT Ath)
	{

		this->init_Qw();

		this->graph_automator();
		std::vector<fT> drarea = this->graph->_get_drainage_area_SFD();

		this->_Qw = (this->_precipitations.size() > 1)
									? this->graph->_accumulate_variable_downstream_SFD(
											this->_precipitations)
									: this->graph->_accumulate_constant_downstream_SFD(
											this->precipitations(0));

		// for(int i = 0;)

		this->_Qw_adder = std::vector<fT>(this->connector->nxy(), 0.);
		this->_water_entry_nodes.clear();
		this->_water_entries.clear();

		int Nentry = 0;
		fT sumprec = 0;

		std::vector<std::uint8_t> done(this->connector->nxy(), false);
		this->water_input_mode = WATER_INPUT::ENTRY_POINTS_H;
		// for(int i= this->connector->nxy()-1; i>=0; --i){
		for (int i = 0; i < this->connector->nxy(); ++i) {
			int node = this->graph->Sstack[i];
			if (drarea[node] > Ath || this->connector->Sreceivers(node) == node) {
				done[node] = true;
				this->_Qw_adder[node] = this->precipitations(node);
				continue;
			}

			if (done[this->connector->Sreceivers(node)]) {
				++Nentry;
				this->_water_entry_nodes.emplace_back(node);
				this->_water_entries.emplace_back(0.);
				this->_Qw_adder[node] = this->_Qw[node];
				sumprec += this->_Qw_adder[node];
			}

			this->_Qw_adder[node] = 0;
		}
		std::cout << "DEBUG::" << Nentry << ":" << sumprec << std::endl;
	}

	void load_up_that_dynastack()
	{

		if (this->hydrostationary) {
			for (size_t i = 0; i < this->_water_entry_nodes.size(); ++i) {

				// this->_Qw[this->_water_entry_nodes[i]] = this->_water_entries[i];
				this->dynastack.emplace(
					dynanode<int, fT>(this->_water_entry_nodes[i],
														this->_surface[this->_water_entry_nodes[i]],
														0.));
			}
		}
	}

	void dynarun()
	{

		if (this->_rec_Qwout.size() == 0) {
			this->enable_Qwout_recording();
			this->_rec_Qwout = std::vector<fT>(this->connector->nxy(), 0.);
		}

		this->tot_Qwin_output = 0;
		this->tot_Qs_output = 0;

		if (this->morphomode != MORPHO::NONE && this->_Qs_adder.size() == 0.)
			this->_Qs_adder = std::vector<fT>(this->connector->nxy(), 0.);

		std::vector<fT> vmot;
		if (this->morphomode != MORPHO::NONE) {
			this->_Qs = std::vector<fT>(this->connector->nxy(), 0.);
			vmot = std::vector<fT>(this->connector->nxy(), 0.);
		}

		// Initialise the water discharge fields according to water input condition
		// and other monitoring features:
		// if(this->hydrostationary == false)
		// 	this->init_Qw();
		this->_Qw = std::vector<fT>(this->connector->nxy(), 0.);
		this->_Qs = std::vector<fT>(this->connector->nxy(), 0.);

		if (this->_hw.size() == 0)
			this->_hw = std::vector<fT>(this->connector->nxy(), 0);

		if (this->dttracker.size() == 0) {
			this->dttracker =
				std::vector<fT>(this->connector->nxy(), this->glob_dynatime - 1);
		}

		// Graph Processing
		// if (this->hydrostationary)
		this->load_up_that_dynastack();

		std::vector<fT> vmot_hw(this->connector->nxy(), 0);

		// Am I in SFD or MFD
		bool SF = (this->hydromode == HYDRO::GRAPH_SFD);

		// To be used if courant dt hydro is selected
		fT tcourant_dt_hydro = std::numeric_limits<fT>::max();

		// Caching neighbours, slopes and weights
		auto receivers = this->connector->get_empty_neighbour();
		std::array<fT, 8> weights, slopes, dxs, dys;
		std::array<int, 16> latnodes;
		std::array<fT, 16> els, dls;
		fT mass_balance = 0;

		// main loop
		int i = -1;
		while (true) {
			// std::cout << "A-1" << std::endl;

			//
			++i;

			bool outs = false;

			dynanode<int, fT> yolo;

			// Getting next node in line
			int node;
			if (this->hydrostationary) {
				if (this->dynastack.empty())
					break;
				yolo = this->dynastack.top();
				node = this->dynastack.top().node;
				while (true) {
					this->dynastack.pop();
					if (this->dynastack.empty())
						break;
					if (this->dynastack.top().node != node)
						break;
					// std::cout << "HAPPENS" << yolo.Qw << std::endl;
					yolo.ingest(this->dynastack.top());
					// std::cout << "AFTER" << yolo.Qw << std::endl;
				}
			} else {
				if (i == this->connector->nxy())
					break;
				node = i;
			}

			bool alreadyDone = this->dttracker[node] == this->glob_dynatime;

			// std::cout << node << std::endl;

			// this->dttracker[node] = this->glob_dynatime;

			if (this->connector->boundaries.can_out(node)) {
				this->dttracker[node] = this->glob_dynatime;

				// std::cout << "happens2::" << node << std::endl;
				continue;
			}

			int nn = this->connector->Neighbours(node, receivers, dxs);
			if (nn == 0) {
				this->dttracker[node] = this->glob_dynatime;
				if (this->connector->boundaries.forcing_io(node))
					std::cout << "happens::" << node << std::endl;
				continue;
			}

			int nrecs = 0;

			int Srecj = -1;
			fT tSS = 0, selSS = 0.;
			;
			fT sumslopes = 0.;
			bool LM = false;
			// std::cout << "A::" << this->_surface[node] << "|" <<
			// this->connector->boundaries.can_out(node) << std::endl;

			while (Srecj == -1) {
				sumslopes = 0.;
				nrecs = 0;
				tSS = 0;
				selSS = 0.;
				for (int j = 0; j < nn; ++j) {
					if (this->_surface[receivers[j]] < this->_surface[node]) {

						receivers[nrecs] = receivers[j];

						if (this->connector->flow_out_model(receivers[nrecs])) {
							dys[0] = this->connector->dy;
							dxs[0] = this->connector->dx;
							nrecs = 1;
							slopes[0] =
								(this->boundhw == BOUNDARY_HW::FIXED_SLOPE)
									? this->bou_fixed_val
									: (this->_surface[node] - this->_surface[receivers[j]]) /
											dxs[j];
							;
							Srecj = 0;
							weights[0] = 1;
							sumslopes = 1;
							outs = true;
							receivers[0] = receivers[j];
							break;
						}

						dxs[nrecs] = dxs[j];

						dys[nrecs] = this->connector->get_travers_dy_from_dx(dxs[j]);

						slopes[nrecs] =
							(this->_surface[node] - this->_surface[receivers[j]]) / dxs[j];

						fT tselSS = (SF) ? slopes[nrecs] * this->connector->randu->get()
														 : slopes[nrecs];

						sumslopes += slopes[nrecs];

						if (slopes[nrecs] > tSS)
							tSS = slopes[nrecs];

						if (tselSS > selSS ||
								this->connector->boundaries.can_out(receivers[j])) {
							selSS = tselSS;
							Srecj = nrecs;
							if (this->connector->boundaries.can_out(receivers[j])) {
								LM = true;
								break;
							}
						}

						++nrecs;
					}
				}
				// std::cout << "A1" << std::endl;

				if (Srecj == -1) {
					LM = true;
					// if(this->hydrostationary){
					fT rdu = this->connector->randu->get() * 1e-4;
					fT tadd = this->dynincr_LM + rdu;
					this->_surface[node] += tadd;
					this->_hw[node] += tadd;
					nn = this->connector->Neighbours(node, receivers, dxs);
					if (nn == 0)
						break;
				}
			}

			if (Srecj == -1 || nrecs == 0) {
				this->dttracker[node] = this->glob_dynatime;
				std::cout << "happens3" << std::endl;
				continue; //  in the case of ddynamic model with LM: I stop here
			}

			fT sumweig = 0;
			if (nrecs == 1) {
				weights[0] = 1.;
				sumweig += 1;
			} else if (sumslopes > 0) {
				for (int j = 0; j < nrecs; ++j) {
					weights[j] = slopes[j] / sumslopes;
					sumweig += weights[j];
				}
			} else {
				for (int j = 0; j < nrecs; ++j) {
					weights[j] = 1. / nrecs;
					sumweig += weights[j];
				}
			}

			if (std::abs(sumweig - 1) > 1e-3 && LM == false)
				std::cout << "sumw" << sumweig << " nrecs " << nrecs << std::endl;

			if ((LM && Srecj > -1) || SF) {
				receivers[0] = receivers[Srecj];
				dxs[0] = dxs[Srecj];
				dys[0] = dys[Srecj];
				slopes[0] = tSS;
				nrecs = 1;
				weights[0] = 1;
				sumslopes = 1;
				Srecj = 0;
			}

			// std::cout << "B" << std::endl;

			// Caching Slope max
			// fT Smax =
			// 	(std::max(this->_surface[node] - this->_surface[receivers[Srecj]],
			// 	 								 1e-6)) /
			// 	dxs[Srecj];

			fT Smax = tSS;

			fT dw0max = dys[Srecj];
			fT dx = dxs[Srecj];
			int recmax = receivers[Srecj];

			// // NOTE:
			// //  No need to calculate the topological number anymore, but keeping it
			// //  for recording its value
			// fT topological_number_v2 = 0.;

			// this->_compute_slopes_weights_et_al(node,
			// 																		SF,
			// 																		Smax,
			// 																		slopes,
			// 																		weights,
			// 																		nrecs,
			// 																		receivers,
			// 																		recmax,
			// 																		dx,
			// 																		dw0max,
			// 																		topological_number_v2);

			// Initialising the total Qout
			fT Qwin = (this->hydrostationary) ? yolo.Qw : this->_Qw[node];

			if (alreadyDone == false)
				Qwin += this->_Qw_adder[node];

			if (this->hydrostationary)
				// if( this->hydrostationary && LM == false)
				this->_Qw[node] = Qwin;

			// precalculating the power
			fT pohw = std::pow(this->_hw[node], TWOTHIRD);
			// std:: cout << "pohw:" << pohw << "|" << this->_hw[node] << std::endl;

			// Squarerooting Smax
			// fT debug_S = Smax;
			auto sqrtSmax = std::sqrt(Smax);

			// Flow Velocity
			fT u_flow = pohw * sqrtSmax / this->mannings(node);
			// Volumetric discahrge
			fT Qwout = dw0max * this->_hw[node] * u_flow;

			// Eventually recording Smax
			if (this->record_Sw)
				this->_rec_Sw[node] = Smax;

			// transfer fluxes
			// Computing the transfer of water and sed
			// this->_compute_transfers(
			// 	nrecs, recmax, node, SF, receivers, weights, Qwin, Qwout);

			// Computing courant based dt
			if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT && u_flow > 0) {
				fT provisional_dt = (this->courant_number * dx) / u_flow;
				provisional_dt = std::min(provisional_dt, this->max_courant_dt_hydro);
				provisional_dt = std::max(provisional_dt, this->min_courant_dt_hydro);
				tcourant_dt_hydro = std::min(provisional_dt, tcourant_dt_hydro);
			}

			this->dttracker[node] = this->glob_dynatime;

			if (this->record_Qw_out)
				this->_rec_Qwout[node] = Qwout;

			if (outs && this->hydrostationary) {
				mass_balance += Qwin;
			}

			fT tQs = 0;

			// MORPHO
			if (alreadyDone == false && this->morphomode != MORPHO::NONE) {

				tQs = yolo.Qs;

				tQs += this->_Qs_adder[node];

				this->_Qs[node] = tQs;

				if (this->connector->boundaries.forcing_io(node) == false) {

					// precaching rates
					fT edot = 0., ddot = 0., totel = 0., totdl = 0.;

					// Lateral dx for lat e/d
					fT dy = dw0max;

					// Calculating sheer stress
					// fT tau = this->rho(node) * this->_hw[node] * GRAVITY * Smax;
					fT slope4shear = 0;
					if (false) {
						slope4shear = 0.;
						for (int j = 0; j < nrecs; ++j)
							slope4shear += slopes[j];

						if (slope4shear <= 0)
							slope4shear = 1e-6;
					} else {
						slope4shear = Smax;
					}

					fT tau = this->rho(node) * this->_hw[node] * GRAVITY * slope4shear;

					if (tau > tau_max)
						this->tau_max = tau;

					// if(tau>150)
					// 	tau = 150;

					// Double checking the orthogonal nodes and if needs be to process
					// them (basically if flow outs model or has boundary exceptions)

					// Calculating erosion if sheer stress is above critical
					if (tau > this->tau_c(node)) {
						// Basal erosion rates
						edot = this->ke(node) *
									 std::pow(tau - this->tau_c(node), this->aexp(node));
						// if recording stuff
						if (this->record_edot)
							this->_rec_edot[node] += edot;
					}

					// Claculating the deposition rates
					ddot = tQs / this->kd(node);

					// std::cout << "A!" << std::endl;

					for (int j = 0; j < nrecs; ++j) {

						int trec = receivers[j];

						if (this->connector->boundaries.can_out(trec) ||
								this->connector->boundaries.forcing_io(trec))
							continue;

						els[j * 2 + 0] = 0.;
						dls[j * 2 + 0] = 0.;
						els[j * 2 + 1] = 0.;
						dls[j * 2 + 1] = 0.;

						// And gathering the orthogonal nodes
						std::pair<int, int> orthonodes =
							this->connector->get_orthogonal_nodes(node, trec);

						latnodes[j * 2 + 0] = orthonodes.first;
						if (this->connector->boundaries.forcing_io(latnodes[j * 2 + 0]) ||
								this->connector->is_in_bound(latnodes[j * 2 + 0]) == false ||
								this->connector->boundaries.no_data(latnodes[j * 2 + 0]) ||
								this->connector->flow_out_model(latnodes[j * 2 + 0]))
							latnodes[j * 2 + 0] = -1;
						latnodes[j * 2 + 1] = orthonodes.second;
						if (this->connector->boundaries.forcing_io(latnodes[j * 2 + 1]) ||
								this->connector->is_in_bound(latnodes[j * 2 + 1]) == false ||
								this->connector->boundaries.no_data(latnodes[j * 2 + 1]) ||
								this->connector->flow_out_model(latnodes[j * 2 + 1]))
							latnodes[j * 2 + 1] = -1;

						// Dealing with lateral deposition if lS > 0 and erosion if lS <0

						if (latnodes[j * 2 + 0] >= 0) {
							fT deltaZ = this->_surface[node] -
													this->_surface[latnodes[j * 2 + 0]] -
													this->_hw[node] + this->_hw[latnodes[j * 2 + 0]];

							// fT tSwl = this->get_Stopo(node, latnodes[j*2+0], dys[j]);
							fT tSwl = deltaZ / dys[j];

							if (tSwl > 0) {
								fT ttdl = std::min(tSwl * this->kd_lateral(node) * ddot,
																	 std::abs(deltaZ) / this->dt_morpho(node));
								totdl += ttdl;
								dls[j * 2 + 0] = ttdl;

							} else {
								fT ttel =
									std::min(std::abs(tSwl) * this->ke_lateral(node) * edot,
													 std::abs(deltaZ) / this->dt_morpho(node));
								totel += ttel;
								els[j * 2 + 0] = ttel;
							}
						}

						if (latnodes[j * 2 + 1] >= 0) {
							fT deltaZ = this->_surface[node] -
													this->_surface[latnodes[j * 2 + 1]] -
													this->_hw[node] + this->_hw[latnodes[j * 2 + 1]];
							// fT tSwl = this->get_Stopo(node, latnodes[j*2+1], dys[j]);
							fT tSwl = deltaZ / dys[j];
							if (tSwl > 0) {
								fT ttdl = std::min(tSwl * this->kd_lateral(node) * ddot,
																	 std::abs(deltaZ) / this->dt_morpho(node));
								totdl += ttdl;
								dls[j * 2 + 1] = ttdl;
							} else {
								fT ttel =
									std::min(std::abs(tSwl) * this->ke_lateral(node) * edot,
													 std::abs(deltaZ) / this->dt_morpho(node));
								totel += ttel;
								els[j * 2 + 1] = ttel;
							}
						}
					}

					// std::cout << "B!" << std::endl;

					// Am I depositing more than I can chew?
					fT totdqs = dx * (totdl + ddot);
					if (totdqs > tQs) {
						// std::cout << " happens??? " << totdqs;
						fT corrqs = tQs / totdqs;
						for (int j = 0; j < nrecs; ++j) {

							if (latnodes[j * 2 + 0] >= 0)
								dls[j * 2 + 0] *= corrqs;
							if (latnodes[j * 2 + 1] >= 0)
								dls[j * 2 + 1] *= corrqs;
						}
						ddot *= corrqs;
						totdl *= corrqs;
					}

					// std::cout << "C!" << std::endl;

					fT sqs = tQs;
					fT fbatch = (ddot + totdl - edot - totel) * dx;
					tQs -= fbatch;

					if (std::isfinite(tQs) == false || std::isfinite(ddot) == false ||
							std::isfinite(edot) == false || std::isfinite(totdl) == false ||
							std::isfinite(totel) == false) {
						std::cout << "QS NAN:" << tQs << " vs " << sqs << " vs " << ddot
											<< " vs " << edot << " ( " << slope4shear << ") "
											<< " vs " << totdl << " vs " << totel << " vs "
											<< std::endl;
						throw std::runtime_error("BITE");
					}

					if (tQs < 0) {
						tQs = 0;
					}

					// std::cout << "D!" << std::endl;
					vmot[node] += ddot - edot;

					for (int j = 0; j < nrecs; ++j) {
						if (this->connector->is_in_bound(latnodes[j * 2 + 0])) {
							vmot[latnodes[j * 2 + 0]] += dls[j * 2 + 0];
							vmot[latnodes[j * 2 + 0]] -= els[j * 2 + 0];
							if (std::isfinite(vmot[latnodes[j * 2 + 0]]) == false) {
								throw std::runtime_error("notfinite LATEREAL");
							}
						}
						if (this->connector->is_in_bound(latnodes[j * 2 + 1])) {
							vmot[latnodes[j * 2 + 1]] += dls[j * 2 + 1];
							vmot[latnodes[j * 2 + 1]] -= els[j * 2 + 1];
							if (std::isfinite(vmot[latnodes[j * 2 + 1]]) == false) {
								throw std::runtime_error("notfinite LATEREAL");
							}
						}
					}

					if (std::isfinite(vmot[node]) == false) {
						std::cout << "edot:" << edot << " ddot" << ddot << std::endl;
						std::cout << "qs:" << sqs << " tau" << tau << std::endl;
						throw std::runtime_error("Non finite vmot gabul	");
					}

					// std::cout << "E!" << std::endl;
				}
			}

			// transfers

			if (SF == false) {

				for (int j = 0; j < nrecs; ++j) {
					fT Qw_next = 0;
					fT Qs_next = 0;

					if (this->hydrostationary) {
						Qw_next = Qwin * weights[j];
						Qs_next = tQs * weights[j];
					} else if (this->hydrostationary == false) {
						if (this->connector->boundaries.forcing_io(node) == false)
							this->_Qw[receivers[j]] += Qwout * weights[j];
						else {
							this->_Qw[receivers[j]] += Qwin * weights[j];
						}
					}

					if (this->hydrostationary) {
						if (this->connector->boundaries.can_out(receivers[j]) == false) {
							this->dynastack.emplace(dynanode<int, fT>(
								receivers[j], this->_surface[receivers[j]], Qw_next, Qs_next));
						} else {
							this->tot_Qwin_output += Qw_next;
							this->tot_Qs_output += Qs_next;
						}
					}
				}

			} else {
				fT Qw_next = 0;
				fT Qs_next = 0;

				if (this->hydrostationary) {
					Qw_next = Qwin;
					Qs_next = tQs;
				} else if (this->hydrostationary == false) {
					this->_Qw[receivers[Srecj]] += Qwout;
				}

				if (this->hydrostationary) {
					if (this->connector->boundaries.can_out(receivers[Srecj]) == false) {
						this->dynastack.emplace(
							dynanode<int, fT>(receivers[Srecj],
																this->_surface[receivers[Srecj]],
																Qw_next,
																Qs_next));
					} else {
						this->tot_Qwin_output += Qw_next;
						this->tot_Qs_output += Qs_next;
					}
				}
			}

			// std::cout << "E2" << std::endl;
		}

		// std::cout << "IOCHECKS::" << mass_balance << std::endl;;
		// if(this->tau_max > 500)
		// 	std::cout << "WARNING::tau_max is " << this->tau_max <<
		// std::endldynanode<int,fT>;

		// END OF MAIN LOOP

		// Computing final courant based dt

		if (this->mode_dt_hydro == PARAM_DT_HYDRO::COURANT) {
			if (this->courant_dt_hydro == -1)
				this->courant_dt_hydro = 1e-3;
			else if (tcourant_dt_hydro > 0 &&
							 tcourant_dt_hydro != std::numeric_limits<fT>::max())
				this->courant_dt_hydro = tcourant_dt_hydro;
		}

		// if(this->hydrostationary == false)
		for (int i = 0; i < this->connector->nxy(); ++i) {
			if (this->connector->flow_out_model(i) == false) {

				// IS THAT NECESSARY????
				if (this->hydrostationary &&
						this->dttracker[i] != this->glob_dynatime && this->_hw[i] > 0) {

					int nn = this->connector->Neighbours(i, receivers, dxs);
					fT tSS = 1e-9;
					fT dw = this->connector->dy;

					for (int j = 0; j < nn; ++j) {
						if (this->_surface[receivers[j]] < this->_surface[i]) {
							fT ttSS =
								(this->_surface[i] - this->_surface[receivers[j]]) / dxs[j];
							if (ttSS > tSS) {
								tSS = ttSS;
								dw = this->connector->get_travers_dy_from_dx(dxs[j]);
							}
						}
					}

					this->_rec_Qwout[i] = std::pow(this->_hw[i], TWOTHIRD) *
																this->_hw[i] / this->mannings(i) * dw *
																std::sqrt(tSS);
				}

				vmot_hw[i] = (this->_Qw[i] - this->_rec_Qwout[i]) /
										 this->connector->get_area_at_node(i);
			}
		}

		for (int i = 0; i < vmot.size(); ++i) {
			if (std::isfinite(vmot[i]) == false) {
				throw std::runtime_error("walkitu");
			}
		}
		// else{
		// 	for (int i = 0; i < this->connector->nxy(); ++i) {
		// 		vmot_hw[i] -= this->_rec_Qwout[i]/
		// 								 this->connector->get_area_at_node(i);
		// 	}
		// }

		// Applying vmots with the right dt
		if (this->morphomode != MORPHO::NONE) {
			this->_compute_vertical_motions(vmot_hw, vmot);
		} else
			this->_compute_vertical_motions_hw(vmot_hw);

		this->glob_dynatime += this->dt_hydro(0);
		// std::cout << "F!" << std::endl;
	}

	void set_out_QS_as_input_QS(fT multiplier)
	{
		if (this->sed_input_mode != SED_INPUT::ENTRY_POINTS_Q) {
			throw std::runtime_error("Cannot set out Qs as input Qs in other mode "
															 "than SED_INPUT::ENTRY_POINTS_Q");
		}

		fT dat_entry =
			this->tot_Qs_output / this->_sed_entry_nodes.size() * multiplier;
		for (size_t i = 0; i < this->_sed_entry_nodes.size(); ++i) {
			this->_Qs_adder[this->_sed_entry_nodes[i]] = dat_entry;
		}
	}

	void set_out_QS_as_input_QS_min(fT multiplier, fT minV)
	{
		if (this->sed_input_mode != SED_INPUT::ENTRY_POINTS_Q) {
			throw std::runtime_error("Cannot set out Qs as input Qs in other mode "
															 "than SED_INPUT::ENTRY_POINTS_Q");
		}

		fT dat_entry = std::max(
			this->tot_Qs_output / this->_sed_entry_nodes.size() * multiplier, minV);

		for (size_t i = 0; i < this->_sed_entry_nodes.size(); ++i) {
			this->_Qs_adder[this->_sed_entry_nodes[i]] = dat_entry;
		}
	}
};
// end of graphflood class

} // namespace DAGGER
