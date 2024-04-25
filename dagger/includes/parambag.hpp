#pragma once

#include "boundary_conditions.hpp"
#include "enumutils.hpp"
#include "graphflood_enums.hpp"
#include "lookup_neighbourer.hpp"
#include "utils.hpp"
#include "wrap_helper.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class ParamBag
{
public:
	// Empty constructor
	ParamBag(){};

	// GRAPHFLOOD PARAMETERS

	// CONSTANTS
	f_t GRAVITY = 9.8;
	f_t rho_water = 1000;
	f_t rho_sed = 2300;

	// Stuff I'll vary later
	f_t MPM_THETHA_C = 0.047;
	f_t tau_c = 2.;
	f_t alpha = 1.5;
	f_t E = 1.;

	// threshold for optimisation (might have consequences)
	f_t minWeightQw = 0.;
	f_t capacityFacQw = 1.;
	f_t capacityFacQs = 1.;

	void set_minWeightQw(f_t val) { this->minWeightQw = val; }
	f_t get_minWeightQw() { return this->minWeightQw; }

	void set_capacityFacQw(f_t val) { this->capacityFacQw = val; }
	f_t get_capacityFacQw() { return this->capacityFacQw; }

	void set_capacityFacQs(f_t val) { this->capacityFacQs = val; }
	f_t get_capacityFacQs() { return this->capacityFacQs; }

	MORPHOMODE morphomode = MORPHOMODE::NONE;
	MORPHOMODE get_morphomode() { return this->morphomode; }
	void set_morphomode(MORPHOMODE mode) { this->morphomode = mode; }

	HYDROMODE hydromode = HYDROMODE::MFD;
	HYDROMODE get_hydromode() { return this->hydromode; }
	void set_hydromode(HYDROMODE mode) { this->hydromode = mode; }

	BOUNDARY_HW gf2Bmode = BOUNDARY_HW::FIXED_SLOPE;
	f_t gf2Bbval = 1e-2;
	void set_gf2Bbval(f_t val) { this->gf2Bbval = val; }
	f_t get_gf2Bbval() { return this->gf2Bbval; }

	bool gf2_morpho = false;
	void enable_gf2_morpho() { this->gf2_morpho = true; }
	void disable_gf2_morpho() { this->gf2_morpho = false; }

	f_t time_dilatation_morpho = 1.;
	void set_time_dilatation_morpho(f_t val)
	{
		this->time_dilatation_morpho = val;
	}
	f_t get_time_dilatation_morpho() { return this->time_dilatation_morpho; }

	// -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
	// -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
	// -~-~-~-~-~-~-~-~-~- ke is the E*xi in MPM in me code -~-~-~-~-~-~-~-~-~
	// -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
	// -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
	f_t _ke = 1e-3;
	std::vector<f_t> _ke_v;
	PARAMTYPE p_ke = PARAMTYPE::CTE;
	void set_ke(f_t val)
	{
		this->p_ke = PARAMTYPE::CTE;
		this->_ke = val;
	}
	f_t get_ke() { return this->_ke; }
	f_t ke(i_t node)
	{
		return (this->p_ke == PARAMTYPE::CTE) ? this->_ke : this->_ke_v[node];
	}

	template<class arrin_t>
	void set_variable_ke(arrin_t& tarr)
	{
		this->p_ke = PARAMTYPE::VAR;
		auto arr = format_input<arrin_t>(tarr);
		this->_ke_v = to_vec(arr);
	}

	void calculate_ke_tau_c_from_MPM(f_t D)
	{
		f_t R = this->rho_sed / this->rho_water - 1;
		this->tau_c = this->rho_water * this->GRAVITY * R * D * this->MPM_THETHA_C;
		this->E = 8 / (std::sqrt(this->rho_water) *
									 (this->rho_sed - this->rho_water) * this->GRAVITY);
		this->set_ke(this->E / this->kd);
		std::cout << "DEBUG:: tau_c = " << this->tau_c
							<< " ke (erosion coeff) : " << E / this->kd << std::endl;
	}

	f_t kd = 10;
	void set_kd(f_t val) { this->kd = val; }
	f_t get_kd() { return this->kd; }

	bool bank_erosion = false;
	void enable_bank_erosion() { this->bank_erosion = true; }
	void disable_bank_erosion() { this->bank_erosion = false; }
	f_t kel = 0.1;
	void set_kel(f_t val) { this->kel = val; }
	f_t get_kel() { return this->kel; }

	f_t kdl = 0.1;
	void set_kdl(f_t val) { this->kdl = val; }
	f_t get_kdl() { return this->kdl; }

	// # Graphflood experimental stuff

	f_t TSG_dist = false;
	void enable_TSG_dist() { this->TSG_dist = true; }
	void disable_TSG_dist() { this->TSG_dist = false; }
	f_t TSG_distmax = 1e9;
	void set_TSG_distmax(f_t val) { this->TSG_distmax = val; }

	bool gf2_diffuse_Qwin = false;
	void enable_gf2_diffuse_Qwin() { this->gf2_diffuse_Qwin = true; }
	void disable_gf2_diffuse_Qwin() { this->gf2_diffuse_Qwin = false; }

	bool transient_flow = false;
	void enable_transient_flow() { this->transient_flow = true; }
	void disable_transient_flow() { this->transient_flow = false; }
};

};
