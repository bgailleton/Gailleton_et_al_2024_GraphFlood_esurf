//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

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

// // Difines the different boundary conditions for a single node

// class BoundaryCondition
// {

// public:
// 	BoundaryCondition();

// 	// Can flow leave the model (if true, flow paths can stop there)
// 	template<class i_t>
// 	virtual bool can_out(i_t i) = 0;
// 	// if true, the flow is stopped there
// 	template<class i_t>
// 	virtual bool force_out(i_t i) = 0;
// 	// if true node is nodata and ignored
// 	template<class i_t>
// 	virtual bool nodata(i_t i) = 0;
// 	// if true, the node is and will always be receiver of its neighbours
// 	template<class i_t>
// 	virtual bool forced_receiver(i_t i) = 0;
// 	// if truem the node is and will always be a donor
// 	template<class i_t>
// 	virtual bool forced_donor(i_t i) = 0;

// };

namespace DAGGER {

enum class BC : std::uint8_t
{
	// Cannot flow at all = nodata
	NO_FLOW = 0,

	// Internal Node (can flow in every directions)
	FLOW = 1,

	// Internal Node (can flow in every directions) BUT neighbours a special flow
	// condition and may need specific care
	FLOW_BUT = 2,

	// flow can out there but can also flow to downstream neighbours
	CAN_OUT = 3,

	// flow can only out from this cell
	OUT = 4,

	// Not only flow HAS to out there: neighbouring flows will be drained there no
	// matter what
	FORCE_OUT = 5,

	// Flows through the cell is possible, but the cell CANNOT out fluxes from
	// this boundary (reserved to model edges, internal boundaries wont give to
	// nodata anyway)
	CANNOT_OUT = 6,

	// Flow can only flow to potential receivers
	IN = 7,

	// Forced INFLOW: flow will flow to all neighbours (except other FORCE_IN)
	FORCE_IN = 8,

	// periodic border
	PERIODIC_BORDER = 9

};

std::string
BC2str(BC tbc)
{
	if (tbc == BC::NO_FLOW)
		return "NO_FLOW";
	if (tbc == BC::FLOW)
		return "FLOW";
	if (tbc == BC::FLOW_BUT)
		return "FLOW_BUT";
	if (tbc == BC::CAN_OUT)
		return "CAN_OUT";
	if (tbc == BC::OUT)
		return "OUT";
	if (tbc == BC::CANNOT_OUT)
		return "CANNOT_OUT";
	if (tbc == BC::FORCE_OUT)
		return "FORCE_OUT";
	if (tbc == BC::FORCE_IN)
		return "FORCE_IN";
	if (tbc == BC::PERIODIC_BORDER)
		return "PERIODIC_BORDER";

	return "UNREGISTERED BC";
}

bool
Fcan_connect(BC A, BC B)
{
	if ((A == BC::FORCE_IN || A == BC::FORCE_OUT) &&
			(B == BC::FORCE_IN || B == BC::FORCE_OUT))
		return false;
	return true;
}

bool
can_out(BC A)
{
	return (A == BC::CAN_OUT || A == BC::OUT || A == BC::FORCE_OUT);
}

bool
can_receive(BC A)
{
	return (A == BC::CAN_OUT || A == BC::OUT || A == BC::FORCE_OUT ||
					A == BC::FLOW || A == BC::PERIODIC_BORDER);
}

bool
nodata(BC A)
{
	return (A == BC::NO_FLOW);
}

// Function testing if a boundary condition object is valid
// template<class TBC>
// void test_TBC(TBC)

// Class complying with the boundary conditions standards
class CodeBC
{
public:
	CodeBC() { ; }

	// vector of BC codes, need to be of node side
	std::vector<BC> codes;

	// return true if the node is an internal "classic" node and no specific
	// conditions are required this will be the case of the huge majority of nodes
	// in most cases, that is why it is important to separate and prioritise their
	// detection for the sake of efficiency
	template<class i_t>
	bool is_normal_node(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::FLOW)
			return true;
		else
			return false;
	}

	// return true if the node is an internal "classic" node but borders a
	// boundary condition requiring specific care] this also is an import test to
	// optimise and prioritise: a receivers of a neigh
	template<class i_t>
	bool is_indirect_bc(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::FLOW_BUT)
			return true;
		else
			return false;
	}

	// Return true if the node is a special boundary, as opposed to a node needind
	// specific care
	template<class i_t>
	bool is_bc(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::FLOW_BUT || tbc == BC::FLOW)
			return false;
		else
			return true;
	}

	// Return true if the node can receive flow
	template<class i_t>
	bool can_receive(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::NO_FLOW || tbc == BC::FORCE_IN || tbc == BC::IN)
			return false;
		return true;
	}

	// Return true if the node can give flow
	template<class i_t>
	bool can_give(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::NO_FLOW || tbc == BC::FORCE_OUT || tbc == BC::OUT)
			return false;
		return true;
	}

	// Return true if the node can give flow
	template<class i_t>
	bool can_out(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::FORCE_OUT || tbc == BC::OUT || tbc == BC::CAN_OUT)
			return true;
		return false;
	}

	// Return true if the node can receive and give flow
	template<class i_t>
	bool can_flow_through(i_t i)
	{
		BC tbc = this->codes[i];
		if (this->can_give(tbc) && this->can_receive(tbc))
			return true;
		return false;
	}

	// Return true if the node is no data
	template<class i_t>
	bool no_data(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::NO_FLOW)
			return true;
		return false;
	}

	// return true if the node can receive or give link
	template<class i_t>
	bool can_create_link(i_t i)
	{
		if (!this->no_data(i) && (this->can_receive(i) || this->can_give(i)))
			return true;

		return false;
	}

	// return true if the node forces a flow condition, either in or out
	template<class i_t>
	bool forcing_io(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::FORCE_IN || tbc == BC::FORCE_OUT)
			return true;
		return false;
	}

	// return true if the node forces the flow to give
	template<class i_t>
	bool force_giving(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::FORCE_IN)
			return true;
		return false;
	}

	// return true if the node forces the flow to receive
	template<class i_t>
	bool force_receiving(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::FORCE_OUT)
			return true;
		return false;
	}

	template<class i_t>
	bool is_periodic(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::PERIODIC_BORDER)
			return true;
		return false;
	}

	template<class i_t>
	bool has_to_out(i_t i)
	{
		BC tbc = this->codes[i];
		if (tbc == BC::OUT || tbc == BC::FORCE_OUT)
			return true;
		return false;
	}

	// Return true if the boundary needs normal neighbouring
	// It does not mean the link can be  created, it just means the n eighbour to
	// get and check are "normallay" following the edges
	template<class i_t>
	bool normal_neighbouring_at_boundary(i_t i)
	{
		return !this->is_periodic(i);
	}

	void set_codes(std::vector<BC>& ncodes)
	{
		if (ncodes.size() == this->codes.size())
			this->codes = std::move(ncodes);
		else
			throw std::runtime_error(
				"Trying to set boundary codes with vector of wrong size");
	}
};

}; // namespace DAGGER

#endif
