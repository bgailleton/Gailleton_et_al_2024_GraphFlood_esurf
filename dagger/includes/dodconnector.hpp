#pragma once

#include "boundary_conditions.hpp"
#include "databag.hpp"
#include "dodcontexts.hpp"
#include "enumutils.hpp"
#include "lookup_neighbourer.hpp"
#include "priority_flood.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class Connector8
{

public:
	// Grid dimensions
	// # total number of nodes
	int _nxy = 0;
	int nxy() const { return this->_nxy; }
	// # Number of columns
	int _nx = 0;
	// # Number of rows
	int _ny = 0;
	// # Spacing in the x dimension
	f_t _dx = 1.;
	// # Spacing in the y dimension
	f_t _dy = 1.;
	// # Spacing in the diagonal dimension
	f_t _dxy = std::sqrt(2.);
	f_t _area = 1.;
	// # Total Length in the x dimension
	f_t _lx = 0.;
	// # Total Length in the y dimension
	f_t _ly = 0.;

	f_t area(i_t i) const { return this->_area; }

	// Data storer
	Hermes<i_t, f_t>* data;

	// Param sheet
	// # Flow topology to compute
	CONFLOWTOPO flowtopo = CONFLOWTOPO::ALL;
	void set_flowtopo(CONFLOWTOPO tfe) { this->flowtopo = tfe; }
	// # Boundary conditions
	CONBOU boutype = CONBOU::EDGES;
	void set_condou(CONBOU tfe) { this->boutype = tfe; }

	// Default Constructor
	Connector8() { ; }

	// Constructor
	Connector8(int nx, int ny, f_t dx, f_t dy, Hermes<i_t, f_t>& data)
	{
		// Initialising the geometry of the grid

		// # Number of col,rows
		this->_nx = nx;
		this->_ny = ny;
		// # Number of nodes
		this->_nxy = nx * ny;
		// # spatial spacing
		this->_dx = dx;
		this->_dy = dy;
		this->_dxy = std::sqrt(dx * dx + dy * dy);
		this->_area = this->_dx * this->_dy;
		// # grid total length
		this->_lx = (nx + 1) * dx;
		this->_ly = (ny + 1) * dy;

		// # Connecting the data bag
		this->data = &data;

		// Initialisatio of the lookup table
		this->data->LK8 =
			lookup8<i_t, f_t>(this->_nx, this->_ny, this->_dx, this->_dy);
	}

	// Initialisation of the model
	// The initialisation has to be run after the construction as it depends on
	// the different params of the parameter sheet
	void init()
	{

		// Initialising the boundary conditions
		// # If boundaries are EDGES, automatically sets them to can out at all
		// edges
		if (this->boutype == CONBOU::EDGES) {
			// # set all boundary codes to normal
			this->data->_boundaries = std::vector<BC>(this->nxy(), BC::FLOW);
			for (int i = 0; i < this->_nx; ++i)
				this->data->_boundaries[i] = BC::OUT;
			for (int i = this->_nxy - this->_nx; i < this->_nxy; ++i)
				this->data->_boundaries[i] = BC::OUT;
			for (int i = 0; i < this->_ny; ++i) {
				this->data->_boundaries[i * this->_nx] = BC::OUT;
				this->data->_boundaries[i * this->_nx + this->_nx - 1] = BC::OUT;
			}

		} else if (this->boutype == CONBOU::PEW) {
			this->data->_boundaries = std::vector<BC>(this->nxy(), BC::FLOW);
			for (int i = 0; i < this->_ny; ++i) {
				this->data->_boundaries[i * this->_nx] = BC::PERIODIC_BORDER;
				this->data->_boundaries[i * this->_nx + this->_nx - 1] =
					BC::PERIODIC_BORDER;
			}
			for (int i = 0; i < this->_nx; ++i)
				this->data->_boundaries[i] = BC::OUT;
			for (int i = this->_nxy - this->_nx; i < this->_nxy; ++i)
				this->data->_boundaries[i] = BC::OUT;

		} else if (this->boutype == CONBOU::CUSTOM) {
			if (this->data->_boundaries.size() != static_cast<size_t>(this->nxy()))
				throw std::runtime_error("cannot init connector: boutype set to custom "
																 "but no boundaires array in the data bag");
		} else {
			throw std::runtime_error("boutype NOT IMPLOEMENTED YET");
		}

		// Computing the neighbour code
		this->_compute_neighbours();

		// Initialising the data concerned by the different topologies
		if (this->flowtopo == CONFLOWTOPO::ALL) {
			this->data->_Sreceivers = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_Sdonors = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_donors = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_receivers = std::vector<std::uint8_t>(this->_nxy, 0);
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			this->data->_donors = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_receivers = std::vector<std::uint8_t>(this->_nxy, 0);
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			this->data->_Sreceivers = std::vector<std::uint8_t>(this->_nxy, 0);
			this->data->_Sdonors = std::vector<std::uint8_t>(this->_nxy, 0);
		}
	}

	// Operation reinitialising the data.
	// More efficient than recalling init at every updates as it skips the memory
	// allocation bits
	void reinit()
	{
		if (this->flowtopo == CONFLOWTOPO::ALL) {
			fillvec(this->data->_Sreceivers, static_cast<std::uint8_t>(0));
			fillvec(this->data->_Sdonors, static_cast<std::uint8_t>(0));
			fillvec(this->data->_donors, static_cast<std::uint8_t>(0));
			fillvec(this->data->_receivers, static_cast<std::uint8_t>(0));
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			fillvec(this->data->_donors, static_cast<std::uint8_t>(0));
			fillvec(this->data->_receivers, static_cast<std::uint8_t>(0));
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			fillvec(this->data->_Sreceivers, static_cast<std::uint8_t>(0));
			fillvec(this->data->_Sdonors, static_cast<std::uint8_t>(0));
		}
	}

	// Access to the the steepest receiver of node i
	i_t Sreceivers(i_t i) const
	{
		return i + this->data->LK8.NeighbourerD8[this->data->LK8.BC2idAdder(
								 i, this->data->_boundaries[i])][this->data->_Sreceivers[i]];
	}

	// Access to the the steepest receiver of node i
	f_t SreceiversDx(i_t i) const
	{
		return this->data->LK8.NeighbourerD8dx[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_Sreceivers[i]];
	}

	// Access to the the steepest receiver of node i
	f_t SreceiversDy(i_t i) const
	{
		f_t tdx = this->data->LK8.NeighbourerD8dx[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_Sreceivers[i]];
		return (tdx == this->_dx) ? this->_dy
															: ((tdx == this->_dy) ? this->_dx : this->_dxy);
	}

	// Access to the the steepest receiver of node i
	std::uint8_t SreceiversBit(i_t i) const { return this->data->_Sreceivers[i]; }

	// Access to the the steepest donors of node i
	i_t Sdonors(i_t i, std::array<i_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourer[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_Sdonors[i]];
		i_t nn = this->nSdonors(i);
		for (int j = 0; j < nn; ++j)
			arr[j] += i;
		return nn;
	}

	// Access to the the number of steepest donors of node i
	i_t nSdonors(i_t i) const
	{
		return this->data->LK8.NeighbourerNN[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_Sdonors[i]];
	}

	// Access to donors of node i
	i_t Donors(i_t i, std::array<i_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourer[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_donors[i]];
		i_t nn = this->nDonors(i);
		for (int j = 0; j < nn; ++j)
			arr[j] += i;
		return nn;
	}

	// Access to donors bitset of node i
	i_t DonorsBits(i_t i, std::array<i_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourer[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_donors[i]];
		i_t nn = this->nDonors(i);
		for (int j = 0; j < nn; ++j)
			arr[j] += i;
		return nn;
	}

	// Access to number of donors of node i
	i_t nDonors(i_t i) const
	{
		return this->data->LK8.NeighbourerNN[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_donors[i]];
	}

	// Access to Receivers of node i
	i_t Receivers(i_t i, std::array<i_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourer[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_receivers[i]];
		i_t nn = this->nReceivers(i);
		for (int j = 0; j < nn; ++j)
			arr[j] += i;
		return nn;
	}

	// Access to Receivers bitcodes of node i
	i_t ReceiversBits(i_t i, std::array<std::uint8_t, 8>& arr) const
	{
		arr = this->data->LK8.NeighbourerBits[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_receivers[i]];
		i_t nn = this->nReceivers(i);
		return nn;
	}

	// Access to number of  Receivers of node i
	i_t nReceivers(i_t i) const
	{
		return this->data->LK8.NeighbourerNN[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_receivers[i]];
	}

	// Access to Neighbours of node i
	i_t Neighbours(i_t i, std::array<i_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourer[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_neighbours[i]];
		i_t nn = this->nNeighbours(i);
		for (int j = 0; j < nn; ++j)
			arr[j] += i;
		return nn;
	}

	// Gets the neighbours from a polar coordinate (similar-ish to D infinity)
	void NeighboursTheta2(i_t i, f_t theta, i_t& n1, i_t& n2, f_t& w1, f_t& w2)
	{
		// correcting if theta > pi or <=pi
		if (std::abs(theta) > NPPI) {
			theta /= NPPI;
			if (std::signbit(theta))
				theta = NPPI - (std::abs(theta) - 1) * NPPI;
			else
				theta = -NPPI + (std::abs(theta) - 1) * NPPI;
		}

		// get a quadran between 0 (right) and 4 (left)
		f_t div = std::abs(theta) / (NPPI / 4.);
		int quarter = std::floor(div);
		div = div - quarter;
		bool isneg = std::signbit(theta);

		std::uint8_t cn1;
		std::uint8_t cn2;

		if (quarter == 4) {
			cn1 = LeftMask8;
			cn2 = TopLeftMask8;
			w1 = 1.;
			w2 = 0.;
		} else if (quarter == 3) {
			if (isneg) {
				cn1 = BottomLeftMask8;
				cn2 = LeftMask8;
				w1 = 1 - div;
				w2 = div;
			} else {
				cn1 = TopLeftMask8;
				cn2 = LeftMask8;
				w1 = 1 - div;
				w2 = div;
			}
		} else if (quarter == 2) {
			if (isneg) {
				cn1 = BottomMask8;
				cn2 = BottomLeftMask8;
				w1 = 1 - div;
				w2 = div;
			} else {
				cn1 = TopMask8;
				cn2 = TopLeftMask8;
				w1 = 1 - div;
				w2 = div;
			}
		} else if (quarter == 1) {
			if (isneg) {
				cn1 = BottomRightMask8;
				cn2 = BottomMask8;
				w1 = 1 - div;
				w2 = div;
			} else {
				cn1 = TopRightMask8;
				cn2 = TopMask8;
				w1 = 1 - div;
				w2 = div;
			}
		} else {
			if (isneg) {
				cn1 = RightMask8;
				cn2 = BottomRightMask8;
				w1 = 1 - div;
				w2 = div;
			} else {
				cn1 = RightMask8;
				cn2 = TopRightMask8;
				w1 = 1 - div;
				w2 = div;
			}
		}

		auto adder = this->data->LK8.BC2idAdder(i, this->data->_boundaries[i]);
		n1 = i + this->data->LK8.NeighbourerD8[adder][cn1];
		n2 = i + this->data->LK8.NeighbourerD8[adder][cn2];
		if (std::abs(1 - w1 - w2) > 1e-4)
			throw std::runtime_error(std::to_string(w1) + "|" + std::to_string(w2));
	}

	// Access to steepest rec of node i but reprocessed on the spot (as opposed to
	// preprocessed from a compute() operation)
	int Sreceivers_raw(i_t i,
										 std::array<i_t, 8>& arr,
										 std::array<f_t, 8>& arrdx,
										 i_t& trec,
										 f_t& tdx) const
	{
		i_t nn = this->Neighbours(i, arr);
		nn = this->NeighboursDx(i, arrdx);
		trec = -1;
		tdx = this->_dx;
		f_t tSS = 0.;
		for (int j = 0; j < nn; ++j) {
			arr[j] += i;
			f_t ttSS =
				(this->data->_surface[i] - this->data->_surface[arr[j]]) / arrdx[j];
			if (ttSS > tSS) {
				tSS = ttSS;
				trec = arr[j];
				tdx = arrdx[j];
			}
		}

		return nn;
	}

	// Access to dx to each neighbours (distance to nodes) of node i
	i_t NeighboursDx(i_t i, std::array<f_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourerdx[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_neighbours[i]];
		i_t nn = this->nNeighbours(i);
		return nn;
	}

	// Access to dx to each neighbours (distance to nodes) of node i
	i_t ReceiversDx(i_t i, std::array<f_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourerdx[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_receivers[i]];
		i_t nn = this->nReceivers(i);
		return nn;
	}

	// Access to dx to each neighbours (distance to nodes) of node i
	i_t DonorsDx(i_t i, std::array<f_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourerdx[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_donors[i]];
		i_t nn = this->nDonors(i);
		return nn;
	}

	// Access to dx to each neighbours (distance to nodes) of node i
	i_t NeighboursDxSign(i_t i, std::array<std::int8_t, 8>& arr) const
	{
		arr = this->data->LK8.NeighbourerdxSign[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_neighbours[i]];
		i_t nn = this->nNeighbours(i);
		return nn;
	}

	// Access to dx to each neighbours (distance to nodes) of node i
	i_t ReceiversDxSign(i_t i, std::array<std::int8_t, 8>& arr) const
	{
		arr = this->data->LK8.NeighbourerdxSign[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_receivers[i]];
		i_t nn = this->nReceivers(i);
		return nn;
	}

	// Access to dx to each neighbours (distance to nodes) of node i
	i_t DonorsDxSign(i_t i, std::array<std::int8_t, 8>& arr) const
	{
		arr = this->data->LK8.NeighbourerdxSign[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_donors[i]];
		i_t nn = this->nDonors(i);
		return nn;
	}

	// Access to dx to each neighbours (distance to nodes) of node i
	i_t NeighboursDySign(i_t i, std::array<std::int8_t, 8>& arr) const
	{
		arr = this->data->LK8.NeighbourerdySign[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_neighbours[i]];
		i_t nn = this->nNeighbours(i);
		return nn;
	}

	// Access to dy to each neighbours (distance to nodes) of node i
	i_t ReceiversDySign(i_t i, std::array<std::int8_t, 8>& arr) const
	{
		arr = this->data->LK8.NeighbourerdySign[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_receivers[i]];
		i_t nn = this->nReceivers(i);
		return nn;
	}

	// Access to dy to each neighbours (distance to nodes) of node i
	i_t DonorsDySign(i_t i, std::array<std::int8_t, 8>& arr) const
	{
		arr = this->data->LK8.NeighbourerdySign[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_donors[i]];
		i_t nn = this->nDonors(i);
		return nn;
	}

	// Access to dy to each neighbours (distance to nodes) of node i
	i_t NeighboursDy(i_t i, std::array<f_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourerdy[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_neighbours[i]];
		i_t nn = this->nNeighbours(i);
		return nn;
	}

	// Access to dy to each neighbours (distance to nodes) of node i
	i_t ReceiversDy(i_t i, std::array<f_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourerdy[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_receivers[i]];
		i_t nn = this->nReceivers(i);
		return nn;
	}

	// Access to dy to each neighbours (distance to nodes) of node i
	i_t DonorsDy(i_t i, std::array<f_t, 8>& arr) const
	{
		arr = this->data->LK8.Neighbourerdy[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_donors[i]];
		i_t nn = this->nDonors(i);
		return nn;
	}

	// Older versions without direct access
	// // Access to dy to each neighbours (distance to orthogonal nodes) of node i
	// i_t NeighboursDy(i_t i, std::array<f_t, 8>& arr) const
	// {
	// 	arr = this->data->LK8.Neighbourerdx[this->data->LK8.BC2idAdder(
	// 		i, this->data->_boundaries[i])][this->data->_neighbours[i]];
	// 	i_t nn = this->nNeighbours(i);
	// 	for (int j = 0; j < nn; ++j)
	// 		arr[j] = (arr[j] == this->_dx)
	// 							 ? this->_dy
	// 							 : ((arr[j] == this->_dy) ? this->_dx : this->_dxy);
	// 	return nn;
	// }

	// // Access to dy to each neighbours (distance to orthogonal nodes) of node i
	// i_t ReceiversDy(i_t i, std::array<f_t, 8>& arr) const
	// {
	// 	arr = this->data->LK8.Neighbourerdx[this->data->LK8.BC2idAdder(
	// 		i, this->data->_boundaries[i])][this->data->_receivers[i]];
	// 	i_t nn = this->nReceivers(i);
	// 	for (int j = 0; j < nn; ++j)
	// 		arr[j] = (arr[j] == this->_dx)
	// 							 ? this->_dy
	// 							 : ((arr[j] == this->_dy) ? this->_dx : this->_dxy);
	// 	return nn;
	// }

	// Access to Neighbours' bitcodes of node i
	i_t NeighboursBits(i_t i, std::array<std::uint8_t, 8>& arr) const
	{
		arr = this->data->LK8.NeighbourerBits[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_neighbours[i]];
		i_t nn = this->nNeighbours(i);
		return nn;
	}

	// Access to number of neighbours of node i
	i_t nNeighbours(i_t i) const
	{
		return this->data->LK8.NeighbourerNN[this->data->LK8.BC2idAdder(
			i, this->data->_boundaries[i])][this->data->_neighbours[i]];
	}

	// one-off computing operation: it computes the neighbours code to loop
	// through the different neighbours of each node
	void _compute_neighbours()
	{

		// Initialising the neighbour array
		this->data->_neighbours = std::vector<std::uint8_t>(this->_nxy, 0);

		// Setting up some refs for code clarity
		std::vector<std::uint8_t>& neighbours = this->data->_neighbours;
		std::vector<BC>& boundaries = this->data->_boundaries;
		auto& LK8 = this->data->LK8;

		// Starting with internal nodes
		// | | | |
		// | |X| |
		// | | | |
		for (int r = 1; r < this->_ny - 1; ++r) {
			for (int c = 1; c < this->_nx - 1; ++c) {
				// Getting the node index
				int node = r * this->_nx + c;

				// Default to nonei
				neighbours[node] = 0;

				// if node is no flow - skip
				if (boundaries[node] == BC::NO_FLOW)
					continue;

				// Going through theoretical neighbours
				for (auto TTN : NeighbourerMask8) {
					// getting the index
					int TN = node + LK8.NeighbourerD8[this->data->LK8.BC2idAdder(
														node, this->data->_boundaries[node])][TTN];
					// if neighbour is no flow ignore it
					if (boundaries[TN] == BC::NO_FLOW)
						continue;
					// Else it's basically a neighbour and that's it
					neighbours[node] |= TTN;
				}
			}
		}

		// Maninging topleft corner
		// |X| | |
		// | | | |
		// | | | |
		if (boundaries[0] != BC::NO_FLOW) {
			if (boundaries[0] != BC::PERIODIC_BORDER) {
				neighbours[0] = LK8.TopLeft_normal_boundary();
			} else {
				neighbours[0] = AllMask8;
			}
		}

		// Maninging topright corner
		// | | |X|
		// | | | |
		// | | | |
		if (boundaries[this->_nx - 1] != BC::NO_FLOW) {
			if (boundaries[this->_nx - 1] != BC::PERIODIC_BORDER) {
				neighbours[this->_nx - 1] = LK8.TopRight_normal_boundary();
			} else {
				neighbours[this->_nx - 1] = AllMask8;
			}
		}

		// Mananging bottomleft corner
		// | | | |
		// | | | |
		// |X| | |
		if (boundaries[this->nxy() - this->_nx] != BC::NO_FLOW) {
			if (boundaries[this->nxy() - this->_nx] != BC::PERIODIC_BORDER) {
				neighbours[this->nxy() - this->_nx] = LK8.BottomLeft_normal_boundary();
			} else {
				neighbours[this->nxy() - this->_nx] = AllMask8;
			}
		}

		// Mananging bottomright corner
		// | | | |
		// | | | |
		// | | |X|
		if (boundaries[this->nxy() - 1] != BC::NO_FLOW) {
			if (boundaries[this->nxy() - 1] != BC::PERIODIC_BORDER) {
				neighbours[this->nxy() - 1] = LK8.BottomRight_normal_boundary();
			} else {
				neighbours[this->nxy() - 1] = AllMask8;
			}
		}

		// Managing first row
		// | |X| |
		// | | | |
		// | | | |
		for (int r = 1; r < this->_ny - 1; ++r) {
			// Getting the node index
			int node = r;
			if (boundaries[node] != BC::NO_FLOW)
				if (boundaries[node] != BC::PERIODIC_BORDER)
					neighbours[node] = LK8.Top_normal_boundary();
				else
					neighbours[node] = AllMask8;
		}

		// Managing last row
		// | | | |
		// | | | |
		// | |X| |
		for (int r = this->nxy() - this->_nx + 1; r < this->nxy() - 1; ++r) {
			// Getting the node index
			int node = r;
			if (boundaries[node] != BC::NO_FLOW)
				if (boundaries[node] != BC::PERIODIC_BORDER)
					neighbours[node] = LK8.Bottom_normal_boundary();
				else
					neighbours[node] = AllMask8;
		}

		// Managing Left and right columns
		// | | | |
		// |X| |X|
		// | | | |
		for (int r = 1; r < this->_ny - 1; ++r) {
			// Getting the node index
			int n1 = r * this->_nx + 0;
			if (boundaries[n1] != BC::NO_FLOW) {
				if (boundaries[n1] != BC::PERIODIC_BORDER)
					neighbours[n1] = LK8.Left_normal_boundary();
				else
					neighbours[n1] = AllMask8;
			}
			n1 = r * this->_nx + this->_nx - 1;
			if (boundaries[n1] != BC::NO_FLOW) {
				if (boundaries[n1] != BC::PERIODIC_BORDER)
					neighbours[n1] = LK8.Right_normal_boundary();
				else
					neighbours[n1] = AllMask8;
			}
		}
	}

	void compute()
	{
		if (this->flowtopo == CONFLOWTOPO::NONE)
			return;

		this->reinit();

		if (this->flowtopo == CONFLOWTOPO::ALL) {
			this->_compute_all();
		} else if (this->flowtopo == CONFLOWTOPO::MFD) {
			this->_compute_mfd_only();
		} else if (this->flowtopo == CONFLOWTOPO::SFD) {
			throw std::runtime_error(
				"compute for sfd only not available yet. Use ALL.");
			// this->_compute_sfd_only();
		}
	}

	void _compute_all()
	{

		// Just double checking
		if (this->data->_surface.size() == 0) {
			throw std::runtime_error("NoTopoError: no topography set in Hermes");
		}

		// Setting up context
		CT_neighbourer_1<i_t, f_t> ctx;

		// Setting up prefetchers

		for (int i = 0; i < this->_nxy; ++i) {
			this->__compute_all_single_node(i, ctx);
		}
	}

	void reset_node(i_t i)
	{
		this->data->_Sreceivers[i] = 0;
		this->data->_receivers[i] = 0;
		this->data->_donors[i] = 0;
		this->data->_Sdonors[i] = 0;
	}

	template<class CTX>
	void __compute_all_single_node(i_t i, CTX& ctx)
	{
		ctx.update(i, *this);
		f_t SS = 0;
		std::uint8_t rcode = 0;
		std::uint8_t srecode = 0;
		std::uint8_t dcode = 0;

		for (int j = 0; j < ctx.nn; ++j) {
			f_t dz = ctx.topo - ctx.neighboursTopo[j];

			if (Fcan_connect(ctx.boundary, ctx.neighboursCode[j]) == false)
				continue;

			// First asserting the connectivity
			if (dz < 0)
				dcode |= ctx.neighboursBits[j];

			else if (dz > 0) {
				rcode |= ctx.neighboursBits[j];
				f_t tSS = dz / ctx.neighboursDx[j];
				if (tSS > SS) {
					SS = tSS;
					srecode = ctx.neighboursBits[j];
				}
			}
		}

		this->data->_Sreceivers[ctx.node] = srecode;
		this->data->_receivers[ctx.node] = rcode;
		this->data->_donors[ctx.node] = dcode;
		this->data->_Sdonors[this->Sreceivers(ctx.node)] |= invBits(srecode);
	}

	template<class CTX>
	void __compute_recs_single_node(i_t i, CTX& ctx)
	{
		ctx.update(i, *this);
		f_t SS = 0;
		std::uint8_t rcode = 0;
		std::uint8_t srecode = 0;

		for (int j = 0; j < ctx.nn; ++j) {
			f_t dz = ctx.topo - ctx.neighboursTopo[j];

			if (Fcan_connect(ctx.boundary, ctx.neighboursCode[j]) == false)
				continue;

			if (dz > 0) {
				rcode |= ctx.neighboursBits[j];
				f_t tSS = dz / ctx.neighboursDx[j];
				if (tSS > SS) {
					SS = tSS;
					srecode = ctx.neighboursBits[j];
				}
			}
		}

		this->data->_Sreceivers[ctx.node] = srecode;
		this->data->_receivers[ctx.node] = rcode;
	}

	template<class CTX>
	void __compute_recs_single_node_mask(i_t i,
																			 CTX& ctx,
																			 std::vector<std::uint8_t>& xtraMask)
	{
		ctx.update(i, *this);
		f_t SS = 0;
		std::uint8_t rcode = 0;
		std::uint8_t srecode = 0;

		for (int j = 0; j < ctx.nn; ++j) {
			f_t dz = ctx.topo - ctx.neighboursTopo[j];

			if (Fcan_connect(ctx.boundary, ctx.neighboursCode[j]) == false ||
					xtraMask[ctx.neighbours[j]] == false)
				continue;

			if (dz > 0) {
				rcode |= ctx.neighboursBits[j];
				f_t tSS = dz / ctx.neighboursDx[j];
				if (tSS > SS) {
					SS = tSS;
					srecode = ctx.neighboursBits[j];
				}
			}
		}

		this->data->_Sreceivers[ctx.node] = srecode;
		this->data->_receivers[ctx.node] = rcode;
	}

	void __invert_recs_at_node(i_t i,
														 std::array<std::uint8_t, 8>& arr,
														 std::array<i_t, 8>& arr2)
	{
		int nr = this->ReceiversBits(i, arr);
		this->Receivers(i, arr2);

		for (int j = 0; j < nr; ++j) {
			int trec = arr2[j];
			this->data->_donors[trec] |= invBits(arr[j]);
		}
		this->data->_Sdonors[this->Sreceivers(i)] |=
			invBits(this->SreceiversBit(i));
	}

	void __compute_all_single_node(i_t i)
	{
		CT_neighbourer_1<i_t, f_t> ctx;
		ctx.update(i, *this);
		f_t SS = 0;
		std::uint8_t rcode = 0;
		std::uint8_t srecode = 0;
		std::uint8_t dcode = 0;

		for (int j = 0; j < ctx.nn; ++j) {
			f_t dz = ctx.topo - ctx.neighboursTopo[j];

			if (Fcan_connect(ctx.boundary, ctx.neighboursCode[j]) == false)
				continue;

			// First asserting the connectivity
			if (dz < 0)
				dcode |= ctx.neighboursBits[j];

			else if (dz > 0) {
				rcode |= ctx.neighboursBits[j];
				f_t tSS = dz / ctx.neighboursDx[j];
				if (tSS > SS) {
					SS = tSS;
					srecode = ctx.neighboursBits[j];
				}
			}
		}

		this->data->_Sreceivers[ctx.node] = srecode;
		this->data->_receivers[ctx.node] = rcode;
		this->data->_donors[ctx.node] = dcode;
		this->data->_Sdonors[this->Sreceivers(ctx.node)] |= invBits(srecode);
	}

#define PFSIZ 512

	void _compute_all_exp1()
	{

		// Just double checking
		if (this->data->_surface.size() == 0) {
			throw std::runtime_error("NoTopoError: no topography set in Hermes");
		}

		std::array<f_t, PFSIZ + 2> topor1, topor2;
		std::array<f_t, PFSIZ> SS1, SS2;
		std::array<std::uint8_t, PFSIZ> rcodes1, rcodes2, dcodes1, dcodes2,
			srcodes1, srcodes2;
		int NN = 0;

		for (int stcol = 1; stcol < this->_nx - 1; stcol += PFSIZ) {

			// Prefetch first batch of data
			int id1 = stcol;
			// int id2 =  stcol + this->_nx;
			for (int ii = -1; ii < PFSIZ + 1; ++ii) {
				topor1[ii + 1] =
					(stcol < this->_nx - 1) ? this->data->_surface[ii + id1] : 0;
				// topor2[ii+1] = this->data->_surface[ii+id2];
			}

			// init the row/cols
			for (int ii = 0; ii < PFSIZ; ++ii) {
				SS1[ii] = 0;
				SS2[ii] = 0;
				rcodes1[ii] = 0;
				rcodes2[ii] = 0;
				dcodes1[ii] = 0;
				dcodes2[ii] = 0;
				srcodes1[ii] = 0;
				srcodes2[ii] = 0;
			}

			std::vector<f_t>& surf = this->data->_surface;

			for (int row = 1; row < this->_ny - 1; ++row) {
				// prefecthing new line
				int idst = row * this->_nx + stcol;

				for (int ii = -1; ii < PFSIZ + 1; ++ii) {
					topor2[ii + 1] = (stcol < this->_nx - 1) ? surf[ii + idst] : 0;
				}

				for (int i = 0; i < PFSIZ; ++i) {
					// corrected index for topo arrays
					int topi = i + 1;

					// top left
					f_t dz = (topor2[topi] - topor1[topi - 1]);
					f_t tSS = dz / this->_dxy;
					if (dz > 0) {
						rcodes2[i] |= TopLeftMask8;
						if (i > 0)
							dcodes1[i - 1] |= BottomRightMask8;
						if (tSS > SS2[i]) {
							SS2[i] = tSS;
							srcodes2[i] = TopLeftMask8;
						}
					} else if (dz < 0) {
						tSS = std::abs(tSS);
						dcodes2[i] |= TopLeftMask8;
						if (i > 0) {
							rcodes1[i - 1] |= BottomRightMask8;
							if (tSS > SS1[i - 1]) {
								SS1[i - 1] = tSS;
								srcodes1[i - 1] = BottomRightMask8;
							}
						}
					}

					// top
					dz = (topor2[topi] - topor1[topi]);
					tSS = dz / this->_dy;
					if (dz > 0) {
						rcodes2[i] |= TopMask8;
						dcodes1[i] |= BottomMask8;
						if (tSS > SS2[i]) {
							SS2[i] = tSS;
							srcodes2[i] = TopMask8;
						}
					} else if (dz < 0) {
						tSS = std::abs(tSS);
						dcodes2[i] |= TopMask8;
						rcodes1[i] |= BottomMask8;
						if (tSS > SS1[i]) {
							SS1[i] = tSS;
							srcodes1[i] = BottomMask8;
						}
					}

					// top right
					dz = (topor2[topi] - topor1[topi + 1]);
					tSS = dz / this->_dxy;
					if (dz > 0) {
						rcodes2[i] |= TopRightMask8;
						if (i < PFSIZ - 1)
							dcodes1[i + 1] |= BottomLeftMask8;
						if (tSS > SS2[i]) {
							SS2[i] = tSS;
							srcodes2[i] = TopRightMask8;
						}
					} else if (dz < 0) {
						tSS = std::abs(tSS);
						dcodes2[i] |= TopRightMask8;
						if (i < PFSIZ - 1) {
							rcodes1[i + 1] |= BottomLeftMask8;
							if (tSS > SS1[i + 1]) {
								SS1[i + 1] = tSS;
								srcodes1[i + 1] = BottomLeftMask8;
							}
						}
					}

					// left
					dz = (topor2[topi] - topor2[topi - 1]);
					tSS = dz / this->_dx;
					if (dz > 0) {
						rcodes2[i] |= LeftMask8;
						if (i > 0)
							dcodes2[i - 1] |= RightMask8;
						if (tSS > SS2[i]) {
							SS2[i] = tSS;
							srcodes2[i] = LeftMask8;
						}
					} else if (dz < 0) {
						tSS = std::abs(tSS);
						dcodes2[i] |= LeftMask8;
						if (i > 0) {
							rcodes2[i - 1] |= RightMask8;
							if (tSS > SS2[i]) {
								SS2[i - 1] = tSS;
								srcodes2[i - 1] = RightMask8;
							}
						}
					}

					// right
					dz = (topor2[topi] - topor2[topi + 1]);
					tSS = dz / this->_dx;
					if (dz > 0) {
						rcodes2[i] |= RightMask8;
						if (i < PFSIZ - 1)
							dcodes2[i + 1] |= LeftMask8;
						if (tSS > SS2[i]) {
							SS2[i] = tSS;
							srcodes2[i] = RightMask8;
						}
					} else if (dz < 0) {
						tSS = std::abs(tSS);
						dcodes2[i] |= RightMask8;
						if (i < PFSIZ - 1) {
							rcodes2[i + 1] |= LeftMask8;
							if (tSS > SS2[i]) {
								SS2[i + 1] = tSS;
								srcodes2[i + 1] = LeftMask8;
							}
						}
					}
				}

				// propagates changes of line 1 (if not first row)
				if (row > 1) {
					for (int ii = 0; ii < PFSIZ; ++ii) {
						if (ii > 0) {
							this->data->_donors[ii + idst] = dcodes1[ii];
							this->data->_receivers[ii + idst] = rcodes1[ii];
							this->data->_Sreceivers[ii + idst] = srcodes1[ii];
						}
					}
				}

				// move row2 into row1
				topor1 = std::move(topor2);
				rcodes1 = std::move(rcodes2);
				dcodes1 = std::move(dcodes2);
				srcodes1 = std::move(srcodes2);
				SS1 = std::move(SS2);

				for (int ii = 0; ii < PFSIZ; ++ii) {
					rcodes2[ii] = 0;
					dcodes2[ii] = 0;
					srcodes2[ii] = 0;
					SS2[ii] = 0;
				}
			}
		}

		// Compute boundaries with normal function
		CT_neighbourer_1<i_t, f_t> ctx;
		for (int i = 0; i < this->_nx; ++i)
			__compute_all_single_node(i, ctx);
		for (int i = this->_nxy - this->_nx; i < this->_nxy; ++i)
			__compute_all_single_node(i, ctx);
		for (int row = 1; row < this->_ny - 1; ++row) {
			int i = row * this->_nx;
			__compute_all_single_node(i, ctx);
			i = (row + 1) * this->_nx - 1;
			__compute_all_single_node(i, ctx);
		}

		for (int i = 0; i < this->_nxy; ++i) {
			int rec = this->Sreceivers(i);
			this->data->_Sdonors[rec] |=
				invBits(this->data->_Sreceivers[i]) * !(rec == i);
		}
	}

	void _compute_all_exp2()
	{

		// Just double checking
		if (this->data->_surface.size() == 0) {
			throw std::runtime_error("NoTopoError: no topography set in Hermes");
		}

		std::vector<f_t> SS(this->_nxy, 0.);

		for (int row = 1; row < this->_ny - 1; ++row) {
			for (int col = 1; col < this->_nx; ++col) {

				int idx = row * this->_nx + col;

				// top left
				int oidx = idx - this->_nx - 1;
				f_t dz = (this->data->_surface[idx] - this->data->_surface[oidx]);
				f_t tSS = dz / this->_dxy;
				if (dz > 0) {
					this->data->_receivers[idx] |= TopLeftMask8;
					this->data->_donors[oidx] |= BottomRightMask8;
					if (tSS > SS[idx]) {
						SS[idx] = tSS;
						this->data->_Sreceivers[idx] = TopLeftMask8;
					}
				} else if (dz < 0) {
					tSS = std::abs(tSS);
					this->data->_donors[idx] |= TopLeftMask8;

					this->data->_receivers[oidx] |= BottomRightMask8;
					if (tSS > SS[oidx]) {
						SS[oidx] = tSS;
						this->data->_Sreceivers[oidx] = BottomRightMask8;
					}
				}

				// top
				oidx += 1;
				dz = (this->data->_surface[idx] - this->data->_surface[oidx]);
				tSS = dz / this->_dy;
				if (dz > 0) {
					this->data->_receivers[idx] |= TopMask8;
					this->data->_donors[oidx] |= BottomMask8;
					if (tSS > SS[idx]) {
						SS[idx] = tSS;
						this->data->_Sreceivers[idx] = TopMask8;
					}
				} else if (dz < 0) {
					tSS = std::abs(tSS);
					this->data->_donors[idx] |= TopMask8;
					this->data->_receivers[oidx] |= BottomMask8;
					if (tSS > SS[oidx]) {
						SS[oidx] = tSS;
						this->data->_Sreceivers[oidx] = BottomMask8;
					}
				}

				// top right
				oidx += 1;
				dz = (this->data->_surface[idx] - this->data->_surface[oidx]);
				tSS = dz / this->_dxy;
				if (dz > 0) {
					this->data->_receivers[idx] |= TopRightMask8;
					this->data->_donors[oidx] |= BottomLeftMask8;
					if (tSS > SS[idx]) {
						SS[idx] = tSS;
						this->data->_Sreceivers[idx] = TopRightMask8;
					}
				} else if (dz < 0) {
					tSS = std::abs(tSS);
					this->data->_donors[idx] |= TopRightMask8;
					this->data->_receivers[oidx] |= BottomLeftMask8;
					if (tSS > SS[oidx]) {
						SS[oidx] = tSS;
						this->data->_Sreceivers[oidx] = BottomLeftMask8;
					}
				}

				// left
				oidx = idx - 1;
				dz = (this->data->_surface[idx] - this->data->_surface[oidx]);
				tSS = dz / this->_dx;
				if (dz > 0) {
					this->data->_receivers[idx] |= LeftMask8;
					this->data->_donors[oidx] |= RightMask8;
					if (tSS > SS[idx]) {
						SS[idx] = tSS;
						this->data->_Sreceivers[idx] = LeftMask8;
					}
				} else if (dz < 0) {
					tSS = std::abs(tSS);
					this->data->_donors[idx] |= LeftMask8;
					this->data->_receivers[oidx] |= RightMask8;
					if (tSS > SS[idx]) {
						SS[oidx] = tSS;
						this->data->_Sreceivers[oidx] = RightMask8;
					}
				}

				// right
				oidx = idx + 1;
				dz = (this->data->_surface[idx] - this->data->_surface[oidx]);
				tSS = dz / this->_dx;
				if (dz > 0) {
					this->data->_receivers[idx] |= RightMask8;
					this->data->_donors[oidx] |= LeftMask8;
					if (tSS > SS[idx]) {
						SS[idx] = tSS;
						this->data->_Sreceivers[idx] = RightMask8;
					}
				} else if (dz < 0) {
					tSS = std::abs(tSS);
					this->data->_donors[idx] |= RightMask8;

					this->data->_receivers[oidx] |= LeftMask8;
					if (tSS > SS[oidx]) {
						SS[oidx] = tSS;
						this->data->_Sreceivers[oidx] = LeftMask8;
					}
				}
			}
		}

		// Compute boundaries with normal function
		CT_neighbourer_1<i_t, f_t> ctx;
		for (int i = 0; i < this->_nx; ++i)
			__compute_all_single_node(i, ctx);
		for (int i = this->_nxy - this->_nx; i < this->_nxy; ++i)
			__compute_all_single_node(i, ctx);
		for (int row = 1; row < this->_ny - 1; ++row) {
			int i = row * this->_nx;
			__compute_all_single_node(i, ctx);
			i = (row + 1) * this->_nx - 1;
			__compute_all_single_node(i, ctx);
		}

		for (int i = 0; i < this->_nxy; ++i) {
			int rec = this->Sreceivers(i);
			this->data->_Sdonors[rec] |=
				invBits(this->data->_Sreceivers[i]) * !(rec == i);
		}
	}

	void _compute_mfd_only()
	{

		// Just double checking
		if (this->data->_surface.size() == 0) {
			throw std::runtime_error("NoTopoError: no topography set in Hermes");
		}

		// Setting up context
		CT_neighbourer_1<i_t, f_t> ctx;

		// Setting up prefetchers

		for (int i = 0; i < this->_nxy; ++i) {
			ctx.update(i, *this);
			std::uint8_t rcode = 0;
			std::uint8_t dcode = 0;

			for (int j = 0; j < ctx.nn; ++j) {
				f_t dz = ctx.topo - ctx.neighboursTopo[j];

				if (Fcan_connect(ctx.boundary, ctx.neighboursCode[j]) == false)
					continue;

				// First asserting the connectivity
				if (dz < 0)
					dcode |= ctx.neighboursBits[j];
				else if (dz > 0) {
					rcode |= ctx.neighboursBits[j];
				}
			}
			this->data->_receivers[ctx.node] = rcode;
			this->data->_donors[ctx.node] = dcode;
		}
	}

	// SLOWER, keep for legacy
	void _compute_mfd_only_exp1()
	{

		// Just double checking
		if (this->data->_surface.size() == 0) {
			throw std::runtime_error("NoTopoError: no topography set in Hermes");
		}

		// Setting up context
		CT_neighbourer_256<i_t, f_t> ctx;

		// Setting up prefetchers

		for (int i = 0; i < this->_nxy - 256; i += 256) {
			ctx.update(i, *this);
			for (int subi = 0; subi < 256; ++subi) {
				std::uint8_t rcode = 0;
				std::uint8_t dcode = 0;

				for (int j = 0; j < ctx.nn[subi]; ++j) {
					f_t dz = ctx.topo[subi] - ctx.neighboursTopo[subi][j];

					if (Fcan_connect(ctx.boundary[subi], ctx.neighboursCode[subi][j]) ==
							false)
						continue;

					// First asserting the connectivity
					if (dz < 0)
						dcode |= ctx.neighboursBits[subi][j];
					else if (dz > 0) {
						rcode |= ctx.neighboursBits[subi][j];
					}
				}
				this->data->_receivers[ctx.node[subi]] = rcode;
				this->data->_donors[ctx.node[subi]] = dcode;
			}
		}
	}

	void _quickSstack()
	{
		// The stack container helper
		this->data->_Sstack = std::vector<int>(this->_nxy, 0);

		std::stack<int, std::vector<int>> stackhelper;
		std::array<i_t, 8> tsdon;
		// std::vector<bool> isdone(this->nnodes,false);
		// going through all the nodes
		int istack = 0;
		for (int i = 0; i < this->_nxy; ++i) {
			// if they are base level I include them in the stack
			if (this->Sreceivers(i) == i) {
				stackhelper.emplace(i);
				// ++istack;
			}

			// While I still have stuff in the stack helper
			while (stackhelper.empty() == false) {
				// I get the next node and pop it from the stack helper
				int nextnode = stackhelper.top();
				stackhelper.pop();
				// std::cout << istack << "/" << this->_nxy << std::endl;
				this->data->_Sstack[istack] = nextnode;
				++istack;

				// as well as all its donors which will be processed next
				i_t nn = this->Sdonors(nextnode, tsdon);
				for (int j = 0; j < nn; ++j) {
					stackhelper.emplace(tsdon[j]);
				}
			}
		}
	}

	void _quickLM()
	{

		// 0: not done
		// 1: connected to the ocean
		// 2: internal pit
		std::vector<std::uint8_t> connected(this->_nxy, 0);

		std::stack<i_t, std::vector<i_t>> helper;

		for (int i = 0; i < this->_nxy; ++i) {
			if (can_out(this->data->_boundaries[i])) {
				helper.emplace(i);
				connected[i] = 1;
			}
		}

		std::array<i_t, 8> tdons;
		while (helper.empty() == false) {
			int next = helper.top();
			helper.pop();
			int nn = this->Donors(next, tdons);
			for (int j = 0; j < nn; ++j) {
				if (connected[tdons[j]] == 0) {
					helper.emplace(tdons[j]);
					connected[tdons[j]] = 1;
				}
			}
		}

		// ocean labelled
		for (int i = 0; i < this->_nxy; ++i) {
			if (connected[i] == 0)
				continue;
			int nn = this->Receivers(i, tdons);
			for (int j = 0; j < nn; ++j) {
				int reci = tdons[j];
				if (connected[reci] == 0 && connected[i] == 1) {
					helper.emplace(i);
				}
			}
		}

		while (helper.empty() == false) {
			int next = helper.top();
			helper.pop();
			int nn = this->Receivers(next, tdons);
			for (int j = 0; j < nn; ++j) {
				int rec = tdons[j];
				if (connected[rec] == 0) {
					this->data->_surface[rec] =
						this->data->_surface[next] + 1e-6 * this->data->randu->get();
					connected[rec] = 1;
					helper.emplace(rec);
				}
			}
		}
	}

	void PFcompute_all(bool affect_topo = true)
	{
		std::vector<f_t> oldtopo;
		if (affect_topo == false) {
			oldtopo = this->data->_surface;
		}

		PQFORPF open;
		std::queue<PQH> pit;
		this->data->_stack = std::vector<i_t>(this->_nxy, 0);

		// Setting up context
		CT_neighbourer_1<i_t, f_t> ctx;

		std::vector<int8_t> closed(this->_nxy, false);
		int istack = 0;

		for (int i = 0; i < this->_nxy; ++i) {

			this->reset_node(i);
			if (nodata(this->data->_boundaries[i])) {
				this->data->_stack[istack] = i;
				closed[i] = true;
				++istack;
				continue;
			}

			if (can_out(this->data->_boundaries[i]) == false)
				continue;

			open.emplace(PQH(i, this->data->_surface[i]));
			closed[i] = true;
		}

		std::array<int, 8> neighbours;
		while (open.size() > 0 || pit.size() > 0) {
			PQH c;
			c = open.top();
			open.pop();

			// this->__compute_all_single_node(c.node, ctx);
			int nn = this->Neighbours(c.node, neighbours);
			this->data->_stack[istack] = c.node;
			++istack;

			f_t ttopo = this->data->_surface[c.node];

			// for (int j = 0; j < ctx.nn; ++j) {
			for (int j = 0; j < nn; ++j) {
				// int n = ctx.neighbours[j];
				int n = neighbours[j];

				if (closed[n])
					continue;

				closed[n] = true;
				ttopo += this->data->randu->get() * 1e-6 + 1e-8;
				this->data->_surface[n] = std::max(ttopo, this->data->_surface[n]);
				open.emplace(n, this->data->_surface[n]);
			}
			this->__compute_all_single_node(c.node, ctx);
		}
		this->_quickSstack();

		if (affect_topo == false) {
			this->data->_surface = oldtopo;
		}
	}
};

} // End of Namespace
