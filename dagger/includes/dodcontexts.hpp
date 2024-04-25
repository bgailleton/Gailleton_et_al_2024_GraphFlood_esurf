#pragma once
#include "boundary_conditions.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class CT_neighbourer_1
{

public:
	i_t node = 0;
	f_t topo;
	BC boundary;

	i_t nn = 0;
	std::uint8_t idAdder = 0;
	std::array<i_t, 8> neighbours;
	std::array<std::uint8_t, 8> neighboursBits;
	std::array<BC, 8> neighboursCode;
	std::array<f_t, 8> neighboursDx;
	std::array<f_t, 8> neighboursTopo;

	template<class CONNECTOR_T>
	void update(i_t node, CONNECTOR_T& con)
	{

		this->node = node;

		this->topo = con.data->_surface[node];
		this->boundary = con.data->_boundaries[node];

		std::uint8_t nid = con.data->_neighbours[node];

		this->idAdder = con.data->LK8.BC2idAdder(node, con.data->_boundaries[node]);

		this->nn = con.data->LK8.NeighbourerNN[this->idAdder][nid];

		this->neighbours = con.data->LK8.Neighbourer[this->idAdder][nid];

		this->neighboursDx = con.data->LK8.Neighbourerdx[this->idAdder][nid];
		this->neighboursBits = con.data->LK8.NeighbourerBits[this->idAdder][nid];

		for (size_t i = 0; i < this->nn; ++i)
			this->neighbours[i] += node;

		for (size_t i = 0; i < nn; ++i)
			this->neighboursTopo[i] = con.data->_surface[this->neighbours[i]];

		for (size_t i = 0; i < nn; ++i)
			neighboursCode[i] = con.data->_boundaries[this->neighbours[i]];
	}
};

template<class i_t, class f_t>
class CT_neighbourer_256
{

public:
	std::array<i_t, 256> node;
	std::array<f_t, 256> topo;
	std::array<BC, 256> boundary;

	std::array<i_t, 256> nn;
	std::array<std::uint8_t, 256> idAdder;
	std::array<std::array<i_t, 8>, 256> neighbours;
	std::array<std::array<std::uint8_t, 8>, 256> neighboursBits;
	std::array<std::array<BC, 8>, 256> neighboursCode;
	std::array<std::array<f_t, 8>, 256> neighboursDx;
	std::array<std::array<f_t, 8>, 256> neighboursTopo;

	template<class CONNECTOR_T>
	void update(i_t node, CONNECTOR_T& con)
	{
		for (int incr = 0; incr < 256; ++incr) {
			this->node[incr] = node + incr;

			this->topo[incr] = con.data->_surface[this->node[incr]];
			this->boundary[incr] = con.data->_boundaries[this->node[incr]];

			std::uint8_t nid = con.data->_neighbours[this->node[incr]];

			this->idAdder[incr] = con.data->LK8.BC2idAdder(
				this->node[incr], con.data->_boundaries[this->node[incr]]);

			this->nn[incr] = con.data->LK8.NeighbourerNN[this->idAdder[incr]][nid];

			this->neighbours[incr] =
				con.data->LK8.Neighbourer[this->idAdder[incr]][nid];

			this->neighboursDx[incr] =
				con.data->LK8.Neighbourerdx[this->idAdder[incr]][nid];
			this->neighboursBits[incr] =
				con.data->LK8.NeighbourerBits[this->idAdder[incr]][nid];

			for (size_t i = 0; i < this->nn[incr]; ++i)
				this->neighbours[incr][i] += this->node[incr];

			for (size_t i = 0; i < nn[incr]; ++i)
				this->neighboursTopo[incr][i] =
					con.data->_surface[this->neighbours[incr][i]];

			for (size_t i = 0; i < nn[incr]; ++i)
				neighboursCode[incr][i] =
					con.data->_boundaries[this->neighbours[incr][i]];
		}
	}
};

template<class i_t, class f_t>
class CT_neighbourer_WaCell
{

public:
	i_t node = 0;
	BC boundary;
	i_t nn = 0;
	i_t nr = 0;
	f_t surface = 0;

	f_t sumslopesdw = 0.;
	i_t SSj = 0;

	std::array<i_t, 8> neighbours;
	std::array<BC, 8> neighboursCode;
	std::array<f_t, 8> neighboursDx;
	std::array<std::int8_t, 8> neighboursDxSign;
	std::array<std::int8_t, 8> neighboursDySign;
	std::array<f_t, 8> neighboursDy;

	std::array<i_t, 8> receivers;
	std::array<BC, 8> receiversCode;
	std::array<f_t, 8> receiversDx;
	std::array<f_t, 8> receiversDy;
	std::array<f_t, 8> receiversSurfaces;
	std::array<f_t, 8> receiversSlopes;
	std::array<f_t, 8> receiversWeights;
	std::array<std::int8_t, 8> receiversDxSign, receiversDySign;

	bool allRecsDone = false;
	bool allNeighsDone = false;
	bool canout = false;

	template<class CONNECTOR_T>
	void update(i_t node, CONNECTOR_T& con)
	{

		// updating to the current node
		this->node = node;
		// fetching the current surface
		this->surface = con.data->_surface[node];
		// node boundary code
		this->boundary = con.data->_boundaries[node];

		// Fetching neighbour infos
		this->nn = con.Neighbours(node, this->neighbours);
		con.NeighboursDx(node, this->neighboursDx);
		con.NeighboursDy(node, this->neighboursDy);
		con.NeighboursDxSign(node, this->neighboursDxSign);
		con.NeighboursDySign(node, this->neighboursDySign);

		// Now dealing with receivers
		this->nr = 0;
		this->receiversSlopes[0] = 0.;
		this->SSj = 0;

		// Checks if all rec/ neighbours have been processed
		this->allRecsDone = true;
		this->allNeighsDone = true;
		this->sumslopesdw = 0;
		this->canout = false;

		// Processing neighbours
		for (size_t i = 0; i < this->nn; ++i) {
			// Assing boundary code
			this->neighboursCode[i] = con.data->_boundaries[this->neighbours[i]];

			if (can_out(this->neighboursCode[i]))
				this->canout = true;

			// checking if node is at the right timing or not
			if (con.data->_timetracker[node] !=
					con.data->_timetracker[this->neighbours[i]])
				this->allNeighsDone = false;

			// Checking if the receivers is done
			if (con.data->_surface[this->neighbours[i]] < this->surface &&
					can_receive(this->neighboursCode[i])) {
				if (con.data->_timetracker[node] !=
						con.data->_timetracker[this->neighbours[i]])
					this->allRecsDone = false;
			}
		}

		// if (this->allRecsDone == false) {
		if (true) {
			bool first = true;
			for (size_t i = 0; i < this->nn; ++i) {
				if (con.data->_surface[this->neighbours[i]] < this->surface &&
						can_receive(
							this->neighboursCode
								[i])) { // &&con.data->_timetracker[node] !=
												// con.data->_timetracker[this->neighbours[i]])
					this->receivers[nr] = this->neighbours[i];
					this->receiversCode[nr] = this->neighboursCode[i];
					this->receiversDx[nr] = this->neighboursDx[i];
					this->receiversDy[nr] = this->neighboursDy[i];
					this->receiversDxSign[nr] = this->neighboursDxSign[i];
					this->receiversDySign[nr] = this->neighboursDySign[i];
					this->receiversSurfaces[nr] = con.data->_surface[this->neighbours[i]];
					this->receiversSlopes[nr] =
						(this->surface - this->receiversSurfaces[nr]) /
						this->receiversDx[nr];
					if (this->receiversSlopes[nr] > this->receiversSlopes[this->SSj] ||
							first) {
						this->SSj = nr;
						first = false;
					}

					this->receiversWeights[nr] =
						this->receiversSlopes[nr] * this->receiversDy[nr];
					this->sumslopesdw += this->receiversWeights[nr];

					++nr;
				}
			}

			for (size_t i = 0; i < this->nr; ++i)
				// this->receiversWeights[i] /= this->sumslopesdw;
				this->receiversWeights[i] = 1. / static_cast<double>(nr);

		} else {
			for (size_t i = 0; i < this->nn; ++i) {
				if (con.data->_surface[this->neighbours[i]] < this->surface &&
						can_receive(this->neighboursCode[i])) {
					this->receivers[nr] = i;
					++nr;
				}
			}

			if (nr == 0)
				return;
			// std::cout << nr << "|";

			this->SSj = 0;
			int tj = this->receivers[std::floor(con.data->randu->get() * this->nr)];
			this->receivers[0] = this->neighbours[tj];
			this->receiversCode[0] = this->neighboursCode[tj];
			this->receiversDx[0] = this->neighboursDx[tj];
			this->receiversDy[0] = this->neighboursDy[tj];
			this->receiversSurfaces[0] = con.data->_surface[this->neighbours[tj]];
			this->receiversSlopes[0] =
				(this->surface - this->receiversSurfaces[0]) / this->neighboursDx[tj];
			this->sumslopesdw = this->receiversSlopes[0] * this->neighboursDy[tj];
			nr = 1;
			this->receiversWeights[0] = 1.;
		}
	}
};

template<class i_t, class f_t>
class CT_neighbours
{

public:
	i_t node = 0;
	BC boundary;
	i_t nn = 0;

	std::array<i_t, 8> neighbours;
	std::array<BC, 8> neighboursCode;
	std::array<f_t, 8> neighboursDx;
	std::array<f_t, 8> neighboursDy;
	std::array<std::int8_t, 8> neighboursDxSign;
	std::array<std::int8_t, 8> neighboursDySign;

	bool canout = false;

	template<class CONNECTOR_T>
	void update(i_t node, CONNECTOR_T& con)
	{

		// updating to the current node
		this->node = node;
		// node boundary code
		this->boundary = con.data->_boundaries[node];

		// Fetching neighbour infos
		this->nn = con.Neighbours(node, this->neighbours);
		con.NeighboursDx(node, this->neighboursDx);
		con.NeighboursDy(node, this->neighboursDy);
		con.NeighboursDxSign(node, this->neighboursDxSign);
		con.NeighboursDySign(node, this->neighboursDySign);
		this->canout = false;

		// Processing neighbours
		for (size_t i = 0; i < this->nn; ++i) {
			// Assing boundary code
			this->neighboursCode[i] = con.data->_boundaries[this->neighbours[i]];

			if (can_out(this->neighboursCode[i]))
				this->canout = true;
		}
	}

	template<class ARRT>
	int idxHighestNeighbour(ARRT& vecin)
	{
		f_t Zmax = vecin[this->node];
		int imax = -1;
		for (int ii = 0; ii < this->nn; ++ii) {
			if (vecin[this->neighbours[ii]] > Zmax) {
				imax = ii;
				Zmax = vecin[this->neighbours[ii]];
			}
		}
		return imax;
	}

	template<class ARRT>
	int idxLowestNeighbour(ARRT& vecin)
	{
		f_t Zmax = vecin[this->node];
		int imax = -1;
		for (int ii = 0; ii < this->nn; ++ii) {
			if (vecin[this->neighbours[ii]] < Zmax) {
				imax = ii;
				Zmax = vecin[this->neighbours[ii]];
			}
		}
		return imax;
	}

	template<class ARRT>
	int idxHighestNeighbour(ARRT& vecadd, ARRT& vecmin)
	{
		f_t Zmax = vecadd[this->node] - vecmin[this->node];
		int imax = -1;
		for (int ii = 0; ii < this->nn; ++ii) {
			if ((vecadd[this->neighbours[ii]] - vecmin[this->neighbours[ii]]) >
					Zmax) {
				imax = ii;
				Zmax = vecadd[this->neighbours[ii]] - vecmin[this->neighbours[ii]];
			}
		}
		return imax;
	}

	template<class ARRT>
	int idxLowestNeighbour(ARRT& vecadd, ARRT& vecmin)
	{
		f_t Zmax = vecadd[this->node] - vecmin[this->node];
		int imax = -1;
		for (int ii = 0; ii < this->nn; ++ii) {
			if ((vecadd[this->neighbours[ii]] - vecmin[this->neighbours[ii]]) <
					Zmax) {
				imax = ii;
				Zmax = vecadd[this->neighbours[ii]] - vecmin[this->neighbours[ii]];
			}
		}
		return imax;
	}

	template<class ARRT>
	int idxHighestNeighbour(ARRT& vecadd, ARRT& vecmin, int notAnOption)
	{
		f_t Zmax = vecadd[this->node] - vecmin[this->node];
		int imax = -1;
		for (int ii = 0; ii < this->nn; ++ii) {
			if (this->neighbours[ii] == notAnOption)
				continue;
			if ((vecadd[this->neighbours[ii]] - vecmin[this->neighbours[ii]]) >
					Zmax) {
				imax = ii;
				Zmax = vecadd[this->neighbours[ii]] - vecmin[this->neighbours[ii]];
			}
		}
		return imax;
	}

	template<class ARRT>
	int idxLowestNeighbour(ARRT& vecadd, ARRT& vecmin, int notAnOption)
	{
		f_t Zmax = vecadd[this->node] - vecmin[this->node];
		int imax = -1;
		for (int ii = 0; ii < this->nn; ++ii) {
			if (this->neighbours[ii] == notAnOption)
				continue;
			if ((vecadd[this->neighbours[ii]] - vecmin[this->neighbours[ii]]) <
					Zmax) {
				imax = ii;
				Zmax = vecadd[this->neighbours[ii]] - vecmin[this->neighbours[ii]];
			}
		}
		return imax;
	}
};

} // end of namespace
