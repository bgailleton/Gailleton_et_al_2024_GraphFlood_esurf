#pragma once

#include "databag.hpp"
#include "dodconnector.hpp"
#include "parambag.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, class f_t, class CONNECTOR_T>
std::vector<f_t>
_compute_min_distance_from_outlets(CONNECTOR_T& con)
{

	// preformatting the output
	std::vector<f_t> out(con.nxy(), 0.);

	std::array<f_t, 8> dxs;
	std::array<i_t, 8> recs;

	for (auto node : con.data->_stack) {

		if (nodata(con.data->_boundaries[node]))
			continue;

		int nd = con.Receivers(node, recs);

		con.ReceiversDx(node, dxs);

		for (int j = 0; j < nd; ++j) {
			out[node] = (out[node] == 0) ? out[recs[j]] + dxs[j]
																	 : std::min(out[node], out[recs[j]] + dxs[j]);
		}
	}

	return out;
}

template<class i_t, class f_t, class CONNECTOR_T>
std::vector<f_t>
_compute_max_distance_from_outlets(CONNECTOR_T& con)
{

	// preformatting the output
	std::vector<f_t> out(con.nxy(), 0.);

	std::array<f_t, 8> dxs;
	std::array<i_t, 8> recs;

	for (auto node : con.data->_stack) {

		if (nodata(con.data->_boundaries[node]))
			continue;

		int nd = con.Receivers(node, recs);

		con.ReceiversDx(node, dxs);

		for (int j = 0; j < nd; ++j) {
			out[node] = (out[node] == 0) ? out[recs[j]] + dxs[j]
																	 : std::max(out[node], out[recs[j]] + dxs[j]);
		}
	}

	return out;
}

template<class i_t, class f_t, class CONNECTOR_T>
std::vector<f_t>
_compute_SFD_distance_from_outlets(CONNECTOR_T& con)
{

	// preformatting the output
	std::vector<f_t> out(con.nxy(), 0.);

	for (auto node : con.data->_stack) {

		if (nodata(con.data->_boundaries[node]))
			continue;

		int rec = con.Sreceivers(node);
		f_t recdx = con.SreceiversDx(node);

		if (rec != node) {
			out[node] = out[rec] + recdx;
		}
	}

	return out;
}

template<class i_t, class f_t, class CONNECTOR_T>
void
compute_min_distance_from_outlets(bool compute_PQ,
																	std::string name,
																	CONNECTOR_T& con)
{

	if (compute_PQ) {
		con.PFcompute_all(false);
	}

	std::vector<f_t> dx2out =
		_compute_min_distance_from_outlets<i_t, f_t, CONNECTOR_T>(con);

	con.data->fbag[name] = dx2out;
}

template<class i_t, class f_t, class CONNECTOR_T>
void
compute_max_distance_from_outlets(bool compute_PQ,
																	std::string name,
																	CONNECTOR_T& con)
{

	if (compute_PQ) {
		con.PFcompute_all(false);
	}

	std::vector<f_t> dx2out =
		_compute_max_distance_from_outlets<i_t, f_t, CONNECTOR_T>(con);

	con.data->fbag[name] = dx2out;
}

template<class i_t, class f_t, class CONNECTOR_T>
void
fill_lake_in_hw(CONNECTOR_T& con)
{
	if (con.data->_hw.size() == 0) {
		con.data->_hw = std::vector<f_t>(con.nxy(), 0.);
	}
	std::vector<f_t> oldZ = con.data->_surface;
	con.PFcompute_all(true);
	for (int i = 0; i < con.nxy(); ++i) {
		con.data->_hw[i] += con.data->_surface[i] - oldZ[i];
	}
}

template<class i_t, class f_t, class CONNECTOR_T>
void
fill_LM(CONNECTOR_T& con)
{
	con.PFcompute_all(true);
}

template<class i_t, class f_t, class CONNECTOR_T>
std::vector<i_t>
_compute_SFD_basin_labels(CONNECTOR_T& con)
{

	std::vector<i_t> basID(con.nxy(), -1);
	int lab = -1;
	for (int i = 0; i < con.nxy(); ++i) {

		int node = con.data->_Sstack[i];

		if (nodata(con.data->_boundaries[node]))
			continue;

		int rec = con.Sreceivers(node);
		if (rec == node) {
			++lab;
		}

		basID[node] = lab;
	}

	return basID;
}

template<class i_t, class f_t, class CONNECTOR_T>
void
compute_SFD_basin_labels(bool compute_PQ, std::string name, CONNECTOR_T& con)
{

	if (compute_PQ) {
		con.PFcompute_all(false);
	}

	std::vector<i_t> out = _compute_SFD_basin_labels<i_t, f_t, CONNECTOR_T>(con);

	con.data->ibag[name] = out;
}

template<class i_t, class f_t, class CONNECTOR_T>
void
compute_SFD_DA(bool compute_PQ, CONNECTOR_T& con)
{

	if (compute_PQ) {
		con.PFcompute_all(false);
	}

	con.data->_DA = std::vector<f_t>(con.nxy(), 0);

	for (int i = con.nxy() - 1; i >= 0; --i) {
		int node = con.data->_Sstack[i];

		if (nodata(con.data->_boundaries[node]))
			continue;

		con.data->_DA[node] += con.area(node);

		int rec = con.Sreceivers(node);

		if (rec == node)
			continue;

		con.data->_DA[rec] += con.data->_DA[node];
	}
}

template<class i_t, class f_t, class CONNECTOR_T>
void
recast_BC_bellow_Z(f_t minelev, CONNECTOR_T& connector)
{

	for (int i = 0; i < connector.nxy(); ++i) {
		if (connector.data->_surface[i] < minelev)
			connector.data->_boundaries[i] = static_cast<BC>(0);
		else
			connector.data->_boundaries[i] = static_cast<BC>(1);
	}

	std::array<i_t, 8> neighs;
	for (int i = 0; i < connector.nxy(); ++i) {
		int nn = connector.Neighbours(i, neighs);
		bool isborder = false;
		for (int j = 0; j < nn; ++j) {
			int ni = neighs[j];
			if (connector.data->_boundaries[ni] == static_cast<BC>(0) &&
					connector.data->_boundaries[i] == static_cast<BC>(1))
				connector.data->_boundaries[i] = static_cast<BC>(3);
		}
	}

	for (int i = 0; i < connector._nx; ++i) {
		if (connector.data->_boundaries[i] == static_cast<BC>(1))
			connector.data->_boundaries[i] = static_cast<BC>(3);
	}

	for (int i = connector.nxy() - connector._nx; i < connector.nxy(); ++i) {
		if (connector.data->_boundaries[i] == static_cast<BC>(1))
			connector.data->_boundaries[i] = static_cast<BC>(3);
	}

	for (int i = 0; i < connector._ny; ++i) {

		if (connector.data->_boundaries[i * connector._nx] == static_cast<BC>(1))
			connector.data->_boundaries[i * connector._nx] = static_cast<BC>(3);
		if (connector.data->_boundaries[(i + 1) * connector._nx - 1] ==
				static_cast<BC>(1))
			connector.data->_boundaries[(i + 1) * connector._nx - 1] =
				static_cast<BC>(3);
	}

	connector.set_condou(CONBOU::CUSTOM);
	connector.init();
}

template<class i_t, class f_t, class CONNECTOR_T>
void
recast_BC_from_outlet(i_t baseNode,
											CONNECTOR_T& connector,
											bool compute_PQ,
											f_t dZ)
{

	if (compute_PQ) {
		connector.PFcompute_all(false);
	}

	std::vector<std::uint8_t> isIn(connector.nxy(), false);
	isIn[baseNode] = true;

	for (int i = 0; i < connector.nxy(); ++i) {
		int node = connector.data->_Sstack[i];

		if (nodata(connector.data->_boundaries[node]) || node == baseNode)
			continue;

		int rec = connector.Sreceivers(node);
		isIn[node] = isIn[rec];
	}

	for (int i = 0; i < connector.nxy(); ++i) {
		if (isIn[i])
			connector.data->_boundaries[i] = BC::FLOW;
		else
			connector.data->_boundaries[i] = BC::NO_FLOW;
	}

	connector.data->_boundaries[baseNode] = BC::CAN_OUT;

	connector.data->_surface[baseNode] += dZ;

	connector.set_condou(CONBOU::CUSTOM);
	connector.init();
}

template<class i_t, class f_t, class CONNECTOR_T>
std::vector<f_t>
_compute_relief(f_t radius, CONNECTOR_T& connector)
{

	std::vector<f_t> relief(connector.nxy(), 0.);
	std::array<i_t, 8> neighbours;
	std::array<f_t, 8> neighboursdx;
	std::queue<i_t> tQ;

	for (int i = 0; i < connector.nxy(); ++i) {

		if (nodata(connector.data->_boundaries[i]))
			continue;

		std::unordered_map<int, f_t> checks;
		checks[i] = 0.;
		tQ.emplace(i);

		f_t tZ = connector.data->_surface[i];
		f_t tR = 0.;

		while (tQ.empty() == false) {
			int node = tQ.front();
			tQ.pop();
			f_t ttR = connector.data->_surface[node] - tZ;
			ttR = std::abs(ttR);

			if (ttR > tR)
				tR = ttR;

			int nn = connector.Neighbours(node, neighbours);
			connector.NeighboursDx(node, neighboursdx);
			f_t tdist = checks[node];
			for (int j = 0; j < nn; ++j) {
				int onode = neighbours[j];
				f_t odx = neighboursdx[j];
				if (checks.find(onode) == checks.end() && odx + tdist <= radius) {
					checks[onode] = odx + tdist;
					tQ.emplace(onode);
				}
			}
		}

		relief[i] = tR;
	}

	return relief;
}

template<class i_t, class f_t, class CONNECTOR_T>
void
compute_relief(f_t radius, std::string name, CONNECTOR_T& con)
{

	std::vector<f_t> relief = _compute_relief<i_t, f_t, CONNECTOR_T>(radius, con);

	con.data->fbag[name] = relief;
}

/// Computes the elevation above the neart connected point on a mask
template<class i_t, class f_t, class CONNECTOR_T>
std::vector<f_t>
_compute_elevation_above_nearest_mask(bool MFD,
																			bool distance,
																			std::vector<std::uint8_t>& mask,
																			CONNECTOR_T& con)
{

	// Computing the full graph
	con.PFcompute_all(false);

	// Init the Zabove vector (Z above the mask) and mindist, min distance to the
	// channel in MFD
	std::vector<f_t> Zabove(con.nxy(), 0), mindist;

	// Only allocate memory if MFD
	if (MFD && distance)
		mindist = std::vector<f_t>(con.nxy(), 0.);

	std::array<i_t, 8> receivers;
	std::array<f_t, 8> receiversDx;

	// Right, let's go through the stack (downstream to upstream)
	for (int i = 0; i < con.nxy(); ++i) {

		int node = (MFD) ? con.data->_stack[i] : con.data->_Sstack[i];
		if (can_out(con.data->_boundaries[node]) ||
				nodata(con.data->_boundaries[node]))
			continue;

		if (mask[node])
			continue;

		if (MFD) {
			int nn = con.Receivers(node, receivers);
			int darec = -1;
			f_t tmindx = 1.;
			if (distance) {
				con.ReceiversDx(node, receiversDx);
				f_t tmindist = std::numeric_limits<f_t>::max();
				for (int j = 0; j < nn; ++j) {
					int trec = receivers[j];
					if (mindist[trec] < tmindist) {
						tmindist = mindist[trec];
						darec = trec;
					}
				}
			} else {
				f_t tmaxZ = std::numeric_limits<f_t>::max();
				for (int j = 0; j < nn; ++j) {
					int trec = receivers[j];
					if (Zabove[trec] < tmaxZ) {
						tmaxZ = Zabove[trec];
						darec = trec;
					}
				}
			}

			if (darec > -1) {
				Zabove[node] =
					Zabove[darec] +
					std::max(con.data->_surface[node] - con.data->_surface[darec], 0.);
				if (distance)
					mindist[node] = mindist[darec] + tmindx;
			}

		} else {

			// Steepest receiver
			i_t Srec = con.Sreceivers(node);

			// Steepest receiver
			Zabove[node] =
				Zabove[Srec] +
				std::max(con.data->_surface[node] - con.data->_surface[Srec], 0.);
		}
	}

	return Zabove;
}

template<class i_t, class f_t, class CONNECTOR_T, class in_t, class out_t>
out_t
compute_elevation_above_nearest_mask(bool MFD,
																		 bool distance,
																		 in_t& mask,
																		 CONNECTOR_T& con)
{

	numvec<std::uint8_t> maski = format_input<in_t>(mask);
	std::vector<std::uint8_t> masko = to_vec(maski);

	auto out = _compute_elevation_above_nearest_mask<i_t, f_t, CONNECTOR_T>(
		MFD, distance, masko, con);

	return format_output<decltype(out), out_t>(out);
}

template<class i_t, class f_t, class CONNECTOR_T>
std::vector<i_t>
_compute_sources_D8(f_t Ath, CONNECTOR_T& con)
{

	con.PFcompute_all(false);

	std::vector<f_t> DA(con.nxy(), 0);
	std::vector<std::uint8_t> isDone(con.nxy(), false);

	for (int i = con.nxy() - 1; i >= 0; --i) {
		int node = con.data->_Sstack[i];

		if (nodata(con.data->_boundaries[node]))
			continue;

		DA[node] += con.area(node);

		int rec = con.Sreceivers(node);

		if (rec == node)
			continue;

		DA[rec] += con.data->_DA[node];
	}

	std::vector<i_t> outi;

	for (int i = con.nxy() - 1; i >= 0; --i) {
		int node = con.data->_Sstack[i];

		if (nodata(con.data->_boundaries[node]))
			continue;
		int rec = con.Sreceivers(node);

		if (rec == node)
			continue;

		if (DA[node] > Ath && isDone[rec] == false) {
			int tnode = node;
			int trec = con.Sreceivers(node);
			outi.emplace_back(node);

			while (tnode != trec) {
				isDone[tnode] = true;
				tnode = trec;
				trec = con.Sreceivers(tnode);
			}
		}
	}

	return outi;
}

template<class i_t, class f_t, class CONNECTOR_T, class out_t>
out_t
compute_sources_D8(f_t Ath, CONNECTOR_T& con)
{

	auto out = _compute_sources_D8<int, double, CONNECTOR_T>(Ath, con);

	return format_output<decltype(out), out_t>(out);
}

/// Calculates the MFD graph, then replace the elevation with a slopped surface.
/// Tries to homogeneify the slope for every links
template<class i_t, class f_t, class CONNECTOR_T>
void
replace_with_slopped_surface(f_t slope, CONNECTOR_T& con)
{

	// preformatting the output
	con.PFcompute_all(false);

	std::array<f_t, 8> dxs;
	std::array<i_t, 8> recs;

	for (auto node : con.data->_stack) {

		if (nodata(con.data->_boundaries[node]))
			continue;

		if (can_out(con.data->_boundaries[node])) {
			con.data->_surface[node] = 0;
			continue;
		}

		int nd = con.Receivers(node, recs);

		con.ReceiversDx(node, dxs);

		for (int j = 0; j < nd; ++j) {
			con.data->_surface[node] = std::max(
				con.data->_surface[node], con.data->_surface[recs[j]] + slope * dxs[j]);
		}
	}
}

}; // end of namespace
