#pragma once

/*

- Quick and dirty tools, will get a proper river objects at some points but I
need that now
- Set of functions to extract rivers quickly
- Only work with the python wrapper
- Ain't got no time ATM


Being slowly deprecated in favour of  RivNet1D.hpp, more generic

*/

#include "D8connector.hpp"
#include "dodconnector.hpp"
#include "graph.hpp"
#include "graphuncs.hpp"
#include "wrap_helper.hpp"

#ifdef DAGGER_FT_PYTHON
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace DAGGER {

template<class fT, class Connector_t, class Graph_t>
py::dict
RiverNetwork(fT threshold,
						 py::array_t<fT, 1>& tAQw,
						 py::array_t<fT, 1>& ttopo,
						 Connector_t& connector,
						 Graph_t& graph)
{

	auto AQw = DAGGER::format_input(tAQw);
	auto topo = DAGGER::format_input(ttopo);

	// output
	py::dict output;

	// basin ID
	std::vector<int> basid;

	// River Source ID
	std::vector<int> sourceid(connector.nxy(), -1);

	// River nodes
	std::vector<int> riverNode;
	riverNode.reserve(std::floor(connector.nxy() / 100));

	// Step 1: calculate basin ID
	basid = graph._get_SFD_basin_labels();

	std::vector<fT> distfromoutlet(connector.nxy(), 0.);
	graph._get_SFD_distance_from_outlets(distfromoutlet);

	int sid = -1;
	std::vector<uint8_t> done(connector.nxy(), false);
	// Step 2: River Nodes
	for (int i = connector.nxy() - 1; i >= 0; --i) {
		int node = graph.Sstack[i];

		if (connector.boundaries.no_data(node))
			continue;

		int rec = connector.Sreceivers(node);

		if (done[node]) {

			if (done[rec] == false) {
				done[rec] = true;
				sourceid[rec] = sourceid[node];
			}

			done[rec] = true;

			riverNode.emplace_back(node);
			continue;
		}

		if (AQw[node] > threshold) {
			++sid;
			sourceid[node] = sid;
			done[node] = true;
			if (done[rec] == false) {
				done[rec] = true;
				sourceid[rec] = sid;
			}

			riverNode.emplace_back(node);
		}
	}

	std::vector<int> riverRecs(riverNode.size());
	std::unordered_map<int, int> arr2riv;

	for (size_t i = 0; i < riverNode.size(); ++i) {
		arr2riv[riverNode[i]] = i;
	}

	for (size_t i = 0; i < riverNode.size(); ++i) {
		int rec = connector.Sreceivers(riverNode[i]);
		riverRecs[i] = (arr2riv.count(rec)) ? arr2riv[rec] : riverNode[i];
	}

	std::vector<int> tempi(riverNode.size(), -1);
	std::vector<fT> tempf(riverNode.size(), 0.);

	// global River nodes
	output["nodes"] = py::array(riverNode.size(), riverNode.data());

	// global River Recs
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempi[i] = connector.Sreceivers(riverNode[i]);
	}
	output["receivers"] = py::array(tempi.size(), tempi.data());

	// local River Recs
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempi[i] = connector.Sreceivers(riverNode[i]);
	}
	output["river_receivers"] = py::array(riverRecs.size(), riverRecs.data());

	// Basin ID
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempi[i] = basid[riverNode[i]];
	}
	output["basinID"] = py::array(tempi.size(), tempi.data());

	// River ID
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempi[i] = sourceid[riverNode[i]];
	}
	output["riverID"] = py::array(tempi.size(), tempi.data());

	// drainge area or discharge
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempf[i] = AQw[riverNode[i]];
	}
	output["A"] = py::array(tempf.size(), tempf.data());

	// Elevation
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempf[i] = topo[riverNode[i]];
	}
	output["Elevation"] = py::array(tempf.size(), tempf.data());

	// Dx
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempf[i] = connector.Sdistance2receivers[riverNode[i]];
	}
	output["dx"] = py::array(tempf.size(), tempf.data());

	// Distance from outlet
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempf[i] = distfromoutlet[riverNode[i]];
	}
	output["flow_distance"] = py::array(tempf.size(), tempf.data());

	return output;
}

template<class fT, class Connector_t, class Graph_t>
py::dict
DrainageDivides(Connector_t& connector,
								Graph_t& graph,
								py::array_t<fT, 1>& ttopo,
								py::array_t<int, 1>& tBID)
{

	// preformatting inputs
	auto topo = DAGGER::format_input(ttopo);
	auto BID = DAGGER::format_input(tBID);

	std::vector<int> nodes;
	nodes.reserve(connector.nxy());
	std::vector<int> basinID;
	basinID.reserve(connector.nxy());
	std::vector<fT> ZZ;
	ZZ.reserve(connector.nxy());

	auto neigh = connector.get_empty_neighbour();
	for (int i = 0; i < connector.nxy(); ++i) {
		// if nodata -> leave
		if (connector.boundaries.no_data(i))
			continue;

		int TB = BID[i];
		// getting neighbours
		int nn = connector.Neighbours(i, neigh);
		// if not the same basin ID -> register, as simple as that
		for (int j = 0; j < nn; ++j) {
			if (BID[neigh[j]] != TB) {
				nodes.emplace_back(i);
				basinID.emplace_back(BID[i]);
				ZZ.emplace_back(topo[i]);
			}
		}
	}

	// Formatting the ouput dictionnary
	py::dict output;
	output["nodes"] = py::array(nodes.size(), nodes.data());
	output["basinID"] = py::array(basinID.size(), basinID.data());
	output["elevation"] = py::array(ZZ.size(), ZZ.data());

	return output;
}

template<class fT, class Connector_t>
py::dict
RiverNetworkC8(fT threshold, Connector_t& connector)
{

	// output
	py::dict output;

	// basin ID
	std::vector<int> basid;

	// River Source ID
	std::vector<int> sourceid(connector.nxy(), -1);

	// River nodes
	std::vector<int> riverNode;
	riverNode.reserve(std::floor(connector.nxy() / 100));

	// Step 1: calculate basin ID
	basid = _compute_SFD_basin_labels<int, fT, Connector_t>(connector);

	std::vector<fT> distfromoutlet =
		_compute_SFD_distance_from_outlets<int, fT, Connector_t>(connector);

	int sid = -1;
	std::vector<uint8_t> done(connector.nxy(), false);
	// Step 2: River Nodes
	for (int i = connector.nxy() - 1; i >= 0; --i) {
		int node = connector.data->_Sstack[i];

		if (nodata(connector.data->_boundaries[node]))
			continue;

		int rec = connector.Sreceivers(node);

		if (done[node]) {

			if (done[rec] == false) {
				done[rec] = true;
				sourceid[rec] = sourceid[node];
			}

			done[rec] = true;

			riverNode.emplace_back(node);
			continue;
		}

		if (connector.data->_DA[node] > threshold) {
			++sid;
			sourceid[node] = sid;
			done[node] = true;
			if (done[rec] == false) {
				done[rec] = true;
				sourceid[rec] = sid;
			}

			riverNode.emplace_back(node);
		}
	}

	std::vector<int> riverRecs(riverNode.size());
	std::unordered_map<int, int> arr2riv;

	for (size_t i = 0; i < riverNode.size(); ++i) {
		arr2riv[riverNode[i]] = i;
	}

	for (size_t i = 0; i < riverNode.size(); ++i) {
		int rec = connector.Sreceivers(riverNode[i]);
		riverRecs[i] = (arr2riv.count(rec)) ? arr2riv[rec] : riverNode[i];
	}

	std::vector<int> tempi(riverNode.size(), -1);
	std::vector<fT> tempf(riverNode.size(), 0.);

	// global River nodes
	output["nodes"] = py::array(riverNode.size(), riverNode.data());

	// global River Recs
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempi[i] = connector.Sreceivers(riverNode[i]);
	}
	output["receivers"] = py::array(tempi.size(), tempi.data());

	// local River Recs
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempi[i] = connector.Sreceivers(riverNode[i]);
	}
	output["river_receivers"] = py::array(riverRecs.size(), riverRecs.data());

	// Basin ID
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempi[i] = basid[riverNode[i]];
	}
	output["basinID"] = py::array(tempi.size(), tempi.data());

	// River ID
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempi[i] = sourceid[riverNode[i]];
	}
	output["riverID"] = py::array(tempi.size(), tempi.data());

	// drainge area or discharge
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempf[i] = connector.data->_DA[riverNode[i]];
	}
	output["A"] = py::array(tempf.size(), tempf.data());

	// Elevation
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempf[i] = connector.data->_surface[riverNode[i]];
	}
	output["Elevation"] = py::array(tempf.size(), tempf.data());

	// Dx
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempf[i] = connector.SreceiversDx(riverNode[i]);
	}
	output["dx"] = py::array(tempf.size(), tempf.data());

	// Distance from outlet
	for (size_t i = 0; i < riverNode.size(); ++i) {
		tempf[i] = distfromoutlet[riverNode[i]];
	}
	output["flow_distance"] = py::array(tempf.size(), tempf.data());

	return output;
}

};

#endif // safeguarding cases where all hpp are included in a non-python context
