#pragma once

#include "enumutils.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t, class f_t>
class Connector8;
template<class i_t, class f_t>
class Hermes;

template<class i_t, class f_t>
class RivNet1D
{

public:
	RivNet1D(Connector8<i_t, f_t>& con, Hermes<i_t, f_t>& dbag)
	{
		this->con = &con;
		this->data = &dbag;
	}

	// Connector adn data bag
	Connector8<i_t, f_t>* con;
	Hermes<i_t, f_t>* data;

	// N nodes size data
	i_t N = 0;

	// # integer data
	std::vector<i_t> nodes;
	std::vector<i_t> recs;
	std::vector<i_t> source_keys;

	// # floating point data
	std::vector<f_t> dXs;
	std::vector<f_t> flow_distance;

	// Specific vectors
	std::vector<i_t> sources;
	std::vector<i_t> outlets;

	void reset()
	{
		this->N = 0;
		this->nodes.clear();
		this->source_keys.clear();
		this->sources.clear();
		this->flow_distance.clear();
		this->dXs.clear();
	}

	template<class int_t>
	void build_from_sources(int_t& tsources)
	{
		auto doodo = format_input<int_t>(tsources);
		std::vector<int> sources = to_vec(doodo);

		if (this->data->_Sstack.size() == 0)
			throw std::runtime_error(
				"Single flow stack needed for river extraction yo");

		this->reset();

		this->sources = sources;

		std::vector<std::uint8_t> isDone(this->con->nxy(), false);
		std::vector<i_t> invind(this->con->nxy(), -1);

		int tsk = -1;
		for (auto source : this->sources) {
			++tsk;
			int node = source;
			while (true) {
				if (isDone[node])
					break;

				this->nodes.emplace_back(node);
				invind[node] = N;
				++this->N;

				this->source_keys.emplace_back(tsk);

				isDone[node] = true;

				if (can_out(this->data->_boundaries[node])) {
					this->outlets.emplace_back(node);
					break;
				}

				node = this->con->Sreceivers(node);
			}
		}

		this->recs = std::vector<i_t>(this->N, 0);
		for (int i = 0; i < this->N; ++i) {
			int node = this->nodes[i];
			this->recs[i] = i;
			if (can_out(this->data->_boundaries[node]))
				continue;
			this->recs[i] = invind[this->con->Sreceivers(node)];
		}

		this->dXs = std::vector<f_t>(this->N, 0.);
		this->flow_distance = std::vector<f_t>(this->N, 0.);

		for (int i = N - 1; i >= 0; --i) {
			int node = this->nodes[i];
			if (can_out(this->data->_boundaries[node])) {
				this->dXs[i] = 0.;
				continue;
			}

			this->dXs[i] = (this->con->SreceiversDx(node));
			this->flow_distance[i] =
				this->flow_distance[this->recs[i]] + this->dXs.back();
		}

		return;
	}

	// Outputting functions
	// Vectors of nodes
	template<class out_t>
	out_t get_nodes()
	{
		return format_output<decltype(this->nodes), out_t>(this->nodes);
	}
	template<class out_t>
	out_t get_recs()
	{
		return format_output<decltype(this->recs), out_t>(this->recs);
	}
	template<class out_t>
	out_t get_source_keys()
	{
		return format_output<decltype(this->source_keys), out_t>(this->source_keys);
	}
	template<class out_t>
	out_t get_dXs()
	{
		return format_output<decltype(this->dXs), out_t>(this->dXs);
	}
	template<class out_t>
	out_t get_flow_distance()
	{
		return format_output<decltype(this->flow_distance), out_t>(
			this->flow_distance);
	}
	template<class out_t>
	out_t get_sources()
	{
		return format_output<decltype(this->sources), out_t>(this->sources);
	}
	template<class out_t>
	out_t get_outlets()
	{
		return format_output<decltype(this->outlets), out_t>(this->outlets);
	}
};

};
