#pragma once

#include "dodconnector.hpp"
#include "dodcontexts.hpp"
#include "utils.hpp"
#include <vector>

namespace DAGGER {

template<class i_t, class f_t, class CON_T, class DATA_T, class PARAM_T>
class TinySubGraph
{
public:
	TinySubGraph(){};
	TinySubGraph(CON_T& con, DATA_T& data, PARAM_T& param)
	{
		this->con = &con;
		this->data = &data;
		this->param = &param;
		this->tempN = std::vector<std::uint8_t>(this->con->nxy(), false);
		this->xtraMask = std::vector<std::uint8_t>(this->con->nxy(), true);
	};

	~TinySubGraph(){};

	CON_T* con;
	DATA_T* data;
	PARAM_T* param;

	std::vector<i_t> stack;
	std::vector<i_t> nodes;
	std::vector<i_t> baseLevels;
	std::vector<std::uint8_t> tempN, xtraMask;

	template<class CONTAINER_INT>
	void build_from_donor_sources(CONTAINER_INT& startingNodes)
	{
		this->label_from_donor_sources(startingNodes);
		this->build_stack();
	}

	template<class CONTAINER_INT>
	void label_from_donor_sources(CONTAINER_INT& startingNodes)
	{

		this->stack.clear();
		this->nodes.clear();
		this->baseLevels.clear();
		fillvec(this->tempN, false);

		// Initialising a node queue
		std::queue<i_t> tQ;

		// Feeding it witht the starting nodes
		for (auto v : startingNodes) {
			nodes.emplace_back(v);
			tempN[v] = true;
			tQ.emplace(v);
		}

		// Setting up context and helper arrays
		CT_neighbourer_1<i_t, f_t> ctx;
		std::array<i_t, 8> recs;
		std::array<f_t, 8> recdxs;
		std::array<std::uint8_t, 8> recbits;

		while (tQ.empty() == false) {

			// Getting next node and popping it out of the Q
			i_t nextnode = tQ.front();
			tQ.pop();

			// saving the node
			this->nodes.emplace_back(nextnode);

			// reset the node recs and donors
			this->con->reset_node(nextnode);

			// compute receivers only at that specific node
			this->con->__compute_recs_single_node_mask(nextnode, ctx, this->xtraMask);

			// gathering receivers
			int nr = this->con->Receivers(nextnode, recs);
			if (this->param->TSG_dist)
				this->con->ReceiversDx(nextnode, recdxs);

			// this->tempN[nextnode] = nr;

			if (nr == 0) {
				this->baseLevels.emplace_back(nextnode);
				continue;
			}

			// Gathering the receivers into the queue
			for (int j = 0; j < nr; ++j) {
				int trec = recs[j];
				if (this->xtraMask[trec] == false)
					continue;

				// double checking they are not already in there/processed
				if (this->tempN[trec] == false) {
					tQ.emplace(trec);
					this->tempN[trec] = true;
				}
			}
		}

		// A that point, all the nodes and baselevels for the graph have been
		// gathered And this->tempN[trec] has the number of receivers per node
		//=================

		// Let's invert the receiver codes into donors
		for (auto v : this->nodes) {
			this->con->__invert_recs_at_node(v, recbits, recs);
			this->tempN[v] = this->con->nReceivers(v);
		}

		this->data->ibag["tsgNodes"] = this->nodes;
		std::vector<int> temp;

		for (auto v : this->xtraMask)
			temp.emplace_back(static_cast<int>(v));

		this->data->ibag["tsgMask"] = temp;
	}

	void build_stack()
	{

		// Initialising a node queue
		std::queue<i_t> tQ;

		std::array<i_t, 8> recs;

		// Build the stack
		this->stack.reserve(this->nodes.size());

		// the stack will start with the base levels
		for (auto v : this->baseLevels)
			tQ.emplace(v);

		while (tQ.empty() == false) {
			// Getting next node and popping it out of the Q
			int nextnode = tQ.front();
			tQ.pop();
			// if the node is inj Q is ready 4 stack
			this->stack.emplace_back(nextnode);

			// Grabbing its donors subgraphically speaking
			int nd = this->con->Donors(nextnode, recs);
			for (int i = 0; i < nd; ++i) {
				// each donors is visited
				int tdon = recs[i];
				// and I keep track of the number of visits
				--this->tempN[tdon];
				// if the donor has been visited by all its receivers, then it is ready
				if (this->tempN[tdon] == 0) {
					tQ.emplace(tdon);
				}
			}
		}

		// std::vector<int> acc(this->con->nxy(),1);
		// for(int i=this->stack.size()-1;i>=0; --i){
		// 	if(this->con->Sreceivers(this->stack[i]) != this->stack[i])
		// 		acc[this->con->Sreceivers(this->stack[i])] += acc[this->stack[i]];
		// }

		// fillvec(this->tempN,false);
		// for(int i=0; i< this->baseLevels.size(); ++i){
		// 	int tbl = this->baseLevels[i];
		// 	tQ.emplace(tbl);
		// 	while(tQ.empty()==false){
		// 		int nextnode = tQ.front();
		// 		acc[nextnode]++;

		// 		tQ.pop();
		// 		int nd = this->con->Donors(nextnode,recs);
		// 		for(int j=0; j<nd; ++j){
		// 			int tdon = recs[j];
		// 			if(this->tempN[tdon]==false){
		// 				this->tempN[tdon] = true;
		// 				tQ.emplace(tdon);
		// 			}
		// 		}
		// 	}
		// }

		// this->data->ibag["acctsg"] = acc;
		// Done
	}
};

}
