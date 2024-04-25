#pragma once

#include "declare_includes.hpp"
using namespace DAGGER;

void
delclare_enums(py::module& m)
{

	py::enum_<DEPRES>(m, "LMR")
		.value("cordonnier_fill", DEPRES::cordonnier_fill)
		.value("cordonnier_carve", DEPRES::cordonnier_carve)
		.value("cordonnier_simple", DEPRES::cordonnier_simple)
		.value("priority_flood", DEPRES::priority_flood)
		.value("priority_flood_opti", DEPRES::priority_flood_opti)
		.value("priority_full_MFD", DEPRES::priority_full_MFD)

		.value("none", DEPRES::none)
		.value("dagger_carve", DEPRES::dagger_carve)
		.value("dagger_fill", DEPRES::dagger_fill);

	py::enum_<MFD_PARTITIONNING>(m, "MFD_PARTITIONNING")
		.value("PROPOSLOPE", MFD_PARTITIONNING::PROPOSLOPE)
		.value("PROPOSLOPE_NODIAG", MFD_PARTITIONNING::PROPOSLOPE_NODIAG)
		.value("SQRTSLOPE", MFD_PARTITIONNING::SQRTSLOPE)
		.value("PROPOREC", MFD_PARTITIONNING::PROPOREC);

	py::enum_<TSC_HILLSLOPE>(m, "TSC_HILLSLOPE")
		.value("NONE", TSC_HILLSLOPE::NONE)
		.value("LINEAR", TSC_HILLSLOPE::LINEAR)
		.value("CIDRE", TSC_HILLSLOPE::CIDRE)
		.value("CIDRE_NOCRIT", TSC_HILLSLOPE::CIDRE_NOCRIT)
		.value("HYLANDS", TSC_HILLSLOPE::HYLANDS);

	py::enum_<TSC_FLUVIAL>(m, "TSC_FLUVIAL")
		.value("NONE", TSC_FLUVIAL::NONE)
		.value("DAVY2009", TSC_FLUVIAL::DAVY2009)
		.value("LATERALDAVY", TSC_FLUVIAL::LATERALDAVY)
		.value("LATERALSPL", TSC_FLUVIAL::LATERALSPL)
		.value("FASTSCAPE", TSC_FLUVIAL::FASTSCAPE);

	py::enum_<TSC_MARINE>(m, "TSC_MARINE")
		.value("NONE", TSC_MARINE::NONE)
		.value("CHARLIE", TSC_MARINE::CHARLIE);

	py::enum_<TSC_FLOW_TOPOLOGY>(m, "TSC_FLOW_TOPOLOGY")
		.value("SFD", TSC_FLOW_TOPOLOGY::SFD)
		.value("MFD", TSC_FLOW_TOPOLOGY::MFD);

	py::enum_<CONVERGENCE>(m, "GRAPHFLOOD_CONVERGENCE")
		.value("NONE", CONVERGENCE::NONE)
		.value("DHW", CONVERGENCE::DHW)
		.value("QWR", CONVERGENCE::QWR)
		.value("ALL", CONVERGENCE::ALL);

	py::enum_<RANDNOISE>(m, "NOISE")
		.value("WHITE", RANDNOISE::WHITE)
		.value("RED", RANDNOISE::RED)
		.value("PERLIN", RANDNOISE::PERLIN);

	py::enum_<CONFLOWTOPO>(m, "CONFLOW")
		.value("ALL", CONFLOWTOPO::ALL)
		.value("SFD", CONFLOWTOPO::SFD)
		.value("MFD", CONFLOWTOPO::MFD)
		.value("NONE", CONFLOWTOPO::NONE);

	py::enum_<CONBOU>(m, "CONBOU")
		.value("EDGES", CONBOU::EDGES)
		.value("PEW", CONBOU::PEW)
		.value("PNS", CONBOU::PNS)
		.value("CUSTOM", CONBOU::CUSTOM);

	py::enum_<MORPHOMODE>(m, "MORPHOMODE")
		.value("NONE", MORPHOMODE::NONE)
		.value("EROS", MORPHOMODE::EROS)
		.value("EROSVEC", MORPHOMODE::EROSVEC)
		.value("MPM", MORPHOMODE::MPM)
		.value("MPMVEC", MORPHOMODE::MPMVEC);

	py::enum_<HYDROMODE>(m, "HYDROMODE")
		.value("MFD", HYDROMODE::MFD)
		.value("SFD", HYDROMODE::SFD)
		.value("VEC", HYDROMODE::VEC);
}
