#ifndef TRACKSCAPE_ENUM_HPP
#define TRACKSCAPE_ENUM_HPP

#pragma once

namespace DAGGER {

enum class TSC_HILLSLOPE { NONE, LINEAR, CIDRE, CIDRE_NOCRIT, HYLANDS };

enum class TSC_FLUVIAL { NONE, DAVY2009, LATERALDAVY, LATERALSPL, FASTSCAPE };

enum class TSC_MARINE { NONE, CHARLIE };

enum class TSC_FLOW_TOPOLOGY { SFD, MFD };

}; // end of namespace DAGGER

#endif
