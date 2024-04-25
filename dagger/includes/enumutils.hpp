#pragma once

namespace DAGGER {

enum class CONFLOWTOPO : std::uint8_t
{
	ALL,
	SFD,
	MFD,
	NONE
};

enum class CONBOU : std::uint8_t
{
	EDGES,
	PEW,
	PNS,
	CUSTOM

};

// Parameter types
enum class PARAMTYPE : std::uint8_t
{
	CTE,
	VAR
};

}
