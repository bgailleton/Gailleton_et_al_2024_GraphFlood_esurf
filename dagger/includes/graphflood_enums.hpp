#pragma once

enum class HYDRO // : std::uint8_t
{
	GRAPH_SFD,
	GRAPH_MFD,
	GRAPH_HYBRID,
};

enum class MORPHO // : std::uint8_t
{
	NONE,
	TL,
};

enum class PARAM_KE // : std::uint8_t
{
	CONSTANT,
	VARIABLE,
	EROSION,
};

enum class PARAM_DT_HYDRO // : std::uint8_t
{
	CONSTANT,
	VARIABLE,
	COURANT,
};

enum class PARAM_DT_MORPHO // : std::uint8_t
{
	CONSTANT,
	VARIABLE,
	COURANT,
	HYDRO,
};

enum class HYDROGRAPH_LM
{
	IGNORE,
	REROUTE,
	FILL,
};

enum class MFD_PARTITIONNING
{
	PROPOSLOPE,
	PROPOSLOPE_NODIAG,
	SQRTSLOPE,
	PROPOREC,
};

enum class WATER_INPUT
{
	PRECIPITATIONS_CONSTANT,
	PRECIPITATIONS_VARIABLE,
	ENTRY_POINTS_H,
	ENTRY_POINTS_QW,
};

enum class SED_INPUT
{
	NONE,
	ENTRY_POINTS_Q,
};

enum class BOUNDARY_HW
{
	FIXED_HW,
	FIXED_SLOPE,
};

enum class CONVERGENCE
{
	NONE,
	DHW,
	QWR,
	ALL
};

// Specific to GF v2
enum class SUBGRAPHMETHOD
{
	V1,
	FILLONLY,
	QWOUTK
};

enum class RUN_GF2
{
	NORMAL,
	COMPUTEqr,
};

enum class MORPHOMODE : std::uint8_t
{
	NONE,
	EROS,
	MPM,
	MPMVEC,
	EROSVEC,

};

enum class HYDROMODE : std::uint8_t
{
	MFD,
	SFD,
	VEC
};
