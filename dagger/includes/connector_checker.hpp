/*
Defines the standardised loose connection of the connectors.

There are 2 different level of standards:

	- the generic standard, defining the functions that will be called in the
different modules. They are the only one allowed in module accepting different
types of connectors.

	- the specifyers standards, defining function specific to build, set and
customise given subsets of connectors. For example a regular grid will need nx
and ny, not needed by a voronoi one. THESE CAN ONLY BE CALLED OUTSIDE OF THE
MODULES. The modules should not need to know how the connector is structured
(loosely connected).

It uses phantom functions, called during the building of a new connector,
checking the said connector has all the required functions to fit the template.


A typical script would construct the connector first, calling generic and
specific functions, then throw the connector inside the module (e.g. graph,
graphflood, trackscape) where it would only use the generic functions.

That way, I can for example build a connector from a regular grid, setting the
nx, ny, dx, dy and boundary condition efficiently and then build graphflood with
that connector.

*/

#ifndef CONNECTOR_CHECKER_HPP
#define CONNECTOR_CHECKER_HPP

/// Generic checker: checks if the connector has the minimal set of functions to
/// be accepted by the other modules
template<class Connector_t,
				 class fT,
				 class CONTAINER_NEIGHBOURS_INT,
				 class CONTAINER_NEIGHBOURS_fT,
				 class VECLIKE>
void
check_connector_template_generic(Connector_t& con)
{

	// Removed the init function from the checker: these are connector specific

	// General accessors:
	// return the number of nodes in the grid
	int (Connector_t::*nxy)() = &Connector_t::nxy;

	// return the empty container of n neighbour max  size
	CONTAINER_NEIGHBOURS_INT(Connector_t::*emptyNeighbour)
	() = &Connector_t::template emptyNeighbour<int>;
	CONTAINER_NEIGHBOURS_fT (Connector_t::*emptyNeighbour)() =
		&Connector_t::template emptyNeighbour<fT>;

	// Generic function precomputing elements in the graph.
	// Its interpretation is free and really depends on the type of
	// grid/connection is needed. This is where, for example, a connector
	// optimised for cpu can compute all the receivers/donors information from
	// topography; an irregular grid could regrid there, or simply, if the
	// neighbouring never caches anything for a memory-saving connector, do
	// nothing.
	void (Connector_t::*compute)() = &Connector_t::compute;

	// Connect a (new) topography, to be used for computing receivers/donors.
	void (Connector_t::*connect_topography)(VECLIKE&) =
		&Connector_t::connect_topography;

	// return the number of neighbours and fill the CONTAINER_NEIGHBOURS_INT with
	// nn neighbours
	int (Connector_t::*Neighbours)(int, CONTAINER_NEIGHBOURS_INT&) =
		&Connector_t::Neighbours;

	// return the number of neighbours and fill the CONTAINER_NEIGHBOURS_INT and
	// with nn neighbours
	int (Connector_t::*Neighbours)(
		int, CONTAINER_NEIGHBOURS_INT&, CONTAINER_NEIGHBOURS_fT&) =
		&Connector_t::neighbouring_nodes_and_distance;

	// // return the number of neighbours and fill the CONTAINER_NEIGHBOURS_INT
	// with nn links int (Connector_t::*Neighbours) (int,
	// CONTAINER_NEIGHBOURS_INT&) = &Connector_t::neighbouring_links;

	// Returns the area at a diven node index
	fT (Connector_t::*Area)(int) = &Connector_t::Area;
}

#endif
