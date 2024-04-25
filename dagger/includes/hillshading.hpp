/*
this file contain the algorithm(s) related to hillshading a dem from a connector
object
*/

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef HILLSHADING_HPP
#define HILLSHADING_HPP

#include "D8connector.hpp"
#include "graph.hpp"
#include "utils.hpp"

#include "wrap_helper.hpp"

namespace DAGGER {

template<class Connector_t, class topo_t, class out_t, class fT>
out_t
hillshade(Connector_t& connector, topo_t& ttopography)
{
	auto topography = format_input(ttopography);
	fT altitude = 45;
	fT azimuth = 315;
	fT z_factor = 1;
	fT pi = 3.1415926;

	// std::vector<fT> hillshade(ptr, ptr + connector.nnodes);
	std::vector<fT> hillshade(connector.nnodes, 0.);

	// convert zenith and azimuth into radians for calculation
	fT zenith_rad = (90 - altitude) * pi / 180.0;
	fT azimuth_math = 360 - azimuth + 90;
	if (azimuth_math >= 360.0)
		azimuth_math = azimuth_math - 360;
	fT azimuth_rad = azimuth_math * pi / 180.0;

	for (int i = 0; i < connector.nnodes; ++i) {
		// Ignoring no data
		if (connector.boundaries.no_data(i) || connector.boundaries.can_out(i))
			continue;

		fT slope_rad = 0;
		fT aspect_rad = 0;
		fT dzdx = 0;
		fT dzdy = 0;

		fT ij = std::abs(topography[i]);
		int j = connector.get_right_idx(i);
		fT ijp1 = (connector.is_in_bound(j)) ? std::abs(topography[j]) : ij;
		j = connector.get_bottom_idx(i);
		fT ip1j = (connector.is_in_bound(j)) ? std::abs(topography[j]) : ij;
		j = connector.get_bottomright_idx(i);
		fT ip1jp1 = (connector.is_in_bound(j)) ? std::abs(topography[j]) : ij;
		j = connector.get_topleft_idx(i);
		fT im1jm1 = (connector.is_in_bound(j)) ? std::abs(topography[j]) : ij;
		j = connector.get_top_idx(i);
		fT im1j = (connector.is_in_bound(j)) ? std::abs(topography[j]) : ij;
		j = connector.get_topright_idx(i);
		fT im1jp1 = (connector.is_in_bound(j)) ? std::abs(topography[j]) : ij;
		j = connector.get_left_idx(i);
		fT ijm1 = (connector.is_in_bound(j)) ? std::abs(topography[j]) : ij;
		j = connector.get_bottomleft_idx(i);
		fT ip1jm1 = (connector.is_in_bound(j)) ? std::abs(topography[j]) : ij;

		if (ij > 0) {
			dzdx = ((ijp1 + 2 * ip1j + ip1jp1) - (im1jm1 + 2 * im1j + im1jp1)) /
						 (8 * connector.dx);
			dzdy = ((im1jp1 + 2 * ijp1 + ip1jp1) - (im1jm1 + 2 * ijm1 + ip1jm1)) /
						 (8 * connector.dy);
			slope_rad = atan(z_factor * sqrt((dzdx * dzdx) + (dzdy * dzdy)));
			if (dzdx != 0) {
				aspect_rad = std::atan2(dzdy, (dzdx * -1));
				if (aspect_rad < 0)
					aspect_rad = 2 * pi + aspect_rad;
			} else {
				if (dzdy > 0)
					aspect_rad = pi / 2;
				else if (dzdy < 0)
					aspect_rad = 2 * pi - pi / 2;
			}

			hillshade[i] = ((std::cos(zenith_rad) * std::cos(slope_rad)) +
											(std::sin(zenith_rad) * std::sin(slope_rad) *
											 std::cos(azimuth_rad - aspect_rad)));
			// std::cout << hillshade[i] << "|";
			if (hillshade[i] < 0)
				hillshade[i] = 0;
		}
	}

	return format_output<decltype(hillshade), out_t>(hillshade);
}

template<class Graph_t, class Connector_t, class topo_t, class out_t, class fT>
out_t
rayshade(Graph_t& graph,
				 Connector_t& connector,
				 topo_t& ttopography,
				 fT ray_slope,	// Slope of the sun
				 fT shadow_mag, // magnitude of the shadow ray
				 fT smooth_r,		// magnitude of final smoothing
				 fT attenuation,
				 fT diagonals)
{
	auto topography = format_input(ttopography);

	std::vector<fT> shade(connector.nnodes, 0.);
	std::vector<fT> tshade(connector.nnodes, 0.);

	for (int i = connector.nnodes - 1; i >= 0; --i) {
		int node = graph.stack[i];
		fT zcasted = topography[node];
		fT oshadowmag = shadow_mag;

		while (true) {
			// Get the link in the right direction
			int link = connector.get_right_idx_links(node);

			if (connector.is_link_valid(link) == false)
				break;

			fT dx = connector.get_dx_from_links_idx(link);

			int onode = connector.get_other_node_from_links(link, node);
			// fT ozcaqs = zcasted;

			zcasted -= ray_slope * dx;

			// std::cout << "zcasted is " << zcasted << " vs " << topography[onode] <<
			// " and was " << ozcaqs << std::endl;
			if (zcasted < topography[onode])
				break;

			node = onode;
			tshade[onode] = std::max(oshadowmag, tshade[onode]);
			oshadowmag *= attenuation;
		}
	}

	for (int i = 0; i < connector.nnodes; ++i) {
		shade[i] += tshade[i];
	}

	if (diagonals > 0) {
		tshade = std::vector<fT>(connector.nnodes, 0.);

		for (int i = connector.nnodes - 1; i >= 0; --i) {
			int node = graph.stack[i];
			fT zcasted = topography[node];
			fT oshadowmag = shadow_mag * diagonals;

			while (true) {
				// Get the link in the right direction
				int link = connector.get_bottomright_idx_links(node);

				if (connector.is_link_valid(link) == false)
					break;

				fT dx = connector.get_dx_from_links_idx(link);

				int onode = connector.get_other_node_from_links(link, node);
				// fT ozcaqs = zcasted;

				zcasted -= ray_slope * dx;

				// std::cout << "zcasted is " << zcasted << " vs " << topography[onode]
				// << " and was " << ozcaqs << std::endl;
				if (zcasted < topography[onode])
					break;

				node = onode;
				tshade[onode] = std::max(oshadowmag, tshade[onode]);
				oshadowmag *= attenuation;
			}
		}

		for (int i = 0; i < connector.nnodes; ++i) {
			shade[i] += tshade[i];
		}

		tshade = std::vector<fT>(connector.nnodes, 0.);

		for (int i = connector.nnodes - 1; i >= 0; --i) {
			int node = graph.stack[i];
			fT zcasted = topography[node];
			fT oshadowmag = shadow_mag * diagonals;

			while (true) {
				// Get the link in the right direction
				int link = connector.get_bottomleft_idx_links(node);

				if (connector.is_link_valid(link) == false)
					break;

				fT dx = connector.get_dx_from_links_idx(link);

				int onode = connector.get_other_node_from_links(link, node);
				// fT ozcaqs = zcasted;

				zcasted -= ray_slope * dx;

				// std::cout << "zcasted is " << zcasted << " vs " << topography[onode]
				// << " and was " << ozcaqs << std::endl;
				if (zcasted < topography[onode])
					break;

				node = onode;
				tshade[onode] = std::max(oshadowmag, tshade[onode]);
				oshadowmag *= attenuation;
			}
		}

		fT tmax = 0;
		for (int i = 0; i < connector.nnodes; ++i) {
			shade[i] += tshade[i];
			tmax = std::max(tmax, shade[i]);
		}

		for (auto& v : shade)
			v /= tmax;
	}

	// REPEAT AND STACK OPERATION IN DIAGONALS
	// visited = std::vector<std::uint8_t>(connector.nnodes,0);

	if (smooth_r >= 1)
		shade = On_gaussian_blur(smooth_r, shade, connector.nx, connector.ny);

	return format_output<decltype(shade), out_t>(shade);
}

// End of namespace DAGGEr
} // namespace DAGGER

// end of definition if
#endif
