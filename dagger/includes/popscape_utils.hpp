#ifndef POPSCAPE_UTILS_HPP
#define POPSCAPE_UTILS_HPP
// -> The connector classes
#include "D4connector.hpp"
#include "D8connector.hpp"
#include "PerlinNoise.hpp"
#include "graph.hpp"

enum class RANDNOISE
{
	WHITE,
	RED,
	PERLIN,
};

template<class T>
DAGGER::D8connector<T>
_create_connector(int nx, int ny, T dx, T dy, T xmin, T ymin)
{
	DAGGER::D8connector<T> con = DAGGER::D8connector<T>();
	con.default_permissive = true;
	con.init_dimensions(nx, ny, dx, dy, xmin, ymin);
	return con;
}

template<class T>
DAGGER::D8connector<T>*
_create_connector_fromptr(int nx, int ny, T dx, T dy, T xmin, T ymin)
{
	DAGGER::D8connector<T>* con = new DAGGER::D8connector<T>();
	con->default_permissive = true;
	con->init_dimensions(nx, ny, dx, dy, xmin, ymin);
	return con;
}

template<class T>
void
_create_connector(int nx,
									int ny,
									T dx,
									T dy,
									T xmin,
									T ymin,
									DAGGER::D4connector<T>& con)
{
	con = DAGGER::D4connector<T>(nx, ny, dx, dy, xmin, ymin);
}

template<class T>
void
_create_connector_fromptr(int nx,
													int ny,
													T dx,
													T dy,
													T xmin,
													T ymin,
													DAGGER::D4connector<T>* con)
{
	auto tcon = DAGGER::D4connector<T>(nx, ny, dx, dy, xmin, ymin);
	con = &tcon;
}

template<class Connector_t, class T>
void
_create_graph(int nxy, Connector_t& con, DAGGER::graph<T, Connector_t>& gr)
{
	gr = DAGGER::graph<T, Connector_t>(con);
}

template<class Connector_t, class T>
void
_create_graph_fromptr(int nxy,
											Connector_t& con,
											DAGGER::graph<T, Connector_t>*& gr)
{
	gr = new DAGGER::graph<T, Connector_t>(con);
	// gr = &tgr;
}
// DAGGER::D4connector _create_connector(int nx, int ny, T dx, T dy, T xmin, T
// ymin){return D4connector(int nx, int ny, T dx, T dy, T xmin, T ymin);}

class Particle
{
public:
	// Construct particle at position _pos
	Particle(int _pos) { pos = _pos; }

	int pos = 0;
	std::pair<double, double> speed = { 0., 0. };

	double volume = 1.0;	 // Total particle volume
	double sediment = 0.0; // Fraction of volume that is sediment!

	std::pair<double, double> get_normalised_speed()
	{
		double tnorm = std::max(std::abs(speed.first), std::abs(speed.second));
		return { speed.first / tnorm, speed.second / tnorm };
	}
	void apply_friction(double frico)
	{
		speed.first *= frico;
		speed.second *= frico;
	}
	void speed_up(std::pair<double, double>& speedin, double weigght)
	{
		this->speed.first += speedin.first * weigght;
		this->speed.second += speedin.second * weigght;
	}
};

template<class fT, class Connector_t>
std::vector<fT>
_generate_perlin_noise_2D(Connector_t& con,
													fT frequency,
													int octaves,
													std::uint32_t seed)
{
	std::vector<fT> noise(con.nx * con.ny, 0.);
	fT fx = frequency / con.nx;
	fT fy = frequency / con.ny;
	siv::PerlinNoise perlin{ seed };

	for (int i = 0; i < con.nnodes; ++i) {
		int xi, yi;
		con.rowcol_from_node_id(i, yi, xi);
		noise[i] = perlin.octave2D_01((xi * fx), (yi * fy), octaves);
	}

	return noise;
}

template<class out_t, class fT, class Connector_t>
out_t
generate_perlin_noise_2D(Connector_t& con,
												 fT frequency,
												 int octaves,
												 std::uint32_t seed)
{
	std::vector<fT> noise =
		_generate_perlin_noise_2D(con, frequency, octaves, seed);

	return DAGGER::format_output<std::vector<fT>, out_t>(noise);
}

#endif
