//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#pragma once

#include "npy.hpp"
#include "vecutils.hpp"

// STL imports
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

const double NPPI = std::acos(-1.0);

namespace DAGGER {

// Enumeration of the possible depression solvers
enum class DEPRES
{
	// No depression resolution
	none,
	// Using Cordonnier et al., 2019 with the filling method (algorithm 4)
	cordonnier_fill,
	// Using Cordonnier et al., 2019 with the carving method (algorithm 3)
	cordonnier_carve,
	// Using Cordonnier et al., 2019 with the carving method (algorithm 3)
	cordonnier_simple,
	// Using the Wei et al 2018 variation of Barnes 2014 priority flood algorithm
	priority_flood,
	priority_full_MFD,
	priority_flood_opti,
	// Using the Barnes et al., 2014 priority flood algorithm
	priority_flood_2014,
	// TO ADD:: Cordonnier simple, Lindsay carving and why not others??
	dagger_carve,
	dagger_fill,
	dagger_simple,

};

// Quick function checking if an element is in a set or not
template<class T>
inline bool
izinset(const std::set<T>& tset, const T& element)
{
	return std::binary_search(tset.begin(), tset.end(), element);
};

// PQ_helper is a simple generic template class aiming to provide
// a structure for priority queue elements made of an ID of type T sorted by
// a score of type U.
// A typical example would be node ID sorted by elevation
template<class T, class U>
class PQ_helper
{
public:
	// empty constructor
	PQ_helper() = default;
	// Constructor by default
	PQ_helper(T node, U score)
	{
		this->node = node;
		this->score = score;
	};
	// Node index
	T node;
	// Score data
	U score;
};

// Custom operator sorting the nodes by scores
template<class T, class U>
inline bool
operator>(const PQ_helper<T, U>& lhs, const PQ_helper<T, U>& rhs)
{
	return lhs.score > rhs.score;
}
template<class T, class U>
inline bool
operator<(const PQ_helper<T, U>& lhs, const PQ_helper<T, U>& rhs)
{
	return lhs.score < rhs.score;
}

// Utility function allowing to access the container behind the priority queue
// Totally unsafe if hte priority queue has to be used again, but quite handy if
// you simply need to iterate through its elements
template<class T, class S, class C>
S&
Container(std::priority_queue<T, S, C>& q)
{
	struct HackedQueue : private std::priority_queue<T, S, C>
	{
		static S& Container(std::priority_queue<T, S, C>& q)
		{
			return q.*&HackedQueue::c;
		}
	};
	return HackedQueue::Container(q);
};

template<class n_t, class dist_t>
class Neighbour
{
public:
	Neighbour() { ; }
	Neighbour(n_t node, dist_t distance)
	{
		this->node = node;
		this->distance = distance;
	}
	n_t node;
	dist_t distance;
};

// Generate a random number between min and max (C style random umber
// generation) Quick and easy option, not the best for true randomness (e.g.
// with emscripten it always return the same series of random number for a given
// context)
template<class T>
T
random_number(T min, T max)
{
	T random = ((T)rand()) / (T)RAND_MAX;
	T range = max - min;
	return (random * range) + min;
}

// Return a vector of sorted indices
// Equivalent of argsort in numpy or other platforms
// Shamelessly ripped from https://stackoverflow.com/a/12399290/7114716
// - credits to Lukasz Wiklendt
template<class U>
std::vector<size_t>
sort_indexes(U& v)
{

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values
	std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
		return v[i1] < v[i2];
	});

	return idx;
}

/*
=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

Implementation of linear time Gaussian blurring.
Adapted from https://blog.ivank.net/fastest-gaussian-blur.html
Credit to Іван Куцкір

The main function is On_gaussian_blur (see bellow)
*/

std::vector<int>
boxesForGauss(double sigma,
							int n) // standard deviation, number of boxes
{
	double wIdeal =
		std::sqrt((12. * sigma * sigma / n) + 1); // Ideal averaging filter width
	int wl = std::floor(wIdeal);
	if (wl % 2 == 0)
		wl--;
	int wu = wl + 2;

	double mIdeal =
		(12. * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4);
	int m = std::round(mIdeal);
	// var sigmaActual = Math.sqrt( (m*wl*wl + (n-m)*wu*wu - n)/12 );

	std::vector<int> sizes;
	sizes.reserve(n);
	for (int i = 0; i < n; ++i) {
		sizes.emplace_back(i < m ? wl : wu);
	}

	return sizes;
}

template<class T>
void
boxBlurH_4(std::vector<T>& scl, std::vector<T>& tcl, int w, int h, double r)
{
	double iarr = 1. / (r + r + 1);
	for (int i = 0; i < h; ++i) {
		int ti = i * w, li = ti, ri = ti + r;
		double fv = scl[ti], lv = scl[ti + w - 1], val = (r + 1) * fv;
		for (int j = 0; j < r; ++j) {
			val += scl[ti + j];
		}
		for (int j = 0; j <= r; ++j) {
			val += scl[ri++] - fv;
			tcl[ti++] = std::round(val * iarr);
		}
		for (int j = r + 1; j < w - r; ++j) {
			val += scl[ri++] - scl[li++];
			tcl[ti++] = val * iarr;
		}
		for (int j = w - r; j < w; ++j) {
			val += lv - scl[li++];
			tcl[ti++] = val * iarr;
		}
	}
}

template<class T>
void
boxBlurT_4(std::vector<T>& scl, std::vector<T>& tcl, int w, int h, double r)
{
	T iarr = 1. / (r + r + 1);
	for (int i = 0; i < w; i++) {
		int ti = i, li = ti, ri = ti + r * w;
		double fv = scl[ti], lv = scl[ti + w * (h - 1.)], val = (r + 1.) * fv;
		for (int j = 0; j < r; ++j)
			val += scl[ti + j * w];
		for (int j = 0; j <= r; ++j) {
			val += scl[ri] - fv;
			tcl[ti] = val * iarr;
			ri += w;
			ti += w;
		}
		for (int j = r + 1; j < h - r; ++j) {
			val += scl[ri] - scl[li];
			tcl[ti] = val * iarr;
			li += w;
			ri += w;
			ti += w;
		}
		for (int j = h - r; j < h; ++j) {
			val += lv - scl[li];
			tcl[ti] = val * iarr;
			li += w;
			ti += w;
		}
	}
}

template<class T>
void
boxBlur_4(std::vector<T>& scl, std::vector<T>& tcl, int w, int h, double r)
{
	for (size_t i = 0; i < scl.size(); ++i)
		tcl[i] = scl[i];
	boxBlurH_4(tcl, scl, w, h, r);
	boxBlurT_4(scl, tcl, w, h, r);
}

template<class T>
void
gaussBlur_4(std::vector<T>& scl, std::vector<T>& tcl, double r, int nx, int ny)
{
	auto bxs = boxesForGauss(r, 3);
	boxBlur_4(scl, tcl, nx, ny, (bxs[0] - 1) / 2);
	boxBlur_4(tcl, scl, nx, ny, (bxs[1] - 1) / 2);
	boxBlur_4(scl, tcl, nx, ny, (bxs[2] - 1) / 2);
}

/// Main function to be called for the Gaussian Blurring
/// r is the radius for the blurring
/// topo is the original array
/// nx,ny are the number of col and rows
template<class T>
std::vector<T>
On_gaussian_blur(double r, std::vector<T>& topo, int nx, int ny)
{
	std::vector<T> newtopo(topo);
	gaussBlur_4(topo, newtopo, r, nx, ny);
	return newtopo;
}

/*
=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
*/

/// Normalises a vector between 0 and 1
template<class T>
void
normalise_vector(std::vector<T>& vec)
{
	T min = std::numeric_limits<T>::max();
	T max = std::numeric_limits<T>::min();

	for (auto v : vec) {
		if (v < min)
			min = v;
		else if (v > max)
			max = v;
	}

	for (auto& v : vec)
		v = (v - min) / (max - min);
}

// Small class emulating some aspect of a vector vector from its pointer
// DAGGER relies heavily on templates to make its internal function usable from
// numpy arrays, R, js, c++... pvector process a vector into a pointer to it for
// this genericity
template<class T>
class pvector
{
public:
	std::shared_ptr<std::vector<T>> data;
	pvector(){};
	pvector(std::vector<T>& dat)
	{
		this->data = std::make_shared<std::vector<T>>(dat);
	}
	pvector(std::vector<T>* dat)
	{
		this->data = std::make_shared<std::vector<T>>(dat);
	}
	size_t size() { return this->data->size(); }
	T& operator[](int i) { return (*this->data)[i]; }
	void emplace_back(T& tin) { this->data->emplace_back(tin); }
	void push_back(T& tin) { this->data->push_back(tin); }
	void shrink_to_fit() { this->data->shrink_to_fit(); }
	void clear() { this->data->clear(); }
};

// template <typename T>
// std::vector<size_t> sort_indexes(T &v) {

//	 // initialize original index locations
//	 int size = int(v.size());
//	 std::vector<size_t> idx(size);
//	 std::iota(idx.begin(), idx.end(), 0);

//	 // sort indexes based on comparing values in v
//	 // using std::stable_sort instead of std::sort
//	 // to avoid unnecessary index re-orderings
//	 // when v contains elements of equal values
//	 std::sort(idx.begin(), idx.end(),
//				[&v](size_t i1, size_t i2) {return v[i1] <
// v[i2];});

//	 return idx;
// }

// Pretty self explanatory, T needs to be a vector-like data structure
// accessible via [] and having a size() attribute
template<class T, class U>
void
add_noise_to_vector(T& vec, U min, U max)
{

	std::random_device rd;	// obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_real_distribution<> distr(min, max); // define the range

	for (size_t i = 0; i < vec.size(); ++i) {
		auto cd = distr(gen);
		vec[i] += cd;
	}
}

/*
Sets of function guaranting compatibility with c++ use of the library
Format input vectors into pvectors
and pvectors into themselves
This is required to make the template generic
Note that I could use c++20 standard with constraints for that but it is a bit
early
*/

template<class T>
pvector<T>
_format_input(std::vector<T>& in)
{
	return pvector<T>(in);
}

template<class T>
pvector<T>
_format_input(pvector<T>& in)
{
	return in;
}

/*
To_vec functions return a copy of any formatted input into vectors
*/

template<class T>
std::vector<T>
to_vec(pvector<T>& in)
{
	return std::vector<T>(*in.data);
}

template<class T>
std::vector<T>
to_vec(std::vector<T>& in)
{
	return std::vector<T>(in);
}

/// Quick and easy lightweigth class to generate random numbers, by default
/// between 0 and 1, and by extension to any range
class easyRand
{
public:
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen;			 // seed the generator
	std::uniform_real_distribution<> distr; // define the range

	easyRand()
	{
		std::random_device rd;
		this->gen = std::mt19937(rd());												// seed the generator
		this->distr = std::uniform_real_distribution<>(0, 1); // define the range
	}

	easyRand(double min, double max)
	{
		std::random_device rd;
		this->gen = std::mt19937(rd()); // seed the generator
		this->distr =
			std::uniform_real_distribution<>(min, max); // define the range
	}

	double get() { return this->distr(this->gen); }
};

/// Quick and easy class to time the code using the std::chrono std library
class ocarina
{
public:
	std::chrono::high_resolution_clock::time_point start, end;

	ocarina(){};

	void tik() { this->start = std::chrono::high_resolution_clock::now(); }
	void tok(std::string message)
	{
		this->end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> ms_double =
			this->end - this->start;
		double timing = ms_double.count();
		std::cout << message << " took " << timing << " ms" << std::endl;
	}
	double tok()
	{
		this->end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> ms_double =
			this->end - this->start;
		double timing = ms_double.count();
		return timing;
	}
};

// simple signum function (returns -1, 0 or 1)
template<typename T>
int
sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}

// thanks https://stackoverflow.com/a/34199939
template<typename E>
constexpr auto
Enum2UnderlyingType(E e)
{
	return static_cast<typename std::underlying_type<E>::type>(e);
}

class DBBTree
{

public:
	// n nodes in the tree
	int nnodes = 0;
	// receivers for each nodes (emty if base)
	std::vector<std::vector<int>> receivers;
	// 0: building - 1: subtree - 2: open
	std::vector<std::uint8_t> label;
	std::vector<int> top;
	std::vector<int> outlet_node;
	std::vector<double> outlet_Z;

	DBBTree() { this->init(200); };

	void init(int resa)
	{
		this->receivers.reserve(resa);
		this->label.reserve(resa);
		this->top.reserve(resa);
		this->outlet_node.reserve(resa);
		this->outlet_Z.reserve(resa);
	}

	int add()
	{
		int tlab = this->nnodes;
		++this->nnodes;
		this->receivers.emplace_back(std::vector<int>());
		this->label.emplace_back(0);
		this->top.emplace_back(tlab);
		this->outlet_node.emplace_back(-1);
		this->outlet_Z.emplace_back(0.);
		return tlab;
	}

	template<class CONTAINER>
	int merge(CONTAINER tomerge)
	{
		int tlab = this->nnodes;
		++this->nnodes;
		this->receivers.emplace_back(std::vector<int>());
		this->label.emplace_back(0);
		this->top.emplace_back(tlab);
		this->outlet_node.emplace_back(-1);
		this->outlet_Z.emplace_back(0.);

		std::queue<int> mergehelper;
		for (auto tm : tomerge)
			mergehelper.emplace(tm);

		while (mergehelper.empty() == false) {
			int next = mergehelper.front();
			mergehelper.pop();
			if (this->top[next] == tlab)
				continue;
			this->top[next] = tlab;
			this->label[next] = 1;
			this->receivers[tlab].emplace_back(next);
			for (auto ttm : this->receivers[next])
				mergehelper.emplace(ttm);
		}

		return tlab;
	}
};

template<class fT>
void
save_vec_to_2Dnpy(std::string fname, int nx, int ny, std::vector<fT>& data)
{

	long unsigned tny = ny;
	long unsigned tnx = nx;

	const std::vector<long unsigned> shape{ tny, tnx };
	const bool fortran_order{ false };
	npy::SaveArrayAsNumpy(fname, fortran_order, shape.size(), shape.data(), data);
}

template<class fT>
void
save_vec_to_1Dnpy(std::string fname, int nx, int ny, std::vector<fT>& data)
{
	long unsigned tny = ny;
	long unsigned tnx = nx;
	const std::vector<long unsigned> shape{ tny * tnx };
	const bool fortran_order{ false };
	npy::SaveArrayAsNumpy(fname, fortran_order, shape.size(), shape.data(), data);
}

template<class fT>
std::vector<fT>
load_npy(std::string path)
{
	npy::npy_data d = npy::read_npy<fT>(path);

	std::vector<fT> data = d.data;

	return data;
}

template<class T>
void
minmax(std::vector<T>& vec, T& tmin, T& tmax)
{
	tmin = std::numeric_limits<T>::max();
	tmax = std::numeric_limits<T>::min();
	for (auto v : vec) {
		if (v < tmin)
			tmin = v;
		if (v > tmax)
			tmax = v;
	}
}

template<class T, class U>
void
fillvec(std::vector<T>& vec, U val)
{
	for (size_t i = 0; i < vec.size(); ++i)
		vec[i] = val;
}

template<class T, class U>
void
fillvec(std::vector<T>& vec, U val, std::vector<int>& nodes)
{
	for (size_t i = 0; i < nodes.size(); ++i)
		vec[nodes[i]] = val;
}

template<class T>
void
fillvecrange(std::vector<T>& vec)
{
	for (size_t i = 0; i < vec.size(); ++i)
		vec[i] = i;
}

template<class f_t>
int
sign(f_t x)
{
	return std::signbit(x) ? -1 : 1;
}

// template<T>
// class V2
// {
// public:

// };

// end of DAGGER namespace
}; // namespace DAGGER
