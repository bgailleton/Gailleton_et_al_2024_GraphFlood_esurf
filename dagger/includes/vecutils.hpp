#pragma once

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

namespace VEC2D {

// Take a cartesian vector and returns its magnitude
template<class f_t>
f_t
magnitude(f_t x, f_t y)
{
	return std::sqrt(std::pow(x, 2) + std::pow(y, 2));
}

// Transform a cartesian vector into a polar one
template<class f_t>
void
cart2pol(f_t x, f_t y, f_t& mag, f_t& theta)
{
	mag = magnitude(x, y);
	theta = std::atan2(y, x);
}

// Transform a unit cartesian vector into a polar one
template<class f_t>
void
cart2pol(f_t x, f_t y, f_t& theta)
{
	theta = std::atan2(y, x);
}

// Transform a cartesian vector into a polar one
template<class f_t>
void
pol2cart(f_t mag, f_t theta, f_t& x, f_t& y)
{
	x = mag * std::cos(theta);
	y = mag * std::sin(theta);
}

// Transform a unit cartesian vector into a polar one
template<class f_t>
void
pol2cart(f_t theta, f_t& x, f_t& y)
{
	x = std::cos(theta);
	y = std::sin(theta);
}

// transform a cartesian vector to a unit one
template<class f_t>
void
unit(f_t x, f_t y, f_t& xu, f_t& yu)
{
	f_t mag = magnitude(x, y);
	xu = x / mag;
	yu = y / mag;
}

template<class f_t>
void
add(f_t x1, f_t y1, f_t x2, f_t y2, f_t& xr, f_t& yr)
{
	xr = x1 + x2;
	yr = y1 + y2;
}

template<class f_t>
void
add_pol(f_t r1, f_t theta1, f_t r2, f_t theta2, f_t& rr, f_t& thetar)
{
	f_t x1, x2, xr, y1, y2, yr;
	pol2cart(r1, theta1, x1, y1);
	pol2cart(r2, theta2, x2, y2);

	xr = x1 + x2;
	yr = y1 + y2;

	cart2pol(xr, yr, rr, thetar);
}

template<class f_t>
void
sub(f_t x1, f_t y1, f_t x2, f_t y2, f_t& xr, f_t& yr)
{
	xr = x1 - x2;
	yr = y1 - y2;
}

// Transform a unit cartesian vector into a polar one
template<class f_t>
void
normal(f_t x, f_t y, f_t& xn, f_t& yn, bool clockwise = true)
{
	if (clockwise) {
		xn = -y;
		yn = x;
	} else {
		xn = y;
		yn = -x;
	}
}

// transform a cartesian vector to a unit one
template<class f_t>
f_t
cross(f_t x1, f_t y1, f_t x2, f_t y2)
{
	return x1 * y2 - y1 * x2;
}

// transform a cartesian vector to a unit one
template<class f_t>
f_t
mean_theta(f_t r1, f_t theta1, f_t r2, f_t theta2)
{
	f_t x1, x2, y1, y2, w1, w2, mtheta, sum = r1 + r2;
	if (sum == 0)
		return 0.;
	w1 = r1 / sum;
	w2 = r2 / sum;
	pol2cart(r1, theta1, x1, y1);
	pol2cart(r2, theta2, x2, y2);
	cart2pol(x1 + x2, y1 + y2, mtheta);
	// std::cout << "x1|y1 " << x1 <<"|" <<y1 << " x2|y2 " << x2 <<"|" <<y2 <<
	// std::endl;

	return mtheta;
}

template<class f_t>
f_t
radius_of_curvature(f_t vxin, f_t vyin, f_t vxout, f_t vyout)
{
	f_t x1, x2, x3, y1, y2, y3;
	x1 = 0;
	y1 = 0;
	x2 = x1 + vxin;
	y2 = y1 + vyin;
	x3 = x2 + vxout;
	y3 = y2 + vyout;
	// Lengths of the sides of the triangle
	f_t a = std::sqrt(std::pow(x2 - x3, 2) + std::pow(y2 - y3, 2));
	f_t b = std::sqrt(std::pow(x1 - x3, 2) + std::pow(y1 - y3, 2));
	f_t c = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));

	// Semiperimeter of the triangle
	f_t s = (a + b + c) / 2;

	// Area of the triangle using Heron's formula
	f_t A = std::sqrt(s * (s - a) * (s - b) * (s - c));

	// Radius of the circumscribed circle
	f_t R = (a * b * c) / (4 * A);
	return (A > 0) ? R : 0.;
}

template<class f_t>
f_t
radius_of_curvature_theta(f_t theta1, f_t theta2)
{
	f_t x1, y1, x2, y2;
	pol2cart(theta1, x1, y1);
	pol2cart(theta2, x2, y2);
	return radius_of_curvature(x1, y1, x2, y2);
}

template<class f_t>
f_t
radius_of_curvature_theta(f_t r1, f_t theta1, f_t r2, f_t theta2)
{
	f_t x1, y1, x2, y2;
	pol2cart(r1, theta1, x1, y1);
	pol2cart(r2, theta2, x2, y2);
	return radius_of_curvature(x1, y1, x2, y2);
}

};
