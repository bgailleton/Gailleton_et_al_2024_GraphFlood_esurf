#pragma once

#include "boundary_conditions.hpp"
#include "utils.hpp"

namespace DAGGER {

constexpr std::uint8_t TopLeftMask8 = 0b10000000;
constexpr std::uint8_t TopMask8 = 0b01000000;
constexpr std::uint8_t TopRightMask8 = 0b00100000;
constexpr std::uint8_t LeftMask8 = 0b00010000;
constexpr std::uint8_t RightMask8 = 0b00001000;
constexpr std::uint8_t BottomLeftMask8 = 0b00000100;
constexpr std::uint8_t BottomMask8 = 0b00000010;
constexpr std::uint8_t BottomRightMask8 = 0b00000001;
constexpr std::uint8_t NopeMask8 = 0b00000000;
constexpr std::uint8_t AllMask8 = 0b11111111;

constexpr std::uint8_t IdAdderNone = 0;
constexpr std::uint8_t IdAdderPerTopLeft = 1;
constexpr std::uint8_t IdAdderPerTop = 2;
constexpr std::uint8_t IdAdderPerTopRight = 3;
constexpr std::uint8_t IdAdderPerLeft = 4;
constexpr std::uint8_t IdAdderPerRight = 5;
constexpr std::uint8_t IdAdderPerBottomLeft = 6;
constexpr std::uint8_t IdAdderPerBottom = 7;
constexpr std::uint8_t IdAdderPerBottomRight = 8;

constexpr std::array<std::uint8_t, 8> NeighbourerMask8 = {
	TopLeftMask8, TopMask8,				 TopRightMask8, LeftMask8,
	RightMask8,		BottomLeftMask8, BottomMask8,		BottomRightMask8
};

std::uint8_t
invBits(std::uint8_t nbit)
{
	switch (nbit) {
		case (TopLeftMask8):
			return BottomRightMask8;
		case (TopMask8):
			return BottomMask8;
		case (TopRightMask8):
			return BottomLeftMask8;
		case (LeftMask8):
			return RightMask8;
		case (RightMask8):
			return LeftMask8;
		case (BottomLeftMask8):
			return TopRightMask8;
		case (BottomMask8):
			return TopMask8;
		case (BottomRightMask8):
			return TopLeftMask8;
		default:
			return NopeMask8;
	}
}

std::string
bits2str(std::uint8_t nbit)
{
	switch (nbit) {
		case (TopLeftMask8):
			return "TopLeftMask8";
		case (TopMask8):
			return "TopMask8";
		case (TopRightMask8):
			return "TopRightMask8";
		case (LeftMask8):
			return "LeftMask8";
		case (RightMask8):
			return "RightMask8";
		case (BottomLeftMask8):
			return "BottomLeftMask8";
		case (BottomMask8):
			return "BottomMask8";
		case (BottomRightMask8):
			return "BottomRightMask8";
		case (NopeMask8):
			return "NopeMask8";
		case (AllMask8):
			return "AllMask8";
		default:
			return "ERROR";
	}
}

const int NADDER = 9;

template<class i_t, class f_t>
class lookup8
{
public:
	i_t nx = 0;
	i_t ny = 0;
	i_t nxy = 0;

	std::array<std::array<i_t, 8>, NADDER> adders;
	std::array<f_t, 8> dxer;
	std::array<f_t, 8> dyer;

	std::array<std::int8_t, 8> dxSigner;
	std::array<std::int8_t, 8> dySigner;

	std::array<std::array<i_t, 256>, NADDER> NeighbourerD8;
	std::array<std::array<f_t, 256>, NADDER> NeighbourerD8dx;
	std::array<std::array<f_t, 256>, NADDER> NeighbourerD8dy;
	std::array<std::array<std::int8_t, 256>, NADDER> NeighbourerD8dxSign;
	std::array<std::array<std::int8_t, 256>, NADDER> NeighbourerD8dySign;

	std::array<std::array<std::uint8_t, 256>, NADDER> NeighbourerNN;
	std::array<std::array<std::array<i_t, 8>, 256>, NADDER> Neighbourer;
	std::array<std::array<std::array<std::uint8_t, 8>, 256>, NADDER>
		NeighbourerBits;
	std::array<std::array<std::array<f_t, 8>, 256>, NADDER> Neighbourerdx;
	std::array<std::array<std::array<f_t, 8>, 256>, NADDER> Neighbourerdy;
	std::array<std::array<std::array<std::int8_t, 8>, 256>, NADDER>
		NeighbourerdxSign;
	std::array<std::array<std::array<std::int8_t, 8>, 256>, NADDER>
		NeighbourerdySign;

	// Helps managing the shenanigans behind periodic conditions or other node
	// changing thingies
	std::array<std::array<i_t, 8>, 256> Corrector;

	lookup8() { ; };
	lookup8(i_t nx, i_t ny, f_t dx, f_t dy)
	{

		f_t dxy = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
		this->nx = nx;
		this->ny = ny;
		this->nxy = nx * ny;

		this->adders[0] = { -this->nx - 1, -this->nx, -this->nx + 1, -1, 1,
												this->nx - 1,	 this->nx,	this->nx + 1 }; // IdAdderNone
		this->adders[1] = { this->ny * this->nx - 1,
												(this->ny - 1) * this->nx,
												(this->ny - 1) * this->nx + 1,
												this->nx - 1,
												1,
												2 * this->nx - 1,
												this->nx,
												this->nx + 1 }; // IdAdderPerTopLeft
		this->adders[2] = { (this->ny - 1) * this->nx - 1,
												(this->ny - 1) * this->nx,
												(this->ny - 1) * this->nx + 1,
												-1,
												1,
												this->nx - 1,
												this->nx,
												this->nx + 1 }; // IdAdderPerTop
		this->adders[3] = { (this->ny - 1) * this->nx - 1,
												(this->ny - 1) * this->nx,
												(this->ny - 2) * this->nx + 1,
												-1,
												-this->nx + 1,
												this->nx - 1,
												this->nx,
												1 }; // IdAdderPerTopRight
		this->adders[4] = {
			-1, -this->nx,				-this->nx + 1, (this->nx - 1),
			1,	2 * this->nx - 1, this->nx,			 this->nx + 1
		}; // IdAdderPerLeft
		this->adders[5] = {
			-this->nx - 1, -this->nx,		 -2 * this->nx + 1, -1,
			-this->nx + 1, this->nx - 1, this->nx,					1
		}; // IdAdderPerRight
		this->adders[6] = { -1,
												-this->nx,
												-this->nx + 1,
												this->nx - 1,
												1,
												-(this->ny - 2) * this->nx - 1,
												-(this->ny - 1) * this->nx,
												-(this->ny - 1) * this->nx +
													1 }; // IdAdderPerBottomLeft
		this->adders[7] = { -this->nx - 1,
												-this->nx,
												-this->nx + 1,
												-1,
												1,
												-(this->ny - 1) * this->nx - 1,
												-(this->ny - 1) * this->nx,
												-(this->ny - 1) * this->nx + 1 }; // IdAdderPerBottom
		this->adders[8] = { -this->nx - 1,
												-this->nx,
												-this->nx + 1,
												-1,
												1 - this->nx + 1,
												-(this->ny - 1) * this->nx - 1,
												-(this->ny - 1) * this->nx,
												-(this->ny) * this->nx + 1 }; // IdAdderPerBottomRight

		this->dxer = { dxy, dy, dxy, dx, dx, dxy, dy, dxy };
		this->dyer = { dxy, dx, dxy, dy, dy, dxy, dx, dxy };

		this->dxSigner = { -1, 0, 1, -1, 1, -1, 0, 1 };
		this->dySigner = { 1, 1, 1, 0, 0, -1, -1, -1 };

		this->_compute_lookup_tables();
	};

	void _local_lookup(uint8_t indices,
										 std::array<int, 8>& arr,
										 std::array<f_t, 8>& arrdx,
										 std::array<f_t, 8>& arrdy,
										 std::array<std::int8_t, 8>& arrdxSign,
										 std::array<std::int8_t, 8>& arrdySign,
										 uint8_t& nn,
										 std::array<std::uint8_t, 8>& arrBits,
										 std::uint8_t idAdder)
	{

		// Retrieve the indices specified by the bits set in the `indices` value
		for (uint8_t i = 0; i < 8; ++i) {
			if (indices & (1 << i)) {
				// Index `i` is set, process the corresponding value
				arr[nn] = this->adders[idAdder][7 - i];
				arrdx[nn] = this->dxer[7 - i];
				arrdy[nn] = this->dyer[7 - i];
				arrBits[nn] = NeighbourerMask8[7 - i];
				arrdxSign[nn] = this->dxSigner[7 - i];
				arrdySign[nn] = this->dySigner[7 - i];
				// std::cout << this->adders[7 - i] << " vs " <<
				// bits2str(NeighbourerMask8[7 - i]) << std::endl;; Just a reminder this
				// checked the validity of Neighbourerbits
				++nn;
			}
		}
	}

	void _compute_lookup_tables()
	{
		for (int idAdder = 0; idAdder < NADDER; ++idAdder)
			for (int i = 0; i < 256; ++i) {
				uint8_t ti = static_cast<uint8_t>(i);
				this->NeighbourerNN[idAdder][ti] = 0;
				this->NeighbourerD8[idAdder][ti] = 0;
				for (auto& v : this->Neighbourer[idAdder][ti])
					v = 0;
				this->_local_lookup(ti,
														this->Neighbourer[idAdder][ti],
														this->Neighbourerdx[idAdder][ti],
														this->Neighbourerdy[idAdder][ti],
														this->NeighbourerdxSign[idAdder][ti],
														this->NeighbourerdySign[idAdder][ti],
														this->NeighbourerNN[idAdder][ti],
														this->NeighbourerBits[idAdder][ti],
														idAdder);

				if (this->NeighbourerNN[idAdder][ti] == 1) {

					this->NeighbourerD8[idAdder][ti] = this->Neighbourer[idAdder][ti][0];

					this->NeighbourerD8dx[idAdder][ti] =
						this->Neighbourerdx[idAdder][ti][0];

					this->NeighbourerD8dy[idAdder][ti] =
						this->Neighbourerdy[idAdder][ti][0];

					this->NeighbourerD8dxSign[idAdder][ti] =
						this->NeighbourerdxSign[idAdder][ti][0];

					this->NeighbourerD8dySign[idAdder][ti] =
						this->NeighbourerdySign[idAdder][ti][0];
				}
			}
	}

	std::uint8_t TopLeft_normal_boundary() const
	{
		std::uint8_t bound = 0;
		bound |= RightMask8;
		bound |= BottomMask8;
		bound |= BottomRightMask8;
		return bound;
	}

	std::uint8_t TopRight_normal_boundary() const
	{
		std::uint8_t bound = 0;
		bound |= LeftMask8;
		bound |= BottomMask8;
		bound |= BottomLeftMask8;
		return bound;
	}

	std::uint8_t Top_normal_boundary() const
	{
		std::uint8_t bound = 0;
		bound |= LeftMask8;
		bound |= RightMask8;
		bound |= BottomMask8;
		bound |= BottomLeftMask8;
		bound |= BottomRightMask8;
		return bound;
	}

	std::uint8_t BottomLeft_normal_boundary() const
	{
		std::uint8_t bound = 0;
		bound |= RightMask8;
		bound |= TopMask8;
		bound |= TopRightMask8;
		return bound;
	}

	std::uint8_t BottomRight_normal_boundary() const
	{
		std::uint8_t bound = 0;
		bound |= LeftMask8;
		bound |= TopMask8;
		bound |= TopLeftMask8;
		return bound;
	}

	std::uint8_t Bottom_normal_boundary() const
	{
		std::uint8_t bound = 0;
		bound |= LeftMask8;
		bound |= RightMask8;
		bound |= TopMask8;
		bound |= TopLeftMask8;
		bound |= TopRightMask8;
		return bound;
	}

	std::uint8_t Left_normal_boundary() const
	{
		std::uint8_t bound = 0;
		bound |= RightMask8;
		bound |= TopMask8;
		bound |= BottomMask8;
		bound |= TopRightMask8;
		bound |= BottomRightMask8;
		return bound;
	}

	std::uint8_t Right_normal_boundary() const
	{
		std::uint8_t bound = 0;
		bound |= LeftMask8;
		bound |= TopMask8;
		bound |= BottomMask8;
		bound |= TopLeftMask8;
		bound |= BottomLeftMask8;
		return bound;
	}

	std::uint8_t BC2idAdder(int i, BC tbc)
	{
		if (tbc != BC::PERIODIC_BORDER)
			return 0;

		else if (i == 0)
			return 1;
		else if (i == this->nxy - 1)
			return 8;
		else if (i == this->nx)
			return 3;
		else if (i == this->nxy - this->nx)
			return 6;
		else if (i < this->nx)
			return 2;
		else if (i > this->nxy - this->nx)
			return 7;
		else if (i % this->nx == 0)
			return 4;
		else if (i % this->nx == this->nx - 1)
			return 5;
		return 0;
	}

	// Keeping legacy things here in case - 03/2024

	// void _local_lookup_archives(uint8_t indices,
	// 									 std::array<int, 8>& arr,
	// 									 std::array<f_t, 8>& arrdx,
	// 									 uint8_t& nn,
	// 									 std::array<std::uint8_t, 8>& arrBits,
	// 									 )
	// {

	// 	// Retrieve the indices specified by the bits set in the `indices` value
	// 	for (uint8_t i = 0; i < 8; ++i) {
	// 		if (indices & (1 << i)) {
	// 			// Index `i` is set, process the corresponding value
	// 			arr[nn] = this->adders[7 - i];
	// 			arrdx[nn] = this->dxer[7 - i];
	// 			arrBits[nn] = NeighbourerMask8[7 - i];
	// 			// std::cout << this->adders[7 - i] << " vs " <<
	// bits2str(NeighbourerMask8[7 - i]) << std::endl;; Just a reminder this
	// checked the validity of Neighbourerbits
	// 			++nn;
	// 		}
	// 	}
	// }

}; // end of class lookup8

} // end of namespace
