from enum import Enum


class Flooder(Enum):
	SFD_STATIC = 0
	MFD_STATIC = 1
	SFD_DYNAMIC = 2
	MFD_DYNAMIC = 3
	CAESAR_LS = 4
	CAESAR_LS_OMP = 5
	FLOODOS = 6


class Topology(Enum):
	D8 = 1
	D4 = 2