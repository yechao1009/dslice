#ifndef DSLICE_BASE_H_
#define DSLICE_BASE_H_

#include <vector>
#include <string>

typedef int ds_int;
typedef unsigned ds_uint;
typedef unsigned long ds_ulong;
typedef float ds_float;
typedef double ds_double;

namespace dslice
{

struct Factor {
	std::vector<ds_uint> factor;
	ds_uint level;
};

struct DPTable {
	ds_double score;
	std::vector<ds_uint> strategy;
};

} // namespace

#endif // DSLICE_BASE_H_
