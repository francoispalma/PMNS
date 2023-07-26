#ifndef PMNS_PARAMS_H_INCLUDED
#define PMNS_PARAMS_H_INCLUDED

#define RHO 53
#define N 5
#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 5")
#define LAMBDA 19
#define M_or_B_is_B

static const int64_t B[N][N] = {
{-2251799813685248, 1, 0, 0, 0}, {0, -2251799813685248, 1, 0, 0}, {0, 0, -2251799813685248, 1, 0}, {0, 0, 0, -2251799813685248, 1}, {-19, 0, 0, 0, 2251799813685248}
		};
static const int64_t B1[N][N] = {
{0u, 0u, 0u, -5825406118003736576u, -8737931403336103397u}, {-1u, 0u, 0u, 0u, 5825406118003736576u}, {-2251799813685248u, -1u, 0u, 0u, 0u}, {0u, -2251799813685248u, -1u, 0u, 0u}, {0u, 0u, -2251799813685248u, -1u, 0u}
		};

static const int64_t __Pi__[N][N] = {
		{19, 0, 67108864, 0, -2251799813685248},
		{19, 0, 0, 2251800082120704, -2251799813685249},
		{2251799813685267, -1, 0, 2251799813685248, -2251798739943425},
		{81604378643, 0, 0, 2251799813685248, -2251799813685249},
		{2251799813685267, 326417514495, 0, 2251799813685248, -2251799813685249}
	};


static _poly __theta__ = { .deg = 5,
	.t = (int64_t[]) { 0x13, 0x2000, 0x0, 0x0, -0x8000000000000 } };
static _mpnum __P__ = { .deg = 4,
		.sign = 1,
		.t = (uint64_t[]) {0xffffffffffffffed, 0xffffffffffffffff, 0xffffffffffffffff, 0x7fffffffffffffff} },
	Gi[] = { { .deg = 1,
		.sign = 1,
		.t = (uint64_t[]) {0x8000000000000} },
	{ .deg = 2,
		.sign = 1,
		.t = (uint64_t[]) {0x0, 0x4000000000} },
	{ .deg = 3,
		.sign = 1,
		.t = (uint64_t[]) {0x0, 0x0, 0x2000000} },
	{ .deg = 4,
		.sign = 1,
		.t = (uint64_t[]) {0x0, 0x0, 0x0, 0x1000} }};
#endif
