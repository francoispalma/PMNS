#ifndef PMNS2048_H
#define PMNS2048_H

#include "structs.h"

extern const mpnum __P2048__;
extern const mpnum Gamma2048;
extern const int8_t N2048;
extern const int8_t RHO2048;

void pmns2048_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif
