#ifndef PMNS1664_H
#define PMNS1664_H

#include "structs.h"

extern const mpnum __P1664__;
extern const mpnum Gamma1664;
extern const int8_t N1664;
extern const int8_t RHO1664;

void pmns1664_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif
