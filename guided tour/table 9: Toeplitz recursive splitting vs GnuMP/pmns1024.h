#ifndef PMNS1024_H
#define PMNS1024_H

#include "structs.h"

extern const mpnum __P1024__;
extern const mpnum Gamma1024;
extern const int8_t N1024;
extern const int8_t RHO1024;

void pmns1024_montg_mult(restrict poly res, const restrict poly A,
	const restrict poly B);

#endif
