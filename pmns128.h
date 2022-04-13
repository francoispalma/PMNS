#ifndef MPPMNS_H
#define MPPMNS_H

#include <stdint.h>
#include <stdlib.h>

#include "structs.h"
#include "params128.h"

extern void amns128_montg_mult(restrict poly128 res, const restrict poly128 A,
	const restrict poly128 B);

#endif
