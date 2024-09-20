#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>

#include "hpmns.h"
#include "bench.h"
// Small hack before a define is applied.
uint8_t ALPHA; uint8_t* pALPHA = &ALPHA;
#include "params.h"

int main(void)
{
	// We set ALPHA (the variable) to ALPHA (the value defined in params.h).
	*pALPHA = ALPHA;
	
	// Initialisation of the randomness.
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	// We do a quick speed test and checks with the parameters to show
	// everything is fine.
	printf("%ld cycles\n", do_bench(pmns_montg_mult, N, RHO, (210/N)*3 + 3));
	printf("Correct: %ld/%d\n", checkPmns(N, RHO, GAMMA, (uint64_t*)__P__, pmns_montg_mult), LOOPCHKNMB);
	return 0;
}
