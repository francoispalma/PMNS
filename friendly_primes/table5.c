#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

#include "structs.h"
#include "eccoptimizedcode.h"
#include "pmns.h"
#include "bench.h"

int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	uint64_t cycles;
	uint64_t gmpcycles[4], optcycles[4], lincycles[4], spacycles[4];
	
	const uint64_t nbrepet = 303; // Must be multiple of 3.
	
	mp_limb_t hp244n4[] = {1, 0, 0x47d060f4c6757810, 0x800003055e858};
	setpmnsparams(4, 1938975295155470336, -1);
	setpmnssparseinv(16507768778554081280u, 16507768778554081280u, 18446744073709551615u);
	cycles = do_bench(pmns_doublesparse_mult, 4, 1938975295155470336, nbrepet);
	printf("DoubleSparse 244: %ld cycles\n", cycles);
	spacycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(4, 1938975295155470336, 1938975295155470336, hp244n4, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp297n5[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0x47514278ffffffff, 0x5a89484e0428b4f8, 0x10000743ce3};
	setpmnsparams(5, 662180435446464512, 3);
	setpmnssparseinv(6369641503052005376u, 6369641503052005376u, 12297829382473034411u);
	cycles = do_bench(pmns_doublesparse_mult, 5, 662180435446464514, nbrepet);
	printf("DoubleSparse 297: %ld cycles\n", cycles);
	spacycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(5, 662180435446464514, 662180435446464512, hp297n5, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp354n6[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0xffffffffffffffff, 0x584f3c6ca5098fff, 0xdca7f94a201803d3, 0x200001120};
	setpmnsparams(6, 513568188978429952, 3);
	setpmnssparseinv(12469018778799177728u, 12469018778799177728u, 12297829382473034411u);
	cycles = do_bench(pmns_doublesparse_mult, 6, 513568188978429954, nbrepet);
	printf("DoubleSparse 354: %ld cycles\n", cycles);
	spacycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(6, 513568188978429954, 513568188978429952, hp354n6, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp480n8[] = {1, 0, 0, 0, 0xc2f7a2b0808c60a1, 0x2c2944875584fbc4, 0xb3e2b86282d7ba82, 0x80000e50};
	setpmnsparams(8, 1057233906744426496, -1);
	setpmnssparseinv(17389510166965125120u, 17389510166965125120u, -1);
	cycles = do_bench(pmns_doublesparse_mult, 8, 1057233906744426496, nbrepet);
	printf("DoubleSparse 480: %ld cycles\n", cycles);
	spacycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(8, 1057233906744426496, 1057233906744426496, hp480n8, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t p244n4[] = {0x80c8c0d7f6ffe511, 0x4492c8088e27e042, 0x96b8c0f8a5d6900c, 0x7ffffffffffff};
	setpmnsparams(4, 1938975120585633030, -1);
	setpmnslastcol((uint64_t[]) {14451576143497112561u, 13568497698566497446u, 13198361674973153764u, 17968173778106601304u});
	cycles = do_bench(pmns_linearred_mult, 4, 1938975120585633030, nbrepet);
	printf("LinearRed 244: %ld cycles\n", cycles);
	lincycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(4, 1938975120585633030, 1938975120585633030, p244n4, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p297n5[] = {0x2327e9d145df11dd, 0x115d3a81bc94c1dc, 0xfa481a70bf496fab, 0x002cfbfcda0339c8, 0x10000000000};
	setpmnsparams(5, 662179517891295902, 3);
	setpmnslastcol((uint64_t[]) {14642749673451933301u, 14583214676991209014u, 8666214240549945684u, 17971678915780893144u, 3853751881226281808u, });
	cycles = do_bench(pmns_linearred_mult, 5, 662179517891295904, nbrepet);
	printf("LinearRed 297: %ld cycles\n", cycles);
	lincycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(5, 662179517891295904, 662179517891295902, p297n5, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p354n6[] = {0xcd2ce1d750150ffd, 0x35c8593e965620ed, 0x80340d8de45e2899, 0x1bc1ad65de681e75, 0x00001ee19e75f62d, 0x200000000};
	setpmnsparams(6, 513568145285335652, 3);
	setpmnslastcol((uint64_t[]) {5525218229414577493u, 4425439584336717620u, 4166206055049746512u, 11237167241322667840u, 8317726432053343488u, 8193909296136762368u, });
	cycles = do_bench(pmns_linearred_mult, 6, 513568145285335654, nbrepet);
	printf("LinearRed 354: %ld cycles\n", cycles);
	lincycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(6, 513568145285335654, 513568145285335652, p354n6, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p480n8[] = {0xb4c1893968210001, 0x64b399583106c723, 0x90ee97e2cd55096c, 0xa8c27fb761fbabc5, 0xb1c9e6cffb98311f, 0xd8e3e48edcf32427, 0x00000f8f595e458e, 0x80000000};
	setpmnsparams(8, 1057233681222091724, -1);
	setpmnslastcol((uint64_t[]) {7055534262214131713u, 6788977602221982668u, 11805067379222014608u, 18284721751823260352u, 490289034654617856u, 8619646048256297984u, 158290366750756864u, 6652374198690168832u, });
	cycles = do_bench(pmns_linearred_mult, 8, 1057233681222091724, nbrepet);
	printf("LinearRed 480: %ld cycles\n", cycles);
	lincycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(8, 1057233681222091724, 1057233681222091724, p480n8, pmns_linearred_mult), LOOPCHKNMB);
	
	cycles = do_pmersbench(multMod244, nbrepet, 5, PoC244convert);
	printf("2^244 - 189: %ld cycles\n", cycles);
	optcycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(4, 0x10000000000000, 189, 5, PoC244convert, multMod244), LOOPCHKNMB);
	
	cycles = do_pmersbench(multMod297, nbrepet, 5, PoC297convert);
	printf("2^297 - 123: %ld cycles\n", cycles);
	optcycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(5, 0x20000000000, 123, 6, PoC297convert, multMod297), LOOPCHKNMB);
	
	cycles = do_pmersbench(multMod354, nbrepet, 6, PoC354convert);
	printf("2^354 - 153: %ld cycles\n", cycles);
	optcycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(6, 0x400000000, 153, 7, PoC354convert, multMod354), LOOPCHKNMB);
	
	cycles = do_pmersbench(multMod480, nbrepet, 9, PoC480convert);
	printf("2^480 - 49: %ld cycles\n", cycles);
	optcycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(8, 0x100000000, 49, 9, PoC480convert, multMod480), LOOPCHKNMB);
	
	cycles = do_gmpbench(gmpmulmod2k, 189, 52, 4, nbrepet);
	printf("gmpmulmod2k 2^244-189: %ld cycles\n", cycles);
	gmpcycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(189, 52, 4), LOOPCHKNMB);
	
	cycles = do_gmpbench(gmpmulmod2k, 123, 41, 5, nbrepet);
	printf("gmpmulmod2k 2^297-123: %ld cycles\n", cycles);
	gmpcycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(123, 41, 5), LOOPCHKNMB);
	
	cycles = do_gmpbench(gmpmulmod2k, 153, 34, 6, nbrepet);
	printf("gmpmulmod2k 2^354-153: %ld cycles\n", cycles);
	gmpcycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(153, 34, 6), LOOPCHKNMB);
	
	cycles = do_gmpbench(gmpmulmod2k, 47, 32, 8, nbrepet);
	printf("gmpmulmod2k 2^480-47: %ld cycles\n", cycles);
	gmpcycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(47, 32, 8), LOOPCHKNMB);
	
	printf("\n===============================================================================\n");
	printf("|   Method \\ Prime size    |    244     |    297     |    354     |    480    |\n");
	printf("===============================================================================\n");
	printf("| Corresponding PseudoMers | 2²⁴⁴ - 189 | 2²⁹⁷ - 123 | 2³⁵⁴ - 153 | 2⁴⁸⁰ - 47 |\n");
	printf("===============================================================================\n");
	printf("| Optimized Multiprecision |     %lu     |     %lu     |    %lu     |    %lu    |\n", optcycles[0], optcycles[1], optcycles[2], optcycles[3]);
	printf("===============================================================================\n");
	printf("|    DoubleSparse PMNS     |     %lu     |     %lu     |     %lu     |    %lu    |\n", spacycles[0], spacycles[1], spacycles[2], spacycles[3]);
	printf("===============================================================================\n");
	printf("|      LinearRed PMNS      |     %lu     |     %lu     |     %lu     |    %lu    |\n", lincycles[0], lincycles[1], lincycles[2], lincycles[3]);
	printf("===============================================================================\n");
	printf("|  Adapted GMP low level   |     %lu     |    %lu     |    %lu     |    %lu    |\n", gmpcycles[0], gmpcycles[1], gmpcycles[2], gmpcycles[3]);
	printf("===============================================================================\n\n");*/
	
	return 0;
}
