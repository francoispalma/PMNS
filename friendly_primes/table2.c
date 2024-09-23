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
	uint64_t gmpcycles[6], optcycles[6], lincycles[6], spacycles[6];
	
	const uint64_t nbrepet = 303; // Must be multiple of 3.
	
	mp_limb_t hp255n5[] = {0xffffffffffffffed, 0xffffffffffffffff, 0xffffffffffffffff, 0x7fffffffffffffff, 0};
	setpmnsparams(5, 2251799813685248, 19);
	setpmnssparseinv(5825406118003736576u, 5825406118003736576u, 9708812670373448219u);
	cycles = do_bench(pmns_doublesparse_mult, 5, 2251799813685266, nbrepet);
	printf("DoubleSparse 255: %ld cycles\n", cycles);
	spacycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(5, 2251799813685266, 2251799813685248, hp255n5, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp383n7[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0xffffffffffffffff, 0xaef6865effffffff, 0x201cb4d65c63dd8c, 0x40006d54c7279d58, 0};
	setpmnsparams(7, 26769392989634560, 3);
	setpmnssparseinv(12306752513469579264u, 12306752513469579264u, 12297829382473034411u);
	cycles = do_bench(pmns_doublesparse_mult, 7, 26769392989634562, nbrepet);
	printf("DoubleSparse 383: %ld cycles\n", cycles);
	spacycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(7, 26769392989634562, 26769392989634560, hp383n7, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp414n7[] = {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x55f9223fffffffff, 0xfe40858d52f68562, 0xe76d37c25fe27a3f, 0x20000194};
	setpmnsparams(7, 636464340236500992, 2);
	setpmnspsi(2);
	setpmnssparseinv(318232170118250496, 636464340236500992, 1);
	cycles = do_bench(pmns_doublesparse_mult, 7, 636464340236500993, nbrepet);
	printf("DoubleSparse 414: %ld cycles\n", cycles);
	spacycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(7, 636464340236500993, 636464340236500992, hp414n7, pmns_doublesparse_mult), LOOPCHKNMB);
	
	/*mp_limb_t hp448n8[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xb9a1e9386a15f060, 0x630c235bfe9221fe, 0x800154a7dd9cb573, 0};
	setpmnsparams(8, 66077440488767488, 3);
	setpmnssparseinv(12319855195969290240u, 12319855195969290240u, 12297829382473034411u);
	cycles = do_bench(pmns_doublesparse_mult, 8, 66077440488767490, nbrepet);
	printf("DoubleSparse 448: %ld cycles\n", cycles);
	spacycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(8, 66077440488767490, 66077440488767488, hp448n8, pmns_doublesparse_mult), LOOPCHKNMB);*/
	
	mp_limb_t hp512n9[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x98fac9ffffffffff,0xee846d77242cfe8f, 0xa545f31d9bad3565, 0x40002ce4a6f2f8a8, 0};
	setpmnsparams(9, 114384436610465792, 3);
	setpmnssparseinv(6187042836773339136u, 6187042836773339136u, 12297829382473034411u);
	cycles = do_bench(pmns_doublesparse_mult9, 9, 114384436610465794, nbrepet);
	printf("DoubleSparse 512: %ld cycles\n", cycles);
	spacycles[4] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 114384436610465794, 114384436610465792, hp512n9, pmns_doublesparse_mult9), LOOPCHKNMB);
	
	mp_limb_t hp521n9[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x554d91ffffffffff,0xccdff61d648cf993, 0x3e27702a6e2f8a3e, 0x0478c452724bdf49, 0x100};
	setpmnsparams(9, 247085628838117376, 3);
	setpmnssparseinv(6231276567515889664u, 6231276567515889664u, 12297829382473034411u);
	cycles = do_bench(pmns_doublesparse_mult9, 9, 247085628838117378, nbrepet);
	printf("DoubleSparse 521: %ld cycles\n", cycles);
	spacycles[5] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 247085628838117378, 247085628838117376, hp521n9, pmns_doublesparse_mult9), LOOPCHKNMB);
	
	mp_limb_t p255n5[] = {0xa80000000003cbfd, 0x0010e0000000000c, 0x0000000b40000000, 0x800000000003c000, 0};
	setpmnsparams(5, 2251799813685260, 3);
	setpmnslastcol((uint64_t[]) {10589719013596916053u, 3328969402090188796u, 12268509515263197136u, 17986819276134612416u,11630609810127578368u,});
	cycles = do_bench(pmns_linearred_mult, 5, 2251799813685262, nbrepet);
	printf("LinearRed 256: %ld cycles\n", cycles);
	lincycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(5, 2251799813685262, 2251799813685260, p255n5, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p383n7[] = {0x827c835db704d27d, 0xfa496d78dfdbb88e, 0xcb892df02209b1d9, 0x5be7c452cf7ab6d1, 0x0d03f447f9a6c8b0, 0x4000000000000915, 0};
	setpmnsparams(7, 26769293307327386, 3);
	setpmnslastcol((uint64_t[]) {8298466069282600661u, 8465212809323277090u, 5194272462830712948u, 4795752030772072904u, 14765568818636884560u, 9515338599022666784u, 12421683519646735168u,});
	cycles = do_bench(pmns_linearred_mult, 7, 26769293307327388, nbrepet);
	printf("LinearRed 383: %ld cycles\n", cycles);
	lincycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(7, 26769293307327388, 26769293307327386, p383n7, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p414n7[] = {0xf8029d845f4d10ab, 0x36400119e6e4c221, 0x06e6000032de98eb, 0x8719700000051986, 0x02d586800000004e, 0x00000e8c00000000, 0x20000000};
	setpmnsparams(7, 576460752303423621, 2);
	setpmnslastcol((uint64_t[]) {2215446049233201667u, 1235801625573266319u, 6987830822409808203u, 13385363970364568567u, 4177830211639852883u, 13201850230578931231u, 2828932912287026715u,});
	cycles = do_bench(pmns_linearred_mult, 7, 576460752303423622, nbrepet);
	printf("LinearRed 414: %ld cycles\n", cycles);
	lincycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(7, 576460752303423622, 576460752303423621, p414n7, pmns_linearred_mult), LOOPCHKNMB);
	
	/*mp_limb_t p448n8[] = {0x624bc8e9c067a0fd, 0xcbf4b2c63b00e317, 0xa0cd9555210a63e0, 0x105354ed7cfe644b, 0xe18007f8a9ad382b, 0x0d9bd8f75e3b4bc5, 0x800000000005776a, 0};
	setpmnsparams(8, 66077105076381050, 3);
	setpmnslastcol((uint64_t[]) {3546674433909423189u, 17272956169590049154u, 4401242359300837876u, 11803910513052176968u, 9760500218265968208u, 13990977856189956640u, 9301119457689969472u, 11609398705935355008u,});
	cycles = do_bench(pmns_linearred_mult, 8, 66077105076381052, nbrepet);
	printf("LinearRed 448: %ld cycles\n", cycles);
	lincycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(8, 66077105076381052, 66077105076381050, p448n8, pmns_linearred_mult), LOOPCHKNMB);*/
	
	mp_limb_t p512n9[] = {0x21cb2f3997fffffd, 0x8dc76bf1f3d1c3ee, 0xb41bfd602771fafd, 0x69772bbf38af18b3, 0x651fe55cbafdeb16, 0x7c8ce4e9c47f9cf8, 0x80d079a56c657e96, 0x8000000000034f27, 0};
	setpmnsparams(9, 123541877815766680, 3);
	setpmnslastcol((uint64_t[]) {1068504517212198229u, 9678330455461638264u, 12298442381156439872u, 9543643910875729408u, 12539807914892349440u, 57035926836051968u, 11747439154570723328u, 15471119491003318272u, 16325372834547236864u,});
	cycles = do_bench(pmns_linearred_mult9, 9, 123541877815766682, nbrepet);
	printf("LinearRed 512: %ld cycles\n", cycles);
	lincycles[4] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 123541877815766682, 123541877815766680, p512n9, pmns_linearred_mult9), LOOPCHKNMB);
	
	mp_limb_t p521n9[] = {0x991e4a2e727472b5, 0x0a363db1d581f52b, 0xcfd4043ed07daf18, 0x3e72292f77fb29cf, 0x7989b43a1b6564cf,0x1dad3dd0c483c688, 0xa47245c62fbe21c4, 0x0000000018682577, 0x100};
	setpmnsparams(9, 247083755631535095, 2);
	setpmnslastcol((uint64_t[]) {9646799981064252317u, 4868595640826006395u, 18224743689174066349u, 978733964884756971u, 721057044947481277u, 16319187962034555739u, 6281386887843620301u, 314023437206896587u, 12739390295029747165u, });
	cycles = do_bench(pmns_linearred_mult9, 9, 247083755631535096, nbrepet);
	printf("LinearRed 521: %ld cycles\n", cycles);
	lincycles[5] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 247083755631535096, 247083755631535095, p521n9, pmns_linearred_mult9), LOOPCHKNMB);
	
	cycles = do_pmersbench(multMod25519, nbrepet, 5, C25519convert);
	printf("C25519: %ld cycles\n", cycles);
	optcycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(4, 0x8000000000000000, 19, 5, C25519convert, multMod25519), LOOPCHKNMB);
	
	cycles = do_pmersbench(multModM383, nbrepet, 7, M383convert);
	printf("M-383: %ld cycles\n", cycles);
	optcycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(6, 0x8000000000000000, 187, 7, M383convert, multModM383), LOOPCHKNMB);
	
	cycles = do_pmersbench(multModC41417, nbrepet, 7, C41417convert);
	printf("Curve41417: %ld cycles\n", cycles);
	optcycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(7, 0x40000000, 17, 7, C41417convert, multModC41417), LOOPCHKNMB);
	
	/*cycles = do_pmersbench(multModEd448, nbrepet, 8, Ed448convert);
	printf("Ed448-Goldilocks: %ld cycles\n", cycles);
	optcycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(0, 0, 0, 0, 0, multModEd448), LOOPCHKNMB);*/
	
	cycles = do_pmersbench(multModM511, nbrepet, 9, M511convert);
	printf("M-511: %ld cycles\n", cycles);
	optcycles[4] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(8, 0x8000000000000000, 187, 9, M511convert, multModM511), LOOPCHKNMB);
	
	cycles = do_pmersbench(multModE521, nbrepet, 9, E521convert);
	printf("E-521: %ld cycles\n", cycles);
	optcycles[5] = cycles;
	printf("Correct: %ld/%d\n", checkOpti(9, 0x200, 1, 9, E521convert, multModE521), LOOPCHKNMB);
	
	cycles = do_gmpbench(gmpmulmod2k, 19, 63, 4, nbrepet);
	printf("gmp 2^255-19: %ld cycles\n", cycles);
	gmpcycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(19, 63, 4), LOOPCHKNMB);
	
	cycles = do_gmpbench(gmpmulmod2k, 187, 63, 6, nbrepet);
	printf("gmp 2^383-187: %ld cycles\n", cycles);
	gmpcycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(187, 63, 6), LOOPCHKNMB);
	
	cycles = do_gmpbench(gmpmulmod2k, 17, 30, 7, nbrepet);
	printf("gmp 2^414-17: %ld cycles\n", cycles);
	gmpcycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(17, 30, 7), LOOPCHKNMB);
	
	cycles = do_gmpbench(gmpmulmod2k, 187, 63, 8, nbrepet);
	printf("gmp 2^511-187: %ld cycles\n", cycles);
	gmpcycles[4] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(187, 63, 8), LOOPCHKNMB);
	
	cycles = do_mersennebench(mersenne521, nbrepet);
	printf("gmp 2^521-1: %ld cycles\n", cycles);
	gmpcycles[5] = cycles;
	printf("Correct: %ld/%d\n", checkmodmul(1, 9, 9), LOOPCHKNMB);
	
	/*printf("\n==============================================================================\n");
	printf("|   Method \\ Prime size    |  255   |  383  |  414   |  448  |  511  |  521  |\n");
	printf("==============================================================================\n");
	printf("|   Corresponding Curve    | C25519 | M-383 | C41417 | Ed448 | M-511 | E-521 |\n");
	printf("==============================================================================\n");
	printf("| Optimized Multiprecision |   %lu   |  %lu  |  %lu   |  %lu  |  %lu  |  %lu  |\n", optcycles[0], optcycles[1], optcycles[2], optcycles[3], optcycles[4], optcycles[5]);
	printf("==============================================================================\n");
	printf("|    DoubleSparse PMNS     |   %lu   |  %lu  |  %lu   |  %lu  |  %lu  |  %lu  |\n", spacycles[0], spacycles[1], spacycles[2], spacycles[3], spacycles[4], spacycles[5]);
	printf("==============================================================================\n");
	printf("|      LinearRed PMNS      |   %lu   |  %lu  |  %lu   |  %lu  |  %lu  |  %lu  |\n", lincycles[0], lincycles[1], lincycles[2], lincycles[3], lincycles[4], lincycles[5]);
	printf("==============================================================================\n");
	printf("|  Adapted GMP low level   |   %lu   |  %lu  |  %lu   |  N/A  |  %lu  |  %lu  |\n", gmpcycles[0], gmpcycles[1], gmpcycles[2], gmpcycles[4], gmpcycles[5]);
	printf("==============================================================================\n\n");*/
	
	printf("\n======================================================================\n");
	printf("|   Method \\ Prime size    |  255   |  383  |  414   |  511  |  521  |\n");
	printf("======================================================================\n");
	printf("|   Corresponding Curve    | C25519 | M-383 | C41417 | M-511 | E-521 |\n");
	printf("======================================================================\n");
	printf("| Optimized Multiprecision |   %lu   |  %lu  |  %lu   |  %lu  |  %lu  |\n", optcycles[0], optcycles[1], optcycles[2], optcycles[4], optcycles[5]);
	printf("======================================================================\n");
	printf("|    DoubleSparse PMNS     |   %lu   |  %lu  |  %lu   |  %lu  |  %lu  |\n", spacycles[0], spacycles[1], spacycles[2], spacycles[4], spacycles[5]);
	printf("======================================================================\n");
	printf("|      LinearRed PMNS      |   %lu   |  %lu  |  %lu   |  %lu  |  %lu  |\n", lincycles[0], lincycles[1], lincycles[2], lincycles[4], lincycles[5]);
	printf("======================================================================\n");
	printf("|  Adapted GMP low level   |   %lu   |  %lu  |  %lu   |  %lu  |  %lu  |\n", gmpcycles[0], gmpcycles[1], gmpcycles[2], gmpcycles[4], gmpcycles[5]);
	printf("======================================================================\n\n");
	
	return 0;
}
