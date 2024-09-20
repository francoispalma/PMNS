#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "gmp.h"
#include "structs.h"
#include "pmns.h"
#include "bench.h"

int main(void)
{
	time_t seed;
	srand((unsigned) (time(&seed)));
	
	uint64_t cycles;
	uint64_t lincycles[6], spacycles[6];
	
	const uint64_t nbrepet = 303; // Must be multiple of 3.
	
	mp_limb_t hp512_X9M2[] = {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x9b343aac0acb4047, 0x000b27db246be8a6, 0x8000c23ba6963b03, 0};
	setpmnsparams(9, 133432951575674880, 2);
	setpmnspsi(2);
	setpmnssparseinv(66716475787837440u, 133432951575674880u, 1);
	cycles = do_bench(pmns_doublesparse_mult9, 9, 133432951575674881, nbrepet);
	printf("DoubleSparse 512 X⁹ - 2: %ld cycles\n", cycles);
	spacycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 133432951575674881, 133432951575674880, hp512_X9M2, pmns_doublesparse_mult9), LOOPCHKNMB);
	
	mp_limb_t hp512_X9M3[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xee263ddaffffffff, 0x452ea3870a328c2b, 0x2840f84fe382d24f, 0x80016fa5bbedcad7, 0};
	setpmnsparams(9, 123542479411609600, 3);
	setpmnssparseinv(12339010208943570944u, 12339010208943570944u, 12297829382473034411u);
	cycles = do_bench(pmns_doublesparse_mult9, 9, 123542479411609602, nbrepet);
	printf("DoubleSparse 512 X⁹ - 3: %ld cycles\n", cycles);
	spacycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 123542479411609602, 123542479411609600, hp512_X9M3, pmns_doublesparse_mult9), LOOPCHKNMB);
	
	mp_limb_t hp512_2X9M1[] = {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xdcb3b885ffffffff, 0x53b08fc89b99a4a3, 0xa878d70f8e09afe3, 0x800156175ff84afc, 0};
	setpmnsparams(9, 114384818862555136, 1);
	setpmnsaux(228769637725110272, 1);
	setpmnssparseinv(228769637725110272u, 114384818862555136u, 1);
	setpmnsalpha(2);
	cycles = do_bench(pmns_doublesparse_multa9, 9, 228769637725110272, nbrepet);
	printf("DoubleSparse 512 2X⁹ - 1: %ld cycles\n", cycles);
	spacycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 228769637725110272, 114384818862555136, hp512_2X9M1, pmns_doublesparse_multa9), LOOPCHKNMB);
	
	mp_limb_t hp512_2X9M3[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xd33a45a5ffffffff, 0xcca8e32b540ac606, 0x20b42a7d01c0f72f, 0x80018373418efc19, 0};
	setpmnsparams(9, 114384887582031872, 3);
	setpmnsaux(228769775164063744, 3);
	setpmnssparseinv(12374085974194388992u, 6187042987097194496u, 12297829382473034411u);
	setpmnsalpha(2);
	cycles = do_bench(pmns_doublesparse_multa9, 9, 228769775164063744, nbrepet);
	printf("DoubleSparse 512 2X⁹ - 3: %ld cycles\n", cycles);
	spacycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 228769775164063744, 114384887582031872, hp512_2X9M3, pmns_doublesparse_multa9), LOOPCHKNMB);
	
	mp_limb_t hp512_3X9M1[] = {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xa8245412ffffffff, 0xc6b02c2ade6cb147, 0xaef020931eca1a8f, 0x80003e214ba8a365, 0};
	setpmnsparams(9, 109345542324092928, 1);
	setpmnsaux(328036626972278784, 1);
	setpmnssparseinv(328036626972278784u, 109345542324092928u, 1);
	setpmnsalpha(3);
	cycles = do_bench(pmns_doublesparse_multa9, 9, 328036626972278784, nbrepet);
	printf("DoubleSparse 512 3X⁹ - 1: %ld cycles\n", cycles);
	spacycles[4] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 328036626972278784, 109345542324092928, hp512_3X9M1, pmns_doublesparse_multa9), LOOPCHKNMB);
	
	mp_limb_t hp512_3X9M2[] = {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x09e2a4ffffffffff, 0x4c96c02d8314a2ae, 0x7c6c958c5ccb4c03, 0x800ad910d5ae0672, 0};
	setpmnsparams(9, 118103964188147712, 2);
	setpmnsaux(177155946282221568, 1);
	setpmnssparseinv(177155946282221568, 118103964188147712, 1);
	setpmnsalpha(3);
	cycles = do_bench(pmns_doublesparse_multa9, 9, 177155946282221568, nbrepet);
	printf("DoubleSparse 512 3X⁹ - 2: %ld cycles\n", cycles);
	spacycles[5] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 177155946282221568, 118103964188147712, hp512_3X9M2, pmns_doublesparse_multa9), LOOPCHKNMB);
	setpmnsalpha(1);
	
	mp_limb_t p512_X9M2[] = {0x05a0aace573cf77b, 0xf2d1cd4f6a93d3c3, 0x3d31400e82ce18dc, 0x60e9ebf72f1ede20, 0xf79302b94bca49d8, 0x266b5985806371a0, 0x7967d16413acb486, 0x8000000000001460, 0};
	setpmnsparams(9, 123541877815766365, 2);
	setpmnslastcol((uint64_t[]) {9975839623651282867u, 2003220533415160583u, 17953570065681390731u, 733836773084401023u, 11528432699608594979u, 715211991917282743u, 13602203582099668603u, 5029592706098270639u, 4152759214268772243u,});
	cycles = do_bench(pmns_linearred_mult9, 9, 123541877815766366, nbrepet);
	printf("LinearRed 512 X⁹ - 2: %ld cycles\n", cycles);
	lincycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 123541877815766366, 123541877815766365, p512_X9M2, pmns_linearred_mult9), LOOPCHKNMB);
	
	mp_limb_t p512_X9M3[] = {0x21cb2f3997fffffd, 0x8dc76bf1f3d1c3ee, 0xb41bfd602771fafd, 0x69772bbf38af18b3, 0x651fe55cbafdeb16, 0x7c8ce4e9c47f9cf8, 0x80d079a56c657e96, 0x8000000000034f27, 0};
	setpmnsparams(9, 123541877815766680, 3);
	setpmnslastcol((uint64_t[]) {1068504517212198229u, 9678330455461638264u, 12298442381156439872u, 9543643910875729408u, 12539807914892349440u, 57035926836051968u, 11747439154570723328u, 15471119491003318272u, 16325372834547236864u, });
	cycles = do_bench(pmns_linearred_mult9, 9, 123541877815766682, nbrepet);
	printf("LinearRed 512 X⁹ - 3: %ld cycles\n", cycles);
	lincycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 123541877815766682, 123541877815766680, p512_X9M3, pmns_linearred_mult9), LOOPCHKNMB);
	
	mp_limb_t p512_2X9M1[] = {0xd7cff2eb73cdbf21, 0xe8de239a3611048f, 0x5711599d95e7611f, 0xd168832df143ba5a, 0x0fba1e30b01ee062, 0x153fd98679015b29, 0x48586b77b39a1965, 0x7ffffffffffdf8c0, 0,};
	setpmnsparams(9, 114384300578104081, 1);
	setpmnsaux(228768601156208162, 1);
	setpmnslastcol((uint64_t[]) {4118677855139693793u, 93625479251121649u, 14728057061835302913u, 1307188972639079185u, 11516234695479693089u, 2595788364181521457u, 4033684542809930305u, 12506448047005286737u, 11529725351497148769u, });
	setpmnsalpha(2);
	cycles = do_bench(pmns_linearred_multa9, 9, 228768601156208162, nbrepet);
	printf("LinearRed 512 2X⁹ - 1: %ld cycles\n", cycles);
	lincycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 228768601156208162, 114384300578104081, p512_2X9M1, pmns_linearred_multa9), LOOPCHKNMB);
	
	mp_limb_t p512_2X9M3[] = {0x1a395d1db007837b, 0x32e16e21e4a983b2, 0x4a124cd12b54650a, 0xdd86089260aa58d1, 0x00027634736ab72a, 0xd841096db7e74e41, 0x9c22a87ec4122de9, 0x7fffffffffffe602, 0};
	setpmnsparams(9, 114384300578104255, 3);
	setpmnsaux(228768601156208510, 3);
	setpmnslastcol((uint64_t[]) {17669114049405379507u, 6942440984207203213u, 4409067010998670643u, 18267986996922550797u, 1705739488380651187u, 11673071576685656205u, 8685528207245359155u, 4657495042385159949u, 14912883622493628851u,});
	setpmnsalpha(2);
	cycles = do_bench(pmns_linearred_multa9, 9, 228768601156208510, nbrepet);
	printf("LinearRed 512 2X⁹ - 3: %ld cycles\n", cycles);
	lincycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 228768601156208510, 114384300578104255, p512_2X9M3, pmns_linearred_multa9), LOOPCHKNMB);
	
	mp_limb_t p512_3X9M1[] = {0xe2907b1738a3ffff, 0xb2e59f910cfc4703, 0xa59d13eb2cb6b81d, 0x1949688a255825f9, 0x553d35ad63982a38, 0x800e88ba9fec6726, 0x51d10289b39ab691, 0x800000000003973b, 0,};
	setpmnsparams(9, 109345452339395084, 1);
	setpmnsaux(328036357018185252, 1);
	setpmnslastcol((uint64_t[]) {9020529565262151679u, 1626805892883252724u, 17779891641117609840u, 16769175410557720896u, 6906438726493253376u, 5863990711463588864u, 8929334579714748416u, 8375488062192173056u, 1722173843283312640u});
	setpmnsalpha(3);
	cycles = do_bench(pmns_linearred_multa9, 9, 328036357018185252, nbrepet);
	printf("LinearRed 512 3X⁹ - 1: %ld cycles\n", cycles);
	lincycles[4] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 328036357018185252, 109345452339395084, p512_3X9M1, pmns_linearred_multa9), LOOPCHKNMB);
	
	mp_limb_t p512_3X9M2[] = {0x174e12654eefc0e7, 0x17c4b445718b7d51, 0x4803fbf80d40f947, 0x094ed44c9ad00817, 0x11d1a950653d8768, 0x89777a6bb60696c7, 0x0592091d091b27e1, 0x800000000005b5e9, 0,};
	setpmnsparams(9, 109345452339395267, 2);
	setpmnsaux(328036357018185801, 2);
	setpmnslastcol((uint64_t[]) {18180413013002900183u, 12454345914356764613u, 5003106998495831311u, 12932702879418884205u, 6808382930289665287u, 3333087395499346517u, 6929706690862407871u, 5614775205741639549u, 12928346647942223415u,});
	setpmnsalpha(3);
	cycles = do_bench(pmns_linearred_multa9, 9, 328036357018185801, nbrepet);
	printf("LinearRed 512 3X⁹ - 2: %ld cycles\n", cycles);
	lincycles[5] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 328036357018185801, 109345452339395267, p512_3X9M2, pmns_linearred_multa9), LOOPCHKNMB);
	
	
	
	
	printf("\n=============================================================================\n");
	printf("|        E        | X⁹ - 2 | X⁹ - 3 | 2X⁹ - 1 | 2X⁹ - 3 | 3X⁹ - 1 | 3X⁹ - 2 |\n");
	printf("=============================================================================\n");
	printf("| DoublSparsePMNS |   %lu  |   %lu  |   %lu   |   %lu   |   %lu   |   %lu   |\n", spacycles[0], spacycles[1], spacycles[2], spacycles[3], spacycles[4], spacycles[5]);
	printf("=============================================================================\n");
	printf("| LinearRed PMNS  |   %lu  |   %lu  |   %lu   |   %lu   |   %lu   |   %lu   |\n", lincycles[0], lincycles[1], lincycles[2], lincycles[3], lincycles[4], lincycles[5]);
	printf("=============================================================================\n");
	
	return 0;
}
