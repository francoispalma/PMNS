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
	
	mp_limb_t hp512n9[] = {0xfffffffffffffffd, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x98fac9ffffffffff,0xee846d77242cfe8f, 0xa545f31d9bad3565, 0x40002ce4a6f2f8a8, 0};
	setpmnsparams(9, 114384436610465792, 3);
	setpmnssparseinv(6187042836773339136u, 6187042836773339136u, 12297829382473034411u);
	cycles = do_bench(pmns_doublesparse_mult, 9, 114384436610465794, nbrepet);
	printf("DoubleSparse 512 n 9: %ld cycles\n", cycles);
	spacycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 114384436610465794, 114384436610465792, hp512n9, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp512n10[] = {0x0000000000000055, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0bea273b1e1f4201, 0x770d03e875d3c28b, 0x80009a7e2192765c, 0, 0};
	setpmnsparams(10, 2413423728001024, -85);
	setpmnssparseinv(14974387384261541888u, 14974387384261541888u, 217020518514230019u);
	cycles = do_bench(pmns_doublesparse_mult, 10, 2413423728001108, nbrepet);
	printf("DoubleSparse 512 n 10: %ld cycles\n", cycles);
	spacycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(10, 2413423728001108, 2413423728001024, hp512n10, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp512n11[] = {0xffffffffffffdac9, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x353fffffffffffff, 0xa6bf492c184c0fc9, 0x8002b883ed2cf311, 0, 0, 0};
	setpmnsparams(11, 96430605729792, 9527);
	setpmnssparseinv(9923329849290653696u, 9923329849290653696u, 1266313700451983495u);
	cycles = do_bench(pmns_doublesparse_mult, 11, 96430605739318, nbrepet);
	printf("DoubleSparse 512 n 11: %ld cycles\n", cycles);
	spacycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(11, 96430605739318, 96430605729792, hp512n11, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp512n12[] = {0x000000000000001f, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x7110c9366947b801, 0x80bc7f3a53591a8d, 0, 0, 0, 0};
	setpmnsparams(12, 6592774799360, -31);
	setpmnssparseinv(2380224829098819584u, 2380224829098819584u, 1190112520884487201u);
	cycles = do_bench(pmns_doublesparse_mult, 12, 6592774799390, nbrepet);
	printf("DoubleSparse 512 n 12: %ld cycles\n", cycles);
	spacycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(12, 6592774799390, 6592774799360, hp512n12, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp512n13[] = {0xffffffffffff8205, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x1204901effffffff, 0x86210c706fd3f85b, 0, 0, 0, 0, 0};
	setpmnsparams(13, 682899800064, 32251);
	setpmnssparseinv(16095929582925905920u, 16095929582925905920u, 2894762077363307827u);
	cycles = do_bench(pmns_doublesparse_mult, 13, 682899832314, nbrepet);
	printf("DoubleSparse 512 n 13: %ld cycles\n", cycles);
	spacycles[4] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(13, 682899832314, 682899800064, hp512n13, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t hp512n14[] = {0xffffffffffff800b, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xa0e2073737609370, 0, 0, 0, 0, 0, 0};
	setpmnsparams(14, 98784247808, 32757);
	setpmnssparseinv(4992790333596106752u, 4992790333596106752u, 11484093534058640477u);
	cycles = do_bench(pmns_doublesparse_mult, 14, 98784280564, nbrepet);
	printf("DoubleSparse 512 n 14: %ld cycles\n", cycles);
	spacycles[5] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(14, 98784280564, 98784247808, hp512n14, pmns_doublesparse_mult), LOOPCHKNMB);
	
	mp_limb_t p512n9[] = {0x21cb2f3997fffffd, 0x8dc76bf1f3d1c3ee, 0xb41bfd602771fafd, 0x69772bbf38af18b3, 0x651fe55cbafdeb16, 0x7c8ce4e9c47f9cf8, 0x80d079a56c657e96, 0x8000000000034f27, 0};
	setpmnsparams(9, 123541877815766680, 3);
	setpmnslastcol((uint64_t[]) {1068504517212198229u, 9678330455461638264u, 12298442381156439872u, 9543643910875729408u, 12539807914892349440u, 57035926836051968u, 11747439154570723328u, 15471119491003318272u, 16325372834547236864u,});
	cycles = do_bench(pmns_linearred_mult, 9, 123541877815766682, nbrepet);
	printf("LinearRed 512 n 9: %ld cycles\n", cycles);
	lincycles[0] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(9, 123541877815766682, 123541877815766680, p512n9, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p512n10[] = {0x25d1a1895134c4cd, 0xd6ef1698ac0ca6bc, 0xaf2a8b0844e4b551, 0x92c0365f510dff23, 0x8cace594059b4443, 0xbad72e65c068d15f, 0x255d6e443567e775, 0x7ffffffffffff4df, 0, 0};
	setpmnsparams(10, 2413419283252018, -205);
	setpmnslastcol((uint64_t[]) {8048274811616020485u, 14548015780920532986u, 14717448781875784916u, 1419881223893630312u, 18380552099849188944u, 13685401514979986336u, 18439981564427248960u, 14167806348610700928u, 15590309945547080960u, 11718697798495468032u,});
	cycles = do_bench(pmns_linearred_mult, 10, 2413419283252222, nbrepet);
	printf("LinearRed 512 n 10: %ld cycles\n", cycles);
	lincycles[1] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(10, 2413419283252222, 2413419283252018, p512n10, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p512n11[] = {0x19f68125383fdad7, 0x6c47d6f365e833c3, 0x87b4fe0f6290751f, 0x9fc45b816eb7d774, 0x570dd00901a2c33d, 0xc1e90cf678116770, 0xf0c39e62f98ca5ea, 0x7ffffffffffaa485, 0, 0, 0};
	setpmnsparams(11, 96429877877380, 9513);
	setpmnslastcol((uint64_t[]) {14724735671038228711u, 24596604350827804u, 17170785449382812272u, 5845011505042620864u, 10747644975828854528u, 15868852375404092416u, 4113434195801567232u, 10598283230708875264u, 5666029066904862720u, 7511305327566651392u, 17614647086092910592u,});
	cycles = do_bench(pmns_linearred_mult, 11, 96429877886892, nbrepet);
	printf("LinearRed 512 n 11: %ld cycles\n", cycles);
	lincycles[2] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(11, 96429877886892, 96429877877380, p512n11, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p512n12[] = {0xdd450f35cc425bed, 0x46735881daf4ba59, 0x0d881cfb7af79be8, 0x5c81fedc4901de86, 0x17f7863415264cff, 0xc115bb6f0b7423ba, 0x2f25032608ea01e2, 0x7fffffffff13c208, 0, 0, 0, 0};
	setpmnsparams(12, 6589624212019, -316);
	setpmnslastcol((uint64_t[]) {16017345202995608037u, 1170012822622800543u, 18149658134610772909u, 6448414547377524087u, 4816006101951513781u, 7838234492792279567u, 4207885199745900285u, 13856534732049047143u, 17502532093746906245u, 11492777003979155583u, 14298226222570977101u, 3334936867563443799u, });
	cycles = do_bench(pmns_linearred_mult, 12, 6589624212334, nbrepet);
	printf("LinearRed 512 n 12: %ld cycles\n", cycles);
	lincycles[3] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(12, 6589624212334, 6589624212019, p512n12, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p512n13[] = {0xff01b498e3ff81c1, 0x46ebf9d1c3a98b5f, 0xcaa7a499c53d5098, 0x5eaf12f8c1f6093c, 0xbb5b2ee7523ad29f, 0x9115ba9b14e3fdcf, 0x18eeadf911c4590b, 0x7ffffffff91ad381, 0, 0, 0, 0, 0};
	setpmnsparams(13, 680447221348, 32319);
	setpmnslastcol((uint64_t[]) {9145749752424926785u, 12916547772262362980u, 17313644160760167184u, 15729382007476253248u, 9806024915076915456u, 3908778838583337984u, 462937427461869568u, 11978120025109839872u, 8992079544097767424u, 9154712279598759936u, 3486524597775368192u, 4426334649127534592u, 9488263963926855680u, });
	cycles = do_bench(pmns_linearred_mult, 13, 680447253666, nbrepet);
	printf("LinearRed 512 n 13: %ld cycles\n", cycles);
	lincycles[4] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(13, 680447253666, 680447221348, p512n13, pmns_linearred_mult), LOOPCHKNMB);
	
	mp_limb_t p512n14[] = {0xe865d1a5247aa50f, 0x5e40f8129da2726b, 0x39eb70a45438f3e6, 0x89f13401fcc46281, 0xd66eea1741fd075f, 0xbfdaa9c9b970d080, 0x42e55687df9ec8bc, 0x7fffffffed802789, 0, 0, 0, 0, 0, 0};
	setpmnsparams(14, 97184015999, 32754);
	setpmnslastcol((uint64_t[]) {11051800913993247215u, 12714427169961707409u, 4342339768435148527u, 6057403313612630673u, 8313579190102352879u, 8584116000186944913u, 6804934225241358575u, 11889436917606445201u, 11495295314222849519u, 7784368896198775697u, 9722824374478469871u, 431508576577402513u, 266525530447730671u, 9313713510985652625u,});
	cycles = do_bench(pmns_linearred_mult, 14, 97184048752, nbrepet);
	printf("LinearRed 512 n 14: %ld cycles\n", cycles);
	lincycles[5] = cycles;
	printf("Correct: %ld/%d\n", checkPmns(14, 97184048752, 97184015999, p512n14, pmns_linearred_mult), LOOPCHKNMB);
	
	
	printf("\n================================================================\n");
	printf("|            n             |  9  | 10  | 11  | 12  | 13  | 14  |\n");
	printf("================================================================\n");
	printf("|    DoubleSparse PMNS     | %lu | %lu | %lu | %lu | %lu | %lu |\n", spacycles[0], spacycles[1], spacycles[2], spacycles[3], spacycles[4], spacycles[5]);
	printf("================================================================\n");
	printf("|      LinearRed PMNS      | %lu | %lu | %lu | %lu | %lu | %lu |\n", lincycles[0], lincycles[1], lincycles[2], lincycles[3], lincycles[4], lincycles[5]);
	printf("================================================================\n");
	
	return 0;
}
