#ifndef PMNS_PARAMS_H_INCLUDED
#define PMNS_PARAMS_H_INCLUDED

#define RHO 57
#define N 32
#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 32")
#define LAMBDA 2
static const int64_t matrM[63] = {1739010899520116, -674848188975296, 1514410973310716, -1511848621879022, -2010807357502540, 4523060390393748, -1863275872830378, 2594909491454056, 1095280543785636, 175998813611876, 4639152212077364, 500450982620158, 849534713996330, 774983099415798, -1267190563080146, 1090862228699608, 693431810589534, 5158412428534522, 2468735548859442, 4746580469325258, 1814100144225008, 2024257337722696, -3095965850556552, -2098882412406078, 8529411884267570, -1602262729721248, -4166894115948778, -3076053510193862, 1230454135973782, 3316288504164370, 13827583825244, 1110168244588365, 869505449760058, -337424094487648, 757205486655358, -755924310939511, -1005403678751270, 2261530195196874, -931637936415189, 1297454745727028, 547640271892818, 87999406805938, 2319576106038682, 250225491310079, 424767356998165, 387491549707899, -633595281540073, 545431114349804, 346715905294767, 2579206214267261, 1234367774429721, 2373290234662629, 907050072112504, 1012128668861348, -1547982925278276, -1049441206203039, 4264705942133785, -801131364860624, -2083447057974389, -1538026755096931, 615227067986891, 1658144252082185, 6913791912622},
	matrM1[63] = {8949124839910702444, -7911302509913314084, -5871549420717787708, -4439991199147255930, 7955690544113509632, 1591716925983414496, 8722894573771328074, 2055191512176694746, 4572975294229151744, 4213986915354815836, -3065358328812684376, -5777888384047912472, -1909640078726749346, 6191591789474107056, -165893077215377568, -5441516223538349202, 538586134159518278, 3325803297402095348, 772284350585546554, -346458637601580378, 8454895749150831384, -3899434260257965798, 8754587829646447308, -7260605670999158822, 4274490378459682384, -3454412523077053014, 1545110820353586414, 8313386676822663768, 5811841165371941266, -8843432686745395626, 9166981748228514764, 1002137203480594193, 4474562419955351222, -3955651254956657042, 6287597326495881954, -2219995599573627965, 3977845272056754816, 795858462991707248, -4861924749969111771, -8195776280766428435, 2286487647114575872, 2106993457677407918, 7690692872448433620, -2888944192023956236, -954820039363374673, -6127576142117722280, -82946538607688784, 6502613925085601207, -8954078969775016669, -7560470388153728134, 386142175292773277, -173229318800790189, -4995924162279360116, -1949717130128982899, 4377293914823223654, -3630302835499579411, 2137245189229841192, -1727206261538526507, 772555410176793207, -5066678698443443924, -6317451454168805175, -4421716343372697813, 4583490874114257382};
static _mpnum __P__ = { .deg = 26,
		.sign = 1,
		.t = (uint64_t[]) {0x6f89b4cf73445b5f, 0x7377846c5ac639d5, 0xd5c2df5e22aaeb11, 0x9ddfd46e0dec0d42, 0x2d368be1b2773a6, 0xacac9c355fdab1b2, 0xc1f566e41a8b16e0, 0x5a86c9bc4701128f, 0x43595e8ecd76fb69, 0x385cbd210c36119c, 0xa873620db115d3fc, 0xc2bd7239ec7c8adc, 0x64d278c4ae1e2eeb, 0x441f6db7c1d635e5, 0x8b2a48141a19e161, 0xd4b15e313958194f, 0xdfb5ec6991957a63, 0x66b26e2b7ffa4891, 0x658dd29b92d631d8, 0xa1068e5c2ee7521, 0x1288cc2edec3becb, 0x7446c03ce5defc9e, 0x8685f5ff6068f6ff, 0xdec53bfc6c41e8a, 0xb6cf3f294d883543, 0x8af3a3481f298816} },
	Gi[] = { { .deg = 26,
		.sign = 1,
		.t = (uint64_t[]) {0x19661bf9229d8766, 0x5732bb9ba3fcd816, 0x1f7791d6b06efdf3, 0xadb548c03cca7128, 0x531781836e94231, 0xad069a621173617e, 0x7ff6f43989072c87, 0xadbc233691d98b66, 0x5629553f657b9fa6, 0xb8c1f42daa65c46e, 0xac023120da7b7c3e, 0x7a493d4788624215, 0x18c0f350a021ecb9, 0x4eb0928b0392bdac, 0x54379146c0974c89, 0x216a25e6a865b4eb, 0x5601593c97dc1090, 0xe2c309fbedb0602c, 0x967f38cf4790ddf0, 0xb205508a40b00f34, 0x73788eac670b0ed, 0x1d5ea453037adf8b, 0xa01282b5b094516f, 0xf6831026373d988e, 0xf21b5fc71f8b523, 0x27a68f58266888a5} }};
#endif
