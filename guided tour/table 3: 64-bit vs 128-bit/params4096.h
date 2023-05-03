#ifndef PMNS_PARAMS_H_INCLUDED
#define PMNS_PARAMS_H_INCLUDED

#define RHO 56
#define N 84
#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 16")
#define LAMBDA 2
#define M_or_B_is_M

static const int64_t matrM[167] = {700640066605402, -222088006664796, -512003595974908, -522046760550590, 566785361601654, -117699209984250, -244361985394042, -207051752986200, 269618769217100, 36296505844490, 221933091089264, -511464644967990, -18680679631802, -711515591108782, 931246133897778, -493732373993232, 854194370119032, -451695989122204, 1164598467558366, 682379142052456, 439534234778212, -213025187168766, -241142319490532, 680790640590940, -385177830969226, -107453544260270, 16468955839494, -1301380384517606, -425849441772444, 629008549402094, 709624658273790, -775208752135628, -684626368996018, 153534716915542, 647759737123094, -333080771821146, 1023731571262530, -1019405253193062, -388236812342370, 406855590912944, -435196723013148, 633034876149618, 484087942789334, 562428066010390, 947763216765432, 420705055474292, -312038651080762, -96156474539328, -191108199789316, -927624207080028, 447798942482786, -130366315831924, 828766481826506, -28182921353664, -602169635366674, -712979223583796, -59863647578606, -26307616661172, 206072992366938, 539578732866348, -79308657633852, -380037498962268, 481770092476796, -604069219921280, 567738653187808, -108146053141910, -37711266592928, 442766762336044, 28746454860538, 212088966300254, -389862293746854, -52453270616764, 476413898228030, -381620449850978, -25976510296744, 235646630527316, 184456431237342, 756707260370198, 233147096588618, 99250369665780, 153547836173872, 204235030572032, -50322521285710, -457287077130029, 350320033302701, -111044003332398, -256001797987454, -261023380275295, 283392680800827, -58849604992125, -122180992697021, -103525876493100, 134809384608550, 18148252922245, 110966545544632, -255732322483995, -9340339815901, -355757795554391, 465623066948889, -246866186996616, 427097185059516, -225847994561102, 582299233779183, 341189571026228, 219767117389106, -106512593584383, -120571159745266, 340395320295470, -192588915484613, -53726772130135, 8234477919747, -650690192258803, -212924720886222, 314504274701047, 354812329136895, -387604376067814, -342313184498009, 76767358457771, 323879868561547, -166540385910573, 511865785631265, -509702626596531, -194118406171185, 203427795456472, -217598361506574, 316517438074809, 242043971394667, 281214033005195, 473881608382716, 210352527737146, -156019325540381, -48078237269664, -95554099894658, -463812103540014, 223899471241393, -65183157915962, 414383240913253, -14091460676832, -301084817683337, -356489611791898, -29931823789303, -13153808330586, 103036496183469, 269789366433174, -39654328816926, -190018749481134, 240885046238398, -302034609960640, 283869326593904, -54073026570955, -18855633296464, 221383381168022, 14373227430269, 106044483150127, -194931146873427, -26226635308382, 238206949114015, -190810224925489, -12988255148372, 117823315263658, 92228215618671, 378353630185099, 116573548294309, 49625184832890, 76773918086936, 102117515286016, -25161260642855},
	matrM1[167] = {8335766512371307838, -1792662776186755542, 4480908124497395854, -6903521565101502260, 8168668568215552420, 3452249563353978950, -4795561953779812790, 4866883045789984828, 1337184985608303680, -6138744571532008554, -1313330125781532184, -6195932585449935712, 962329761421277570, -2223124231731562104, -580962037014515642, -1862059307901801252, 7639224810761173686, 5505022026240967156, 4794637990099902326, 6591730730561909370, -3454203167195277720, -5286065386995494124, -7146227418688532542, 802617681991533780, 3854554337387371764, -2607680370296427460, -8541889600463442324, 194909735100862632, -7115651108653577438, 5000672693263845356, -1003312231814024068, -3757384273773851410, 1081113484758758456, -1428209149948484998, 5299149313851297904, -6344692687272164196, 7855108935401878552, -3718587113718175458, 3696659388905020996, -348427355305096602, -1416903025790807126, -6693177763908390748, 4957897772067394604, 6577969139228039840, -8182463118354777286, -4379331983420075506, 3566462340584493948, 1584033527420141148, -5537996681132793832, -3928849273563394242, 1951835593294121776, -8627048703170332692, 6017243997896465356, 898260441522857312, 2286613885136399826, -7503985463638346318, 2858145576492390002, 4868046674856750054, -1116357477971897006, 6394388416477972530, 3290254183582786644, -4665649680202855122, 8687084158633393372, -428010237432844482, -6301453206345275710, 4411577167852773964, 3296570690005474592, 8321412329270620172, -8676009397540174982, 391376424097861050, -1976545055626969878, 1562755626182323844, 3389146782252697646, -820070092770986564, -7927841331949413406, 4383611642392222640, 9170710649355004874, -8837472401501379022, 3434095087275702564, 8075020453512712050, 5492878967475716854, -6626643465058840524, 2627641254441328356, 1742975883062342877, -5055488780669121889, 8327040648761398037, 2240454062248697927, -3451760782550751130, 4084334284107776210, 1726124781676989475, 6825591059964869413, -6789930513959783394, -8554779544050623968, 6153999751088771531, 8566706973964009716, -3097966292724967856, 481164880710638785, 8111809920988994756, 8932891018347517987, 8292342382903875182, -5403759631474188965, 2752511013120483578, -6826053041804824645, -5927506671573821123, -1727101583597638860, 6580339343357028746, 5650258327510509537, 401308840995766890, 1927277168693685882, -1303840185148213730, 4952427236623054646, -9125917169304344492, 5665546482527987089, 2500336346631922678, -501656115907012034, 7344679899967850103, 540556742379379228, 8509267461880533309, -6573797379929126856, -3172346343636082098, -5295817569153836532, 7364078479995688079, 1848329694452510498, 9049158359202227507, 8514920523959372245, -3346588881954195374, -6744423150821078506, 3288984569614019920, 5132140477677387165, 7033706045144738055, -7440140866562528834, 792016763710070574, 6454373696288378892, 7258947400073078687, 975917796647060888, 4909847685269609462, 3008621998948232678, -8774241816093347152, 1143306942568199913, 5471379305035602649, -7794299248608580807, 2434023337428375027, 8665193297868827305, -6026177828615789543, -7578244945063382486, 6890547196753348247, 4343542079316696686, 9009366918138353567, 6072645433682137953, 2205788583926386982, -7575086691852038512, -5062665872219465722, -4338004698770087491, 195688212048930525, 8235099509041290869, 781377813091161922, -7528798645728426985, 8813336990469282526, 5259451370880069105, 2191805821196111320, -4638016712177273371, 4804635836104086297, -7506324493216924526, 4037510226756356025, 2746439483737858427, -3313321732529420262, 1313820627220664178};


#endif
