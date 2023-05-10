#ifndef PMNS_PARAMS_H_INCLUDED
#define PMNS_PARAMS_H_INCLUDED

#define RHO 55
#define N 189
#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 16")
#define LAMBDA 2

static const int64_t matrM[377] = {-74835604927954, -140753560851948, -96749279868414, -22452966589650, -82385225707050, 45305026920656, -39923394668990, 43369589992726, -53912200785528, -21652048977494, -22373041228154, 17647119305530, -40532864308352, 23506614663924, 13128383145810, -13689079475010, 14986399157214, -34826920040820, -70226087758372, 153023692807534, -68156564455580, -59549987098354, 33580765893838, -20268623633672, -100865113365838, 37821386085946, -14840123595630, 31847365063054, -17043500213138, 61592994411222, -55578836358914, -52340117770360, 10703435169754, -66646697870262, -96860660901054, 91206032964430, 47138394215506, -34324005024512, 35803587251208, -57394054448230, 5116667629040, -11539067171050, 85008298705240, 192774950423918, 95880606804552, -70165088294832, -61528430043140, 66036816265540, -3645318830098, -66608449848398, 60841053485538, -98659092089250, -11951055156222, -441468149298, 25914888476214, -82699863510466, 87862392068102, -120068981287740, 96178234501452, 59004181670900, -165968254220516, -118331452866108, 6125179753380, -21285114607718, 72515901352694, 112586886534204, -72349495520364, 6011776827138, 83318716615560, 63582750160880, 94604597191066, 83829836576096, -52737770328806, -54763028408460, -162586612534772, -60877109719624, 63108628272386, 12359207126838, -16521778922188, -31197521861238, 112252709414122, -93910904797138, -100505063280426, -8987768085514, 59363267044402, -34716468140698, -120627267627488, -56814806997910, 125717061267076, -27287275800506, -38715120007346, 123456761925524, 90120752249086, -66653485325348, 18954597965140, 67795176247848, -27203508982494, -3931723513014, -63885282857808, 109565245725306, -106906487103952, -27820291274706, -124522460937072, 130919304412450, 39373995609602, 126864064292258, -1041987590198, 122110014212870, -16519761505156, -19674379306072, -13960459459886, -61938575431436, 14251357788396, 52758537779456, 13788004114436, 125178308017548, 37640800976318, 43632518100626, -127179494101742, 106404606768320, 7709963174410, -88297191710706, 131040011740574, 96456340130180, -105504462705618, -30174960041700, -86981964056906, 14414233646558, -148031599122680, 37463080638038, 27359370333492, -18545230203348, 28569214353168, 86359307034588, 63601906866132, 133106354686944, 112701861455146, -30401438390458, 56553775062666, -66982666214166, 14935731282934, 67616482347802, 18924486007172, 31514854092670, 101465348598932, 101179687968750, -92716792342040, 101764990024508, 78047203251456, -99539001061816, 166155149098406, -4724943382228, -55142570483468, -113261186270216, 19015779850972, -87190789904118, -18844930633556, -63027966674348, -157434176102674, 77180830289890, -5942947268252, -42917098828740, 39992259308794, -71905051715644, 71307115317278, -53774922793672, 19111782636970, -40183352514398, -54414691991622, 18314062906732, -100415846284536, 57346013784922, -86298112333276, 15997553373022, -99209744211012, -21901307596878, 31671976647082, -44675055482198, 63671746394366, 70581920594710, -50573407105698, 12643465564648, 15745141607930, 19832114980962, -15893670863988, -227026241686, -2858716038692, -1807100220, -25575841249679, -37417802463977, -70376780425974, -48374639934207, -11226483294825, -41192612853525, 22652513460328, -19961697334495, 21684794996363, -26956100392764, -10826024488747, -11186520614077, 8823559652765, -20266432154176, 11753307331962, 6564191572905, -6844539737505, 7493199578607, -17413460020410, -35113043879186, 76511846403767, -34078282227790, -29774993549177, 16790382946919, -10134311816836, -50432556682919, 18910693042973, -7420061797815, 15923682531527, -8521750106569, 30796497205611, -27789418179457, -26170058885180, 5351717584877, -33323348935131, -48430330450527, 45603016482215, 23569197107753, -17162002512256, 17901793625604, -28697027224115, 2558333814520, -5769533585525, 42504149352620, 96387475211959, 47940303402276, -35082544147416, -30764215021570, 33018408132770, -1822659415049, -33304224924199, 30420526742769, -49329546044625, -5975527578111, -220734074649, 12957444238107, -41349931755233, 43931196034051, -60034490643870, 48089117250726, 29502090835450, -82984127110258, -59165726433054, 3062589876690, -10642557303859, 36257950676347, 56293443267102, -36174747760182, 3005888413569, 41659358307780, 31791375080440, 47302298595533, 41914918288048, -26368885164403, -27381514204230, -81293306267386, -30438554859812, 31554314136193, 6179603563419, -8260889461094, -15598760930619, 56126354707061, -46955452398569, -50252531640213, -4493884042757, 29681633522201, -17358234070349, -60313633813744, -28407403498955, 62858530633538, -13643637900253, -19357560003673, 61728380962762, 45060376124543, -33326742662674, 9477298982570, 33897588123924, -13601754491247, -1965861756507, -31942641428904, 54782622862653, -53453243551976, -13910145637353, -62261230468536, 65459652206225, 19686997804801, 63432032146129, -520993795099, 61055007106435, -8259880752578, -9837189653036, -6980229729943, -30969287715718, 7125678894198, 26379268889728, 6894002057218, 62589154008774, 18820400488159, 21816259050313, -63589747050871, 53202303384160, 3854981587205, -44148595855353, 65520005870287, 48228170065090, -52752231352809, -15087480020850, -43490982028453, 7207116823279, -74015799561340, 18731540319019, 13679685166746, -9272615101674, 14284607176584, 43179653517294, 31800953433066, 66553177343472, 56350930727573, -15200719195229, 28276887531333, -33491333107083, 7467865641467, 33808241173901, 9462243003586, 15757427046335, 50732674299466, 50589843984375, -46358396171020, 50882495012254, 39023601625728, -49769500530908, 83077574549203, -2362471691114, -27571285241734, -56630593135108, 9507889925486, -43595394952059, -9422465316778, -31513983337174, -78717088051337, 38590415144945, -2971473634126, -21458549414370, 19996129654397, -35952525857822, 35653557658639, -26887461396836, 9555891318485, -20091676257199, -27207345995811, 9157031453366, -50207923142268, 28673006892461, -43149056166638, 7998776686511, -49604872105506, -10950653798439, 15835988323541, -22337527741099, 31835873197183, 35290960297355, -25286703552849, 6321732782324, 7872570803965, 9916057490481, -7946835431994, -113513120843, -1429358019346, -903550110},
	matrM1[377] = {1073890529184970794, -261048136680857178, -1242995272374954896, 7251512640639377312, -2126306029951346378, 5848316238302271854, 6648263751765544330, 5194326501046165528, -8380246138171717360, 7503180577893123382, -6045780688690967974, 5741037116902739760, 3920161917413348730, -7444209102613586944, -4562127960898053618, -8287295613025968694, -2658809744823763326, -4339625448124964202, 6434059182244071336, 4547467726642974538, -3270507004156614300, 1052306189488882824, 4713018099419449864, -136149417393030666, 1518103432246098624, 5059910293790221160, -3780563989290745600, 5719313154293661056, -5828445959493663586, 4016793945242201798, 2581349746848417198, 8444566828760615518, -9148077788785647556, 626313868520374576, -442162065774150642, 6533153722339249304, 9140260683658171186, 5185591553050485568, 8205166344541595366, 1133293146525444396, -3117886660632024518, -7557906122393156700, 6624923388417453368, -4897378404832804928, -5282643430322060836, 2520015257039040406, -9209420661174869610, -4184358435498716716, -4660669365297370294, 9020770690252883436, -5890456144463783120, 4620986804328261186, -8670299656847069006, 6737987097996223904, 6525561292965203720, -5108019265678780302, 8126933007216076, -5651742220516228178, -2609409060205658878, 7050423063508378164, -914584049969331340, 3060521559113984034, 7817534328530710022, 7109798426463143612, -741997656599913108, -4639020514583104192, 4708467551571266068, 4311188326796953392, 4059339883690762670, -7539702073567677018, -1782436237605485298, 2940338433512431518, 4488475182461174432, 1474974666966450542, 3048433669693393486, 7936926633526126882, 4780889388898187660, -6550655424305303950, -8579937363298201132, -7891544701796006054, -4990560993222180584, -3890366958500613814, 4916960078635125422, 3283634917223561338, -2749797670938810780, -2960875585793943686, -3282293106057411242, -524237588834578266, -2629797125829868490, -4827347228856174316, 7741415272880031910, -3996247631887791464, -6534291820509776676, -9048763187705695544, 3802877068587173528, -7368227662440062816, -8050473706719025334, 5442566432272588664, -6430304082847404344, 13163247448394530, 6469095830841769192, 4921077679612306854, 4896706852676989170, 8133801386635846408, -7546416466393426188, -6881458976563766754, -2952380147771676294, -262571260500081218, -1402482706411703300, -5357580094965339274, 4502441622219303638, -2628283023246599798, 736093293738640362, 2118549304825722280, 1303763581386589616, 1051217148903866982, -3283846189201139960, -7868097585662523094, 6331770774609625462, 2620940424779361774, 3107441351002955088, 8933106849768966216, 4425431337897735356, 6530237680528745078, -1853951713766117318, 9187071417688402646, 5528330762155289774, -8221436824244766532, -3364792893436990086, -2324554452125868964, -8605555377140631930, -2683226158623346226, -4686488952446553248, 4132947827663035360, -2644626724402694176, -8111118389841741336, 8548410708174584960, 4233808461872167300, -5668971403342208642, 5234602925206639782, 9126797344004516036, 1880233687197233390, -656285522982497004, 1011233662338988784, 4989855614796124580, 8005218747224500272, -7487276849805032168, -1762313838034748982, 8763899481343804378, 7431859907850516718, -7976143953154301760, 9010300109212615584, 979703013575455472, 1342382433752865956, -8980813227868130028, -1743900573113025864, -4904790364010718710, -4217193070916090776, -4520058218089983752, -7824644370943555574, 4821478370971393272, -9094569428066275642, -6438080841803750664, -4331456458563597422, -7105759202509993402, -1610197553606605846, 3311093769893110624, 7814212375498446474, -1274995804477942138, -6780268493687085332, -7308626551318628318, -8440533428451912924, -9208512224879584964, 2368850643103328736, -5214424264142963102, 5047351242342825296, -506535475139478152, 6266998030693190736, -6533339641904815058, 7438425511748948104, 7980709017612912846, -5563813451839813118, 7471742125929322688, 4259208483353032568, 4467329527276739108, 1417440082664658198, -7116829354105099092, -2079037093648219392, 753863518897951137, -8686426772262290411, -130524068340428589, -621497636187477448, -5597615716535087152, 8160219021879102619, 2924158119151135927, -5899240160972003643, 2597163250523082764, 5033248967768917128, -5471781747908214117, 6200481692509291821, -6352853478403405928, -7263291078148101443, -3722104551306793472, 6942308056405748999, 5079724230341791461, -1329404872411881663, -2169812724062482101, 3217029591122035668, 2273733863321487269, -1635253502078307150, 526153094744441412, 2356509049709724932, -68074708696515333, -8464320320731726496, -6693416889959665228, -1890281994645372800, -6363715459707945280, 6309149057107944015, -7214975064233674909, -7932697163430567209, 4222283414380307759, 4649333142461952030, -8910215102594588520, -221081032887075321, 3266576861169624652, 4570130341829085593, -6630576260329533024, -5120788864583978125, -8656725463592053610, -1558943330316012259, -3778953061196578350, -5910910342646049124, -2448689202416402464, 6582050321693745390, -7963364408335255605, -4604710330587434805, 7131192819105417450, 6893037354206090661, 4510385345126441718, -2945228072231891560, 2310493402164130593, 4888222208431241305, -5854378487856663856, 3262780646482601860, 6669362404015385657, 4063466503608038, 6397500926596661719, 7918667506751946369, 3525211531754189082, 8766080011870110138, 1530260779556992017, -5314604872589420797, 3554899213231571806, -370998828299956554, -2319510257291552096, 2354233775785633034, 2155594163398476696, 2029669941845381335, 5453521000070937299, -891218118802742649, 1470169216756215759, 2244237591230587216, -8485884703371550537, -7699155202008079065, -5254908720091712367, 2390444694449093830, -3275327712152651975, -4289968681649100566, 5277599685956772781, -2495280496611090292, 7278188557604468901, 2458480039317562711, -7581554578242995139, -1374898835469405390, -1480437792896971843, 7582225483826070187, 8961253242437486675, 7908473473939841563, -2413673614428087158, 3870707636440015955, 7225248220910880076, 5956226126599887470, -4524381593852847772, -7321933502561189044, -3684113831220031408, 5198135183495263141, -6502088820718481476, -3215152041423702172, -9216790413130578543, 3234547915420884596, -6762833197048622381, -6775018610516281223, -5156471343536852604, 5450163803658062714, -3440729488281883377, 7747181962968937661, 9092086406604735199, -701241353205851650, -2678790047482669637, -6972151225745123989, 7909230525231475909, -8855325389985455627, -8164097384441914668, -8571490246161481000, 525608574451933491, -1641923094600569980, -3934048792831261547, -6057486649549963077, -7912901824465094921, 1553720675501477544, -4756818611970292700, -7010656367905908130, 3265118840264372539, -926975856883058659, -4629836328010574485, -6459206655777130921, -4110718412122383266, 7540975590136280765, 8061094810791841326, 4920594348284459843, 7881758957543102695, -2343244476223276624, -7156898123023258128, -1322313362201347088, -4055559194920870668, -4949166682767483328, 2116904230936083650, 6388886335183671487, -6606070574251455917, 4563398672002258018, 940116843598616695, -328142761491248502, 505616831169494392, 2494927807398062290, -5220762663242525672, -3743638424902516084, -881156919017374491, -4841422296182873619, 3715929953925258359, 5235300060277624928, 4505150054606307792, -8733520530067048072, 671191216876432978, -4490406613934065014, 8351421750298262876, 6770976854849416453, -2108596535458045388, 6963342927809783932, 5311049851382998021, -6812632851369079172, 4676087322821637987, 6004331615952900476, 7057643807572977097, 5670492435599779107, -805098776803302923, -7567825151908220496, -5316265849105552571, -637497902238971069, -3390134246843542666, 5569058761195461649, -4220266714225956462, -4604256112439792482, 1184425321551664368, -2607212132071481551, 2523675621171412648, -253267737569739076, 3133499015346595368, -3266669820952407529, 3719212755874474052, 3990354508806456423, 6441465310934869249, -5487500973890114464, 2129604241676516284, 2233664763638369554, -8514651995522446709, 5664957359802226262, -1039518546824109696};


#endif
