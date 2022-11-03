from commonpmns import pmnsWBdicts
from commonpmns import primesdict

(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsWBdicts[1024][primesdict[1024][0]]
#lam = [lam]
phi = 2**64
#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = (151557439827487301966960191905248420233409327748281249278640976537066266993988257219581147225805470889269899198735325674900327895538259785923482269207281672809731126255742686220820238394778878387647139546578384166978770043318856199095456814458300821186859906675966709623666377298322314372989786653793067293709, 20, 55529167333924931813547044165534606520774105315930203155652996896434240435165352314683910969818036125406366662120227324560427778281969046616316196428192805442528130670612588752931583815932158159842209292990151378076459303907140486715922697961694130441659567081245216208128077425413577600806343206289162522446, [1, 1], 56, [(1414837684904611, -414081814840150, -413391471858640, -127766233006937, 829527898116747, -298878292916479, -614181044643504, 348528569938809, -1783390822846956, 22139998511170, 987484838653877, 79149968742888, 1242595136257866, -1380621251883242, -1433749953981254, 216696762668728, -309109397772540, 368176042980452, 711627591568893, -100161670475483), (-143753864769633, 812080157219963, 827430876516685, -954074522151485, 366286059797800, -554358386681699, 714203351854180, -811895690112258, 891188900921903, -1193555433723812, 204877959326596, -728900321503521, -1731992288841426, -348647649968071, 1410568557490639, 126936871450257, 1208722967733439, 870076036142482, 691554705512250, -207221682432660),  (-61150188182945, 1939674980932783, 143301419411364, -892413519346329, -528391081009533, 264700617288401, -363225539588150, -760497156106728, 520401252442662, -770471714887050, 252664862033965, -762772490027948, 518704999184298, 1776657009525433, 986650112389251, 292292404453164, 1796380946742583, 30045299100678, -33955000690841, -430101930860897),(1918221492529040, -984293404301456, -208284062241008, -590052083814689, 1159377758095734, -1182284543558250, 316931735335602, 469002718437132, -399684066407809, -170418856802797, -810559392735317, 552577167708725, -474040278500291, -1138654547104253, 716210849554552, 1631025413739676, -557612679908466, 806075736350963, 295407775342194, 161730060245292),(923023294928398, -864991228024481, -334292311717835, 158308945529420, 512191436629487, 908063455843986, 496058064852867, -910466961765452, -383212754803161, 315856113757401, 1590584140061603, 27922238117419, -1019660753919940, -375581405129647, 148724272880352, -98362087699258, 321429539425876, 2117572926765171, 489683773810828, -371291413869482),(-1191957436476374, 554358386681699, -714203351854180, 811895690112258, -891188900921903, 1193555433723812, -204877959326596, 728900321503521, 1731992288841426, 348647649968071, -1410568557490639, -126936871450257, -1208722967733439, -870076036142482, -691554705512250, 207221682432660, 143753864769633, -955834021989596, 128403145472911, 825671376678574),(1016984210172541, -1112244104416098, -1097375941562849, 1258699579144129, -25967305672166, 449502734565779, -448670150524108, 1651686057028631, -1713956686166474, 975349674213646, -981565183537486, -969219798813478, -867352649152369, -366088452034794, -859713240938994, 916430563280275, -926304910600101, 661509406411572, -173266681741819, 286348066091264),  (430101930860897, 491252119043842, -1939674980932783, -143301419411364, 892413519346329, 528391081009533, -264700617288401, 363225539588150, 760497156106728, -520401252442662, 770471714887050, -252664862033965, 762772490027948, -518704999184298, -1776657009525433, -986650112389251, -292292404453164, -1796380946742583, -30045299100678, 33955000690841),   (-691554705512250, -484333023079590, 350975547202293, -812080157219963, -827430876516685, 954074522151485, -366286059797800, 554358386681699, -714203351854180, 811895690112258, -891188900921903, 1193555433723812, -204877959326596, 728900321503521, 1731992288841426, 348647649968071, -1410568557490639, -126936871450257, -1208722967733439, -870076036142482),   (714390629684028, 1511457756402999, -845308107285489, 153733538679103, -1279030632682526, 747548443440587, -1037505012385127, 1365428116227665, 808041148633310, 959425185026316, -18265039840399, 788202680409481, -876506684223072, 2240334492822236, 517319390700979, 709608147931373, -352400008639032, -194909361238633, -997975657660157, -916822539697058),  (665863795322505, -765014951860338, 1893920912378385, -368383081397862, -100837842553307, 266263752587124, -90991266508250, 1606405273779233, 359519745100384, -301173897058570, -543805080963179, 694275111910156, 1000368580347792, -1708810892350811, -692281284087056, 510954987529740, -1111315207672677, -194459329347704, 1177459389050625, -579529114622728),   (489683773810828, 118392359941346, 551731881058916, -864991228024481, -334292311717835, 158308945529420, 512191436629487, 908063455843986, 496058064852867, -910466961765452, -383212754803161, 315856113757401, 1590584140061603, 27922238117419, -1019660753919940, -375581405129647, 148724272880352, -98362087699258, 321429539425876, 2117572926765171), (-417968992505942, -514976353720999, -547118674421739, -1168359902965904, 242209132336297, -599527636101266, -102758630515858, 920729430565078, -140676737524931, 80119171375937, -2715799934211302, 57348464412538, -1601824182585025, -1326639535527095, 322981074928405, -100865613828076, 1341609769590704, 847762469770137, -188092512709125, 759527195300626), (28604249959603, 479022047487689, 400624848002596, 564827280828346, 64347246671671, 146316111463224, -171872682503853, -1012919107959906, -230524863522795, 1750257328681825, -439555030441410, -534061873267567, -2367271364272493, -1726042358434418, -1579684184073855, -339154696873218, 402131043671293, 1141729522429790, -39011482292538, -585987484211117), (-355250351591250, -928107216479319, -83216674767979, -105688236157591, -937482705174451, 902060996054216, -84160773291743, -211990250186326, 1174097758140074, 138452327648848, -1365903836806632, 511065590970923, 494137994210364, 1053241782050358, 547213561321867, 1043342717884301, 405206639420986, -1224205892605201, 1985474449819006, -76018351036194), (371291413869482, -551731881058916, 864991228024481, 334292311717835, -158308945529420, -512191436629487, -908063455843986, -496058064852867, 910466961765452, 383212754803161, -315856113757401, -1590584140061603, -27922238117419, 1019660753919940, 375581405129647, -148724272880352, 98362087699258, -321429539425876, -2117572926765171, -489683773810828),  (-615960131762736, -187751695401183, 1037151925406692, -1668366925217626, -304987005411793, -715153077416732, -1146675309315184, 929947667744189, -1183268639918869, 638653723364350, -110581445890482, 1052224309035314, 259333699462016, 197287673715763, -1567117352778337, 905673205115046, 108853308846751, -462005104164028, -1235282432111300, -161046599377972),(1788014522952879, 334292311717835, -158308945529420, -512191436629487, -908063455843986, -496058064852867, 910466961765452, 383212754803161, -315856113757401, -1590584140061603, -27922238117419, 1019660753919940, 375581405129647, -148724272880352, 98362087699258, -321429539425876, -2117572926765171, -489683773810828, 371291413869482, -923023294928398), (-589845444286311, 1196283654487782, -965813695899066, 451599756165841, 206526078710898, 671218952587327, -811069729545966, -1522244500487490, -147529494914058, -872923861081504, 405352753314331, 671628724896476, -1511434171318715, 1214672898140447, -360960497963302, -1058168548851607, 67972489788376, -210747310073282, 46746503554576, -1405945335196278),(-692281284087056, -181326296557316, -600360220142937, -1305774537020381, 983000059702921, 597930274427897, 86334680699777, -765014951860338, 1893920912378385, -368383081397862, -100837842553307, 266263752587124, -90991266508250, 1606405273779233, 359519745100384, -301173897058570, -543805080963179, 694275111910156, 1000368580347792, -1708810892350811)], [(16715620595547032752, 4148838620841847864, 3597172859447728654, 7147879807400772910, 15325302466994821025, 10149355799312004575, 5482139343293310071, 1298904718642239327, 12905173153395464177, 6400998728209654478, 3236558598097896140, 5569380232627906982, 7360377178497902058, 8909629826417828477, 12192879859356672649, 3512108077211541495, 1266459952607530582, 8667482064756412904, 17679466101863801584, 4264366641295374686), (1781570040899904858, 8608321285403524393, 6879056493083033924, 15092654548600467429, 18381514882328226489, 4724721066949717149, 10779285392892366556, 557551249852544815, 5834545757700524698, 11878229124787744453, 13075864468269153839, 5040736555707075768, 1266459952607530582, 2846007644579156125, 4264366641295374686, 11748278737820163231, 7335792920767557195, 3085266210791345201, 335957382899063361, 7360377178497902058),   (2403106819189706829, 14452671095549271205, 10720815483223933944, 18084850450966957272, 5044335911969660924, 10634481471468245531, 17679466101863801584, 10418671131542534042, 454426738221826118, 17316977488817039325, 5555540301308744290, 1641953695102114398, 7335792920767557195, 15599414267750807455, 7360377178497902058, 8667482064756412904, 793809201993448930, 14436798631075612475, 3346358466480214464, 1266459952607530582),(8247975941210233816, 16630612927756125862, 1673402968830445353, 11898011804299759070, 3616227687747099791, 7315492956444172328, 335957382899063361, 5429530003340369883, 8227713309106259750, 16287912564109625432, 2484114021836561781, 5223353487533353490, 793809201993448930, 3000126284188067307, 1266459952607530582, 3085266210791345201, 17259217116829263336, 8584028860982832243, 13374823289938260868, 7335792920767557195),(16601594295384204202, 1941385415785603542, 9789805908450394318, 8504542623640766464, 11679548077419337353, 741308967241306945, 3346358466480214464, 13617426274693425323, 16615845839983654695, 2067560962198384488, 9050681998341014909, 2706525069560306952, 17259217116829263336, 10928190095091924231, 7335792920767557195, 14436798631075612475, 3004359821572906260, 7833995554219420400, 3914283183600742831, 793809201993448930),  (31187280893934331, 3656687220497275097, 13398389166539047998, 4907911911931560475, 9079685037776833603, 10625819475727315181, 13374823289938260868, 4157069519486330892, 8101375074099592836, 623140119685461042, 11504668907828833771, 17078622692242827096, 3004359821572906260, 15891798723648029881, 793809201993448930, 8584028860982832243, 13909320770423347743, 18432678704676272675, 8681547433789691845, 17259217116829263336),  (11802338157196571016, 16014834199636685995, 16751759292972963326, 5275535788655107628, 4946138521313392424, 4722502932641115502, 3914283183600742831, 12635699379828881459, 15029213391463063339, 295700184353542415, 15836102007023484706, 394383667041817235, 13909320770423347743, 9850815842933532590, 17259217116829263336, 7833995554219420400, 9469496455566317356, 6568311592386295223, 14963698692114601594, 3004359821572906260),(225998561258944820, 8790711838187215639, 7449573515611732117, 659495256922999283, 17536014587601268237, 514903175527400760, 8681547433789691845, 2020562924890528604, 8715690113173075528, 7294959730883255534, 2720590438593585893, 14068290374420966093, 9469496455566317356, 11700620206191985117, 3004359821572906260, 18432678704676272675, 14634829724421246104, 1514542470033625206, 12922256161041987269, 13909320770423347743),   (8791488289532500321, 17689874589783434536, 2211988467190293567, 5625965064525758528, 7767686993139730333, 16689219295288724193, 14963698692114601594, 7765890147840213225, 2781029580167956091, 10877088173915454912, 10510311099856531873, 13861235695372168996, 14634829724421246104, 11699980710954125441, 13909320770423347743, 6568311592386295223, 13094143812461661051, 14373272219451087372, 12657500715241296487, 9469496455566317356),   (8397267462272954994, 5935777016860589864, 5684261551882823438, 6441781517725304896, 6087409551515912390, 11697583283570644384, 12922256161041987269, 10234715615096455370, 3507402826161678053, 16869449473672468034, 17326585270717743645, 12763755467355110790, 13094143812461661051, 3994427672161937656, 9469496455566317356, 1514542470033625206, 7366537958505338711, 12063259043341311630, 14217480584587660386, 14634829724421246104), (16442754257803327316, 16445131314702071671, 9766671021013350805, 16250293053578471326, 4493670347680145960, 14743732578154023082, 12657500715241296487, 2773834719998673959, 10886499868720601451, 13485616130315323967, 18141762228679430337, 17077147063899284713, 7366537958505338711, 5375369507645642164, 14634829724421246104, 14373272219451087372, 9190391432516118879, 5521340523419126199, 12911338181257146476, 13094143812461661051),  (17462369527157700779, 274898463937787746, 8667279602191769562, 1867045436290342427, 15925147787627861034, 14463933150766407379, 14217480584587660386, 5853080360366665078, 9007407257510882435, 11726053087631088402, 1797976652030857366, 14364944669396656100, 9190391432516118879, 14517823060815261424, 13094143812461661051, 12063259043341311630, 12256096878530868354, 8996685186032507425, 16018370480063339414, 7366537958505338711), (14284434601243591425, 3784008174642219629, 7436835748190176343, 12009382890761067480, 16348618354322138880, 7104272720772390163, 12911338181257146476, 4145827718544188302, 6176508010862735048, 7785911005531985867, 7242414943935984591, 2066915106105540317, 12256096878530868354, 13418454044487977768, 7366537958505338711, 5521340523419126199, 13089206661823777951, 12007151494848318052, 7954135048894858770, 9190391432516118879),   (8250478224579813817, 2821406846460062271, 14023789498329257497, 7386083203769911896, 5191648071732654630, 9360795802310314802, 16018370480063339414, 6042035231096234280, 13228968657042219977, 6383016120057498115, 8080461877866777288, 194447287102511547, 13089206661823777951, 12498044946445157399, 9190391432516118879, 8996685186032507425, 4109592904757194731, 13075067477507829737, 1298204741943257635, 12256096878530868354), (4733095049699878925, 9995951453487399506, 14307654109545657736, 13387290829571656413, 4069046797189073424, 1233698467748780916, 7954135048894858770, 15823949675809147662, 12624057485613719889, 4114017961826585551, 2357793174548338048, 5433217365266022651, 4109592904757194731, 16705460832419197537, 12256096878530868354, 12007151494848318052, 13034521369939708430, 11981024149913641043, 82234893427670587, 13089206661823777951),   (9608222910374234283, 13842023417793262622, 14663424894651952430, 8401450452188776199, 8244337340373764625, 16323397365593370855, 1298204741943257635, 2010807428432378561, 979019200958829833, 11173669059934183791, 7438591702307262196, 8479944440352100659, 13034521369939708430, 15279030726980074526, 13089206661823777951, 13075067477507829737, 7190840120347697238, 11321427686601071434, 12646656910527606381, 4109592904757194731),  (5834453813465634313, 18112898720006446055, 8438590095119888391, 15726247099103979923, 6815231268062906594, 11221482623080581061, 82234893427670587, 7795087880655489563, 3158606220103390406, 13301012402575216061, 6660167210898422120, 5148103793160894687, 7190840120347697238, 197545839361771388, 4109592904757194731, 11981024149913641043, 10945818834316280790, 2757636617536683885, 6239048568495823540, 13034521369939708430),   (13860720141614291918, 18187294031468305943, 5652786109905890511, 5524761983453645386, 6198684643722092624, 1267088669909423107, 12646656910527606381, 11304145774792531969, 6304065553772928217, 14197193078271232357, 12558533752374502833, 13126065807689203089, 10945818834316280790, 10060249242477294115, 13034521369939708430, 11321427686601071434, 16457246500652047335, 2848521101057816818, 9158630393790501234, 7190840120347697238),   (3483782409369791413, 5523719300869400139, 13024482248111133866, 10686077523335335370, 6786725145677955645, 5774217338169227064, 6239048568495823540, 10868234026306140110, 10780514555840803672, 16110113933835325254, 5722307822815416774, 2153677424735371032, 16457246500652047335, 2935988197966182613, 7190840120347697238, 2757636617536683885, 11624743819793276744, 15260386815031704726, 16261424736185676627, 10945818834316280790), (2602581213976375193, 8450535902980837088, 7584785046893675639, 13376434868400867366, 4991194135867004658, 15839931264513069157, 9158630393790501234, 13144796210964435862, 9947911398358600337, 12005226280127993965, 2299582692103077869, 1345658524983512560, 11624743819793276744, 14837614487286209406, 10945818834316280790, 2848521101057816818, 8626837131105432640, 1969016728867024519, 10012007421046616524, 16457246500652047335)])
