#ifndef PMNS_PARAMS128_H_INCLUDED
#define PMNS_PARAMS128_H_INCLUDED

#define RHO 118
#define N 18
#define LAMBDA 2

#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 18")

static const uint64_t Blo[N][N] = {
{3470414713307816357u, 7903536351458640381u, 5764977005191454068u, 5933600014869213668u, 16701521397656342453u, 17139294038958903764u, 10054533042546738942u, 1150185687112210211u, 2155481186135059958u, 15406447097124422997u, 11919030303352447879u, 1430881475012654698u, 2890391294705074332u, 7729399942530821681u, 10268746376625683920u, 11453682010611315102u, 15808255107448654012u, 2664399730523436155u}, {13302847384567125314u, 751265108697769604u, 12180204935537748761u, 767093895220179054u, 16679329493376635099u, 6202597489297903926u, 5271976952202843285u, 13022663182262102126u, 12011299516645638745u, 8867025309758780250u, 17465738731556679889u, 4133144466152863675u, 138632888775908338u, 7028189882574966185u, 12385810038165215221u, 4233032251590503182u, 7212353433323584936u, 8761341902355894872u}, {13049191079335965028u, 7125003677735759136u, 14555528085633161783u, 14525507471846706323u, 5163198935710851286u, 12914419158432525034u, 12623135365790469160u, 15602971656705112049u, 18061969092547674153u, 10554250467082746404u, 15971674295351981248u, 6277026923808910564u, 6666359356795999174u, 1761693480656552144u, 7365067243773894949u, 12637354918384717077u, 17905949686861771462u, 14609745198198036036u}, {14610686590622480978u, 13016564572244021363u, 3826382830012989673u, 17604232215136407379u, 13382301755956933906u, 17312781658607287602u, 4887949710442307515u, 7267899172467232720u, 2397257458298336763u, 13621074272760579003u, 8756095045666888985u, 13366498880575089380u, 13562639125427913612u, 1153145345174171704u, 11482954189268099547u, 16343950516494743640u, 14551819128265915860u, 2369532115651207564u}, {8226421791646689884u, 17552047572554203043u, 6154296740515555579u, 4864549065635525981u, 5988502586750379582u, 6277859104244690316u, 12807529254569474497u, 17445384039414121751u, 15282392720191035799u, 6593084381281889130u, 3488051538744758705u, 18255012314641170567u, 10308185839729087494u, 5282299211889741167u, 16543508648392223435u, 3083798899367679703u, 18208524553000323697u, 15938579069282435714u}, {12180204935537748761u, 767093895220179054u, 16679329493376635099u, 6202597489297903926u, 5271976952202843285u, 13022663182262102126u, 12011299516645638745u, 8867025309758780250u, 17465738731556679889u, 4133144466152863675u, 138632888775908338u, 7028189882574966185u, 12385810038165215221u, 4233032251590503182u, 7212353433323584936u, 8761341902355894872u, 6651423692283562657u, 9599004591203660610u}, {2059104314753125268u, 5998760325374448904u, 10860015863404971582u, 15069723739583311808u, 2905569751878202733u, 16382161121512226961u, 12638635940501566356u, 6256535506850769339u, 17699727325121029884u, 3773043620848927026u, 13870244853972624387u, 9371987303198134043u, 2005452045517419843u, 6510765500764420740u, 3330176402693387508u, 876527624518273778u, 16100807597298225796u, 259613905578875349u}, {13758413576626527269u, 6930613516043767803u, 6711821669134621108u, 7351158413587625875u, 2580308474443009923u, 6050669575902035408u, 1687225157323966154u, 16952679637504852975u, 2143882457656046889u, 6527726468020090836u, 13180247671791137575u, 13426001279318231344u, 8404322666794213895u, 10693596253538186526u, 5848403295842141164u, 9096543884239195665u, 3186869284518994766u, 1187661575047706511u}, {17029606393058198837u, 2777967145943033214u, 2322299207985314407u, 5259038421725444647u, 17517015045803992206u, 106761464208719094u, 10145128318354488652u, 3289166055726423156u, 5226630607526746499u, 12731680892952172648u, 10137449928285087783u, 10166634526199972644u, 5123779428703024877u, 7169394987068918529u, 5568731913585132104u, 12874692959100137209u, 189212297406783382u, 3297475481624254691u}, {14376766508233572188u, 8358269960016871788u, 8530287277877039313u, 9578726314234642503u, 10010331966175180299u, 3930754150488474748u, 11385967387240342354u, 5982572617179164051u, 10595224395498764375u, 4674883413351445145u, 1438054576111920770u, 16336787379973607917u, 13504929292713716426u, 13579600092683583708u, 11003216614271921771u, 5585683770358505497u, 8647465585990731739u, 6539057402515675421u}, {6679059096969625331u, 3413862871215167197u, 13423275448701061128u, 1726440887588789295u, 9383743281720594381u, 8678973550228421465u, 4013006513686260745u, 4819695501675988889u, 18381097535035361480u, 9802746926493950848u, 17877351561637826473u, 13993398748723730417u, 17025044861294578334u, 518843363517857608u, 11592074901386382850u, 14194134867265499690u, 16307543493448916234u, 12016492220143089022u}, {5998760325374448904u, 10860015863404971582u, 15069723739583311808u, 2905569751878202733u, 16382161121512226961u, 12638635940501566356u, 6256535506850769339u, 17699727325121029884u, 3773043620848927026u, 13870244853972624387u, 9371987303198134043u, 2005452045517419843u, 6510765500764420740u, 3330176402693387508u, 876527624518273778u, 16100807597298225796u, 259613905578875349u, 1029552157376562634u}, {7275376532992540622u, 1151812901728088138u, 5368629268678200774u, 4069977565475979428u, 10088474113692679828u, 9916456795832512303u, 8868017759474909113u, 8436412107534371317u, 14515989923221076868u, 7060776686469209262u, 12464171456530387565u, 7851519678210787241u, 13771860660358106471u, 17008689497597630846u, 2109956693735943699u, 4941814780995835190u, 4867143981025967908u, 7443527459437629845u}, {17019786662283625904u, 11224287893149809581u, 13671995505027074633u, 12753298536215910619u, 2908577370391326890u, 115586283461283903u, 3299058579288291336u, 8482950143338025070u, 6699758460220888109u, 13102238573996533005u, 4484955004026473163u, 8280868235191433958u, 12318167549944625621u, 14279226895262359995u, 18113109202597698336u, 16312732142577337318u, 11981707636726854126u, 12362820593348417186u}, {2826022184124566635u, 14403316165118598773u, 15564177650154992038u, 1201649570962552240u, 2972926826169434605u, 57240905433076801u, 932883134736301844u, 16423840307419773461u, 16339563636356840005u, 3447724863152445923u, 13561874365868826281u, 4302558937388293731u, 11615848834193783037u, 10737698823321288073u, 15475493516755247188u, 538610570592255722u, 6178627310091105084u, 9913142066710160115u}, {10415605116315741716u, 17068649228700736610u, 5133546540064216239u, 13005702215006432030u, 2161228247547888678u, 14764295728668374300u, 6715479588351140174u, 12145906191903104114u, 13766525824019881424u, 10154398197877770780u, 5405576648493930438u, 1356860509479676767u, 1785931520949422668u, 13487115705090124451u, 391458791112259473u, 9051693137835232342u, 4426439856253108553u, 17691478043793086952u}, {5133546540064216239u, 13005702215006432030u, 2161228247547888678u, 14764295728668374300u, 6715479588351140174u, 12145906191903104114u, 13766525824019881424u, 10154398197877770780u, 5405576648493930438u, 1356860509479676767u, 1785931520949422668u, 13487115705090124451u, 391458791112259473u, 9051693137835232342u, 4426439856253108553u, 17691478043793086952u, 14431174595012646666u, 8534324614350368305u}, {10802110090528933755u, 5900640435284395293u, 11238287634673515513u, 8363825736845792604u, 1589528607161665969u, 1291398696903690684u, 5710461634839191243u, 4186807060069110058u, 9173392855724899053u, 9538721114646794113u, 1495493398255585105u, 8814121403524388853u, 7426181669545788056u, 4624491042702762655u, 16264046571158817278u, 13187781758609003425u, 5896157662367097993u, 5583435112506755660u}
		},
	B1lo[N][N] = {
{7812344729830787328u, 16460181580158049629u, 16366222863686211595u, 17891067109724719178u, 10229618157069374760u, 4478302601324552700u, 7880618252948522152u, 4902332784916359185u, 17320646795513107470u, 7980761880335200181u, 11744682912095446933u, 5554963389411039413u, 895796881558069816u, 8278155644529319555u, 2424013458343575806u, 13780869086957075575u, 246446421144431782u, 679010219931519137u}, {7724970028543815705u, 5167196668253977955u, 825365234343611059u, 9442764098519105983u, 7279888078040751643u, 13036297019280773264u, 3560167718423233661u, 2837533374638937822u, 3309215015760116778u, 3824665645718876740u, 13988811577703914284u, 2163902276137636470u, 8322460282980804193u, 1907157874088570567u, 8452443521456379137u, 4120896187056503405u, 7812344729830787328u, 246446421144431782u}, {14588502076493490645u, 13430893497949434028u, 8681945219712297854u, 17823481253326337247u, 1512864510673484300u, 3248126113579990833u, 1449580604984506855u, 12675584659677253117u, 16831680894293576558u, 14970506698044573999u, 14975776470549588569u, 7745970522550436829u, 13700059513197898878u, 12026572638552404381u, 7980761880335200181u, 2591823907687008453u, 7724970028543815705u, 7812344729830787328u}, {14096063199118396195u, 2637433306916003451u, 6317972766414656477u, 7390409076963804611u, 6700194207798763248u, 7561427602657588624u, 7774770457521139596u, 9372067076188343157u, 16672970097572197514u, 7891715052823109393u, 7964768259680814696u, 1745003238412056669u, 4580876341000356810u, 15378732243382679106u, 3824665645718876740u, 5462831396623179501u, 14588502076493490645u, 7724970028543815705u}, {11651962624734638320u, 7172694645022574142u, 7994432061988281656u, 17422353137154634694u, 260530736704241260u, 4804699671343179892u, 11280935272782324774u, 7542131862014987016u, 14620923507236451622u, 6838518577632698583u, 9815167132115340745u, 8205997100311947355u, 11571420893607960775u, 13575517387900354089u, 14970506698044573999u, 14259220740814487741u, 14096063199118396195u, 14588502076493490645u}, {17946680416075449129u, 3903035313403071326u, 1305936856040338793u, 3518491066472923691u, 11499263110504817135u, 4705733568972490243u, 2866668279838480158u, 15878785737320088330u, 16931982212510499140u, 6136604706616836835u, 6002050447675964422u, 168417150504412853u, 5197484830307709735u, 13509789980106083141u, 7891715052823109393u, 3589302623622775997u, 11651962624734638320u, 14096063199118396195u}, {8067234080930513105u, 6521355258235600137u, 14921455342909233134u, 13451638015808704339u, 13187001250098363024u, 17160147657401486141u, 4067831922782266328u, 9597585718399059039u, 10299499589414920787u, 11992221230012411032u, 17054984912942639526u, 14546439690884612651u, 4238867068401440245u, 7260396042103376890u, 6838518577632698583u, 11288571941799071682u, 17946680416075449129u, 11651962624734638320u}, {15097576868117044588u, 6732913656851006618u, 16226104797409906921u, 171530286477465621u, 11078572443177220981u, 96745484552685049u, 5519854005817127928u, 2475677653333857192u, 15066805203647834611u, 7169080709038813554u, 1548688972135624559u, 382004621551083354u, 15814136609006981955u, 15819608714800932685u, 6136604706616836835u, 3922557387149924367u, 8067234080930513105u, 17946680416075449129u}, {497188416646485181u, 14705609002489230386u, 10325768932024122907u, 14265955786535969819u, 5636043959223353247u, 16508236299506893916u, 8074644629438272943u, 4838090431883066058u, 9661230952856663536u, 18064522632325444845u, 14240551259715146682u, 12242395927414264814u, 16266157483685265933u, 13259159033201506664u, 11992221230012411032u, 11813448040205610701u, 15097576868117044588u, 8067234080930513105u}, {5011640524814135666u, 15513334720536213021u, 2038115878878510119u, 2037398851144762419u, 15467996330898195346u, 8033168803031326980u, 599714682699790265u, 4215889287866062675u, 11767191438288983000u, 7733901325284040396u, 15890059664001108320u, 13633300026607611879u, 2546987935836785555u, 3089056023386567692u, 7169080709038813554u, 6034341298927733709u, 497188416646485181u, 15097576868117044588u}, {3234411327692083510u, 2305209738918995113u, 5699387672091059038u, 932061486369959472u, 348758974113916518u, 9945510568729811945u, 17044160361927837271u, 16744654924382302943u, 8628532157714155776u, 8253417819245073274u, 6033859886294620750u, 6636668364250939322u, 16671657746450622835u, 17243082658157124127u, 18064522632325444845u, 11514406566311044560u, 5011640524814135666u, 497188416646485181u}, {14560244007462724878u, 6564254699678602821u, 10180952547380499783u, 13033640954585775951u, 5332925630693348171u, 3443194462893089288u, 14453330573426244628u, 8221748569828681626u, 7702232075867648748u, 1734454092810870110u, 11041721550444605160u, 3924987149081897927u, 11596158140851195696u, 7199344906972272438u, 7733901325284040396u, 17497506925909515364u, 3234411327692083510u, 5011640524814135666u}, {6385066765628262387u, 9139233931296733087u, 5835812799588724481u, 6754327172753076689u, 7377669830385517823u, 15967902014753154001u, 6459093804954700187u, 17403753943935234160u, 5475973596659885708u, 5147532354232685940u, 4262007149246902846u, 10518247914631320582u, 4186005503792361066u, 1146920126900987105u, 8253417819245073274u, 9695648002509907072u, 14560244007462724878u, 3234411327692083510u}, {12026542267688242450u, 2980902569702409581u, 11706723300144416426u, 2332165060730801213u, 14191943613183567680u, 15743811680173112026u, 17525588199405649022u, 1377546719156708366u, 16411212531444444152u, 4689724972800769660u, 12960674967336639789u, 17567334215678959664u, 18293259311678852231u, 7531154592416217208u, 1734454092810870110u, 13909413594239376612u, 6385066765628262387u, 14560244007462724878u}, {689059069414593449u, 4850482650718025037u, 4625257478502243008u, 12026159523424607819u, 15422670260774041144u, 6802952823438340935u, 16595753313438107573u, 4795251645421000727u, 5253059026852917007u, 5016863621689866389u, 5215275793911135907u, 2722260800469904730u, 13922285740750326809u, 269657872433240871u, 5147532354232685940u, 16895258750983641901u, 12026542267688242450u, 6385066765628262387u}, {8305431726585820478u, 7052939948514154814u, 1959327341484165664u, 1694199486812421031u, 8911802092814120898u, 17858679870387059769u, 1382995016340309625u, 261538079073770137u, 18381143503894225024u, 5350628697977142520u, 15957691863599849981u, 5019006491552989764u, 15458982410834347207u, 18410431923649245553u, 4689724972800769660u, 7614577019937872627u, 689059069414593449u, 12026542267688242450u}, {1358020439863038274u, 11829769760378740849u, 8882513673059100369u, 235859762310110630u, 11545827641122433199u, 10108423251826223992u, 10067354325683304950u, 7782867720278301754u, 7449102643991025623u, 4848026916687151612u, 9657858622109665129u, 5620954159057696848u, 3495069941431286368u, 5416732151382834530u, 5016863621689866389u, 17569438260412322593u, 8305431726585820478u, 689059069414593449u}, {492892842288863564u, 7267691267516784505u, 13578198133730624292u, 2500419580498744722u, 18347905052174778742u, 10837157722228707648u, 2593707842121113719u, 7217585963917033909u, 7720187745207024817u, 16904887042912758274u, 1807862990279351364u, 17209209662313975169u, 18210180571113052830u, 11782391143718931985u, 5350628697977142520u, 295005248660502759u, 1358020439863038274u, 8305431726585820478u}
		};
static const int64_t Bhi[N][N] = {
{120504028279714, -77781460621766, -529289864238405, 246495106126142, -317628950176215, -324601843700929, 102138683600342, -278105361430707, 141083402911286, -74682024192910, -167600286745817, -75974573036454, -154221974672283, 419146881980658, -290013519940122, -90322987879816, 699091835821457, -7494239115492}, {-92148030197068, 565304274839905, 213102248073222, 305124040935491, -37310002752963, -827705404150772, -355908806327331, 151756892690627, 192168991282597, 4421314917929, 59124649232222, -317792503755251, -169651972319849, 8729619712191, -134205596790396, 195619835976466, -332183889656473, -59117862922795}, {-100827504463973, -4311996683369, 64043269197422, -212600996406453, 67877888722166, -908844647757512, 309247832667714, -176193426733082, 161682401926832, -30602368389831, -505633548852325, -197007276880040, -94863423267333, -229648881817252, -682452888111462, 177490854174517, 199717981380975, 33885251056230}, {-86606882414865, 83669584912573, -56905399020535, 238305057170824, 186199416322663, -42175968209275, -447781734549379, 233515969781439, 332843658129445, -228568101821560, 946850582658149, 272452263183867, -349825346964293, -40276673763716, -149690190781006, -52357113214329, -514411199273876, 116724074224439}, {69147642990481, 184741608668216, -334334272932398, 426062722142120, -67963690477075, 134323998406341, -117522540290527, -446618217854663, -637967699064937, 265878104574521, -119145178507378, 83456543143462, 198068454273664, -151892317518882, 145268875863075, -6767536017895, 832203703029125, 52927898095409}, {213102248073222, 305124040935491, -37310002752963, -827705404150772, -355908806327331, 151756892690627, 192168991282597, 4421314917929, 59124649232222, -317792503755251, -169651972319849, 8729619712191, -134205596790396, 195619835976466, -332183889656473, -59117862922795, -46074015098534, 282652137419952}, {570845422622108, -268532441854109, 35116393841733, -104128986517629, -604195985075146, 429620629614166, 59883964468579, 273928068373409, 145095981764776, -173864767507266, 569933429670675, 420592794619269, -171443754932252, -183211890266303, 180135241985856, -580160838847267, -241345172540198, 129767922048700}, {-517725037341944, 105187891475128, -81139243606740, 665156638995045, -327950319423709, -30486589355765, -35023683307760, -564758198084548, 120785226875211, 74788549052516, -238378501529444, -548247291321067, -18128981801949, 531901871037448, 93003113979025, -4339737133453, -284808135761637, -74529489437900}, {178120904546075, -114330335471438, 57976353928807, -341758717540072, -149944428430997, 431176872861746, -99218316515195, 127736526733943, -158250255213444, -223358154112982, -472009041093530, 272321404841510, -516297981363998, -197617654082613, -57147091744697, 219537312370364, 651863467003538, 353138232436229}, {-151391046215180, 173246085172339, -521063080452434, -101835391378707, 324776004904027, -138233321816375, -183199176284594, -56553163927743, -239530331407047, 476221836280946, -833525214799995, -38500410045668, 564681244991330, -91824907645474, -458790417279016, -117776643254805, 170859077523920, -112277250285127}, {-561702978625799, 291695331532041, -358514315296708, 582023734372333, 203505701569164, -295378629791213, 377137706882024, -236463040035177, 3473623295687, 884252310130238, -294007543191121, 113834168546676, -162840462022585, 147355868031974, -395332817222769, 213105507605365, -37263570458133, 40034144349233}, {-268532441854109, 35116393841733, -104128986517629, -604195985075146, 429620629614166, 59883964468579, 273928068373409, 145095981764776, -173864767507266, 569933429670675, 420592794619269, -171443754932252, -183211890266303, 180135241985856, -580160838847267, -241345172540198, 129767922048700, 285422711311054}, {235553286509609, -341718155047841, 224554500570253, 151391046215179, -173246085172340, 521063080452433, 101835391378706, -324776004904028, 138233321816374, 183199176284593, 56553163927742, 239530331407046, -476221836280947, 833525214799994, 38500410045667, -564681244991331, 91824907645473, 458790417279015}, {458850666495527, -107839137306062, -252886288599194, 387259117477179, -38501019198599, -263943905661166, -178873885247068, 358378255303435, -123483183526613, -94786470174117, -258609201046990, -298451019884963, 533645998958765, 577106195178121, 164291483143305, -506620100161846, -181782154906880, 296892911958380}, {-125401264874111, -282801683653417, 120516386666663, -325108064287420, -781474181042133, 431848295583331, 133837714344667, 56053750009527, -634454534129661, -707289504276542, 387320869913518, 133819565912710, 276165829880738, 86595985106226, 610188715249806, -345207089014555, -571217037601889, -210167900818367}, {567647110225473, 157645588553044, -648137788134794, -106243546628192, 610682734797897, -27492232608272, -164091541506026, -719926452744000, 22767625012179, -121196864081411, 168160735892108, 218084832155552, -375419363523074, 103098506146724, 24437588939126, 87037852109134, 343074274630854, 80872931391994}, {-648137788134794, -106243546628192, 610682734797897, -27492232608272, -164091541506026, -719926452744000, 22767625012179, -121196864081411, 168160735892108, 218084832155552, -375419363523074, 103098506146724, 24437588939126, 87037852109134, 343074274630854, 80872931391994, 283823555112736, 78822794276522}, {370747836626267, -343013747199303, -143553549381154, -217022669352875, -383401038935602, -12334648815398, -527757461461402, 27188939930109, -62072214849188, -149631767863143, 48432859835703, -366689743810883, -31107090643671, 220057424915592, -245146037547339, 283956411708059, 34798916293461, 566475692532689}
		},
	B1hi[N][N] = {
{-4168718777810231416, 7057639375154698997, 1272802723724519386, 3390525983468263436, 3561475714090653878, -397764861216007485, 2035985941170081424, 1467465920171330608, -6939690633089322140, -658848411167338976, -8361043319694276461, 6913020394224356622, -6928164415455031534, -3073387142370463717, -1400478575611675825, 8143374428719719571, 6382016063637001435, -5180590204937266133}, {4931245107932535813, 4524380729525479919, -304827776628204544, -772827655811459478, 7771559307531860244, -2485201232273471841, 856641160671590865, 3926255626340509436, -6854962058552726282, 3580603343091630000, -5603207417914228884, 3972368859906885723, -8408748962441555709, -6476519397177342029, 7006677977648976922, 563349453307889425, -4168718777810231416, 6382016063637001435}, {9004235129870973429, 4649506764444838414, 7393116646156475991, -3612986864262832666, -2976467504525120942, 36449853104780304, 900396790109179046, -5738432866665873323, 8468859452784709317, -7507867063161313973, -2742378702195664414, -560422599686866488, -4362725086933784998, -6312459679243906375, -658848411167338976, 354248438291148947, 4931245107932535813, -4168718777810231416}, {7299745233054286884, 8813814065436046838, -6641892446206056866, -877325020685656753, 109637269918304379, 3992413093653477207, 999366459468964444, 2283038577262285745, -5075235090194321899, 7299969685012678362, -5913554353236280948, 3780111629325590308, 6372894478169659615, 9097382091014771058, 3580603343091630000, 6868642792926728936, 9004235129870973429, 4931245107932535813}, {6451662732926523698, 8994854235659525700, 4383764162418763038, -8120592004404873446, 5388765034991576229, -6926268118269578084, -7477702680984314931, -820057400208351209, 2818849925327208072, -112700059349308877, 7394527329723866672, 624661914011975487, -660915576279234517, 770552846197538895, -7507867063161313973, -6909182659985564115, 7299745233054286884, 9004235129870973429}, {3849697553567728273, 7199465895350087470, 7437062114121245406, 4021218153680227712, 8855803917505346472, 4281441606328268610, -631651875533182746, -245749637545414952, 7166920106816572287, -2288672990366134492, 586557158541393278, 3509870366642892440, -7735334850719845000, -5322644187998130388, 7299969685012678362, -3701266341202101317, 6451662732926523698, 7299745233054286884}, {-1135850700955022356, 4626596998988424272, 2898624138610497531, -7868398319323370301, -7856862250058867197, 8953301255447766362, -8793115546416427548, 5475470820979688147, 7651443751996513762, -8076387084160064789, 4445314854822610569, -1350641764840905814, 7939401239365859620, 916402678139486852, -112700059349308877, -2154285889641375180, 3849697553567728273, 6451662732926523698}, {-4556546210056417529, -7366011024264362933, -1121821176201840287, -699612829499088646, -8638877827473594737, -1030506974649689872, 2484014814540305648, -334700043128520741, 5453014050665409555, -8077159923027954683, -4860563274661344390, 1934641880677264214, 749948462967254560, -8606810713026121758, -2288672990366134492, 1200342696746984347, -1135850700955022356, 3849697553567728273}, {3503445836784229309, 7022455472814929481, 5420946936217936576, -3266506113204245202, -9202081407249975460, -2794604076890835356, -2240765668900113470, -7798242253295885452, -7797469414427995558, -6751529716124361246, 3924101611587376392, 224451958391477, 8734703916369853448, 1073162329866103430, -8076387084160064789, -8898410949676077192, -4556546210056417529, -1135850700955022356}, {5583424934039900662, 340674455725645727, 374030922165477168, -7281250622344189177, 7350249309945350216, -5906456670516664082, 4763121366757046176, 7959052508714232542, 6633422301810639104, -1005000872572813191, 7248320932523345600, -6564362792275832576, -684730979115402920, -8517350428134572540, -8077159923027954683, 4377219128769278381, 3503445836784229309, -4556546210056417529}, {5556743528941383104, -5843522960148736417, 4054277966181010245, 609537015344623128, 8549604509264010767, -4260565057344849012, -4394307880865695705, 3734422390348749811, -2012106453202798243, -1418741803384101067, 7489949100481584883, -6138370543933862766, -25948349644128681, 7376197659589478896, -6751529716124361246, -216150174552259784, 5583424934039900662, 3503445836784229309}, {-156928703221016795, 2026757779348670310, -838699603528266373, -1777947970469469560, 1962257547695085277, 1620810253101549655, -1120597865358665564, -1730270253528623755, -1316529322717335879, -7691257685040186888, 2992995823945608341, -6940536383205042433, -8752674670898569471, -1144464893546971379, -1005000872572813191, 1698116992310160825, 5556743528941383104, 5583424934039900662}, {5560757284258862849, 4871829990908509723, 1790193118524720777, 6385432178527799563, 6223581236936318548, -5054009669720903214, 6688353902158498898, -6099868606105857621, 172647275550228201, 7517056651271754450, -7368069723555873616, -3520613712971537155, -7971801567993428722, -8512684958021037617, -1418741803384101067, -6200913587155280736, -156928703221016795, 5556743528941383104}, {796987555019198791, 4602718965749820021, -3537830603201967250, -8588285637582940616, -2220133780867289299, 1499637283335588445, -8270923176962891640, 3385610052303368883, 6624039789700979160, -3823715730241639580, 8933544567341409390, 8191768520800961061, -4951084896088904659, -7272077940684328410, -7691257685040186888, 6979889754789473610, 5560757284258862849, -156928703221016795}, {3525787966115907121, -5882776954552224966, -6770760124191533344, 5274710423488003234, 4514777263334758271, -770793658385039328, -68763661457236672, 6163203437153038650, -942768255043118936, -6976218387779873049, -8641238884030054152, -6588425806612713853, -6455509120464393890, 4235375339597104591, 7517056651271754450, -8197493932671682188, 796987555019198791, 5560757284258862849}, {-3374481343709267547, -1037003684411430883, -663366331305465255, -4218256018824814602, -4094376278692388927, 3477132947611256846, -6026801555218990485, -5440218140491587605, -2287715482953354137, -4495326543083000522, -7629709927515293833, -6975485332325484171, 3264422878327383972, 1250354385007374298, -3823715730241639580, -5417129167049032294, 3525787966115907121, 796987555019198791}, {8085563663835019350, 4765978493075089257, -7632446146653117362, -7443055904166015232, -3785641442968922471, -1141294225814307783, 2155586899598297009, 3247883603234357207, 766991758537484680, -2800957151223351650, 7239312079181944399, -7534328981819170094, -2857485754901007775, -1236890523791381152, -6976218387779873049, -3502934665846898089, -3374481343709267547, 3525787966115907121}, {-5682711946435548746, 4547668850597985080, -1781759160640056639, -5323642507812490090, -1556359041941438886, 6586616334666149583, 3974648110049603658, -7763904556865925326, 8988470124983977418, -4433388118411597772, 9196739556246403156, 1956299367012891600, -8672147078562422231, 4886505635593499759, -4495326543083000522, -7681885119216656864, 8085563663835019350, -3374481343709267547}
		};


#endif