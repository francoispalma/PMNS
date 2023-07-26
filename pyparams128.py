from commonpmns import pmnsWBdicts
from commonpmns import primesdict
from generated.BKZpmns8192128 import pmnsdict
#from generated.pmns1024128 import pmnsdict
#from tmp8192128 import pmns128dict as pmnsdict
from findm import findm

siz = 1024
#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsWBdicts[siz * 1000 + 128][primesdict[siz][0]]
phi = 2**128

primes = list(pmnsdict.keys())

#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdict[primes[0]] # n = 73 XnX1 in tmp

(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdict[primes[4]] #n = 72, in BKZ

#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdict[primes[1]] # n = 36
#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdict[primes[16]] # n = 36 XnX1

#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdict[primes[1]] # n = 18
#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdict[primes[291]] # n = 18 XnX1

#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = pmnsdict[primes[0]] # n = 9

#M_or_B, M1_or_B1 = findm(p, n, gamma, lam, rho, M_or_B, M1_or_B1, phi)

#lam = [lam]

#1664
#(p, n, gamma, lam, rho, M_or_B, M1_or_B1) = [474553065352780189310452755182331644834508592750533842937874345731334915146234261501522574808132104494743319377192885948653558902010649427117192895702514187551482120028413325852756890745393760250298693189517277523145025312193179015452155475616637521046961313888177203793142395465988380733731832164398194308005629753432606988166651587042593773270097777354576357463115463309405321999803927219065702539470856455802056793761608582648983106081845067729494912560380898415586731176634727418923713142246780559,
#15,
#456581101466483487340221467617318557370791639125190250901460472859083453743725524505226769201894828388155467190121163318356141080154347479465155815013845109741093381120157700457359391376764025652215928124963446744563289074609819669825359927509760216535825784192871712417312216689684225655952623341268112043907048130668098209432174048057156305146355930839569818053477906086917277120180510508722410682006783288823405525353740313823907868230388152909688336736753691643170808515553725985998234576389943107,
#12,
#118,
#[-382837356738320151909192039331849, 821424127318151282128136314178416, -296085021691196253854211699827900, -56134193052630029476248545690273, -25454695038184416145035795749880, 421785273796699401935889999841574, 992059024646600838549144519686158, -454988792596567170655074255220896, 410769225221808071486676045468691, -1043295540990889811472974679248161, 1071703967161051114602388095452925, 1246049059249217086583998957174801, 600326163866251507977387169980393, -883838602903953404006036653137916, -4148185930555168281019388013775],
#[121897930843180443110771309865212641297, 140326010033178531395067048013662004160, -156660515535317672153258034069106820640, -74627135765825397005971553736138170869, -100517335577786143078824260935665106412, 116857790044374803722346039063303791994, 144114839416902779323053745319253837371, -151304674028743981846742656054021503048, 115754034077543296821875974097219339395, -160678587815357282067570484220500032974, 79416448328561064018971769343617896505, -92567131706565730670469153883587783947, -156522387086451086681305740181514693770, 100516607113792905755127540083311556314, -19217165342075806096135706397138348536]]



