#ifndef PMNS_PARAMS_H_INCLUDED
#define PMNS_PARAMS_H_INCLUDED

#define RHO 54
#define N 10
#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 10")
#define LAMBDA 5

static const int64_t M[N] = {469799818280136, -152080239661648, 598994790531987, 516377690304149, 1142053037191459, -784482490732716, -552314955453948, 566486891798520, 1487472002814664, -1196753499178744},
	M1[N] = {-6029253850808398197, 56231205061316773, 6792934487449497806, -2765138884635369035, -4028868174327776185, -6846182431872386498, 4695272658421107979, -1792797856598094539, -8261918006478193905, -8203057425321911138},
	matrM[19] = {-760401198308240, 2994973952659935, 2581888451520745, 5710265185957295, -3922412453663580, -2761574777269740, 2832434458992600, 7437360014073320, -5983767495893720, 469799818280136, -152080239661648, 598994790531987, 516377690304149, 1142053037191459, -784482490732716, -552314955453948, 566486891798520, 1487472002814664, -1196753499178744},
	matrM1[19] = {281156025306583865, -2928815710171614202, 4621049650532706441, -1697596797929329309, 2662575988057170742, 5029619218395988279, -8963989282990472695, -4416101884971866293, -4121798979190452458, -6029253850808398197, 56231205061316773, 6792934487449497806, -2765138884635369035, -4028868174327776185, -6846182431872386498, 4695272658421107979, -1792797856598094539, -8261918006478193905, -8203057425321911138};


#endif
