#ifndef PMNS_PARAMS_H_INCLUDED
#define PMNS_PARAMS_H_INCLUDED

#define RHO 54
#define N 5
#define _PRAGMAGCCUNROLLLOOP_ _Pragma("GCC unroll 5")
#define LAMBDA 2

static const int64_t B[N][N] = {
{-1377551199037903, -955712434975478, 719787268038937, -941057892625767, -337855473576492}, {675710947152984, 1377551199037903, 955712434975478, -719787268038937, 941057892625767}, {-504564586213631, 280001487822494, -2097338467076840, -14654542349711, 1057642741615429}, {1439574536077874, -1882115785251534, -675710947152984, -1377551199037903, -955712434975478}, {-29309084699422, 2115285483230858, -504564586213631, 280001487822494, -2097338467076840}
		},
	B1[N][N] = {
{9093953327563851615u, 3196985667173371090u, 3556254531690257572u, -3108658914698264879u, 5091040640737529957u}, {3876436699603423945u, -5537698795873594043u, -6305644581871635969u, -1534786109047272385u, -3108658914698264879u}, {-2214405142570208100u, 8264662792234491702u, 4002912686826321658u, -3196985667173371090u, -1534786109047272385u}, {3915523396207144054u, 6217317829396529758u, 6985095614301688824u, 5537698795873594043u, -3196985667173371090u}, {-7073590367869677895u, 3069572218094544770u, -679619033522935715u, -8264662792234491702u, 5537698795873594043u}
		};

static const int64_t M[N] = {-1377551199037903, -955712434975478, 719787268038937, -941057892625767, -337855473576492},
	M1[N] = {5537698795873594043, -3196985667173371090, -1534786109047272385, -3108658914698264879, 5091040640737529957},
	matrM[9] = {-1911424869950956, 1439574536077874, -1882115785251534, -675710947152984, -1377551199037903, -955712434975478, 719787268038937, -941057892625767, -337855473576492},
	matrM1[9] = {-6393971334346742180, -3069572218094544770, -6217317829396529758, -8264662792234491702, 5537698795873594043, -3196985667173371090, -1534786109047272385, -3108658914698264879, 5091040640737529957};


#endif