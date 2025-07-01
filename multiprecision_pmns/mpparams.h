#define N 10
#define NBCHUNKS 2
#define LAMBDA 2
#define RHOver2 54

static const uint64_t RHOhi = 873421328803747u;

static const int64_t M[2][19] = { {
		6422187908595974, 17117383092504324, 3670584704703816, 4671739493311360, 3474872017372778, 7295686071961204, 14459207993696308, 13411012672650260, 2716106671599424, 12121510176586993, 12218293209038979, 17565890800993154, 1835292352351908, 2335869746655680, 10744635263427381, 3647843035980602, 7229603996848154, 6705506336325130, 10365252590540704
	}, {
-50345473046631, -283382866594249, -339409645786896, -6717742189740, -322972412790991, -10368688049592, -95199516250230, -472798579662628, -143266483518599, -22381249717943, -25172736523316, -141691433297125, -169704822893448, -3358871094870, -161486206395496, -5184344024796, -47599758125115, -236399289831314, -71633241759300
} };

static const int64_t M1[] = {
7351057839742562, 9397544450031902, 15803276564548878, 4412649804587610, 17625867920583496, 11804776505780734, 8196505759245976, 12502075737265114, 8417859779598044, 7977340740016611, 12682728174612273, 4698772225015951, 16908837537015431, 11213524157034797, 8812933960291748, 5902388252890367, 4098252879622988, 6251037868632557, 4208929889799022
};

static inline void mtoep10(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N/2)

static inline void m1toep10(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1SCHOOLBOOK(N/2)


static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP22TOP(N, mtoep10)

static inline void ptoeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	pTOEP22TOP(N, mtoep10)

static inline void m1toep20(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1TOEP22TOP(N, m1toep10)


static inline void mtoeplitz_vm0(__int128* restrict rop, const int64_t* restrict vect)
{
	toeplitz_vm(rop, vect, M[0]);
}



static inline void mtoeplitz_vm1(__int128* restrict rop, const int64_t* restrict vect)
{
	toeplitz_vm(rop, vect, M[1]);
}



static inline void m1toeplitz_vm(int64_t* restrict rop, const __int128* restrict vect)
{
	int64_t v[N];
	for(int i = 0; i < N; i++)
		v[i] = vect[i];
	m1toep20(rop, v, M1);
}


void (* const multbym[2]) (__int128* restrict rop, const int64_t* restrict vect) = { mtoeplitz_vm0, mtoeplitz_vm1 };

