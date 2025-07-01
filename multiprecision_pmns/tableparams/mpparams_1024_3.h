#define N 6
#define NBCHUNKS 3
#define LAMBDA 2
#define RHOver2 59

static const uint64_t RHOhi = 10565647096983743u;

static const int64_t M[3][11] = { {
		96445245815240584, 316565494012253068, 119627053769411638, 163975446349048938, 456817811811266578, 251722168629509213, 48222622907620292, 158282747006126534, 348043903036417563, 370218099326236213, 516639282057345033
	}, {
		566408934989478172, 69800799140443740, 32814931887751129, 317699452827019499, 320071741768797523, 14301428553985666, 283204467494739086, 323130775721933614, 16407465943875564, 158849726413509749, 448266247036110505
	}, {
6897749540641076, 10294432740438585, 1254163267683738, 1786660513642644, -226774358821603, 671513772739839, 3448874770320538, 5147216370219292, 627081633841869, 893330256821322, -113387179410802
} };

static const int64_t M1[] = {
356743099159279152, 253301482161708832, 355309363641048794, 283610169609023234, 276939852823379266, 274816468548000833, 466601925731351320, 414881117232566160, 465885057972236141, 141805084804511617, 426700302563401377
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



static inline void mtoeplitz_vm2(__int128* restrict rop, const int64_t* restrict vect)
{
	toeplitz_vm(rop, vect, M[2]);
}



static inline void m1toeplitz_vm(int64_t* restrict rop, const __int128* restrict vect)
{
	int64_t v[N];
	for(int i = 0; i < N; i++)
		v[i] = vect[i];
	m1toep20(rop, v, M1);
}


void (* const multbym[3]) (__int128* restrict rop, const int64_t* restrict vect) = { mtoeplitz_vm0, mtoeplitz_vm1, mtoeplitz_vm2 };

