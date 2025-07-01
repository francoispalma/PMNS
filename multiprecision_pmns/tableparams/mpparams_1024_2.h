#define N 9
#define NBCHUNKS 2
#define LAMBDA 2
#define RHOver2 59

static const uint64_t RHOhi = 56658347597845894u;

static const int64_t M[2][17] = { {
		2057075501242702, 218606003891118846, 547281604355546076, 527102666073639754, 437351122802424174, 555426076915698772, 358359893266728720, 221123485733968416, 94749002638713627, 289258913902333095, 397533378097271167, 561871178329484782, 551781709188531621, 218675561401212087, 565943414609561130, 179179946633364360, 110561742866984208
	}, {
-972388893042483, 13007816906647319, 1954655172407619, 18918932767587661, 9499580449486904, -31049310289814293, 911131010172082, 11996362221285386, -25006517485248040, -486194446521242, 6503908453323659, 977327586203809, 9459466383793830, 4749790224743452, -15524655144907147, 455565505086041, 5998181110642693
} };

static const int64_t M1[] = {
111733446451438278, 130586174345112656, 163025183127864450, 187643182845796512, 447543211848197562, 236144087179884320, 80750179380517052, 324340445475130454, 141811040085303793, 55866723225719139, 65293087172556328, 369742967715643969, 93821591422898256, 223771605924098781, 406302419741653904, 40375089690258526, 450400598889276971
};

static inline void mtoep10(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N/3)

static inline void m1toep10(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1SCHOOLBOOK(N/3)


static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP33TOP(N, mtoep10)

static inline void ptoeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	pTOEP33TOP(N, mtoep10)

static inline void m1toep20(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1TOEP33TOP(N, m1toep10)


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

