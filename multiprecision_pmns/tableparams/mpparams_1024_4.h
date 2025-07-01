#define N 5
#define NBCHUNKS 4
#define LAMBDA 2
#define RHOver2 53

static const uint64_t RHOhi = 82963814466964u;

static const int64_t M[4][9] = { {
		366719323402696, 7646832282300524, 5796106324703996, 5300204987831944, 8185043839211171, 183359661701348, 3823416141150262, 7401652789722494, 7153702121286468
	}, {
		4210182939788582, 3659659364230074, 5412251772586437, 6522101678724401, 3483005597633375, 6608691097264787, 1829829682115037, 7209725513663714, 3261050839362200
	}, {
		4941167828394579, 7971707500426092, 5935424314443747, 5627311745587440, 1016847567857192, 6974183541567785, 3985853750213046, 7471311784592369, 7317255500164216
	}, {
123041560673823, 8373904759096, 1832824802999, 15636509616639, 17042829081370, 61520780336911, 4186952379548, 916412401499, 7818254808319
} };

static const int64_t M1[] = {
5324185585587448, 91185611302988, 98608890784476, 6644188667831152, 6088080244074277, 7165692420164220, 45592805651494, 49304445392238, 3322094333915576
};

static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N)

static inline void ptoeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	pSCHOOLBOOK(N)

static inline void m1toep20(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1SCHOOLBOOK(N)


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



static inline void mtoeplitz_vm3(__int128* restrict rop, const int64_t* restrict vect)
{
	toeplitz_vm(rop, vect, M[3]);
}



static inline void m1toeplitz_vm(int64_t* restrict rop, const __int128* restrict vect)
{
	int64_t v[N];
	for(int i = 0; i < N; i++)
		v[i] = vect[i];
	m1toep20(rop, v, M1);
}


void (* const multbym[4]) (__int128* restrict rop, const int64_t* restrict vect) = { mtoeplitz_vm0, mtoeplitz_vm1, mtoeplitz_vm2, mtoeplitz_vm3 };

