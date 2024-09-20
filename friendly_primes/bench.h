#ifndef BENCH_H_INCLUDED
#define BENCH_H_INCLUDED

#define LOOPCHKNMB 1000000
#include "structs.h"

uint64_t do_bench(void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly), const uint16_t N, const uint64_t RHO, const uint64_t W);
uint64_t do_gmpbench(void (*gmp_mult)(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMBMASK, const mp_limb_t MPLIMB), const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMBMASK, const mp_limb_t MPLIMB, const uint64_t W);
uint64_t do_mersennebench(void (*gmp_mult)(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[]), const uint64_t W);
uint64_t do_pmersbench(void (*pmersmult)(uint64_t c[], uint64_t a[], uint64_t b[]), const uint64_t W, const uint64_t NBLIMB, void (*convertfunc)(uint64_t a[], uint64_t b[]));
void gmpmulmod2k(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMB, const mp_limb_t MPLIMB);
void gmpmulmod2kbig(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[], const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMB, const mp_limb_t MPLIMB);
void mersenne521(mp_limb_t c[], mp_limb_t a[], mp_limb_t b[]);

void PoC244convert(uint64_t a[], uint64_t alang[]);
void PoC244convertv2(uint64_t a[], uint64_t alang[]);
void PoC297convert(uint64_t a[], uint64_t alang[]);
void PoC354convert(uint64_t a[], uint64_t alang[]);
void PoC480convert(uint64_t a[], uint64_t alang[]);
void C25519convert(uint64_t a[], uint64_t alang[]);
void M383convert(uint64_t a[], uint64_t alang[]);
void C41417convert(uint64_t a[], uint64_t b[]);
void Ed448convert(uint64_t a[], uint64_t b[]);
void M511convert(uint64_t a[], uint64_t b[]);
void E521convert(uint64_t a[], uint64_t b[]);

uint64_t checkmodmul(const mp_limb_t KPRIMEC, const mp_limb_t LASTLIMB, const mp_limb_t MPLIMB);
int64_t checkOpti(const uint16_t SIZEP, const uint64_t PMASK, const uint64_t PSUB, const uint16_t NSIZE, void (*convertfunc)(uint64_t a[], uint64_t b[]), void (*pmersmult)(uint64_t c[], uint64_t a[], uint64_t b[]));
int64_t checkPmns(const uint16_t N, const uint64_t RHO, const uint64_t GAMMA, mp_limb_t prime[], void (*pmns_mult)(restrict poly, const restrict poly, const restrict poly));

#endif
