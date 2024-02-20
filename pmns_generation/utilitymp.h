#ifndef UTILITYMP_CORE_H
#define UTILITYMP_CORE_H

void __print128(__int128);
void convert_string_to_binary(mpnum*, const char*);
_Bool mp_iszero(const mpnum);
int8_t mp_comp(mpnum, mpnum);
void mp_leftshift(mpnum*);
void mp_rightshift(mpnum);
void mp_alignleft(mpnum*, uint16_t);
void mp_add(mpnum*, const mpnum, const mpnum);
void mp_sub(mpnum*, const mpnum, const mpnum);
void mp_uadd(mpnum*, const mpnum, const mpnum);
void mp_usub(mpnum*, const mpnum, const mpnum);
void mp_mult(mpnum*, const mpnum, const mpnum);
void mp_mod(mpnum*, const mpnum, const mpnum);
void mp_utmod(mpnum*, const mpnum, const mpnum);

#endif

