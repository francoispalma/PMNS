#ifndef UTILITYMP_CORE_H
#define UTILITYMP_CORE_H

void __print128(__int128);
void convert_string_to_multipre(restrict mpnum*, const char*);
_Bool mp_iszero(const restrict mpnum);
int8_t mp_comp(restrict mpnum, restrict mpnum);
void mp_leftshift(restrict mpnum*);
void mp_rightshift(restrict mpnum);
void mp_alignleft(restrict mpnum*, uint16_t);
void mp_add(restrict mpnum*, restrict const mpnum, restrict const mpnum);
void mp_sub(restrict mpnum*, restrict const mpnum, restrict const mpnum);
void mp_uadd(restrict mpnum*, restrict const mpnum, restrict const mpnum);
void mp_usub(restrict mpnum*, restrict const mpnum, restrict const mpnum);
void mp_mult(restrict mpnum*, restrict const mpnum, restrict const mpnum);
void mp_mod(restrict mpnum*, restrict const mpnum, restrict const mpnum);
void mp_utmod(restrict mpnum*, restrict const mpnum, restrict const mpnum);

#endif
