#ifndef UTILITYMP_CORE_H
#define UTILITYMP_CORE_H

void __print128(__int128);
void mp_print(const restrict poly);
void convert_string_to_multipre(restrict poly*, const char*);
void mp_copy(restrict poly*, restrict const poly);
int8_t mp_comp(restrict poly, restrict poly);
void mp_leftshift(restrict poly*);
void mp_rightshift(restrict poly);
void mp_alignleft(restrict poly*, uint16_t);
void mp_add(restrict poly*, restrict const poly, restrict const poly);
void mp_sub(restrict poly*, restrict const poly, restrict const poly);
void mp_mult(restrict poly*, restrict const poly, restrict const poly);
void mp_mod(restrict poly*, restrict const poly, restrict const poly);

#endif
