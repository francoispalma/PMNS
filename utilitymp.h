#ifndef UTILITYMP_H
#define UTILITYMP_H

void __print128(__int128);
void mp_print(const restrict poly);
void convert_string_to_amns(restrict poly, const char*);
void convert_amns_to_poly(restrict poly*, const restrict poly);

#endif
