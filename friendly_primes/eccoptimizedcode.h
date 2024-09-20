#ifndef ECCOPTIMIZEDCODE_H
#define ECCOPTIMIZEDCODE_H

void multMod25519(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModM383(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModC41417(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModEd448(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModM511(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModE521(uint64_t* output, uint64_t * in, uint64_t * in2);

void multMod244(uint64_t* output, uint64_t * in, uint64_t * in2);
void multMod244v2(uint64_t* output, uint64_t * in, uint64_t * in2);
void multMod297(uint64_t* output, uint64_t * in, uint64_t * in2);
void multMod354(uint64_t* output, uint64_t * in, uint64_t * in2);
void multMod480(uint64_t* output, uint64_t * in, uint64_t * in2);

#endif
