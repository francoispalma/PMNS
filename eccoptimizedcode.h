#ifndef ECCOPTIMIZEDCODE_H
#define ECCOPTIMIZEDCODE_H

void multMod25519(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModM383(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModC41417(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModEd448(uint64_t* output, uint64_t * in, uint64_t * in2);
void multModCE521(uint64_t* output, uint64_t * in, uint64_t * in2);

#endif
