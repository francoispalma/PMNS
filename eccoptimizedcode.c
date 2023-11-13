#include <stdint.h>

#include "eccoptimizedcode.h"

void multMod25519(uint64_t* output, uint64_t * in, uint64_t * in2) {
	/* Copyright 2008, Google Inc.
	 * All rights reserved.
	 *
	 * Code released into the public domain by Adam Langley <agl@imperialviolet.org>
	 *
	 * See https://code.google.com/archive/p/curve25519-donna/
	 *
	 * Derived from public domain C code by Daniel J. Bernstein <djb@cr.yp.to>
	 */
	unsigned __int128 t[5];
	uint64_t r0,r1,r2,r3,r4,s0,s1,s2,s3,s4,c;
	
	r0 = in[0];
	r1 = in[1];
	r2 = in[2];
	r3 = in[3];
	r4 = in[4];
	
	s0 = in2[0];
	s1 = in2[1];
	s2 = in2[2];
	s3 = in2[3];
	s4 = in2[4];
	
	t[0] = ((unsigned __int128) r0) * s0;
	t[1] = ((unsigned __int128) r0) * s1 + ((unsigned __int128) r1) * s0;
	t[2] = ((unsigned __int128) r0) * s2 + ((unsigned __int128) r2) * s0 + ((unsigned __int128) r1) * s1;
	t[3] = ((unsigned __int128) r0) * s3 + ((unsigned __int128) r3) * s0 + ((unsigned __int128) r1) * s2 + ((unsigned __int128) r2) * s1;
	t[4] = ((unsigned __int128) r0) * s4 + ((unsigned __int128) r4) * s0 + ((unsigned __int128) r3) * s1 + ((unsigned __int128) r1) * s3 + ((unsigned __int128) r2) * s2;
	
	r4 *= 19;
	r1 *= 19;
	r2 *= 19;
	r3 *= 19;
	
	t[0] += ((unsigned __int128) r4) * s1 + ((unsigned __int128) r1) * s4 + ((unsigned __int128) r2) * s3 + ((unsigned __int128) r3) * s2;
	t[1] += ((unsigned __int128) r4) * s2 + ((unsigned __int128) r2) * s4 + ((unsigned __int128) r3) * s3;
	t[2] += ((unsigned __int128) r4) * s3 + ((unsigned __int128) r3) * s4;
	t[3] += ((unsigned __int128) r4) * s4;
	
							r0 = (uint64_t)t[0] & 0x7ffffffffffff; c = (uint64_t)(t[0] >> 51);
	t[1] += c;	r1 = (uint64_t)t[1] & 0x7ffffffffffff; c = (uint64_t)(t[1] >> 51);
	t[2] += c;	r2 = (uint64_t)t[2] & 0x7ffffffffffff; c = (uint64_t)(t[2] >> 51);
	t[3] += c;	r3 = (uint64_t)t[3] & 0x7ffffffffffff; c = (uint64_t)(t[3] >> 51);
	t[4] += c;	r4 = (uint64_t)t[4] & 0x7ffffffffffff; c = (uint64_t)(t[4] >> 51);
	r0 += c * 19;	c = r0 >> 51; r0 = r0 & 0x7ffffffffffff;
	r1 += c;			c = r1 >> 51; r1 = r1 & 0x7ffffffffffff;
	r2 += c;
	
	output[0] = r0;
	output[1] = r1;
	output[2] = r2;
	output[3] = r3;
	output[4] = r4;
}

void multModM383(uint64_t* output, uint64_t * in, uint64_t * in2)
{
	/* Adapted from multMod25519 for M-383
	 */
	unsigned __int128 t[7];
	uint64_t r0,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,c;
	unsigned __int128 c2;
	
	r0 = in[0]*187;
	r1 = in[1]*187;
	r2 = in[2]*187;
	r3 = in[3]*187;
	r4 = in[4]*187;
	r5 = in[5]*187;
	r6 = in[6]*187;
	
	s1 = in2[1]*4;
	s2 = in2[2]*4;
	s3 = in2[3]*4;
	s4 = in2[4]*4;
	s5 = in2[5]*4;
	s6 = in2[6]*4;
	
	t[0] = ((((unsigned __int128) r6) * s1 + ((unsigned __int128) r1) * s6 + ((unsigned __int128) r2) * s5 + ((unsigned __int128) r5) * s2 + ((unsigned __int128) r4) * s3 + ((unsigned __int128) r3) * s4));
	t[1] = ((((unsigned __int128) r6) * s2 + ((unsigned __int128) r2) * s6 + ((unsigned __int128) r3) * s5 + ((unsigned __int128) r5) * s3 + ((unsigned __int128) r4) * s4));
	t[2] = ((((unsigned __int128) r6) * s3 + ((unsigned __int128) r3) * s6 + ((unsigned __int128) r5) * s4 + ((unsigned __int128) r4) * s5));
	t[3] = ((((unsigned __int128) r6) * s4 + ((unsigned __int128) r4) * s6 + ((unsigned __int128) r5) * s5));
	t[4] = ((((unsigned __int128) r6) * s5 + ((unsigned __int128) r5) * s6));
	t[5] = ((((unsigned __int128) r6) * s6));
	
	t[0] += ((unsigned __int128) in[0]) * in2[0];
	t[1] += ((unsigned __int128) in[0]) * in2[1] + ((unsigned __int128) in[1]) * in2[0];
	t[2] += ((unsigned __int128) in[0]) * in2[2] + ((unsigned __int128) in[2]) * in2[0] + ((unsigned __int128) in[1]) * in2[1];
	t[3] += ((unsigned __int128) in[0]) * in2[3] + ((unsigned __int128) in[3]) * in2[0] + ((unsigned __int128) in[1]) * in2[2] + ((unsigned __int128) in[2]) * in2[1];
	t[4] += ((unsigned __int128) in[0]) * in2[4] + ((unsigned __int128) in[4]) * in2[0] + ((unsigned __int128) in[3]) * in2[1] + ((unsigned __int128) in[1]) * in2[3] + ((unsigned __int128) in[2]) * in2[2];
	t[5] += ((unsigned __int128) in[0]) * in2[5] + ((unsigned __int128) in[1]) * in2[4] + ((unsigned __int128) in[2]) * in2[3] + ((unsigned __int128) in[3]) * in2[2] + ((unsigned __int128) in[4]) * in2[1] + ((unsigned __int128) in[5]) * in2[0];
	t[6] = ((unsigned __int128) in[0]) * in2[6] + ((unsigned __int128) in[1]) * in2[5] + ((unsigned __int128) in[2]) * in2[4] + ((unsigned __int128) in[3]) * in2[3] + ((unsigned __int128) in[4]) * in2[2] + ((unsigned __int128) in[5]) * in2[1] + ((unsigned __int128) in[6]) * in2[0];
	
/*	r1 *= 187; //s1 *= 4;*/
/*	r2 *= 187; //s2 *= 4;*/
/*	r3 *= 187; //s3 *= 4;*/
/*	r4 *= 187; //s4 *= 4;*/
/*	r5 *= 187; //s5 *= 4;*/
/*	r6 *= 187; //s6 *= 4;*/
/*	*/
/*	t[0] += 4*(((unsigned __int128) r6) * s1 + ((unsigned __int128) r1) * s6 + ((unsigned __int128) r2) * s5 + ((unsigned __int128) r5) * s2 + ((unsigned __int128) r4) * s3 + ((unsigned __int128) r3) * s4);*/
/*	t[1] += 4*(((unsigned __int128) r6) * s2 + ((unsigned __int128) r2) * s6 + ((unsigned __int128) r3) * s5 + ((unsigned __int128) r5) * s3 + ((unsigned __int128) r4) * s4);*/
/*	t[2] += 4*(((unsigned __int128) r6) * s3 + ((unsigned __int128) r3) * s6 + ((unsigned __int128) r5) * s4 + ((unsigned __int128) r4) * s5);*/
/*	t[3] += 4*(((unsigned __int128) r6) * s4 + ((unsigned __int128) r4) * s6 + ((unsigned __int128) r5) * s5);*/
/*	t[4] += 4*(((unsigned __int128) r6) * s5 + ((unsigned __int128) r5) * s6);*/
/*	t[5] += 4*(((unsigned __int128) r6) * s6);*/
	
	//printf("[0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx]\n", (uint64_t)(t[0]>>64), (uint64_t)t[0], (uint64_t)(t[1]>>64), (uint64_t)t[1], (uint64_t)(t[2]>>64), (uint64_t)t[2], (uint64_t)(t[3]>>64), (uint64_t)t[3], (uint64_t)(t[4]>>64), (uint64_t)t[4], (uint64_t)(t[5]>>64), (uint64_t)t[5], (uint64_t)(t[6]>>64), (uint64_t)t[6]);
	
							r0 = (uint64_t)t[0] & 0x7fffffffffffff; c2 = (t[0] >> 55);
	t[1] += c2;	r1 = (uint64_t)t[1] & 0x7fffffffffffff; c2 = (t[1] >> 55);
	t[2] += c2;	r2 = (uint64_t)t[2] & 0x7fffffffffffff; c2 = (t[2] >> 55);
	t[3] += c2;	r3 = (uint64_t)t[3] & 0x7fffffffffffff; c2 = (t[3] >> 55);
	t[4] += c2;	r4 = (uint64_t)t[4] & 0x7fffffffffffff; c = (uint64_t)(t[4] >> 55);
	t[5] += c;	r5 = (uint64_t)t[5] & 0x7fffffffffffff; c = (uint64_t)(t[5] >> 55);
	t[6] += c;	r6 = (uint64_t)t[6] & 0x1fffffffffffff; c2 =(t[6] >> 53)*187+r0;
	r0 = (uint64_t)c2 & 0x7fffffffffffff;	c = (c2>>55) + r1;
	r1 = (uint64_t)c & 0x7fffffffffffff;	c = c>>55;
	r2 += c;
	
	output[0] = r0;
	output[1] = r1;
	output[2] = r2;
	output[3] = r3;
	output[4] = r4;
	output[5] = r5;
	output[6] = r6;
}

void multModC41417(uint64_t* output, uint64_t * in, uint64_t * in2)
{
	/* Adapted from multMod25519 for Curve41417
	 */
	unsigned __int128 t[7];
	uint64_t r0,r1,r2,r3,r4,r5,r6,s0,s1,s2,s3,s4,s5,s6,c;
	unsigned __int128 c2;
	
	r0 = in[0];
	r1 = in[1];
	r2 = in[2];
	r3 = in[3];
	r4 = in[4];
	r5 = in[5];
	r6 = in[6];
	
	s0 = in2[0];
	s1 = in2[1];
	s2 = in2[2];
	s3 = in2[3];
	s4 = in2[4];
	s5 = in2[5];
	s6 = in2[6];
	
	t[0] = ((unsigned __int128) r0) * s0;
	t[1] = ((unsigned __int128) r0) * s1 + ((unsigned __int128) r1) * s0;
	t[2] = ((unsigned __int128) r0) * s2 + ((unsigned __int128) r2) * s0 + 2*((unsigned __int128) r1) * s1;
	t[3] = ((unsigned __int128) r0) * s3 + ((unsigned __int128) r3) * s0 + 2*(((unsigned __int128) r1) * s2 + ((unsigned __int128) r2) * s1);
	t[4] = ((unsigned __int128) r0) * s4 + ((unsigned __int128) r4) * s0 + 2*(((unsigned __int128) r3) * s1 + ((unsigned __int128) r1) * s3 + ((unsigned __int128) r2) * s2);
	t[5] = ((unsigned __int128) r0) * s5 + 2*(((unsigned __int128) r1) * s4 + ((unsigned __int128) r2) * s3 + ((unsigned __int128) r3) * s2 + ((unsigned __int128) r4) * s1) + ((unsigned __int128) r5) * s0;
	t[6] = ((unsigned __int128) r0) * s6 + 2*(((unsigned __int128) r1) * s5 + ((unsigned __int128) r2) * s4 + ((unsigned __int128) r3) * s3 + ((unsigned __int128) r4) * s2 + ((unsigned __int128) r5) * s1) + ((unsigned __int128) r6) * s0;
	
	r1 *= 17;
	r2 *= 17;
	r3 *= 17;
	r4 *= 17;
	r5 *= 17;
	r6 *= 17;
	
	t[0] += (((unsigned __int128) r6) * s1 + ((unsigned __int128) r1) * s6 + ((unsigned __int128) r2) * s5 + ((unsigned __int128) r5) * s2 + ((unsigned __int128) r4) * s3 + ((unsigned __int128) r3) * s4)*2;
	t[1] += ((unsigned __int128) r6) * s2 + ((unsigned __int128) r2) * s6 + ((unsigned __int128) r3) * s5 + ((unsigned __int128) r5) * s3 + ((unsigned __int128) r4) * s4;
	t[2] += ((unsigned __int128) r6) * s3 + ((unsigned __int128) r3) * s6 + ((unsigned __int128) r5) * s4 + ((unsigned __int128) r4) * s5;
	t[3] += ((unsigned __int128) r6) * s4 + ((unsigned __int128) r4) * s6 + ((unsigned __int128) r5) * s5;
	t[4] += ((unsigned __int128) r6) * s5 + ((unsigned __int128) r5) * s6;
	t[5] += ((unsigned __int128) r6) * s6;
	
	//printf("[0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx]\n", (uint64_t)(t[0]>>64), (uint64_t)t[0], (uint64_t)(t[1]>>64), (uint64_t)t[1], (uint64_t)(t[2]>>64), (uint64_t)t[2], (uint64_t)(t[3]>>64), (uint64_t)t[3], (uint64_t)(t[4]>>64), (uint64_t)t[4], (uint64_t)(t[5]>>64), (uint64_t)t[5], (uint64_t)(t[6]>>64), (uint64_t)t[6]);
	
							r0 = (uint64_t)t[0] & 0xfffffffffffffff; c2 = (t[0] >> 60);
	t[1] += c2;	r1 = (uint64_t)t[1] & 0x7ffffffffffffff; c2 = (t[1] >> 59);
	t[2] += c2;	r2 = (uint64_t)t[2] & 0x7ffffffffffffff; c2 = (t[2] >> 59);
	t[3] += c2;	r3 = (uint64_t)t[3] & 0x7ffffffffffffff; c2 = (t[3] >> 59);
	t[4] += c2;	r4 = (uint64_t)t[4] & 0x7ffffffffffffff; c2 = (t[4] >> 59);
	t[5] += c2;	r5 = (uint64_t)t[5] & 0x7ffffffffffffff; c = (uint64_t)(t[5] >> 59);
	t[6] += c;	r6 = (uint64_t)t[6] & 0x7ffffffffffffff; c2 =(t[6] >> 59)*17+r0;
	r0 = (uint64_t)c2 & 0xfffffffffffffff;	c = (c2>>60) + r1;
	r1 = (uint64_t)c & 0x7ffffffffffffff;	c = c>>59;
	r2 += c;
	
	output[0] = r0;
	output[1] = r1;
	output[2] = r2;
	output[3] = r3;
	output[4] = r4;
	output[5] = r5;
	output[6] = r6;
}

void multModEd448(uint64_t* output, uint64_t * in, uint64_t * in2)
{
	/* Adapted from multMod25519 for Ed448-Goldilocks
	 */
	unsigned __int128 t[8], t2[7];
	uint64_t r0,r1,r2,r3,r4,r5,r6,r7,s0,s1,s2,s3,s4,s5,s6,s7,c,c2;
	
	r0 = in[0];
	r1 = in[1];
	r2 = in[2];
	r3 = in[3];
	r4 = in[4];
	r5 = in[5];
	r6 = in[6];
	r7 = in[7];
	
	s0 = in2[0];
	s1 = in2[1];
	s2 = in2[2];
	s3 = in2[3];
	s4 = in2[4];
	s5 = in2[5];
	s6 = in2[6];
	s7 = in2[7];
	
	
	
	t2[0] = ((unsigned __int128) r7) * s1 + ((unsigned __int128) r1) * s7 + ((unsigned __int128) r2) * s6 + ((unsigned __int128) r6) * s2 + ((unsigned __int128) r5) * s3 + ((unsigned __int128) r3) * s5 + ((unsigned __int128) r4) * s4;
	t2[1] = ((unsigned __int128) r7) * s2 + ((unsigned __int128) r2) * s7 + ((unsigned __int128) r6) * s3 + ((unsigned __int128) r3) * s6 + ((unsigned __int128) r5) * s4 + ((unsigned __int128) r4) * s5;
	t2[2] = ((unsigned __int128) r7) * s3 + ((unsigned __int128) r3) * s7 + ((unsigned __int128) r6) * s4 + ((unsigned __int128) r4) * s6 + ((unsigned __int128) r5) * s5;
	t2[3] = ((unsigned __int128) r7) * s4 + ((unsigned __int128) r4) * s7 + ((unsigned __int128) r6) * s5 + ((unsigned __int128) r5) * s6;
	t2[4] = ((unsigned __int128) r7) * s5 + ((unsigned __int128) r5) * s7 + ((unsigned __int128) r6) * s6;
	t2[0] += t2[4];
	t2[5] = ((unsigned __int128) r7) * s6 + ((unsigned __int128) r6) * s7;
	t2[1] += t2[5];
	t2[6] = ((unsigned __int128) r7) * s7;
	t2[2] += t2[6];
	
	t[0] = ((unsigned __int128) r0) * s0 + t2[0];
	t[1] = ((unsigned __int128) r0) * s1 + ((unsigned __int128) r1) * s0 + t2[1];
	t[2] = ((unsigned __int128) r0) * s2 + ((unsigned __int128) r2) * s0 + ((unsigned __int128) r1) * s1 + t2[2];
	t[3] = ((unsigned __int128) r0) * s3 + ((unsigned __int128) r3) * s0 + ((unsigned __int128) r1) * s2 + ((unsigned __int128) r2) * s1 + t2[3];
	t[4] = ((unsigned __int128) r0) * s4 + ((unsigned __int128) r4) * s0 + ((unsigned __int128) r3) * s1 + ((unsigned __int128) r1) * s3 + ((unsigned __int128) r2) * s2 + t2[0] + t2[4];
	t[5] = ((unsigned __int128) r0) * s5 + ((unsigned __int128) r1) * s4 + ((unsigned __int128) r2) * s3 + ((unsigned __int128) r3) * s2 + ((unsigned __int128) r4) * s1 + ((unsigned __int128) r5) * s0 + t2[1] + t2[5];
	t[6] = ((unsigned __int128) r0) * s6 + ((unsigned __int128) r1) * s5 + ((unsigned __int128) r2) * s4 + ((unsigned __int128) r3) * s3 + ((unsigned __int128) r4) * s2 + ((unsigned __int128) r5) * s1 + ((unsigned __int128) r6) * s0 + t2[2] + t2[6];
	t[7] = ((unsigned __int128) r0) * s7 + ((unsigned __int128) r1) * s6 + ((unsigned __int128) r2) * s5 + ((unsigned __int128) r3) * s4 + ((unsigned __int128) r4) * s3 + ((unsigned __int128) r5) * s2 + ((unsigned __int128) r6) * s1 + ((unsigned __int128) r7) * s0 + t2[3];
	
	//printf("[0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx]\n", (uint64_t)(t[0]>>64), (uint64_t)t[0], (uint64_t)(t[1]>>64), (uint64_t)t[1], (uint64_t)(t[2]>>64), (uint64_t)t[2], (uint64_t)(t[3]>>64), (uint64_t)t[3], (uint64_t)(t[4]>>64), (uint64_t)t[4], (uint64_t)(t[5]>>64), (uint64_t)t[5], (uint64_t)(t[6]>>64), (uint64_t)t[6]);
	
							r0 = (uint64_t)t[0] & 0xffffffffffffff; c = (uint64_t)(t[0] >> 56);
	t[1] += c;	r1 = (uint64_t)t[1] & 0xffffffffffffff; c = (uint64_t)(t[1] >> 56);
	t[2] += c;	r2 = (uint64_t)t[2] & 0xffffffffffffff; c = (uint64_t)(t[2] >> 56);
	t[3] += c;	r3 = (uint64_t)t[3] & 0xffffffffffffff; c = (uint64_t)(t[3] >> 56);
	t[4] += c;	r4 = (uint64_t)t[4] & 0xffffffffffffff; c = (uint64_t)(t[4] >> 56);
	t[5] += c;	r5 = (uint64_t)t[5] & 0xffffffffffffff; c = (uint64_t)(t[5] >> 56);
	t[6] += c;	r6 = (uint64_t)t[6] & 0xffffffffffffff; c = (uint64_t)(t[6] >> 56);
	t[7] += c;	r7 = (uint64_t)t[7] & 0xffffffffffffff; c = (uint64_t)(t[7] >> 56);
	r0 += c;	r4 += c; c = (r0>>56); r0 &= 0xffffffffffffff; c2 = (r4>>56); r4 &= 0xffffffffffffff;
	r1 += c;	r5 += c2; 
	
	output[0] = r0;
	output[1] = r1;
	output[2] = r2;
	output[3] = r3;
	output[4] = r4;
	output[5] = r5;
	output[6] = r6;
	output[7] = r7;
}

void multModM511(uint64_t* output, uint64_t * in, uint64_t * in2)
{
	/* Adapted from multMod25519 for M-511
	 */
	unsigned __int128 t[9];
	uint64_t rr0,rr1,rr2,rr3,rr4,rr5,rr6,rr7,rr8,s1,s2,s3,s4,s5,s6,s7,s8,c;
	unsigned __int128 c2;
	
	rr0 = in[0]*17;
	rr1 = in[1]*17;
	rr2 = in[2]*17;
	rr3 = in[3]*17;
	rr4 = in[4]*17;
	rr5 = in[5]*17;
	rr6 = in[6]*17;
	rr7 = in[7]*17;
	rr8 = in[8]*17;
	
	s1 = in2[1]*44;
	s2 = in2[2]*44;
	s3 = in2[3]*44;
	s4 = in2[4]*44;
	s5 = in2[5]*44;
	s6 = in2[6]*44;
	s7 = in2[7]*44;
	s8 = in2[8]*44;
	
	t[0] = ((unsigned __int128) rr8) * s1 + ((unsigned __int128) rr7) * s2 + ((unsigned __int128) rr6) * s3 + ((unsigned __int128) rr5) * s4 + ((unsigned __int128) rr4) * s5 + ((unsigned __int128) rr3) * s6 + ((unsigned __int128) rr2) * s7 + ((unsigned __int128) rr1) * s8;
	t[1] = ((unsigned __int128) rr8) * s2 + ((unsigned __int128) rr7) * s3 + ((unsigned __int128) rr6) * s4 + ((unsigned __int128) rr5) * s5 + ((unsigned __int128) rr4) * s6 + ((unsigned __int128) rr3) * s7 + ((unsigned __int128) rr2) * s8;
	t[2] = ((unsigned __int128) rr8) * s3 + ((unsigned __int128) rr7) * s4 + ((unsigned __int128) rr6) * s5 + ((unsigned __int128) rr5) * s6 + ((unsigned __int128) rr4) * s7 + ((unsigned __int128) rr3) * s8;
	t[3] = ((unsigned __int128) rr8) * s4 + ((unsigned __int128) rr7) * s5 + ((unsigned __int128) rr6) * s6 + ((unsigned __int128) rr5) * s7 + ((unsigned __int128) rr4) * s8;
	t[4] = ((unsigned __int128) rr8) * s5 + ((unsigned __int128) rr7) * s6 + ((unsigned __int128) rr6) * s7 + ((unsigned __int128) rr5) * s8;
	t[5] = ((unsigned __int128) rr8) * s6 + ((unsigned __int128) rr7) * s7 + ((unsigned __int128) rr6) * s8;
	t[6] = ((unsigned __int128) rr8) * s7 + ((unsigned __int128) rr7) * s8;
	t[7] = ((unsigned __int128) rr8) * s8;
	
	t[0] += ((unsigned __int128) in[0]) * in2[0];
	t[1] += ((unsigned __int128) in[0]) * in2[1] + ((unsigned __int128) in[1]) * in2[0];
	t[2] += ((unsigned __int128) in[0]) * in2[2] + ((unsigned __int128) in[2]) * in2[0] + ((unsigned __int128) in[1]) * in2[1];
	t[3] += ((unsigned __int128) in[0]) * in2[3] + ((unsigned __int128) in[3]) * in2[0] + ((unsigned __int128) in[1]) * in2[2] + ((unsigned __int128) in[2]) * in2[1];
	t[4] += ((unsigned __int128) in[0]) * in2[4] + ((unsigned __int128) in[4]) * in2[0] + ((unsigned __int128) in[3]) * in2[1] + ((unsigned __int128) in[1]) * in2[3] + ((unsigned __int128) in[2]) * in2[2];
	t[5] += ((unsigned __int128) in[0]) * in2[5] + ((unsigned __int128) in[1]) * in2[4] + ((unsigned __int128) in[2]) * in2[3] + ((unsigned __int128) in[3]) * in2[2] + ((unsigned __int128) in[4]) * in2[1] + ((unsigned __int128) in[5]) * in2[0];
	t[6] += ((unsigned __int128) in[0]) * in2[6] + ((unsigned __int128) in[1]) * in2[5] + ((unsigned __int128) in[2]) * in2[4] + ((unsigned __int128) in[3]) * in2[3] + ((unsigned __int128) in[4]) * in2[2] + ((unsigned __int128) in[5]) * in2[1] + ((unsigned __int128) in[6]) * in2[0];
	t[7] += ((unsigned __int128) in[0]) * in2[7] + ((unsigned __int128) in[1]) * in2[6] + ((unsigned __int128) in[2]) * in2[5] + ((unsigned __int128) in[3]) * in2[4] + ((unsigned __int128) in[4]) * in2[3] + ((unsigned __int128) in[5]) * in2[2] + ((unsigned __int128) in[6]) * in2[1] + ((unsigned __int128) in[7]) * in2[0];
	t[8] = ((unsigned __int128) in[0]) * in2[8] + ((unsigned __int128) in[1]) * in2[7] + ((unsigned __int128) in[2]) * in2[6] + ((unsigned __int128) in[3]) * in2[5] + ((unsigned __int128) in[4]) * in2[4] + ((unsigned __int128) in[5]) * in2[3] + ((unsigned __int128) in[6]) * in2[2] + ((unsigned __int128) in[7]) * in2[1] + ((unsigned __int128) in[8]) * in2[0];
	
	
	//printf("[0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx]\n", (uint64_t)(t[0]>>64), (uint64_t)t[0], (uint64_t)(t[1]>>64), (uint64_t)t[1], (uint64_t)(t[2]>>64), (uint64_t)t[2], (uint64_t)(t[3]>>64), (uint64_t)t[3], (uint64_t)(t[4]>>64), (uint64_t)t[4], (uint64_t)(t[5]>>64), (uint64_t)t[5], (uint64_t)(t[6]>>64), (uint64_t)t[6]);
	
							rr0 = (uint64_t)t[0] & 0x1ffffffffffffff; c2 = (t[0] >> 57);
	t[1] += c2;	rr1 = (uint64_t)t[1] & 0x1ffffffffffffff; c2 = (t[1] >> 57);
	t[2] += c2;	rr2 = (uint64_t)t[2] & 0x1ffffffffffffff; c2 = (t[2] >> 57);
	t[3] += c2;	rr3 = (uint64_t)t[3] & 0x1ffffffffffffff; c2 = (t[3] >> 57);
	t[4] += c2;	rr4 = (uint64_t)t[4] & 0x1ffffffffffffff; c2 = (t[4] >> 57);
	t[5] += c2;	rr5 = (uint64_t)t[5] & 0x1ffffffffffffff; c2 = (t[5] >> 57);
	t[6] += c2;	rr6 = (uint64_t)t[6] & 0x1ffffffffffffff; c2 = (t[6] >> 57);
	t[7] += c2;	rr7 = (uint64_t)t[7] & 0x1ffffffffffffff; c = (uint64_t)(t[7] >> 57);
	t[8] += c;	rr8 = (uint64_t)t[8] & 0x7fffffffffffff; c2 =(t[8] >> 55)*187+rr0;
	rr0 = (uint64_t)c2 & 0x1ffffffffffffff;	c = (c2>>57) + rr1;
	rr1 = (uint64_t)c & 0x1ffffffffffffff;	c = c>>57;
	rr2 += c;
	
	output[0] = rr0;
	output[1] = rr1;
	output[2] = rr2;
	output[3] = rr3;
	output[4] = rr4;
	output[5] = rr5;
	output[6] = rr6;
	output[7] = rr7;
	output[8] = rr8;
}

void multModE521(uint64_t* output, uint64_t * in, uint64_t * in2)
{
	/* Adapted from multMod25519 for E-521
	 */
	unsigned __int128 t[9];
	uint64_t r0,r1,r2,r3,r4,r5,r6,r7,r8,s0,s1,s2,s3,s4,s5,s6,s7,s8,c,rr1,rr2,rr3,rr4,rr5,rr6,rr7,rr8;
	
	r0 = in[0];
	r1 = in[1];
	r2 = in[2];
	r3 = in[3];
	r4 = in[4];
	r5 = in[5];
	r6 = in[6];
	r7 = in[7];
	r8 = in[8];
	
	rr1 = in[1]*2;
	rr2 = in[2]*2;
	rr3 = in[3]*2;
	rr4 = in[4]*2;
	rr5 = in[5]*2;
	rr6 = in[6]*2;
	rr7 = in[7]*2;
	rr8 = in[8]*2;
	
	s0 = in2[0];
	s1 = in2[1];
	s2 = in2[2];
	s3 = in2[3];
	s4 = in2[4];
	s5 = in2[5];
	s6 = in2[6];
	s7 = in2[7];
	s8 = in2[8];
	
	t[0] = ((unsigned __int128) r0) * s0 + ((unsigned __int128) rr8) * s1 + ((unsigned __int128) rr7) * s2 + ((unsigned __int128) rr6) * s3 + ((unsigned __int128) rr5) * s4 + ((unsigned __int128) rr4) * s5 + ((unsigned __int128) rr3) * s6 + ((unsigned __int128) rr2) * s7 + ((unsigned __int128) rr1) * s8;
	t[1] = ((unsigned __int128) r0) * s1 + ((unsigned __int128) r1) * s0 + ((unsigned __int128) rr8) * s2 + ((unsigned __int128) rr7) * s3 + ((unsigned __int128) rr6) * s4 + ((unsigned __int128) rr5) * s5 + ((unsigned __int128) rr4) * s6 + ((unsigned __int128) rr3) * s7 + ((unsigned __int128) rr2) * s8;
	t[2] = ((unsigned __int128) r0) * s2 + ((unsigned __int128) r2) * s0 + ((unsigned __int128) r1) * s1 + ((unsigned __int128) rr8) * s3 + ((unsigned __int128) rr7) * s4 + ((unsigned __int128) rr6) * s5 + ((unsigned __int128) rr5) * s6 + ((unsigned __int128) rr4) * s7 + ((unsigned __int128) rr3) * s8;
	t[3] = ((unsigned __int128) r0) * s3 + ((unsigned __int128) r3) * s0 + ((unsigned __int128) r1) * s2 + ((unsigned __int128) r2) * s1 + ((unsigned __int128) rr8) * s4 + ((unsigned __int128) rr7) * s5 + ((unsigned __int128) rr6) * s6 + ((unsigned __int128) rr5) * s7 + ((unsigned __int128) rr4) * s8;
	t[4] = ((unsigned __int128) r0) * s4 + ((unsigned __int128) r4) * s0 + ((unsigned __int128) r3) * s1 + ((unsigned __int128) r1) * s3 + ((unsigned __int128) r2) * s2 + ((unsigned __int128) rr8) * s5 + ((unsigned __int128) rr7) * s6 + ((unsigned __int128) rr6) * s7 + ((unsigned __int128) rr5) * s8;
	t[5] = ((unsigned __int128) r0) * s5 + ((unsigned __int128) r1) * s4 + ((unsigned __int128) r2) * s3 + ((unsigned __int128) r3) * s2 + ((unsigned __int128) r4) * s1 + ((unsigned __int128) r5) * s0 + ((unsigned __int128) rr8) * s6 + ((unsigned __int128) rr7) * s7 + ((unsigned __int128) rr6) * s8;
	t[6] = ((unsigned __int128) r0) * s6 + ((unsigned __int128) r1) * s5 + ((unsigned __int128) r2) * s4 + ((unsigned __int128) r3) * s3 + ((unsigned __int128) r4) * s2 + ((unsigned __int128) r5) * s1 + ((unsigned __int128) r6) * s0 + ((unsigned __int128) rr8) * s7 + ((unsigned __int128) rr7) * s8;
	t[7] = ((unsigned __int128) r0) * s7 + ((unsigned __int128) r1) * s6 + ((unsigned __int128) r2) * s5 + ((unsigned __int128) r3) * s4 + ((unsigned __int128) r4) * s3 + ((unsigned __int128) r5) * s2 + ((unsigned __int128) r6) * s1 + ((unsigned __int128) r7) * s0 + ((unsigned __int128) rr8) * s8;
	t[8] = ((unsigned __int128) r0) * s8 + ((unsigned __int128) r1) * s7 + ((unsigned __int128) r2) * s6 + ((unsigned __int128) r3) * s5 + ((unsigned __int128) r4) * s4 + ((unsigned __int128) r5) * s3 + ((unsigned __int128) r6) * s2 + ((unsigned __int128) r7) * s1 + ((unsigned __int128) r8) * s0;
	
	//printf("[0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx, 0x%lx%016lx]\n", (uint64_t)(t[0]>>64), (uint64_t)t[0], (uint64_t)(t[1]>>64), (uint64_t)t[1], (uint64_t)(t[2]>>64), (uint64_t)t[2], (uint64_t)(t[3]>>64), (uint64_t)t[3], (uint64_t)(t[4]>>64), (uint64_t)t[4], (uint64_t)(t[5]>>64), (uint64_t)t[5], (uint64_t)(t[6]>>64), (uint64_t)t[6]);
	
							r0 = (uint64_t)t[0] & 0x3ffffffffffffff; c = (uint64_t)(t[0] >> 58);
	t[1] += c;	r1 = (uint64_t)t[1] & 0x3ffffffffffffff; c = (uint64_t)(t[1] >> 58);
	t[2] += c;	r2 = (uint64_t)t[2] & 0x3ffffffffffffff; c = (uint64_t)(t[2] >> 58);
	t[3] += c;	r3 = (uint64_t)t[3] & 0x3ffffffffffffff; c = (uint64_t)(t[3] >> 58);
	t[4] += c;	r4 = (uint64_t)t[4] & 0x3ffffffffffffff; c = (uint64_t)(t[4] >> 58);
	t[5] += c;	r5 = (uint64_t)t[5] & 0x3ffffffffffffff; c = (uint64_t)(t[5] >> 58);
	t[6] += c;	r6 = (uint64_t)t[6] & 0x3ffffffffffffff; c = (uint64_t)(t[6] >> 58);
	t[7] += c;	r7 = (uint64_t)t[7] & 0x3ffffffffffffff; c = (uint64_t)(t[7] >> 58);
	t[8] += c;	r8 = (uint64_t)t[8] & 0x1ffffffffffffff; c = (uint64_t)(t[8] >> 57);
	r0 += c;	c = (r0>>58); r0 &= 0x3ffffffffffffff;
	r1 += c;
	
	output[0] = r0;
	output[1] = r1;
	output[2] = r2;
	output[3] = r3;
	output[4] = r4;
	output[5] = r5;
	output[6] = r6;
	output[7] = r7;
	output[8] = r8;
}
