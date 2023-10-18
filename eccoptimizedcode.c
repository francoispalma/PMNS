#include <stdint.h>

#include "eccoptimizedcode.h"

void multMod25519(uint64_t* output, uint64_t * in, uint64_t * in2) {
	/* Copyright 2008, Google Inc.
	 * All rights reserved.
	 *
	 * Code released into the public domain by Adam Langley <agl@imperialviolet.org>
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
	uint64_t r0,r1,r2,r3,r4,r5,r6,s0,s1,s2,s3,s4,s5,s6,c;
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
	/* Adapted from multMod25519 for Curve41419
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
