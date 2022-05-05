a = [0x1ee7aa9b53fc36df4825109007c98e, 0x3b8a66f0255bd7c849a77f98b98286, -0x24ec71dc0a10f82cd775a39480cb41, 0x90f1b3b82e0ff68ba97b2524cbd7f, 0x187fa971c3c416a5ce5aea856d491, 0x2cf90f813012dad56da8867522691f, -0x35501d90b0b20b001f05a9834c02fb, 0x102e3caf4d3f2430041a7a99a1da70, 0x3cccbf6cb320df6611363e50fa55c0]
b = [0x3eb80b8e806591760d4e7d901c64f2, 0x17791fb1d51d3dc6354ca14d7cd508, -0x2c5914274470d0519127d6afcccc37, -0x3d8031784b57b76e7ccabb40133064, -0x3f9af08f2bb796621d8ae90dd40e0b, 0x36d325944f54c8c1dbd52d3a24ad6c, -0x8944c793da129366e53346e56c5d5, -0x2553e15f3d65cd7577d1847d9f2b33, -0x1bc3370eb82d8d8e434997a3226436]
c = [-0x37029ab909ae89129b3c81d2252d0, -0x10cd975f7155c26c28c32b3c95392, -0x13a7856540b1581e5309754b2aaec, 0x7a75b80e3425a635269da8bfc241, 0x1b4b7bcebe8bd5e2a06d1141774e5, 0x84e85502012a8bade07f4cd4846e, -0xbefbdc4a2f5627c5b8ccb36a9259, -0xd0f49ce41e0fcb1db2e3e3c946e4, -0x1d98072a7d1169e6e996c3dd7cff6]
cc = [-0x37029ab909ae7b4ab1eb48d92a178, -0x10cd975f7155b6306ac41b678a6aa, -0x13a7856540b15052f34a6a05d5c34, 0x7a75b80e34261300a46d92f20b97, 0x1b4b7bcebe8be48cf021750daffaf, 0x84e855020129a7fced4075a8f1be, -0xbefbdc4a2f569d944f3eb50691ba, -0xd0f49ce41e0f913b86657f70feaf, -0x1d98072a7d1164d09caa9b54f7490]

LOW = lambda x: x % 2**64
HIGH = lambda x: x >> 64
HI = lambda x: (x >> 64) % 2**64

def multadd128k(R, A, B):
	A0B0 = LOW(A) * LOW(B)
	A1B1 = HIGH(A) * HIGH(B)
	A1B0_A0B1 = A1B1 + A0B0 - ((LOW(A) - HIGH(A)) * (LOW(B) - HIGH(B)))
	aux3 = HI(A0B0) + LOW(A1B0_A0B1)
	aux2 = (HIGH(aux3) + HIGH(A1B0_A0B1) + LOW(A1B1)) % 2**128
	aux1 = HIGH(A1B1)
	tmplo = (R % 2**128)
	Rlo = ((R % 2**128) + LOW(A0B0) + (aux3 << 64)) % 2**128
	Rhi = ((R >> 128) + aux2 + (aux1 << 64) + (tmplo > Rlo)) % 2**128
	
	if (R + A * B) % 2**256 != Rlo + (Rhi << 128):
		print(hex(A * B % 2**256))
		print(hex(Rlo + (Rhi << 128)))

for i in range(len(a)):
	multadd128k(0, a[i], b[i])
	


/*		a->lo[0] = 0xdf4825109007c98e;*/
/*		a->hi[0] = 0x1ee7aa9b53fc36;*/
/*		b->lo[0] = 0x760d4e7d901c64f2;*/
/*		b->hi[0] = 0x3eb80b8e806591;*/
/*		a->lo[1] = 0xc849a77f98b98286;*/
/*		a->hi[1] = 0x3b8a66f0255bd7;*/
/*		b->lo[1] = 0xc6354ca14d7cd508;*/
/*		b->hi[1] = 0x17791fb1d51d3d;*/
/*		a->lo[2] = 0xd3288a5c6b7f34bf;*/
/*		a->hi[2] = 0xffdb138e23f5ef07;*/
/*		b->lo[2] = 0xae6ed829503333c9;*/
/*		b->hi[2] = 0xffd3a6ebd8bb8f2f;*/
/*		a->lo[3] = 0x68ba97b2524cbd7f;*/
/*		a->hi[3] = 0x90f1b3b82e0ff;*/
/*		b->lo[3] = 0x91833544bfeccf9c;*/
/*		b->hi[3] = 0xffc27fce87b4a848;*/
/*		a->lo[4] = 0x6a5ce5aea856d491;*/
/*		a->hi[4] = 0x187fa971c3c41;*/
/*		b->lo[4] = 0x9de27516f22bf1f5;*/
/*		b->hi[4] = 0xffc0650f70d44869;*/
/*		a->lo[5] = 0xd56da8867522691f;*/
/*		a->hi[5] = 0x2cf90f813012da;*/
/*		b->lo[5] = 0xc1dbd52d3a24ad6c;*/
/*		b->hi[5] = 0x36d325944f54c8;*/
/*		a->lo[6] = 0xffe0fa567cb3fd05;*/
/*		a->hi[6] = 0xffcaafe26f4f4df4;*/
/*		b->lo[6] = 0xc991accb91a93a2b;*/
/*		b->hi[6] = 0xfff76bb386c25ed6;*/
/*		a->lo[7] = 0x30041a7a99a1da70;*/
/*		a->hi[7] = 0x102e3caf4d3f24;*/
/*		b->lo[7] = 0x8a882e7b8260d4cd;*/
/*		b->hi[7] = 0xffdaac1ea0c29a32;*/
/*		a->lo[8] = 0x6611363e50fa55c0;*/
/*		a->hi[8] = 0x3cccbf6cb320df;*/
/*		b->lo[8] = 0x71bcb6685cdd9bca;*/
/*		b->hi[8] = 0xffe43cc8f147d272;*/
