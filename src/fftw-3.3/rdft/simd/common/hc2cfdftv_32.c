/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Thu Dec 10 07:06:57 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2cdft_c.native -fma -simd -compact -variables 4 -pipeline-latency 8 -trivial-stores -variables 32 -no-generate-bytw -n 32 -dit -name hc2cfdftv_32 -include rdft/simd/hc2cfv.h */

/*
 * This function contains 249 FP additions, 224 FP multiplications,
 * (or, 119 additions, 94 multiplications, 130 fused multiply/add),
 * 154 stack variables, 8 constants, and 64 memory accesses
 */
#include "rdft/simd/hc2cfv.h"

static void hc2cfdftv_32(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DVK(KP668178637, +0.668178637919298919997757686523080761552472251);
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DVK(KP198912367, +0.198912367379658006911597622644676228597850501);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DVK(KP414213562, +0.414213562373095048801688724209698078569671875);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * ((TWVL / VL) * 62)); m < me; m = m + VL, Rp = Rp + (VL * ms), Ip = Ip + (VL * ms), Rm = Rm - (VL * ms), Im = Im - (VL * ms), W = W + (TWVL * 62), MAKE_VOLATILE_STRIDE(128, rs)) {
	       V T47, T48, T4l, T3w, T3F, T3B, T41, Ts, T2y, T1Q, T2B, T27, T2J, T3a, T40;
	       V T1X, T2C, T43, T44, T4a, T4b, T4m, T3p, T3E, T15, T2K, T1u, T2F, T3h, T3C;
	       V T1n, T2E, T2a, T2z, T1a, T18, TU, T3m, T3f, T1r, T1p, T13, T3n, T3e, TB;
	       V T3k, T1l, T3c, TK, T3j, T1g, T3b, T3l, T3o, TL, T14, T1s, T1t, T3d, T3g;
	       V T1b, T1m, T28, T29, T3Q, T3W, T3T, T3X, T3O, T3P, T3R, T3S, T3U, T3Z, T3V;
	       V T3Y;
	       {
		    V T1U, T1S, T3, T3u, T7, T1z, T1D, T3t, T24, T22, Tc, Tg, Th, T3q, T1J;
		    V Tl, Tp, Tq, T3r, T1O, T3s, T3v, T3z, T3A, T8, Tr, T1E, T1P, T25, T26;
		    V T38, T39, T1V, T1W;
		    {
			 V T1, T2, T5, T6, T1T, T1R, T4, T1x, T1y, T1B, T1C, T1w, T1A, T23, T21;
			 V T1I, T1G, Ta, Tb, T9, T1H, Te, Tf, Td, T1F, T1N, T1L, Tj, Tk, Ti;
			 V T1M, Tn, To, Tm, T1K;
			 T1 = LD(&(Rp[0]), ms, &(Rp[0]));
			 T2 = LD(&(Rm[0]), -ms, &(Rm[0]));
			 T1T = LDW(&(W[0]));
			 T1U = VZMULIJ(T1T, VFNMSCONJ(T2, T1));
			 T5 = LD(&(Rp[WS(rs, 8)]), ms, &(Rp[0]));
			 T6 = LD(&(Rm[WS(rs, 8)]), -ms, &(Rm[0]));
			 T1R = LDW(&(W[TWVL * 32]));
			 T1S = VZMULIJ(T1R, VFNMSCONJ(T6, T5));
			 T3 = VFMACONJ(T2, T1);
			 T3u = VADD(T1U, T1S);
			 T4 = LDW(&(W[TWVL * 30]));
			 T7 = VZMULJ(T4, VFMACONJ(T6, T5));
			 T1x = LD(&(Rp[WS(rs, 12)]), ms, &(Rp[0]));
			 T1y = LD(&(Rm[WS(rs, 12)]), -ms, &(Rm[0]));
			 T1w = LDW(&(W[TWVL * 48]));
			 T1z = VZMULIJ(T1w, VFNMSCONJ(T1y, T1x));
			 T1B = LD(&(Rp[WS(rs, 4)]), ms, &(Rp[0]));
			 T1C = LD(&(Rm[WS(rs, 4)]), -ms, &(Rm[0]));
			 T1A = LDW(&(W[TWVL * 16]));
			 T1D = VZMULIJ(T1A, VFNMSCONJ(T1C, T1B));
			 T3t = VADD(T1D, T1z);
			 T23 = LDW(&(W[TWVL * 46]));
			 T24 = VZMULJ(T23, VFMACONJ(T1y, T1x));
			 T21 = LDW(&(W[TWVL * 14]));
			 T22 = VZMULJ(T21, VFMACONJ(T1C, T1B));
			 Ta = LD(&(Rp[WS(rs, 2)]), ms, &(Rp[0]));
			 Tb = LD(&(Rm[WS(rs, 2)]), -ms, &(Rm[0]));
			 T9 = LDW(&(W[TWVL * 6]));
			 Tc = VZMULJ(T9, VFMACONJ(Tb, Ta));
			 T1H = LDW(&(W[TWVL * 8]));
			 T1I = VZMULIJ(T1H, VFNMSCONJ(Tb, Ta));
			 Te = LD(&(Rp[WS(rs, 10)]), ms, &(Rp[0]));
			 Tf = LD(&(Rm[WS(rs, 10)]), -ms, &(Rm[0]));
			 Td = LDW(&(W[TWVL * 38]));
			 Tg = VZMULJ(Td, VFMACONJ(Tf, Te));
			 T1F = LDW(&(W[TWVL * 40]));
			 T1G = VZMULIJ(T1F, VFNMSCONJ(Tf, Te));
			 Th = VSUB(Tc, Tg);
			 T3q = VADD(T1I, T1G);
			 T1J = VSUB(T1G, T1I);
			 Tj = LD(&(Rp[WS(rs, 14)]), ms, &(Rp[0]));
			 Tk = LD(&(Rm[WS(rs, 14)]), -ms, &(Rm[0]));
			 Ti = LDW(&(W[TWVL * 54]));
			 Tl = VZMULJ(Ti, VFMACONJ(Tk, Tj));
			 T1M = LDW(&(W[TWVL * 56]));
			 T1N = VZMULIJ(T1M, VFNMSCONJ(Tk, Tj));
			 Tn = LD(&(Rp[WS(rs, 6)]), ms, &(Rp[0]));
			 To = LD(&(Rm[WS(rs, 6)]), -ms, &(Rm[0]));
			 Tm = LDW(&(W[TWVL * 22]));
			 Tp = VZMULJ(Tm, VFMACONJ(To, Tn));
			 T1K = LDW(&(W[TWVL * 24]));
			 T1L = VZMULIJ(T1K, VFNMSCONJ(To, Tn));
			 Tq = VSUB(Tl, Tp);
			 T3r = VADD(T1N, T1L);
			 T1O = VSUB(T1L, T1N);
		    }
		    T47 = VADD(T3u, T3t);
		    T48 = VADD(T3q, T3r);
		    T4l = VSUB(T48, T47);
		    T3s = VSUB(T3q, T3r);
		    T3v = VSUB(T3t, T3u);
		    T3w = VFNMS(LDK(KP414213562), T3v, T3s);
		    T3F = VFMA(LDK(KP414213562), T3s, T3v);
		    T3z = VADD(Tl, Tp);
		    T3A = VADD(Tc, Tg);
		    T3B = VSUB(T3z, T3A);
		    T41 = VADD(T3A, T3z);
		    T8 = VSUB(T3, T7);
		    Tr = VADD(Th, Tq);
		    Ts = VFNMS(LDK(KP707106781), Tr, T8);
		    T2y = VFMA(LDK(KP707106781), Tr, T8);
		    T1E = VSUB(T1z, T1D);
		    T1P = VSUB(T1J, T1O);
		    T1Q = VFNMS(LDK(KP707106781), T1P, T1E);
		    T2B = VFMA(LDK(KP707106781), T1P, T1E);
		    T25 = VSUB(T22, T24);
		    T26 = VSUB(Tq, Th);
		    T27 = VFMA(LDK(KP707106781), T26, T25);
		    T2J = VFNMS(LDK(KP707106781), T26, T25);
		    T38 = VADD(T3, T7);
		    T39 = VADD(T22, T24);
		    T3a = VSUB(T38, T39);
		    T40 = VADD(T38, T39);
		    T1V = VSUB(T1S, T1U);
		    T1W = VADD(T1J, T1O);
		    T1X = VFNMS(LDK(KP707106781), T1W, T1V);
		    T2C = VFMA(LDK(KP707106781), T1W, T1V);
	       }
	       {
		    V TP, TT, TN, TO, TM, T19, TR, TS, TQ, T17, TY, T12, TW, TX, TV;
		    V T1q, T10, T11, TZ, T1o, Tw, T1i, TA, T1k, Tu, Tv, Tt, T1h, Ty, Tz;
		    V Tx, T1j, TF, T1f, TJ, T1d, TD, TE, TC, T1e, TH, TI, TG, T1c;
		    TN = LD(&(Rp[WS(rs, 3)]), ms, &(Rp[WS(rs, 1)]));
		    TO = LD(&(Rm[WS(rs, 3)]), -ms, &(Rm[WS(rs, 1)]));
		    TM = LDW(&(W[TWVL * 10]));
		    TP = VZMULJ(TM, VFMACONJ(TO, TN));
		    T19 = LDW(&(W[TWVL * 12]));
		    T1a = VZMULIJ(T19, VFNMSCONJ(TO, TN));
		    TR = LD(&(Rp[WS(rs, 11)]), ms, &(Rp[WS(rs, 1)]));
		    TS = LD(&(Rm[WS(rs, 11)]), -ms, &(Rm[WS(rs, 1)]));
		    TQ = LDW(&(W[TWVL * 42]));
		    TT = VZMULJ(TQ, VFMACONJ(TS, TR));
		    T17 = LDW(&(W[TWVL * 44]));
		    T18 = VZMULIJ(T17, VFNMSCONJ(TS, TR));
		    TU = VSUB(TP, TT);
		    T3m = VADD(T1a, T18);
		    T3f = VADD(TP, TT);
		    TW = LD(&(Rp[WS(rs, 15)]), ms, &(Rp[WS(rs, 1)]));
		    TX = LD(&(Rm[WS(rs, 15)]), -ms, &(Rm[WS(rs, 1)]));
		    TV = LDW(&(W[TWVL * 58]));
		    TY = VZMULJ(TV, VFMACONJ(TX, TW));
		    T1q = LDW(&(W[TWVL * 60]));
		    T1r = VZMULIJ(T1q, VFNMSCONJ(TX, TW));
		    T10 = LD(&(Rp[WS(rs, 7)]), ms, &(Rp[WS(rs, 1)]));
		    T11 = LD(&(Rm[WS(rs, 7)]), -ms, &(Rm[WS(rs, 1)]));
		    TZ = LDW(&(W[TWVL * 26]));
		    T12 = VZMULJ(TZ, VFMACONJ(T11, T10));
		    T1o = LDW(&(W[TWVL * 28]));
		    T1p = VZMULIJ(T1o, VFNMSCONJ(T11, T10));
		    T13 = VSUB(TY, T12);
		    T3n = VADD(T1r, T1p);
		    T3e = VADD(TY, T12);
		    Tu = LD(&(Rp[WS(rs, 5)]), ms, &(Rp[WS(rs, 1)]));
		    Tv = LD(&(Rm[WS(rs, 5)]), -ms, &(Rm[WS(rs, 1)]));
		    Tt = LDW(&(W[TWVL * 18]));
		    Tw = VZMULJ(Tt, VFMACONJ(Tv, Tu));
		    T1h = LDW(&(W[TWVL * 20]));
		    T1i = VZMULIJ(T1h, VFNMSCONJ(Tv, Tu));
		    Ty = LD(&(Rp[WS(rs, 13)]), ms, &(Rp[WS(rs, 1)]));
		    Tz = LD(&(Rm[WS(rs, 13)]), -ms, &(Rm[WS(rs, 1)]));
		    Tx = LDW(&(W[TWVL * 50]));
		    TA = VZMULJ(Tx, VFMACONJ(Tz, Ty));
		    T1j = LDW(&(W[TWVL * 52]));
		    T1k = VZMULIJ(T1j, VFNMSCONJ(Tz, Ty));
		    TB = VSUB(Tw, TA);
		    T3k = VADD(T1i, T1k);
		    T1l = VSUB(T1i, T1k);
		    T3c = VADD(Tw, TA);
		    TD = LD(&(Rp[WS(rs, 1)]), ms, &(Rp[WS(rs, 1)]));
		    TE = LD(&(Rm[WS(rs, 1)]), -ms, &(Rm[WS(rs, 1)]));
		    TC = LDW(&(W[TWVL * 2]));
		    TF = VZMULJ(TC, VFMACONJ(TE, TD));
		    T1e = LDW(&(W[TWVL * 4]));
		    T1f = VZMULIJ(T1e, VFNMSCONJ(TE, TD));
		    TH = LD(&(Rp[WS(rs, 9)]), ms, &(Rp[WS(rs, 1)]));
		    TI = LD(&(Rm[WS(rs, 9)]), -ms, &(Rm[WS(rs, 1)]));
		    TG = LDW(&(W[TWVL * 34]));
		    TJ = VZMULJ(TG, VFMACONJ(TI, TH));
		    T1c = LDW(&(W[TWVL * 36]));
		    T1d = VZMULIJ(T1c, VFNMSCONJ(TI, TH));
		    TK = VSUB(TF, TJ);
		    T3j = VADD(T1f, T1d);
		    T1g = VSUB(T1d, T1f);
		    T3b = VADD(TF, TJ);
	       }
	       T43 = VADD(T3b, T3c);
	       T44 = VADD(T3e, T3f);
	       T4a = VADD(T3j, T3k);
	       T4b = VADD(T3n, T3m);
	       T4m = VSUB(T4a, T4b);
	       T3l = VSUB(T3j, T3k);
	       T3o = VSUB(T3m, T3n);
	       T3p = VFMA(LDK(KP414213562), T3o, T3l);
	       T3E = VFNMS(LDK(KP414213562), T3l, T3o);
	       TL = VFMA(LDK(KP414213562), TK, TB);
	       T14 = VFNMS(LDK(KP414213562), T13, TU);
	       T15 = VSUB(TL, T14);
	       T2K = VADD(TL, T14);
	       T1s = VSUB(T1p, T1r);
	       T1t = VADD(T1g, T1l);
	       T1u = VFNMS(LDK(KP707106781), T1t, T1s);
	       T2F = VFMA(LDK(KP707106781), T1t, T1s);
	       T3d = VSUB(T3b, T3c);
	       T3g = VSUB(T3e, T3f);
	       T3h = VADD(T3d, T3g);
	       T3C = VSUB(T3g, T3d);
	       T1b = VSUB(T18, T1a);
	       T1m = VSUB(T1g, T1l);
	       T1n = VFNMS(LDK(KP707106781), T1m, T1b);
	       T2E = VFMA(LDK(KP707106781), T1m, T1b);
	       T28 = VFMA(LDK(KP414213562), TU, T13);
	       T29 = VFNMS(LDK(KP414213562), TB, TK);
	       T2a = VSUB(T28, T29);
	       T2z = VADD(T29, T28);
	       {
		    V T4o, T4u, T4r, T4v, T4k, T4n, T4p, T4q, T4s, T4x, T4t, T4w, T3y, T3K, T3H;
		    V T3L, T3i, T3x, T3D, T3G, T3I, T3N, T3J, T3M, T46, T4g, T4d, T4h, T42, T45;
		    V T49, T4c, T4e, T4j, T4f, T4i;
		    T4k = VSUB(T40, T41);
		    T4n = VADD(T4l, T4m);
		    T4o = VFMA(LDK(KP707106781), T4n, T4k);
		    T4u = VFNMS(LDK(KP707106781), T4n, T4k);
		    T4p = VSUB(T44, T43);
		    T4q = VSUB(T4m, T4l);
		    T4r = VFMA(LDK(KP707106781), T4q, T4p);
		    T4v = VFNMS(LDK(KP707106781), T4q, T4p);
		    T4s = VCONJ(VMUL(LDK(KP500000000), VFNMSI(T4r, T4o)));
		    ST(&(Rm[WS(rs, 3)]), T4s, -ms, &(Rm[WS(rs, 1)]));
		    T4x = VCONJ(VMUL(LDK(KP500000000), VFMAI(T4v, T4u)));
		    ST(&(Rm[WS(rs, 11)]), T4x, -ms, &(Rm[WS(rs, 1)]));
		    T4t = VMUL(LDK(KP500000000), VFMAI(T4r, T4o));
		    ST(&(Rp[WS(rs, 4)]), T4t, ms, &(Rp[0]));
		    T4w = VMUL(LDK(KP500000000), VFNMSI(T4v, T4u));
		    ST(&(Rp[WS(rs, 12)]), T4w, ms, &(Rp[0]));
		    T3i = VFNMS(LDK(KP707106781), T3h, T3a);
		    T3x = VSUB(T3p, T3w);
		    T3y = VFMA(LDK(KP923879532), T3x, T3i);
		    T3K = VFNMS(LDK(KP923879532), T3x, T3i);
		    T3D = VFNMS(LDK(KP707106781), T3C, T3B);
		    T3G = VSUB(T3E, T3F);
		    T3H = VFNMS(LDK(KP923879532), T3G, T3D);
		    T3L = VFMA(LDK(KP923879532), T3G, T3D);
		    T3I = VMUL(LDK(KP500000000), VFNMSI(T3H, T3y));
		    ST(&(Rp[WS(rs, 6)]), T3I, ms, &(Rp[0]));
		    T3N = VMUL(LDK(KP500000000), VFMAI(T3L, T3K));
		    ST(&(Rp[WS(rs, 10)]), T3N, ms, &(Rp[0]));
		    T3J = VCONJ(VMUL(LDK(KP500000000), VFMAI(T3H, T3y)));
		    ST(&(Rm[WS(rs, 5)]), T3J, -ms, &(Rm[WS(rs, 1)]));
		    T3M = VCONJ(VMUL(LDK(KP500000000), VFNMSI(T3L, T3K)));
		    ST(&(Rm[WS(rs, 9)]), T3M, -ms, &(Rm[WS(rs, 1)]));
		    T42 = VADD(T40, T41);
		    T45 = VADD(T43, T44);
		    T46 = VSUB(T42, T45);
		    T4g = VADD(T42, T45);
		    T49 = VADD(T47, T48);
		    T4c = VADD(T4a, T4b);
		    T4d = VSUB(T49, T4c);
		    T4h = VADD(T49, T4c);
		    T4e = VMUL(LDK(KP500000000), VFMAI(T4d, T46));
		    ST(&(Rp[WS(rs, 8)]), T4e, ms, &(Rp[0]));
		    T4j = VCONJ(VMUL(LDK(KP500000000), VADD(T4h, T4g)));
		    ST(&(Rm[WS(rs, 15)]), T4j, -ms, &(Rm[WS(rs, 1)]));
		    T4f = VCONJ(VMUL(LDK(KP500000000), VFNMSI(T4d, T46)));
		    ST(&(Rm[WS(rs, 7)]), T4f, -ms, &(Rm[WS(rs, 1)]));
		    T4i = VMUL(LDK(KP500000000), VSUB(T4g, T4h));
		    ST(&(Rp[0]), T4i, ms, &(Rp[0]));
	       }
	       T3O = VFMA(LDK(KP707106781), T3h, T3a);
	       T3P = VADD(T3F, T3E);
	       T3Q = VFMA(LDK(KP923879532), T3P, T3O);
	       T3W = VFNMS(LDK(KP923879532), T3P, T3O);
	       T3R = VFMA(LDK(KP707106781), T3C, T3B);
	       T3S = VADD(T3w, T3p);
	       T3T = VFMA(LDK(KP923879532), T3S, T3R);
	       T3X = VFNMS(LDK(KP923879532), T3S, T3R);
	       T3U = VCONJ(VMUL(LDK(KP500000000), VFNMSI(T3T, T3Q)));
	       ST(&(Rm[WS(rs, 1)]), T3U, -ms, &(Rm[WS(rs, 1)]));
	       T3Z = VCONJ(VMUL(LDK(KP500000000), VFMAI(T3X, T3W)));
	       ST(&(Rm[WS(rs, 13)]), T3Z, -ms, &(Rm[WS(rs, 1)]));
	       T3V = VMUL(LDK(KP500000000), VFMAI(T3T, T3Q));
	       ST(&(Rp[WS(rs, 2)]), T3V, ms, &(Rp[0]));
	       T3Y = VMUL(LDK(KP500000000), VFNMSI(T3X, T3W));
	       ST(&(Rp[WS(rs, 14)]), T3Y, ms, &(Rp[0]));
	       {
		    V T2I, T35, T2S, T31, T2P, T34, T2T, T2Y, T2A, T2Z, T2H, T30, T2D, T2G, T2L;
		    V T2W, T2O, T2X, T2M, T2N, T2Q, T36, T37, T2R, T2U, T32, T33, T2V, T20, T2v;
		    V T2i, T2r, T2f, T2u, T2j, T2o, T16, T2p, T1Z, T2q, T1v, T1Y, T2b, T2m, T2e;
		    V T2n, T2c, T2d, T2g, T2w, T2x, T2h, T2k, T2s, T2t, T2l;
		    T2A = VFNMS(LDK(KP923879532), T2z, T2y);
		    T2Z = VFMA(LDK(KP923879532), T2K, T2J);
		    T2D = VFMA(LDK(KP198912367), T2C, T2B);
		    T2G = VFNMS(LDK(KP198912367), T2F, T2E);
		    T2H = VSUB(T2D, T2G);
		    T30 = VADD(T2D, T2G);
		    T2I = VFMA(LDK(KP980785280), T2H, T2A);
		    T35 = VFNMS(LDK(KP980785280), T30, T2Z);
		    T2S = VFNMS(LDK(KP980785280), T2H, T2A);
		    T31 = VFMA(LDK(KP980785280), T30, T2Z);
		    T2L = VFNMS(LDK(KP923879532), T2K, T2J);
		    T2W = VFMA(LDK(KP923879532), T2z, T2y);
		    T2M = VFMA(LDK(KP198912367), T2E, T2F);
		    T2N = VFNMS(LDK(KP198912367), T2B, T2C);
		    T2O = VSUB(T2M, T2N);
		    T2X = VADD(T2N, T2M);
		    T2P = VFMA(LDK(KP980785280), T2O, T2L);
		    T34 = VFNMS(LDK(KP980785280), T2X, T2W);
		    T2T = VFNMS(LDK(KP980785280), T2O, T2L);
		    T2Y = VFMA(LDK(KP980785280), T2X, T2W);
		    T2Q = VCONJ(VMUL(LDK(KP500000000), VFNMSI(T2P, T2I)));
		    ST(&(Rm[WS(rs, 6)]), T2Q, -ms, &(Rm[0]));
		    T36 = VCONJ(VMUL(LDK(KP500000000), VFNMSI(T35, T34)));
		    ST(&(Rm[WS(rs, 14)]), T36, -ms, &(Rm[0]));
		    T37 = VMUL(LDK(KP500000000), VFMAI(T35, T34));
		    ST(&(Rp[WS(rs, 15)]), T37, ms, &(Rp[WS(rs, 1)]));
		    T2R = VMUL(LDK(KP500000000), VFMAI(T2P, T2I));
		    ST(&(Rp[WS(rs, 7)]), T2R, ms, &(Rp[WS(rs, 1)]));
		    T2U = VMUL(LDK(KP500000000), VFNMSI(T2T, T2S));
		    ST(&(Rp[WS(rs, 9)]), T2U, ms, &(Rp[WS(rs, 1)]));
		    T32 = VMUL(LDK(KP500000000), VFNMSI(T31, T2Y));
		    ST(&(Rp[WS(rs, 1)]), T32, ms, &(Rp[WS(rs, 1)]));
		    T33 = VCONJ(VMUL(LDK(KP500000000), VFMAI(T31, T2Y)));
		    ST(&(Rm[0]), T33, -ms, &(Rm[0]));
		    T2V = VCONJ(VMUL(LDK(KP500000000), VFMAI(T2T, T2S)));
		    ST(&(Rm[WS(rs, 8)]), T2V, -ms, &(Rm[0]));
		    T16 = VFNMS(LDK(KP923879532), T15, Ts);
		    T2p = VFMA(LDK(KP923879532), T2a, T27);
		    T1v = VFMA(LDK(KP668178637), T1u, T1n);
		    T1Y = VFNMS(LDK(KP668178637), T1X, T1Q);
		    T1Z = VSUB(T1v, T1Y);
		    T2q = VADD(T1Y, T1v);
		    T20 = VFMA(LDK(KP831469612), T1Z, T16);
		    T2v = VFNMS(LDK(KP831469612), T2q, T2p);
		    T2i = VFNMS(LDK(KP831469612), T1Z, T16);
		    T2r = VFMA(LDK(KP831469612), T2q, T2p);
		    T2b = VFNMS(LDK(KP923879532), T2a, T27);
		    T2m = VFMA(LDK(KP923879532), T15, Ts);
		    T2c = VFNMS(LDK(KP668178637), T1n, T1u);
		    T2d = VFMA(LDK(KP668178637), T1Q, T1X);
		    T2e = VSUB(T2c, T2d);
		    T2n = VADD(T2d, T2c);
		    T2f = VFNMS(LDK(KP831469612), T2e, T2b);
		    T2u = VFNMS(LDK(KP831469612), T2n, T2m);
		    T2j = VFMA(LDK(KP831469612), T2e, T2b);
		    T2o = VFMA(LDK(KP831469612), T2n, T2m);
		    T2g = VMUL(LDK(KP500000000), VFNMSI(T2f, T20));
		    ST(&(Rp[WS(rs, 5)]), T2g, ms, &(Rp[WS(rs, 1)]));
		    T2w = VMUL(LDK(KP500000000), VFNMSI(T2v, T2u));
		    ST(&(Rp[WS(rs, 13)]), T2w, ms, &(Rp[WS(rs, 1)]));
		    T2x = VCONJ(VMUL(LDK(KP500000000), VFMAI(T2v, T2u)));
		    ST(&(Rm[WS(rs, 12)]), T2x, -ms, &(Rm[0]));
		    T2h = VCONJ(VMUL(LDK(KP500000000), VFMAI(T2f, T20)));
		    ST(&(Rm[WS(rs, 4)]), T2h, -ms, &(Rm[0]));
		    T2k = VMUL(LDK(KP500000000), VFMAI(T2j, T2i));
		    ST(&(Rp[WS(rs, 11)]), T2k, ms, &(Rp[WS(rs, 1)]));
		    T2s = VMUL(LDK(KP500000000), VFMAI(T2r, T2o));
		    ST(&(Rp[WS(rs, 3)]), T2s, ms, &(Rp[WS(rs, 1)]));
		    T2t = VCONJ(VMUL(LDK(KP500000000), VFNMSI(T2r, T2o)));
		    ST(&(Rm[WS(rs, 2)]), T2t, -ms, &(Rm[0]));
		    T2l = VCONJ(VMUL(LDK(KP500000000), VFNMSI(T2j, T2i)));
		    ST(&(Rm[WS(rs, 10)]), T2l, -ms, &(Rm[0]));
	       }
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(1, 1),
     VTW(1, 2),
     VTW(1, 3),
     VTW(1, 4),
     VTW(1, 5),
     VTW(1, 6),
     VTW(1, 7),
     VTW(1, 8),
     VTW(1, 9),
     VTW(1, 10),
     VTW(1, 11),
     VTW(1, 12),
     VTW(1, 13),
     VTW(1, 14),
     VTW(1, 15),
     VTW(1, 16),
     VTW(1, 17),
     VTW(1, 18),
     VTW(1, 19),
     VTW(1, 20),
     VTW(1, 21),
     VTW(1, 22),
     VTW(1, 23),
     VTW(1, 24),
     VTW(1, 25),
     VTW(1, 26),
     VTW(1, 27),
     VTW(1, 28),
     VTW(1, 29),
     VTW(1, 30),
     VTW(1, 31),
     { TW_NEXT, VL, 0 }
};

static const hc2c_desc desc = { 32, XSIMD_STRING("hc2cfdftv_32"), twinstr, &GENUS, { 119, 94, 130, 0 } };

void XSIMD(codelet_hc2cfdftv_32) (planner *p) {
     X(khc2c_register) (p, hc2cfdftv_32, &desc, HC2C_VIA_DFT);
}
#else

/* Generated by: ../../../genfft/gen_hc2cdft_c.native -simd -compact -variables 4 -pipeline-latency 8 -trivial-stores -variables 32 -no-generate-bytw -n 32 -dit -name hc2cfdftv_32 -include rdft/simd/hc2cfv.h */

/*
 * This function contains 249 FP additions, 133 FP multiplications,
 * (or, 233 additions, 117 multiplications, 16 fused multiply/add),
 * 130 stack variables, 9 constants, and 64 memory accesses
 */
#include "rdft/simd/hc2cfv.h"

static void hc2cfdftv_32(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP555570233, +0.555570233019602224742830813948532874374937191);
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DVK(KP195090322, +0.195090322016128267848284868477022240927691618);
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DVK(KP353553390, +0.353553390593273762200422181052424519642417969);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * ((TWVL / VL) * 62)); m < me; m = m + VL, Rp = Rp + (VL * ms), Ip = Ip + (VL * ms), Rm = Rm - (VL * ms), Im = Im - (VL * ms), W = W + (TWVL * 62), MAKE_VOLATILE_STRIDE(128, rs)) {
	       V Ta, T2m, Tx, T2h, T3R, T4h, T3q, T4g, T3B, T4n, T3E, T4o, T1B, T2S, T1O;
	       V T2R, TV, T2p, T1i, T2o, T3L, T4q, T3I, T4r, T3w, T4k, T3t, T4j, T26, T2V;
	       V T2d, T2U;
	       {
		    V T4, T1m, T1H, T2j, T1M, T2l, T9, T1o, Tf, T1r, Tq, T1w, Tv, T1y, Tk;
		    V T1t, Tl, Tw, T3P, T3Q, T3o, T3p, T3z, T3A, T3C, T3D, T1p, T1N, T1A, T1C;
		    V T1u, T1z;
		    {
			 V T1, T3, T2, T1l, T1G, T1F, T1E, T1D, T2i, T1L, T1K, T1J, T1I, T2k, T6;
			 V T8, T7, T5, T1n, Tc, Te, Td, Tb, T1q, Tn, Tp, To, Tm, T1v, Ts;
			 V Tu, Tt, Tr, T1x, Th, Tj, Ti, Tg, T1s;
			 T1 = LD(&(Rp[0]), ms, &(Rp[0]));
			 T2 = LD(&(Rm[0]), -ms, &(Rm[0]));
			 T3 = VCONJ(T2);
			 T4 = VADD(T1, T3);
			 T1l = LDW(&(W[0]));
			 T1m = VZMULIJ(T1l, VSUB(T3, T1));
			 T1G = LD(&(Rp[WS(rs, 4)]), ms, &(Rp[0]));
			 T1E = LD(&(Rm[WS(rs, 4)]), -ms, &(Rm[0]));
			 T1F = VCONJ(T1E);
			 T1D = LDW(&(W[TWVL * 16]));
			 T1H = VZMULIJ(T1D, VSUB(T1F, T1G));
			 T2i = LDW(&(W[TWVL * 14]));
			 T2j = VZMULJ(T2i, VADD(T1G, T1F));
			 T1L = LD(&(Rp[WS(rs, 12)]), ms, &(Rp[0]));
			 T1J = LD(&(Rm[WS(rs, 12)]), -ms, &(Rm[0]));
			 T1K = VCONJ(T1J);
			 T1I = LDW(&(W[TWVL * 48]));
			 T1M = VZMULIJ(T1I, VSUB(T1K, T1L));
			 T2k = LDW(&(W[TWVL * 46]));
			 T2l = VZMULJ(T2k, VADD(T1L, T1K));
			 T6 = LD(&(Rp[WS(rs, 8)]), ms, &(Rp[0]));
			 T7 = LD(&(Rm[WS(rs, 8)]), -ms, &(Rm[0]));
			 T8 = VCONJ(T7);
			 T5 = LDW(&(W[TWVL * 30]));
			 T9 = VZMULJ(T5, VADD(T6, T8));
			 T1n = LDW(&(W[TWVL * 32]));
			 T1o = VZMULIJ(T1n, VSUB(T8, T6));
			 Tc = LD(&(Rp[WS(rs, 2)]), ms, &(Rp[0]));
			 Td = LD(&(Rm[WS(rs, 2)]), -ms, &(Rm[0]));
			 Te = VCONJ(Td);
			 Tb = LDW(&(W[TWVL * 6]));
			 Tf = VZMULJ(Tb, VADD(Tc, Te));
			 T1q = LDW(&(W[TWVL * 8]));
			 T1r = VZMULIJ(T1q, VSUB(Te, Tc));
			 Tn = LD(&(Rp[WS(rs, 14)]), ms, &(Rp[0]));
			 To = LD(&(Rm[WS(rs, 14)]), -ms, &(Rm[0]));
			 Tp = VCONJ(To);
			 Tm = LDW(&(W[TWVL * 54]));
			 Tq = VZMULJ(Tm, VADD(Tn, Tp));
			 T1v = LDW(&(W[TWVL * 56]));
			 T1w = VZMULIJ(T1v, VSUB(Tp, Tn));
			 Ts = LD(&(Rp[WS(rs, 6)]), ms, &(Rp[0]));
			 Tt = LD(&(Rm[WS(rs, 6)]), -ms, &(Rm[0]));
			 Tu = VCONJ(Tt);
			 Tr = LDW(&(W[TWVL * 22]));
			 Tv = VZMULJ(Tr, VADD(Ts, Tu));
			 T1x = LDW(&(W[TWVL * 24]));
			 T1y = VZMULIJ(T1x, VSUB(Tu, Ts));
			 Th = LD(&(Rp[WS(rs, 10)]), ms, &(Rp[0]));
			 Ti = LD(&(Rm[WS(rs, 10)]), -ms, &(Rm[0]));
			 Tj = VCONJ(Ti);
			 Tg = LDW(&(W[TWVL * 38]));
			 Tk = VZMULJ(Tg, VADD(Th, Tj));
			 T1s = LDW(&(W[TWVL * 40]));
			 T1t = VZMULIJ(T1s, VSUB(Tj, Th));
		    }
		    Ta = VMUL(LDK(KP500000000), VSUB(T4, T9));
		    T2m = VSUB(T2j, T2l);
		    Tl = VSUB(Tf, Tk);
		    Tw = VSUB(Tq, Tv);
		    Tx = VMUL(LDK(KP353553390), VADD(Tl, Tw));
		    T2h = VMUL(LDK(KP707106781), VSUB(Tw, Tl));
		    T3P = VADD(Tq, Tv);
		    T3Q = VADD(Tf, Tk);
		    T3R = VSUB(T3P, T3Q);
		    T4h = VADD(T3Q, T3P);
		    T3o = VADD(T4, T9);
		    T3p = VADD(T2j, T2l);
		    T3q = VMUL(LDK(KP500000000), VSUB(T3o, T3p));
		    T4g = VADD(T3o, T3p);
		    T3z = VADD(T1m, T1o);
		    T3A = VADD(T1H, T1M);
		    T3B = VSUB(T3z, T3A);
		    T4n = VADD(T3z, T3A);
		    T3C = VADD(T1w, T1y);
		    T3D = VADD(T1r, T1t);
		    T3E = VSUB(T3C, T3D);
		    T4o = VADD(T3D, T3C);
		    T1p = VSUB(T1m, T1o);
		    T1N = VSUB(T1H, T1M);
		    T1u = VSUB(T1r, T1t);
		    T1z = VSUB(T1w, T1y);
		    T1A = VMUL(LDK(KP707106781), VADD(T1u, T1z));
		    T1C = VMUL(LDK(KP707106781), VSUB(T1z, T1u));
		    T1B = VADD(T1p, T1A);
		    T2S = VADD(T1N, T1C);
		    T1O = VSUB(T1C, T1N);
		    T2R = VSUB(T1p, T1A);
	       }
	       {
		    V TD, T1R, T1b, T29, T1g, T2b, TI, T1T, TO, T1Y, T10, T22, T15, T24, TT;
		    V T1W, TJ, TU, T16, T1h, T3J, T3K, T3G, T3H, T3u, T3v, T3r, T3s, T25, T2c;
		    V T20, T27, T1U, T1Z;
		    {
			 V TA, TC, TB, Tz, T1Q, T18, T1a, T19, T17, T28, T1d, T1f, T1e, T1c, T2a;
			 V TF, TH, TG, TE, T1S, TL, TN, TM, TK, T1X, TX, TZ, TY, TW, T21;
			 V T12, T14, T13, T11, T23, TQ, TS, TR, TP, T1V;
			 TA = LD(&(Rp[WS(rs, 1)]), ms, &(Rp[WS(rs, 1)]));
			 TB = LD(&(Rm[WS(rs, 1)]), -ms, &(Rm[WS(rs, 1)]));
			 TC = VCONJ(TB);
			 Tz = LDW(&(W[TWVL * 2]));
			 TD = VZMULJ(Tz, VADD(TA, TC));
			 T1Q = LDW(&(W[TWVL * 4]));
			 T1R = VZMULIJ(T1Q, VSUB(TC, TA));
			 T18 = LD(&(Rp[WS(rs, 3)]), ms, &(Rp[WS(rs, 1)]));
			 T19 = LD(&(Rm[WS(rs, 3)]), -ms, &(Rm[WS(rs, 1)]));
			 T1a = VCONJ(T19);
			 T17 = LDW(&(W[TWVL * 10]));
			 T1b = VZMULJ(T17, VADD(T18, T1a));
			 T28 = LDW(&(W[TWVL * 12]));
			 T29 = VZMULIJ(T28, VSUB(T1a, T18));
			 T1d = LD(&(Rp[WS(rs, 11)]), ms, &(Rp[WS(rs, 1)]));
			 T1e = LD(&(Rm[WS(rs, 11)]), -ms, &(Rm[WS(rs, 1)]));
			 T1f = VCONJ(T1e);
			 T1c = LDW(&(W[TWVL * 42]));
			 T1g = VZMULJ(T1c, VADD(T1d, T1f));
			 T2a = LDW(&(W[TWVL * 44]));
			 T2b = VZMULIJ(T2a, VSUB(T1f, T1d));
			 TF = LD(&(Rp[WS(rs, 9)]), ms, &(Rp[WS(rs, 1)]));
			 TG = LD(&(Rm[WS(rs, 9)]), -ms, &(Rm[WS(rs, 1)]));
			 TH = VCONJ(TG);
			 TE = LDW(&(W[TWVL * 34]));
			 TI = VZMULJ(TE, VADD(TF, TH));
			 T1S = LDW(&(W[TWVL * 36]));
			 T1T = VZMULIJ(T1S, VSUB(TH, TF));
			 TL = LD(&(Rp[WS(rs, 5)]), ms, &(Rp[WS(rs, 1)]));
			 TM = LD(&(Rm[WS(rs, 5)]), -ms, &(Rm[WS(rs, 1)]));
			 TN = VCONJ(TM);
			 TK = LDW(&(W[TWVL * 18]));
			 TO = VZMULJ(TK, VADD(TL, TN));
			 T1X = LDW(&(W[TWVL * 20]));
			 T1Y = VZMULIJ(T1X, VSUB(TN, TL));
			 TX = LD(&(Rp[WS(rs, 15)]), ms, &(Rp[WS(rs, 1)]));
			 TY = LD(&(Rm[WS(rs, 15)]), -ms, &(Rm[WS(rs, 1)]));
			 TZ = VCONJ(TY);
			 TW = LDW(&(W[TWVL * 58]));
			 T10 = VZMULJ(TW, VADD(TX, TZ));
			 T21 = LDW(&(W[TWVL * 60]));
			 T22 = VZMULIJ(T21, VSUB(TZ, TX));
			 T12 = LD(&(Rp[WS(rs, 7)]), ms, &(Rp[WS(rs, 1)]));
			 T13 = LD(&(Rm[WS(rs, 7)]), -ms, &(Rm[WS(rs, 1)]));
			 T14 = VCONJ(T13);
			 T11 = LDW(&(W[TWVL * 26]));
			 T15 = VZMULJ(T11, VADD(T12, T14));
			 T23 = LDW(&(W[TWVL * 28]));
			 T24 = VZMULIJ(T23, VSUB(T14, T12));
			 TQ = LD(&(Rp[WS(rs, 13)]), ms, &(Rp[WS(rs, 1)]));
			 TR = LD(&(Rm[WS(rs, 13)]), -ms, &(Rm[WS(rs, 1)]));
			 TS = VCONJ(TR);
			 TP = LDW(&(W[TWVL * 50]));
			 TT = VZMULJ(TP, VADD(TQ, TS));
			 T1V = LDW(&(W[TWVL * 52]));
			 T1W = VZMULIJ(T1V, VSUB(TS, TQ));
		    }
		    TJ = VSUB(TD, TI);
		    TU = VSUB(TO, TT);
		    TV = VFNMS(LDK(KP382683432), TU, VMUL(LDK(KP923879532), TJ));
		    T2p = VFMA(LDK(KP382683432), TJ, VMUL(LDK(KP923879532), TU));
		    T16 = VSUB(T10, T15);
		    T1h = VSUB(T1b, T1g);
		    T1i = VFMA(LDK(KP923879532), T16, VMUL(LDK(KP382683432), T1h));
		    T2o = VFNMS(LDK(KP923879532), T1h, VMUL(LDK(KP382683432), T16));
		    T3J = VADD(T1Y, T1W);
		    T3K = VADD(T1R, T1T);
		    T3L = VSUB(T3J, T3K);
		    T4q = VADD(T3K, T3J);
		    T3G = VADD(T22, T24);
		    T3H = VADD(T29, T2b);
		    T3I = VSUB(T3G, T3H);
		    T4r = VADD(T3G, T3H);
		    T3u = VADD(T10, T15);
		    T3v = VADD(T1b, T1g);
		    T3w = VSUB(T3u, T3v);
		    T4k = VADD(T3u, T3v);
		    T3r = VADD(TD, TI);
		    T3s = VADD(TO, TT);
		    T3t = VSUB(T3r, T3s);
		    T4j = VADD(T3r, T3s);
		    T25 = VSUB(T22, T24);
		    T2c = VSUB(T29, T2b);
		    T1U = VSUB(T1R, T1T);
		    T1Z = VSUB(T1W, T1Y);
		    T20 = VMUL(LDK(KP707106781), VADD(T1U, T1Z));
		    T27 = VMUL(LDK(KP707106781), VSUB(T1Z, T1U));
		    T26 = VADD(T20, T25);
		    T2V = VADD(T27, T2c);
		    T2d = VSUB(T27, T2c);
		    T2U = VSUB(T25, T20);
	       }
	       {
		    V T4m, T4w, T4t, T4x, T4i, T4l, T4p, T4s, T4u, T4z, T4v, T4y, T4E, T4L, T4H;
		    V T4K, T4A, T4F, T4D, T4G, T4B, T4C, T4I, T4N, T4J, T4M, T3O, T4c, T4d, T3X;
		    V T40, T46, T49, T41, T3y, T47, T3T, T45, T3N, T44, T3W, T48, T3x, T3S, T3F;
		    V T3M, T3U, T3V, T3Y, T4e, T4f, T3Z, T42, T4a, T4b, T43;
		    T4i = VADD(T4g, T4h);
		    T4l = VADD(T4j, T4k);
		    T4m = VADD(T4i, T4l);
		    T4w = VSUB(T4i, T4l);
		    T4p = VADD(T4n, T4o);
		    T4s = VADD(T4q, T4r);
		    T4t = VADD(T4p, T4s);
		    T4x = VBYI(VSUB(T4s, T4p));
		    T4u = VCONJ(VMUL(LDK(KP500000000), VSUB(T4m, T4t)));
		    ST(&(Rm[WS(rs, 15)]), T4u, -ms, &(Rm[WS(rs, 1)]));
		    T4z = VMUL(LDK(KP500000000), VADD(T4w, T4x));
		    ST(&(Rp[WS(rs, 8)]), T4z, ms, &(Rp[0]));
		    T4v = VMUL(LDK(KP500000000), VADD(T4m, T4t));
		    ST(&(Rp[0]), T4v, ms, &(Rp[0]));
		    T4y = VCONJ(VMUL(LDK(KP500000000), VSUB(T4w, T4x)));
		    ST(&(Rm[WS(rs, 7)]), T4y, -ms, &(Rm[WS(rs, 1)]));
		    T4A = VMUL(LDK(KP500000000), VSUB(T4g, T4h));
		    T4F = VSUB(T4k, T4j);
		    T4B = VSUB(T4n, T4o);
		    T4C = VSUB(T4r, T4q);
		    T4D = VMUL(LDK(KP353553390), VADD(T4B, T4C));
		    T4G = VMUL(LDK(KP707106781), VSUB(T4C, T4B));
		    T4E = VADD(T4A, T4D);
		    T4L = VMUL(LDK(KP500000000), VBYI(VSUB(T4G, T4F)));
		    T4H = VMUL(LDK(KP500000000), VBYI(VADD(T4F, T4G)));
		    T4K = VSUB(T4A, T4D);
		    T4I = VCONJ(VSUB(T4E, T4H));
		    ST(&(Rm[WS(rs, 3)]), T4I, -ms, &(Rm[WS(rs, 1)]));
		    T4N = VADD(T4K, T4L);
		    ST(&(Rp[WS(rs, 12)]), T4N, ms, &(Rp[0]));
		    T4J = VADD(T4E, T4H);
		    ST(&(Rp[WS(rs, 4)]), T4J, ms, &(Rp[0]));
		    T4M = VCONJ(VSUB(T4K, T4L));
		    ST(&(Rm[WS(rs, 11)]), T4M, -ms, &(Rm[WS(rs, 1)]));
		    T3x = VMUL(LDK(KP353553390), VADD(T3t, T3w));
		    T3y = VADD(T3q, T3x);
		    T47 = VSUB(T3q, T3x);
		    T3S = VMUL(LDK(KP707106781), VSUB(T3w, T3t));
		    T3T = VADD(T3R, T3S);
		    T45 = VSUB(T3S, T3R);
		    T3F = VFMA(LDK(KP923879532), T3B, VMUL(LDK(KP382683432), T3E));
		    T3M = VFNMS(LDK(KP382683432), T3L, VMUL(LDK(KP923879532), T3I));
		    T3N = VMUL(LDK(KP500000000), VADD(T3F, T3M));
		    T44 = VSUB(T3M, T3F);
		    T3U = VFNMS(LDK(KP382683432), T3B, VMUL(LDK(KP923879532), T3E));
		    T3V = VFMA(LDK(KP923879532), T3L, VMUL(LDK(KP382683432), T3I));
		    T3W = VADD(T3U, T3V);
		    T48 = VMUL(LDK(KP500000000), VSUB(T3V, T3U));
		    T3O = VADD(T3y, T3N);
		    T4c = VMUL(LDK(KP500000000), VBYI(VADD(T45, T44)));
		    T4d = VADD(T47, T48);
		    T3X = VMUL(LDK(KP500000000), VBYI(VADD(T3T, T3W)));
		    T40 = VSUB(T3y, T3N);
		    T46 = VMUL(LDK(KP500000000), VBYI(VSUB(T44, T45)));
		    T49 = VSUB(T47, T48);
		    T41 = VMUL(LDK(KP500000000), VBYI(VSUB(T3W, T3T)));
		    T3Y = VCONJ(VSUB(T3O, T3X));
		    ST(&(Rm[WS(rs, 1)]), T3Y, -ms, &(Rm[WS(rs, 1)]));
		    T4e = VADD(T4c, T4d);
		    ST(&(Rp[WS(rs, 6)]), T4e, ms, &(Rp[0]));
		    T4f = VCONJ(VSUB(T4d, T4c));
		    ST(&(Rm[WS(rs, 5)]), T4f, -ms, &(Rm[WS(rs, 1)]));
		    T3Z = VADD(T3O, T3X);
		    ST(&(Rp[WS(rs, 2)]), T3Z, ms, &(Rp[0]));
		    T42 = VCONJ(VSUB(T40, T41));
		    ST(&(Rm[WS(rs, 13)]), T42, -ms, &(Rm[WS(rs, 1)]));
		    T4a = VADD(T46, T49);
		    ST(&(Rp[WS(rs, 10)]), T4a, ms, &(Rp[0]));
		    T4b = VCONJ(VSUB(T49, T46));
		    ST(&(Rm[WS(rs, 9)]), T4b, -ms, &(Rm[WS(rs, 1)]));
		    T43 = VADD(T40, T41);
		    ST(&(Rp[WS(rs, 14)]), T43, ms, &(Rp[0]));
		    {
			 V T2g, T2K, T2L, T2v, T2y, T2E, T2H, T2z, T1k, T2F, T2u, T2G, T2f, T2C, T2r;
			 V T2D, Ty, T1j, T2s, T2t, T1P, T2e, T2n, T2q, T2w, T2M, T2N, T2x, T2A, T2I;
			 V T2J, T2B;
			 Ty = VADD(Ta, Tx);
			 T1j = VMUL(LDK(KP500000000), VADD(TV, T1i));
			 T1k = VADD(Ty, T1j);
			 T2F = VSUB(Ty, T1j);
			 T2s = VFNMS(LDK(KP195090322), T1B, VMUL(LDK(KP980785280), T1O));
			 T2t = VFMA(LDK(KP195090322), T26, VMUL(LDK(KP980785280), T2d));
			 T2u = VADD(T2s, T2t);
			 T2G = VMUL(LDK(KP500000000), VSUB(T2t, T2s));
			 T1P = VFMA(LDK(KP980785280), T1B, VMUL(LDK(KP195090322), T1O));
			 T2e = VFNMS(LDK(KP195090322), T2d, VMUL(LDK(KP980785280), T26));
			 T2f = VMUL(LDK(KP500000000), VADD(T1P, T2e));
			 T2C = VSUB(T2e, T1P);
			 T2n = VSUB(T2h, T2m);
			 T2q = VSUB(T2o, T2p);
			 T2r = VADD(T2n, T2q);
			 T2D = VSUB(T2q, T2n);
			 T2g = VADD(T1k, T2f);
			 T2K = VMUL(LDK(KP500000000), VBYI(VADD(T2D, T2C)));
			 T2L = VADD(T2F, T2G);
			 T2v = VMUL(LDK(KP500000000), VBYI(VADD(T2r, T2u)));
			 T2y = VSUB(T1k, T2f);
			 T2E = VMUL(LDK(KP500000000), VBYI(VSUB(T2C, T2D)));
			 T2H = VSUB(T2F, T2G);
			 T2z = VMUL(LDK(KP500000000), VBYI(VSUB(T2u, T2r)));
			 T2w = VCONJ(VSUB(T2g, T2v));
			 ST(&(Rm[0]), T2w, -ms, &(Rm[0]));
			 T2M = VADD(T2K, T2L);
			 ST(&(Rp[WS(rs, 7)]), T2M, ms, &(Rp[WS(rs, 1)]));
			 T2N = VCONJ(VSUB(T2L, T2K));
			 ST(&(Rm[WS(rs, 6)]), T2N, -ms, &(Rm[0]));
			 T2x = VADD(T2g, T2v);
			 ST(&(Rp[WS(rs, 1)]), T2x, ms, &(Rp[WS(rs, 1)]));
			 T2A = VCONJ(VSUB(T2y, T2z));
			 ST(&(Rm[WS(rs, 14)]), T2A, -ms, &(Rm[0]));
			 T2I = VADD(T2E, T2H);
			 ST(&(Rp[WS(rs, 9)]), T2I, ms, &(Rp[WS(rs, 1)]));
			 T2J = VCONJ(VSUB(T2H, T2E));
			 ST(&(Rm[WS(rs, 8)]), T2J, -ms, &(Rm[0]));
			 T2B = VADD(T2y, T2z);
			 ST(&(Rp[WS(rs, 15)]), T2B, ms, &(Rp[WS(rs, 1)]));
		    }
		    {
			 V T2Y, T3k, T3l, T35, T38, T3e, T3h, T39, T2Q, T3f, T34, T3g, T2X, T3c, T31;
			 V T3d, T2O, T2P, T32, T33, T2T, T2W, T2Z, T30, T36, T3m, T3n, T37, T3a, T3i;
			 V T3j, T3b;
			 T2O = VSUB(Ta, Tx);
			 T2P = VMUL(LDK(KP500000000), VADD(T2p, T2o));
			 T2Q = VADD(T2O, T2P);
			 T3f = VSUB(T2O, T2P);
			 T32 = VFNMS(LDK(KP555570233), T2R, VMUL(LDK(KP831469612), T2S));
			 T33 = VFMA(LDK(KP555570233), T2U, VMUL(LDK(KP831469612), T2V));
			 T34 = VADD(T32, T33);
			 T3g = VMUL(LDK(KP500000000), VSUB(T33, T32));
			 T2T = VFMA(LDK(KP831469612), T2R, VMUL(LDK(KP555570233), T2S));
			 T2W = VFNMS(LDK(KP555570233), T2V, VMUL(LDK(KP831469612), T2U));
			 T2X = VMUL(LDK(KP500000000), VADD(T2T, T2W));
			 T3c = VSUB(T2W, T2T);
			 T2Z = VADD(T2m, T2h);
			 T30 = VSUB(T1i, TV);
			 T31 = VADD(T2Z, T30);
			 T3d = VSUB(T30, T2Z);
			 T2Y = VADD(T2Q, T2X);
			 T3k = VMUL(LDK(KP500000000), VBYI(VADD(T3d, T3c)));
			 T3l = VADD(T3f, T3g);
			 T35 = VMUL(LDK(KP500000000), VBYI(VADD(T31, T34)));
			 T38 = VSUB(T2Q, T2X);
			 T3e = VMUL(LDK(KP500000000), VBYI(VSUB(T3c, T3d)));
			 T3h = VSUB(T3f, T3g);
			 T39 = VMUL(LDK(KP500000000), VBYI(VSUB(T34, T31)));
			 T36 = VCONJ(VSUB(T2Y, T35));
			 ST(&(Rm[WS(rs, 2)]), T36, -ms, &(Rm[0]));
			 T3m = VADD(T3k, T3l);
			 ST(&(Rp[WS(rs, 5)]), T3m, ms, &(Rp[WS(rs, 1)]));
			 T3n = VCONJ(VSUB(T3l, T3k));
			 ST(&(Rm[WS(rs, 4)]), T3n, -ms, &(Rm[0]));
			 T37 = VADD(T2Y, T35);
			 ST(&(Rp[WS(rs, 3)]), T37, ms, &(Rp[WS(rs, 1)]));
			 T3a = VCONJ(VSUB(T38, T39));
			 ST(&(Rm[WS(rs, 12)]), T3a, -ms, &(Rm[0]));
			 T3i = VADD(T3e, T3h);
			 ST(&(Rp[WS(rs, 11)]), T3i, ms, &(Rp[WS(rs, 1)]));
			 T3j = VCONJ(VSUB(T3h, T3e));
			 ST(&(Rm[WS(rs, 10)]), T3j, -ms, &(Rm[0]));
			 T3b = VADD(T38, T39);
			 ST(&(Rp[WS(rs, 13)]), T3b, ms, &(Rp[WS(rs, 1)]));
		    }
	       }
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(1, 1),
     VTW(1, 2),
     VTW(1, 3),
     VTW(1, 4),
     VTW(1, 5),
     VTW(1, 6),
     VTW(1, 7),
     VTW(1, 8),
     VTW(1, 9),
     VTW(1, 10),
     VTW(1, 11),
     VTW(1, 12),
     VTW(1, 13),
     VTW(1, 14),
     VTW(1, 15),
     VTW(1, 16),
     VTW(1, 17),
     VTW(1, 18),
     VTW(1, 19),
     VTW(1, 20),
     VTW(1, 21),
     VTW(1, 22),
     VTW(1, 23),
     VTW(1, 24),
     VTW(1, 25),
     VTW(1, 26),
     VTW(1, 27),
     VTW(1, 28),
     VTW(1, 29),
     VTW(1, 30),
     VTW(1, 31),
     { TW_NEXT, VL, 0 }
};

static const hc2c_desc desc = { 32, XSIMD_STRING("hc2cfdftv_32"), twinstr, &GENUS, { 233, 117, 16, 0 } };

void XSIMD(codelet_hc2cfdftv_32) (planner *p) {
     X(khc2c_register) (p, hc2cfdftv_32, &desc, HC2C_VIA_DFT);
}
#endif
