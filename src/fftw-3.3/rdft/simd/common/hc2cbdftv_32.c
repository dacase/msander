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
/* Generated on Thu Dec 10 07:06:58 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2cdft_c.native -fma -simd -compact -variables 4 -pipeline-latency 8 -trivial-stores -variables 32 -no-generate-bytw -n 32 -dif -sign 1 -name hc2cbdftv_32 -include rdft/simd/hc2cbv.h */

/*
 * This function contains 249 FP additions, 192 FP multiplications,
 * (or, 119 additions, 62 multiplications, 130 fused multiply/add),
 * 143 stack variables, 7 constants, and 64 memory accesses
 */
#include "rdft/simd/hc2cbv.h"

static void hc2cbdftv_32(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DVK(KP198912367, +0.198912367379658006911597622644676228597850501);
     DVK(KP668178637, +0.668178637919298919997757686523080761552472251);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DVK(KP414213562, +0.414213562373095048801688724209698078569671875);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * ((TWVL / VL) * 62)); m < me; m = m + VL, Rp = Rp + (VL * ms), Ip = Ip + (VL * ms), Rm = Rm - (VL * ms), Im = Im - (VL * ms), W = W + (TWVL * 62), MAKE_VOLATILE_STRIDE(128, rs)) {
	       V Ts, T1S, T3p, T45, T3A, T48, T1b, T1V, T1o, T2G, T2o, T2Y, T2z, T31, T1L;
	       V T2H, T2J, T2K, TJ, T1c, T3D, T46, T10, T1d, T2r, T2A, T3w, T49, T1D, T1M;
	       V T2u, T2B;
	       {
		    V T4, T1i, T15, T1j, Tb, T1m, T16, T1l, T1G, T1F, Tj, T3m, T18, T1J, T1I;
		    V Tq, T3n, T19, T2, T3, T13, T14, T5, T6, T7, T8, T9, Ta, Tf, Ti;
		    V Td, Te, Tg, Th, Tm, Tp, Tk, Tl, Tn, To, Tc, Tr, T3l, T3o, T3y;
		    V T3z, T17, T1a, T1k, T1n, T2m, T2n, T2x, T2y, T1H, T1K;
		    T2 = LD(&(Rp[0]), ms, &(Rp[0]));
		    T3 = LD(&(Rm[WS(rs, 15)]), -ms, &(Rm[WS(rs, 1)]));
		    T4 = VFNMSCONJ(T3, T2);
		    T1i = VFMACONJ(T3, T2);
		    T13 = LD(&(Rp[WS(rs, 8)]), ms, &(Rp[0]));
		    T14 = LD(&(Rm[WS(rs, 7)]), -ms, &(Rm[WS(rs, 1)]));
		    T15 = VFNMSCONJ(T14, T13);
		    T1j = VFMACONJ(T14, T13);
		    T5 = LD(&(Rp[WS(rs, 4)]), ms, &(Rp[0]));
		    T6 = LD(&(Rm[WS(rs, 11)]), -ms, &(Rm[WS(rs, 1)]));
		    T7 = VFNMSCONJ(T6, T5);
		    T8 = LD(&(Rp[WS(rs, 12)]), ms, &(Rp[0]));
		    T9 = LD(&(Rm[WS(rs, 3)]), -ms, &(Rm[WS(rs, 1)]));
		    Ta = VFMSCONJ(T9, T8);
		    Tb = VADD(T7, Ta);
		    T1m = VFMACONJ(T9, T8);
		    T16 = VSUB(T7, Ta);
		    T1l = VFMACONJ(T6, T5);
		    Td = LD(&(Rp[WS(rs, 10)]), ms, &(Rp[0]));
		    Te = LD(&(Rm[WS(rs, 5)]), -ms, &(Rm[WS(rs, 1)]));
		    Tf = VFNMSCONJ(Te, Td);
		    T1G = VFMACONJ(Te, Td);
		    Tg = LD(&(Rp[WS(rs, 2)]), ms, &(Rp[0]));
		    Th = LD(&(Rm[WS(rs, 13)]), -ms, &(Rm[WS(rs, 1)]));
		    Ti = VFNMSCONJ(Th, Tg);
		    T1F = VFMACONJ(Th, Tg);
		    Tj = VFMA(LDK(KP414213562), Ti, Tf);
		    T3m = VSUB(T1F, T1G);
		    T18 = VFNMS(LDK(KP414213562), Tf, Ti);
		    Tk = LD(&(Rp[WS(rs, 6)]), ms, &(Rp[0]));
		    Tl = LD(&(Rm[WS(rs, 9)]), -ms, &(Rm[WS(rs, 1)]));
		    Tm = VFNMSCONJ(Tl, Tk);
		    T1J = VFMACONJ(Tl, Tk);
		    Tn = LD(&(Rp[WS(rs, 14)]), ms, &(Rp[0]));
		    To = LD(&(Rm[WS(rs, 1)]), -ms, &(Rm[WS(rs, 1)]));
		    Tp = VFMSCONJ(To, Tn);
		    T1I = VFMACONJ(To, Tn);
		    Tq = VFNMS(LDK(KP414213562), Tp, Tm);
		    T3n = VSUB(T1I, T1J);
		    T19 = VFMA(LDK(KP414213562), Tm, Tp);
		    Tc = VFNMS(LDK(KP707106781), Tb, T4);
		    Tr = VSUB(Tj, Tq);
		    Ts = VFMA(LDK(KP923879532), Tr, Tc);
		    T1S = VFNMS(LDK(KP923879532), Tr, Tc);
		    T3l = VSUB(T1i, T1j);
		    T3o = VADD(T3m, T3n);
		    T3p = VFMA(LDK(KP707106781), T3o, T3l);
		    T45 = VFNMS(LDK(KP707106781), T3o, T3l);
		    T3y = VSUB(T1l, T1m);
		    T3z = VSUB(T3m, T3n);
		    T3A = VFMA(LDK(KP707106781), T3z, T3y);
		    T48 = VFNMS(LDK(KP707106781), T3z, T3y);
		    T17 = VFNMS(LDK(KP707106781), T16, T15);
		    T1a = VSUB(T18, T19);
		    T1b = VFNMS(LDK(KP923879532), T1a, T17);
		    T1V = VFMA(LDK(KP923879532), T1a, T17);
		    T1k = VADD(T1i, T1j);
		    T1n = VADD(T1l, T1m);
		    T1o = VSUB(T1k, T1n);
		    T2G = VADD(T1k, T1n);
		    T2m = VFMA(LDK(KP707106781), Tb, T4);
		    T2n = VADD(T18, T19);
		    T2o = VFNMS(LDK(KP923879532), T2n, T2m);
		    T2Y = VFMA(LDK(KP923879532), T2n, T2m);
		    T2x = VFMA(LDK(KP707106781), T16, T15);
		    T2y = VADD(Tj, Tq);
		    T2z = VFNMS(LDK(KP923879532), T2y, T2x);
		    T31 = VFMA(LDK(KP923879532), T2y, T2x);
		    T1H = VADD(T1F, T1G);
		    T1K = VADD(T1I, T1J);
		    T1L = VSUB(T1H, T1K);
		    T2H = VADD(T1H, T1K);
	       }
	       {
		    V Tv, T3q, TG, T1r, TM, T3t, TX, T1y, TC, T3r, TH, T1u, TT, T3u, TY;
		    V T1B, Tt, Tu, T1p, TE, TF, T1q, TK, TL, T1w, TV, TW, T1x, Ty, T1s;
		    V TB, T1t, Tw, Tx, Tz, TA, TP, T1z, TS, T1A, TN, TO, TQ, TR, TD;
		    V TI, T3B, T3C, TU, TZ, T2p, T2q, T3s, T3v, T1v, T1C, T2s, T2t;
		    Tt = LD(&(Rp[WS(rs, 1)]), ms, &(Rp[WS(rs, 1)]));
		    Tu = LD(&(Rm[WS(rs, 14)]), -ms, &(Rm[0]));
		    T1p = VFMACONJ(Tu, Tt);
		    TE = LD(&(Rp[WS(rs, 9)]), ms, &(Rp[WS(rs, 1)]));
		    TF = LD(&(Rm[WS(rs, 6)]), -ms, &(Rm[0]));
		    T1q = VFMACONJ(TF, TE);
		    Tv = VFNMSCONJ(Tu, Tt);
		    T3q = VSUB(T1p, T1q);
		    TG = VFNMSCONJ(TF, TE);
		    T1r = VADD(T1p, T1q);
		    TK = LD(&(Rp[WS(rs, 15)]), ms, &(Rp[WS(rs, 1)]));
		    TL = LD(&(Rm[0]), -ms, &(Rm[0]));
		    T1w = VFMACONJ(TL, TK);
		    TV = LD(&(Rp[WS(rs, 7)]), ms, &(Rp[WS(rs, 1)]));
		    TW = LD(&(Rm[WS(rs, 8)]), -ms, &(Rm[0]));
		    T1x = VFMACONJ(TW, TV);
		    TM = VFMSCONJ(TL, TK);
		    T3t = VSUB(T1w, T1x);
		    TX = VFNMSCONJ(TW, TV);
		    T1y = VADD(T1w, T1x);
		    Tw = LD(&(Rp[WS(rs, 5)]), ms, &(Rp[WS(rs, 1)]));
		    Tx = LD(&(Rm[WS(rs, 10)]), -ms, &(Rm[0]));
		    Ty = VFNMSCONJ(Tx, Tw);
		    T1s = VFMACONJ(Tx, Tw);
		    Tz = LD(&(Rp[WS(rs, 13)]), ms, &(Rp[WS(rs, 1)]));
		    TA = LD(&(Rm[WS(rs, 2)]), -ms, &(Rm[0]));
		    TB = VFMSCONJ(TA, Tz);
		    T1t = VFMACONJ(TA, Tz);
		    TC = VADD(Ty, TB);
		    T3r = VSUB(T1s, T1t);
		    TH = VSUB(Ty, TB);
		    T1u = VADD(T1s, T1t);
		    TN = LD(&(Rp[WS(rs, 3)]), ms, &(Rp[WS(rs, 1)]));
		    TO = LD(&(Rm[WS(rs, 12)]), -ms, &(Rm[0]));
		    TP = VFNMSCONJ(TO, TN);
		    T1z = VFMACONJ(TO, TN);
		    TQ = LD(&(Rp[WS(rs, 11)]), ms, &(Rp[WS(rs, 1)]));
		    TR = LD(&(Rm[WS(rs, 4)]), -ms, &(Rm[0]));
		    TS = VFMSCONJ(TR, TQ);
		    T1A = VFMACONJ(TR, TQ);
		    TT = VADD(TP, TS);
		    T3u = VSUB(T1A, T1z);
		    TY = VSUB(TS, TP);
		    T1B = VADD(T1z, T1A);
		    T2J = VADD(T1r, T1u);
		    T2K = VADD(T1y, T1B);
		    TD = VFNMS(LDK(KP707106781), TC, Tv);
		    TI = VFNMS(LDK(KP707106781), TH, TG);
		    TJ = VFMA(LDK(KP668178637), TI, TD);
		    T1c = VFNMS(LDK(KP668178637), TD, TI);
		    T3B = VFMA(LDK(KP414213562), T3q, T3r);
		    T3C = VFMA(LDK(KP414213562), T3t, T3u);
		    T3D = VSUB(T3B, T3C);
		    T46 = VADD(T3B, T3C);
		    TU = VFNMS(LDK(KP707106781), TT, TM);
		    TZ = VFMA(LDK(KP707106781), TY, TX);
		    T10 = VFNMS(LDK(KP668178637), TZ, TU);
		    T1d = VFMA(LDK(KP668178637), TU, TZ);
		    T2p = VFMA(LDK(KP707106781), TH, TG);
		    T2q = VFMA(LDK(KP707106781), TC, Tv);
		    T2r = VFMA(LDK(KP198912367), T2q, T2p);
		    T2A = VFNMS(LDK(KP198912367), T2p, T2q);
		    T3s = VFNMS(LDK(KP414213562), T3r, T3q);
		    T3v = VFNMS(LDK(KP414213562), T3u, T3t);
		    T3w = VADD(T3s, T3v);
		    T49 = VSUB(T3s, T3v);
		    T1v = VSUB(T1r, T1u);
		    T1C = VSUB(T1y, T1B);
		    T1D = VADD(T1v, T1C);
		    T1M = VSUB(T1v, T1C);
		    T2s = VFNMS(LDK(KP707106781), TY, TX);
		    T2t = VFMA(LDK(KP707106781), TT, TM);
		    T2u = VFNMS(LDK(KP198912367), T2t, T2s);
		    T2B = VFMA(LDK(KP198912367), T2s, T2t);
	       }
	       {
		    V T3f, T38, T4p, T4v, T3T, T3Z, T2a, T2i, T4b, T4h, T1O, T20, T2M, T2U, T3F;
		    V T3L, T1g, T3X, T2g, T3J, T2E, T4l, T2S, T4f, T1Y, T4t, T26, T43, T34, T3P;
		    V T3e, T3j, T36, T37, T35, T4n, T4o, T4m, T4u, T3R, T3S, T3Q, T3Y, T28, T29;
		    V T27, T2h, T47, T4a, T44, T4g, T1E, T1N, T1h, T1Z;
		    T36 = VADD(T2G, T2H);
		    T37 = VADD(T2J, T2K);
		    T3f = VADD(T36, T37);
		    T35 = LDW(&(W[TWVL * 30]));
		    T38 = VZMUL(T35, VSUB(T36, T37));
		    T4n = VFMA(LDK(KP923879532), T46, T45);
		    T4o = VFNMS(LDK(KP923879532), T49, T48);
		    T4m = LDW(&(W[TWVL * 10]));
		    T4p = VZMUL(T4m, VFNMSI(T4o, T4n));
		    T4u = LDW(&(W[TWVL * 50]));
		    T4v = VZMUL(T4u, VFMAI(T4o, T4n));
		    T3R = VFMA(LDK(KP923879532), T3w, T3p);
		    T3S = VFMA(LDK(KP923879532), T3D, T3A);
		    T3Q = LDW(&(W[TWVL * 58]));
		    T3T = VZMUL(T3Q, VFNMSI(T3S, T3R));
		    T3Y = LDW(&(W[TWVL * 2]));
		    T3Z = VZMUL(T3Y, VFMAI(T3S, T3R));
		    T28 = VFMA(LDK(KP707106781), T1D, T1o);
		    T29 = VFMA(LDK(KP707106781), T1M, T1L);
		    T27 = LDW(&(W[TWVL * 6]));
		    T2a = VZMUL(T27, VFMAI(T29, T28));
		    T2h = LDW(&(W[TWVL * 54]));
		    T2i = VZMUL(T2h, VFNMSI(T29, T28));
		    T47 = VFNMS(LDK(KP923879532), T46, T45);
		    T4a = VFMA(LDK(KP923879532), T49, T48);
		    T44 = LDW(&(W[TWVL * 18]));
		    T4b = VZMUL(T44, VFMAI(T4a, T47));
		    T4g = LDW(&(W[TWVL * 42]));
		    T4h = VZMUL(T4g, VFNMSI(T4a, T47));
		    T1E = VFNMS(LDK(KP707106781), T1D, T1o);
		    T1N = VFNMS(LDK(KP707106781), T1M, T1L);
		    T1h = LDW(&(W[TWVL * 22]));
		    T1O = VZMUL(T1h, VFNMSI(T1N, T1E));
		    T1Z = LDW(&(W[TWVL * 38]));
		    T20 = VZMUL(T1Z, VFMAI(T1N, T1E));
		    {
			 V T2I, T2L, T2F, T2T, T3x, T3E, T3k, T3K, T12, T2e, T1f, T2f, T11, T1e, T1;
			 V T3W, T2d, T3I, T2w, T2Q, T2D, T2R, T2v, T2C, T2l, T4k, T2P, T4e, T1U, T24;
			 V T1X, T25, T1T, T1W, T1R, T4s, T23, T42, T30, T3c, T33, T3d, T2Z, T32, T2X;
			 V T3O, T3b, T3i;
			 T2I = VSUB(T2G, T2H);
			 T2L = VSUB(T2J, T2K);
			 T2F = LDW(&(W[TWVL * 46]));
			 T2M = VZMUL(T2F, VFNMSI(T2L, T2I));
			 T2T = LDW(&(W[TWVL * 14]));
			 T2U = VZMUL(T2T, VFMAI(T2L, T2I));
			 T3x = VFNMS(LDK(KP923879532), T3w, T3p);
			 T3E = VFNMS(LDK(KP923879532), T3D, T3A);
			 T3k = LDW(&(W[TWVL * 26]));
			 T3F = VZMUL(T3k, VFNMSI(T3E, T3x));
			 T3K = LDW(&(W[TWVL * 34]));
			 T3L = VZMUL(T3K, VFMAI(T3E, T3x));
			 T11 = VADD(TJ, T10);
			 T12 = VFNMS(LDK(KP831469612), T11, Ts);
			 T2e = VFMA(LDK(KP831469612), T11, Ts);
			 T1e = VADD(T1c, T1d);
			 T1f = VFNMS(LDK(KP831469612), T1e, T1b);
			 T2f = VFMA(LDK(KP831469612), T1e, T1b);
			 T1 = LDW(&(W[TWVL * 24]));
			 T1g = VZMULI(T1, VFMAI(T1f, T12));
			 T3W = LDW(&(W[TWVL * 4]));
			 T3X = VZMULI(T3W, VFNMSI(T2f, T2e));
			 T2d = LDW(&(W[TWVL * 56]));
			 T2g = VZMULI(T2d, VFMAI(T2f, T2e));
			 T3I = LDW(&(W[TWVL * 36]));
			 T3J = VZMULI(T3I, VFNMSI(T1f, T12));
			 T2v = VSUB(T2r, T2u);
			 T2w = VFMA(LDK(KP980785280), T2v, T2o);
			 T2Q = VFNMS(LDK(KP980785280), T2v, T2o);
			 T2C = VSUB(T2A, T2B);
			 T2D = VFNMS(LDK(KP980785280), T2C, T2z);
			 T2R = VFMA(LDK(KP980785280), T2C, T2z);
			 T2l = LDW(&(W[TWVL * 48]));
			 T2E = VZMULI(T2l, VFMAI(T2D, T2w));
			 T4k = LDW(&(W[TWVL * 12]));
			 T4l = VZMULI(T4k, VFNMSI(T2D, T2w));
			 T2P = LDW(&(W[TWVL * 16]));
			 T2S = VZMULI(T2P, VFMAI(T2R, T2Q));
			 T4e = LDW(&(W[TWVL * 44]));
			 T4f = VZMULI(T4e, VFNMSI(T2R, T2Q));
			 T1T = VSUB(T1d, T1c);
			 T1U = VFNMS(LDK(KP831469612), T1T, T1S);
			 T24 = VFMA(LDK(KP831469612), T1T, T1S);
			 T1W = VSUB(TJ, T10);
			 T1X = VFNMS(LDK(KP831469612), T1W, T1V);
			 T25 = VFMA(LDK(KP831469612), T1W, T1V);
			 T1R = LDW(&(W[TWVL * 40]));
			 T1Y = VZMULI(T1R, VFMAI(T1X, T1U));
			 T4s = LDW(&(W[TWVL * 52]));
			 T4t = VZMULI(T4s, VFNMSI(T25, T24));
			 T23 = LDW(&(W[TWVL * 8]));
			 T26 = VZMULI(T23, VFMAI(T25, T24));
			 T42 = LDW(&(W[TWVL * 20]));
			 T43 = VZMULI(T42, VFNMSI(T1X, T1U));
			 T2Z = VADD(T2A, T2B);
			 T30 = VFNMS(LDK(KP980785280), T2Z, T2Y);
			 T3c = VFMA(LDK(KP980785280), T2Z, T2Y);
			 T32 = VADD(T2r, T2u);
			 T33 = VFNMS(LDK(KP980785280), T32, T31);
			 T3d = VFMA(LDK(KP980785280), T32, T31);
			 T2X = LDW(&(W[TWVL * 32]));
			 T34 = VZMULI(T2X, VFMAI(T33, T30));
			 T3O = LDW(&(W[TWVL * 60]));
			 T3P = VZMULI(T3O, VFNMSI(T3d, T3c));
			 T3b = LDW(&(W[0]));
			 T3e = VZMULI(T3b, VFMAI(T3d, T3c));
			 T3i = LDW(&(W[TWVL * 28]));
			 T3j = VZMULI(T3i, VFNMSI(T33, T30));
		    }
		    {
			 V T1P, T4w, T2j, T4c, T4x, T1Q, T4d, T2k, T21, T4q, T2b, T4i, T4r, T22, T4j;
			 V T2c, T2N, T40, T3g, T3G, T41, T2O, T3H, T3h, T2V, T3U, T39, T3M, T3V, T2W;
			 V T3N, T3a;
			 T1P = VADD(T1g, T1O);
			 ST(&(Rp[WS(rs, 6)]), T1P, ms, &(Rp[0]));
			 T4w = VADD(T4t, T4v);
			 ST(&(Rp[WS(rs, 13)]), T4w, ms, &(Rp[WS(rs, 1)]));
			 T2j = VADD(T2g, T2i);
			 ST(&(Rp[WS(rs, 14)]), T2j, ms, &(Rp[0]));
			 T4c = VADD(T43, T4b);
			 ST(&(Rp[WS(rs, 5)]), T4c, ms, &(Rp[WS(rs, 1)]));
			 T4x = VCONJ(VSUB(T4v, T4t));
			 ST(&(Rm[WS(rs, 13)]), T4x, -ms, &(Rm[WS(rs, 1)]));
			 T1Q = VCONJ(VSUB(T1O, T1g));
			 ST(&(Rm[WS(rs, 6)]), T1Q, -ms, &(Rm[0]));
			 T4d = VCONJ(VSUB(T4b, T43));
			 ST(&(Rm[WS(rs, 5)]), T4d, -ms, &(Rm[WS(rs, 1)]));
			 T2k = VCONJ(VSUB(T2i, T2g));
			 ST(&(Rm[WS(rs, 14)]), T2k, -ms, &(Rm[0]));
			 T21 = VADD(T1Y, T20);
			 ST(&(Rp[WS(rs, 10)]), T21, ms, &(Rp[0]));
			 T4q = VADD(T4l, T4p);
			 ST(&(Rp[WS(rs, 3)]), T4q, ms, &(Rp[WS(rs, 1)]));
			 T2b = VADD(T26, T2a);
			 ST(&(Rp[WS(rs, 2)]), T2b, ms, &(Rp[0]));
			 T4i = VADD(T4f, T4h);
			 ST(&(Rp[WS(rs, 11)]), T4i, ms, &(Rp[WS(rs, 1)]));
			 T4r = VCONJ(VSUB(T4p, T4l));
			 ST(&(Rm[WS(rs, 3)]), T4r, -ms, &(Rm[WS(rs, 1)]));
			 T22 = VCONJ(VSUB(T20, T1Y));
			 ST(&(Rm[WS(rs, 10)]), T22, -ms, &(Rm[0]));
			 T4j = VCONJ(VSUB(T4h, T4f));
			 ST(&(Rm[WS(rs, 11)]), T4j, -ms, &(Rm[WS(rs, 1)]));
			 T2c = VCONJ(VSUB(T2a, T26));
			 ST(&(Rm[WS(rs, 2)]), T2c, -ms, &(Rm[0]));
			 T2N = VADD(T2E, T2M);
			 ST(&(Rp[WS(rs, 12)]), T2N, ms, &(Rp[0]));
			 T40 = VADD(T3X, T3Z);
			 ST(&(Rp[WS(rs, 1)]), T40, ms, &(Rp[WS(rs, 1)]));
			 T3g = VADD(T3e, T3f);
			 ST(&(Rp[0]), T3g, ms, &(Rp[0]));
			 T3G = VADD(T3j, T3F);
			 ST(&(Rp[WS(rs, 7)]), T3G, ms, &(Rp[WS(rs, 1)]));
			 T41 = VCONJ(VSUB(T3Z, T3X));
			 ST(&(Rm[WS(rs, 1)]), T41, -ms, &(Rm[WS(rs, 1)]));
			 T2O = VCONJ(VSUB(T2M, T2E));
			 ST(&(Rm[WS(rs, 12)]), T2O, -ms, &(Rm[0]));
			 T3H = VCONJ(VSUB(T3F, T3j));
			 ST(&(Rm[WS(rs, 7)]), T3H, -ms, &(Rm[WS(rs, 1)]));
			 T3h = VCONJ(VSUB(T3f, T3e));
			 ST(&(Rm[0]), T3h, -ms, &(Rm[0]));
			 T2V = VADD(T2S, T2U);
			 ST(&(Rp[WS(rs, 4)]), T2V, ms, &(Rp[0]));
			 T3U = VADD(T3P, T3T);
			 ST(&(Rp[WS(rs, 15)]), T3U, ms, &(Rp[WS(rs, 1)]));
			 T39 = VADD(T34, T38);
			 ST(&(Rp[WS(rs, 8)]), T39, ms, &(Rp[0]));
			 T3M = VADD(T3J, T3L);
			 ST(&(Rp[WS(rs, 9)]), T3M, ms, &(Rp[WS(rs, 1)]));
			 T3V = VCONJ(VSUB(T3T, T3P));
			 ST(&(Rm[WS(rs, 15)]), T3V, -ms, &(Rm[WS(rs, 1)]));
			 T2W = VCONJ(VSUB(T2U, T2S));
			 ST(&(Rm[WS(rs, 4)]), T2W, -ms, &(Rm[0]));
			 T3N = VCONJ(VSUB(T3L, T3J));
			 ST(&(Rm[WS(rs, 9)]), T3N, -ms, &(Rm[WS(rs, 1)]));
			 T3a = VCONJ(VSUB(T38, T34));
			 ST(&(Rm[WS(rs, 8)]), T3a, -ms, &(Rm[0]));
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

static const hc2c_desc desc = { 32, XSIMD_STRING("hc2cbdftv_32"), twinstr, &GENUS, { 119, 62, 130, 0 } };

void XSIMD(codelet_hc2cbdftv_32) (planner *p) {
     X(khc2c_register) (p, hc2cbdftv_32, &desc, HC2C_VIA_DFT);
}
#else

/* Generated by: ../../../genfft/gen_hc2cdft_c.native -simd -compact -variables 4 -pipeline-latency 8 -trivial-stores -variables 32 -no-generate-bytw -n 32 -dif -sign 1 -name hc2cbdftv_32 -include rdft/simd/hc2cbv.h */

/*
 * This function contains 249 FP additions, 104 FP multiplications,
 * (or, 233 additions, 88 multiplications, 16 fused multiply/add),
 * 161 stack variables, 7 constants, and 64 memory accesses
 */
#include "rdft/simd/hc2cbv.h"

static void hc2cbdftv_32(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP195090322, +0.195090322016128267848284868477022240927691618);
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DVK(KP555570233, +0.555570233019602224742830813948532874374937191);
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * ((TWVL / VL) * 62)); m < me; m = m + VL, Rp = Rp + (VL * ms), Ip = Ip + (VL * ms), Rm = Rm - (VL * ms), Im = Im - (VL * ms), W = W + (TWVL * 62), MAKE_VOLATILE_STRIDE(128, rs)) {
	       V T1W, T21, Tf, T2c, T1t, T2r, T3T, T4m, Ty, T2q, T3P, T4n, T1n, T2d, T1T;
	       V T22, T1E, T24, T3I, T4p, TU, T2n, T1i, T2h, T1L, T25, T3L, T4q, T1f, T2o;
	       V T1j, T2k;
	       {
		    V T2, T4, T1Z, T1p, T1r, T20, T9, T1U, Td, T1V, T3, T1q, T6, T8, T7;
		    V Tc, Tb, Ta, T5, Te, T1o, T1s, T3R, T3S, Tj, T1N, Tw, T1Q, Tn, T1O;
		    V Ts, T1R, Tg, Ti, Th, Tv, Tu, Tt, Tk, Tm, Tl, Tp, Tr, Tq, To;
		    V Tx, T3N, T3O, T1l, T1m, T1P, T1S;
		    T2 = LD(&(Rp[0]), ms, &(Rp[0]));
		    T3 = LD(&(Rm[WS(rs, 15)]), -ms, &(Rm[WS(rs, 1)]));
		    T4 = VCONJ(T3);
		    T1Z = VADD(T2, T4);
		    T1p = LD(&(Rp[WS(rs, 8)]), ms, &(Rp[0]));
		    T1q = LD(&(Rm[WS(rs, 7)]), -ms, &(Rm[WS(rs, 1)]));
		    T1r = VCONJ(T1q);
		    T20 = VADD(T1p, T1r);
		    T6 = LD(&(Rp[WS(rs, 4)]), ms, &(Rp[0]));
		    T7 = LD(&(Rm[WS(rs, 11)]), -ms, &(Rm[WS(rs, 1)]));
		    T8 = VCONJ(T7);
		    T9 = VSUB(T6, T8);
		    T1U = VADD(T6, T8);
		    Tc = LD(&(Rp[WS(rs, 12)]), ms, &(Rp[0]));
		    Ta = LD(&(Rm[WS(rs, 3)]), -ms, &(Rm[WS(rs, 1)]));
		    Tb = VCONJ(Ta);
		    Td = VSUB(Tb, Tc);
		    T1V = VADD(Tb, Tc);
		    T1W = VSUB(T1U, T1V);
		    T21 = VSUB(T1Z, T20);
		    T5 = VSUB(T2, T4);
		    Te = VMUL(LDK(KP707106781), VADD(T9, Td));
		    Tf = VSUB(T5, Te);
		    T2c = VADD(T5, Te);
		    T1o = VMUL(LDK(KP707106781), VSUB(T9, Td));
		    T1s = VSUB(T1p, T1r);
		    T1t = VSUB(T1o, T1s);
		    T2r = VADD(T1s, T1o);
		    T3R = VADD(T1Z, T20);
		    T3S = VADD(T1U, T1V);
		    T3T = VSUB(T3R, T3S);
		    T4m = VADD(T3R, T3S);
		    Tg = LD(&(Rp[WS(rs, 2)]), ms, &(Rp[0]));
		    Th = LD(&(Rm[WS(rs, 13)]), -ms, &(Rm[WS(rs, 1)]));
		    Ti = VCONJ(Th);
		    Tj = VSUB(Tg, Ti);
		    T1N = VADD(Tg, Ti);
		    Tv = LD(&(Rp[WS(rs, 14)]), ms, &(Rp[0]));
		    Tt = LD(&(Rm[WS(rs, 1)]), -ms, &(Rm[WS(rs, 1)]));
		    Tu = VCONJ(Tt);
		    Tw = VSUB(Tu, Tv);
		    T1Q = VADD(Tu, Tv);
		    Tk = LD(&(Rp[WS(rs, 10)]), ms, &(Rp[0]));
		    Tl = LD(&(Rm[WS(rs, 5)]), -ms, &(Rm[WS(rs, 1)]));
		    Tm = VCONJ(Tl);
		    Tn = VSUB(Tk, Tm);
		    T1O = VADD(Tk, Tm);
		    Tp = LD(&(Rp[WS(rs, 6)]), ms, &(Rp[0]));
		    Tq = LD(&(Rm[WS(rs, 9)]), -ms, &(Rm[WS(rs, 1)]));
		    Tr = VCONJ(Tq);
		    Ts = VSUB(Tp, Tr);
		    T1R = VADD(Tp, Tr);
		    To = VFMA(LDK(KP382683432), Tj, VMUL(LDK(KP923879532), Tn));
		    Tx = VFNMS(LDK(KP382683432), Tw, VMUL(LDK(KP923879532), Ts));
		    Ty = VSUB(To, Tx);
		    T2q = VADD(To, Tx);
		    T3N = VADD(T1N, T1O);
		    T3O = VADD(T1Q, T1R);
		    T3P = VSUB(T3N, T3O);
		    T4n = VADD(T3N, T3O);
		    T1l = VFNMS(LDK(KP382683432), Tn, VMUL(LDK(KP923879532), Tj));
		    T1m = VFMA(LDK(KP923879532), Tw, VMUL(LDK(KP382683432), Ts));
		    T1n = VSUB(T1l, T1m);
		    T2d = VADD(T1l, T1m);
		    T1P = VSUB(T1N, T1O);
		    T1S = VSUB(T1Q, T1R);
		    T1T = VMUL(LDK(KP707106781), VSUB(T1P, T1S));
		    T22 = VMUL(LDK(KP707106781), VADD(T1P, T1S));
	       }
	       {
		    V TD, T1B, TR, T1y, TH, T1C, TM, T1z, TA, TC, TB, TO, TQ, TP, TG;
		    V TF, TE, TJ, TL, TK, T1A, T1D, T3G, T3H, TN, T2f, TT, T2g, TI, TS;
		    V TY, T1I, T1c, T1F, T12, T1J, T17, T1G, TV, TX, TW, T1b, T1a, T19, T11;
		    V T10, TZ, T14, T16, T15, T1H, T1K, T3J, T3K, T18, T2i, T1e, T2j, T13, T1d;
		    TA = LD(&(Rp[WS(rs, 5)]), ms, &(Rp[WS(rs, 1)]));
		    TB = LD(&(Rm[WS(rs, 10)]), -ms, &(Rm[0]));
		    TC = VCONJ(TB);
		    TD = VSUB(TA, TC);
		    T1B = VADD(TA, TC);
		    TO = LD(&(Rp[WS(rs, 1)]), ms, &(Rp[WS(rs, 1)]));
		    TP = LD(&(Rm[WS(rs, 14)]), -ms, &(Rm[0]));
		    TQ = VCONJ(TP);
		    TR = VSUB(TO, TQ);
		    T1y = VADD(TO, TQ);
		    TG = LD(&(Rp[WS(rs, 13)]), ms, &(Rp[WS(rs, 1)]));
		    TE = LD(&(Rm[WS(rs, 2)]), -ms, &(Rm[0]));
		    TF = VCONJ(TE);
		    TH = VSUB(TF, TG);
		    T1C = VADD(TF, TG);
		    TJ = LD(&(Rp[WS(rs, 9)]), ms, &(Rp[WS(rs, 1)]));
		    TK = LD(&(Rm[WS(rs, 6)]), -ms, &(Rm[0]));
		    TL = VCONJ(TK);
		    TM = VSUB(TJ, TL);
		    T1z = VADD(TJ, TL);
		    T1A = VSUB(T1y, T1z);
		    T1D = VSUB(T1B, T1C);
		    T1E = VFNMS(LDK(KP382683432), T1D, VMUL(LDK(KP923879532), T1A));
		    T24 = VFMA(LDK(KP382683432), T1A, VMUL(LDK(KP923879532), T1D));
		    T3G = VADD(T1y, T1z);
		    T3H = VADD(T1B, T1C);
		    T3I = VSUB(T3G, T3H);
		    T4p = VADD(T3G, T3H);
		    TI = VMUL(LDK(KP707106781), VSUB(TD, TH));
		    TN = VSUB(TI, TM);
		    T2f = VADD(TM, TI);
		    TS = VMUL(LDK(KP707106781), VADD(TD, TH));
		    TT = VSUB(TR, TS);
		    T2g = VADD(TR, TS);
		    TU = VFMA(LDK(KP831469612), TN, VMUL(LDK(KP555570233), TT));
		    T2n = VFNMS(LDK(KP195090322), T2f, VMUL(LDK(KP980785280), T2g));
		    T1i = VFNMS(LDK(KP555570233), TN, VMUL(LDK(KP831469612), TT));
		    T2h = VFMA(LDK(KP980785280), T2f, VMUL(LDK(KP195090322), T2g));
		    TV = LD(&(Rp[WS(rs, 3)]), ms, &(Rp[WS(rs, 1)]));
		    TW = LD(&(Rm[WS(rs, 12)]), -ms, &(Rm[0]));
		    TX = VCONJ(TW);
		    TY = VSUB(TV, TX);
		    T1I = VADD(TV, TX);
		    T1b = LD(&(Rp[WS(rs, 15)]), ms, &(Rp[WS(rs, 1)]));
		    T19 = LD(&(Rm[0]), -ms, &(Rm[0]));
		    T1a = VCONJ(T19);
		    T1c = VSUB(T1a, T1b);
		    T1F = VADD(T1a, T1b);
		    T11 = LD(&(Rp[WS(rs, 11)]), ms, &(Rp[WS(rs, 1)]));
		    TZ = LD(&(Rm[WS(rs, 4)]), -ms, &(Rm[0]));
		    T10 = VCONJ(TZ);
		    T12 = VSUB(T10, T11);
		    T1J = VADD(T10, T11);
		    T14 = LD(&(Rp[WS(rs, 7)]), ms, &(Rp[WS(rs, 1)]));
		    T15 = LD(&(Rm[WS(rs, 8)]), -ms, &(Rm[0]));
		    T16 = VCONJ(T15);
		    T17 = VSUB(T14, T16);
		    T1G = VADD(T14, T16);
		    T1H = VSUB(T1F, T1G);
		    T1K = VSUB(T1I, T1J);
		    T1L = VFMA(LDK(KP923879532), T1H, VMUL(LDK(KP382683432), T1K));
		    T25 = VFNMS(LDK(KP382683432), T1H, VMUL(LDK(KP923879532), T1K));
		    T3J = VADD(T1F, T1G);
		    T3K = VADD(T1I, T1J);
		    T3L = VSUB(T3J, T3K);
		    T4q = VADD(T3J, T3K);
		    T13 = VMUL(LDK(KP707106781), VSUB(TY, T12));
		    T18 = VSUB(T13, T17);
		    T2i = VADD(T17, T13);
		    T1d = VMUL(LDK(KP707106781), VADD(TY, T12));
		    T1e = VSUB(T1c, T1d);
		    T2j = VADD(T1c, T1d);
		    T1f = VFNMS(LDK(KP555570233), T1e, VMUL(LDK(KP831469612), T18));
		    T2o = VFMA(LDK(KP195090322), T2i, VMUL(LDK(KP980785280), T2j));
		    T1j = VFMA(LDK(KP555570233), T18, VMUL(LDK(KP831469612), T1e));
		    T2k = VFNMS(LDK(KP195090322), T2j, VMUL(LDK(KP980785280), T2i));
	       }
	       {
		    V T4L, T4G, T4s, T4y, T3W, T4g, T42, T4a, T3g, T4e, T3o, T3E, T1w, T46, T2M;
		    V T40, T2u, T4w, T2C, T4k, T36, T3A, T3i, T3s, T28, T2O, T2w, T2G, T2Y, T4K;
		    V T3y, T4C;
		    {
			 V T4E, T4F, T4D, T4o, T4r, T4l, T4x, T3Q, T48, T3V, T49, T3M, T3U, T3F, T4f;
			 V T41, T47, T3c, T3n, T3f, T3m, T3a, T3b, T3d, T3e, T39, T4d, T3l, T3D, T1h;
			 V T2K, T1v, T2L, Tz, T1g, T1k, T1u, T1, T45, T2J, T3Z, T2m, T2A, T2t, T2B;
			 V T2e, T2l, T2p, T2s, T2b, T4v, T2z, T4j;
			 T4E = VADD(T4m, T4n);
			 T4F = VADD(T4p, T4q);
			 T4L = VADD(T4E, T4F);
			 T4D = LDW(&(W[TWVL * 30]));
			 T4G = VZMUL(T4D, VSUB(T4E, T4F));
			 T4o = VSUB(T4m, T4n);
			 T4r = VBYI(VSUB(T4p, T4q));
			 T4l = LDW(&(W[TWVL * 46]));
			 T4s = VZMUL(T4l, VSUB(T4o, T4r));
			 T4x = LDW(&(W[TWVL * 14]));
			 T4y = VZMUL(T4x, VADD(T4o, T4r));
			 T3M = VMUL(LDK(KP707106781), VSUB(T3I, T3L));
			 T3Q = VBYI(VSUB(T3M, T3P));
			 T48 = VBYI(VADD(T3P, T3M));
			 T3U = VMUL(LDK(KP707106781), VADD(T3I, T3L));
			 T3V = VSUB(T3T, T3U);
			 T49 = VADD(T3T, T3U);
			 T3F = LDW(&(W[TWVL * 22]));
			 T3W = VZMUL(T3F, VADD(T3Q, T3V));
			 T4f = LDW(&(W[TWVL * 54]));
			 T4g = VZMUL(T4f, VSUB(T49, T48));
			 T41 = LDW(&(W[TWVL * 38]));
			 T42 = VZMUL(T41, VSUB(T3V, T3Q));
			 T47 = LDW(&(W[TWVL * 6]));
			 T4a = VZMUL(T47, VADD(T48, T49));
			 T3a = VADD(T1t, T1n);
			 T3b = VADD(TU, T1f);
			 T3c = VBYI(VADD(T3a, T3b));
			 T3n = VBYI(VSUB(T3b, T3a));
			 T3d = VADD(Tf, Ty);
			 T3e = VADD(T1i, T1j);
			 T3f = VADD(T3d, T3e);
			 T3m = VSUB(T3d, T3e);
			 T39 = LDW(&(W[TWVL * 4]));
			 T3g = VZMULI(T39, VADD(T3c, T3f));
			 T4d = LDW(&(W[TWVL * 56]));
			 T4e = VZMULI(T4d, VSUB(T3f, T3c));
			 T3l = LDW(&(W[TWVL * 36]));
			 T3o = VZMULI(T3l, VSUB(T3m, T3n));
			 T3D = LDW(&(W[TWVL * 24]));
			 T3E = VZMULI(T3D, VADD(T3n, T3m));
			 Tz = VSUB(Tf, Ty);
			 T1g = VSUB(TU, T1f);
			 T1h = VSUB(Tz, T1g);
			 T2K = VADD(Tz, T1g);
			 T1k = VSUB(T1i, T1j);
			 T1u = VSUB(T1n, T1t);
			 T1v = VBYI(VSUB(T1k, T1u));
			 T2L = VBYI(VADD(T1u, T1k));
			 T1 = LDW(&(W[TWVL * 20]));
			 T1w = VZMULI(T1, VADD(T1h, T1v));
			 T45 = LDW(&(W[TWVL * 8]));
			 T46 = VZMULI(T45, VADD(T2K, T2L));
			 T2J = LDW(&(W[TWVL * 52]));
			 T2M = VZMULI(T2J, VSUB(T2K, T2L));
			 T3Z = LDW(&(W[TWVL * 40]));
			 T40 = VZMULI(T3Z, VSUB(T1h, T1v));
			 T2e = VSUB(T2c, T2d);
			 T2l = VSUB(T2h, T2k);
			 T2m = VSUB(T2e, T2l);
			 T2A = VADD(T2e, T2l);
			 T2p = VSUB(T2n, T2o);
			 T2s = VSUB(T2q, T2r);
			 T2t = VBYI(VSUB(T2p, T2s));
			 T2B = VBYI(VADD(T2s, T2p));
			 T2b = LDW(&(W[TWVL * 44]));
			 T2u = VZMULI(T2b, VSUB(T2m, T2t));
			 T4v = LDW(&(W[TWVL * 16]));
			 T4w = VZMULI(T4v, VADD(T2m, T2t));
			 T2z = LDW(&(W[TWVL * 12]));
			 T2C = VZMULI(T2z, VADD(T2A, T2B));
			 T4j = LDW(&(W[TWVL * 48]));
			 T4k = VZMULI(T4j, VSUB(T2A, T2B));
			 {
			      V T32, T3q, T35, T3r, T30, T31, T33, T34, T2Z, T3z, T3h, T3p, T1Y, T2E, T27;
			      V T2F, T1M, T1X, T23, T26, T1x, T2N, T2v, T2D, T2U, T3x, T2X, T3w, T2S, T2T;
			      V T2V, T2W, T2R, T4J, T3v, T4B;
			      T30 = VADD(T21, T22);
			      T31 = VADD(T1E, T1L);
			      T32 = VADD(T30, T31);
			      T3q = VSUB(T30, T31);
			      T33 = VADD(T1W, T1T);
			      T34 = VADD(T24, T25);
			      T35 = VBYI(VADD(T33, T34));
			      T3r = VBYI(VSUB(T34, T33));
			      T2Z = LDW(&(W[TWVL * 58]));
			      T36 = VZMUL(T2Z, VSUB(T32, T35));
			      T3z = LDW(&(W[TWVL * 26]));
			      T3A = VZMUL(T3z, VADD(T3q, T3r));
			      T3h = LDW(&(W[TWVL * 2]));
			      T3i = VZMUL(T3h, VADD(T32, T35));
			      T3p = LDW(&(W[TWVL * 34]));
			      T3s = VZMUL(T3p, VSUB(T3q, T3r));
			      T1M = VSUB(T1E, T1L);
			      T1X = VSUB(T1T, T1W);
			      T1Y = VBYI(VSUB(T1M, T1X));
			      T2E = VBYI(VADD(T1X, T1M));
			      T23 = VSUB(T21, T22);
			      T26 = VSUB(T24, T25);
			      T27 = VSUB(T23, T26);
			      T2F = VADD(T23, T26);
			      T1x = LDW(&(W[TWVL * 18]));
			      T28 = VZMUL(T1x, VADD(T1Y, T27));
			      T2N = LDW(&(W[TWVL * 50]));
			      T2O = VZMUL(T2N, VSUB(T2F, T2E));
			      T2v = LDW(&(W[TWVL * 42]));
			      T2w = VZMUL(T2v, VSUB(T27, T1Y));
			      T2D = LDW(&(W[TWVL * 10]));
			      T2G = VZMUL(T2D, VADD(T2E, T2F));
			      T2S = VADD(T2c, T2d);
			      T2T = VADD(T2n, T2o);
			      T2U = VADD(T2S, T2T);
			      T3x = VSUB(T2S, T2T);
			      T2V = VADD(T2r, T2q);
			      T2W = VADD(T2h, T2k);
			      T2X = VBYI(VADD(T2V, T2W));
			      T3w = VBYI(VSUB(T2W, T2V));
			      T2R = LDW(&(W[TWVL * 60]));
			      T2Y = VZMULI(T2R, VSUB(T2U, T2X));
			      T4J = LDW(&(W[0]));
			      T4K = VZMULI(T4J, VADD(T2X, T2U));
			      T3v = LDW(&(W[TWVL * 28]));
			      T3y = VZMULI(T3v, VADD(T3w, T3x));
			      T4B = LDW(&(W[TWVL * 32]));
			      T4C = VZMULI(T4B, VSUB(T3x, T3w));
			 }
		    }
		    {
			 V T29, T4M, T2P, T4t, T4N, T2a, T4u, T2Q, T2x, T4H, T2H, T4z, T4I, T2y, T4A;
			 V T2I, T37, T4h, T3B, T3X, T4i, T38, T3Y, T3C, T3j, T4b, T3t, T43, T4c, T3k;
			 V T44, T3u;
			 T29 = VADD(T1w, T28);
			 ST(&(Rp[WS(rs, 5)]), T29, ms, &(Rp[WS(rs, 1)]));
			 T4M = VADD(T4K, T4L);
			 ST(&(Rp[0]), T4M, ms, &(Rp[0]));
			 T2P = VADD(T2M, T2O);
			 ST(&(Rp[WS(rs, 13)]), T2P, ms, &(Rp[WS(rs, 1)]));
			 T4t = VADD(T4k, T4s);
			 ST(&(Rp[WS(rs, 12)]), T4t, ms, &(Rp[0]));
			 T4N = VCONJ(VSUB(T4L, T4K));
			 ST(&(Rm[0]), T4N, -ms, &(Rm[0]));
			 T2a = VCONJ(VSUB(T28, T1w));
			 ST(&(Rm[WS(rs, 5)]), T2a, -ms, &(Rm[WS(rs, 1)]));
			 T4u = VCONJ(VSUB(T4s, T4k));
			 ST(&(Rm[WS(rs, 12)]), T4u, -ms, &(Rm[0]));
			 T2Q = VCONJ(VSUB(T2O, T2M));
			 ST(&(Rm[WS(rs, 13)]), T2Q, -ms, &(Rm[WS(rs, 1)]));
			 T2x = VADD(T2u, T2w);
			 ST(&(Rp[WS(rs, 11)]), T2x, ms, &(Rp[WS(rs, 1)]));
			 T4H = VADD(T4C, T4G);
			 ST(&(Rp[WS(rs, 8)]), T4H, ms, &(Rp[0]));
			 T2H = VADD(T2C, T2G);
			 ST(&(Rp[WS(rs, 3)]), T2H, ms, &(Rp[WS(rs, 1)]));
			 T4z = VADD(T4w, T4y);
			 ST(&(Rp[WS(rs, 4)]), T4z, ms, &(Rp[0]));
			 T4I = VCONJ(VSUB(T4G, T4C));
			 ST(&(Rm[WS(rs, 8)]), T4I, -ms, &(Rm[0]));
			 T2y = VCONJ(VSUB(T2w, T2u));
			 ST(&(Rm[WS(rs, 11)]), T2y, -ms, &(Rm[WS(rs, 1)]));
			 T4A = VCONJ(VSUB(T4y, T4w));
			 ST(&(Rm[WS(rs, 4)]), T4A, -ms, &(Rm[0]));
			 T2I = VCONJ(VSUB(T2G, T2C));
			 ST(&(Rm[WS(rs, 3)]), T2I, -ms, &(Rm[WS(rs, 1)]));
			 T37 = VADD(T2Y, T36);
			 ST(&(Rp[WS(rs, 15)]), T37, ms, &(Rp[WS(rs, 1)]));
			 T4h = VADD(T4e, T4g);
			 ST(&(Rp[WS(rs, 14)]), T4h, ms, &(Rp[0]));
			 T3B = VADD(T3y, T3A);
			 ST(&(Rp[WS(rs, 7)]), T3B, ms, &(Rp[WS(rs, 1)]));
			 T3X = VADD(T3E, T3W);
			 ST(&(Rp[WS(rs, 6)]), T3X, ms, &(Rp[0]));
			 T4i = VCONJ(VSUB(T4g, T4e));
			 ST(&(Rm[WS(rs, 14)]), T4i, -ms, &(Rm[0]));
			 T38 = VCONJ(VSUB(T36, T2Y));
			 ST(&(Rm[WS(rs, 15)]), T38, -ms, &(Rm[WS(rs, 1)]));
			 T3Y = VCONJ(VSUB(T3W, T3E));
			 ST(&(Rm[WS(rs, 6)]), T3Y, -ms, &(Rm[0]));
			 T3C = VCONJ(VSUB(T3A, T3y));
			 ST(&(Rm[WS(rs, 7)]), T3C, -ms, &(Rm[WS(rs, 1)]));
			 T3j = VADD(T3g, T3i);
			 ST(&(Rp[WS(rs, 1)]), T3j, ms, &(Rp[WS(rs, 1)]));
			 T4b = VADD(T46, T4a);
			 ST(&(Rp[WS(rs, 2)]), T4b, ms, &(Rp[0]));
			 T3t = VADD(T3o, T3s);
			 ST(&(Rp[WS(rs, 9)]), T3t, ms, &(Rp[WS(rs, 1)]));
			 T43 = VADD(T40, T42);
			 ST(&(Rp[WS(rs, 10)]), T43, ms, &(Rp[0]));
			 T4c = VCONJ(VSUB(T4a, T46));
			 ST(&(Rm[WS(rs, 2)]), T4c, -ms, &(Rm[0]));
			 T3k = VCONJ(VSUB(T3i, T3g));
			 ST(&(Rm[WS(rs, 1)]), T3k, -ms, &(Rm[WS(rs, 1)]));
			 T44 = VCONJ(VSUB(T42, T40));
			 ST(&(Rm[WS(rs, 10)]), T44, -ms, &(Rm[0]));
			 T3u = VCONJ(VSUB(T3s, T3o));
			 ST(&(Rm[WS(rs, 9)]), T3u, -ms, &(Rm[WS(rs, 1)]));
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

static const hc2c_desc desc = { 32, XSIMD_STRING("hc2cbdftv_32"), twinstr, &GENUS, { 233, 88, 16, 0 } };

void XSIMD(codelet_hc2cbdftv_32) (planner *p) {
     X(khc2c_register) (p, hc2cbdftv_32, &desc, HC2C_VIA_DFT);
}
#endif
