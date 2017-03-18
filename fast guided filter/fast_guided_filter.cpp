#include "fast_guided_filter.h"

#include <stdlib.h>
#include <memory.h>

//****************************************************************
//resample sse
#if use_sse


#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <emmintrin.h> // SSE2:<e*.h>, SSE3:<p*.h>, SSE4:<s*.h>

#define RETf inline __m128
#define RETi inline __m128i

// set, load and store values
RETf SET(const float &x) { return _mm_set1_ps(x); }
RETf SET(float x, float y, float z, float w) { return _mm_set_ps(x, y, z, w); }
RETi SET(const int &x) { return _mm_set1_epi32(x); }
RETf LD(const float &x) { return _mm_load_ps(&x); }
RETf LDu(const float &x) { return _mm_loadu_ps(&x); }
RETf STR(float &x, const __m128 y) { _mm_store_ps(&x, y); return y; }
RETf STR1(float &x, const __m128 y) { _mm_store_ss(&x, y); return y; }
RETf STRu(float &x, const __m128 y) { _mm_storeu_ps(&x, y); return y; }
RETf STR(float &x, const float y) { return STR(x, SET(y)); }

// arithmetic operators
RETi ADD(const __m128i x, const __m128i y) { return _mm_add_epi32(x, y); }
RETf ADD(const __m128 x, const __m128 y) { return _mm_add_ps(x, y); }
RETf ADD(const __m128 x, const __m128 y, const __m128 z) {
	return ADD(ADD(x, y), z);
}
RETf ADD(const __m128 a, const __m128 b, const __m128 c, const __m128 &d) {
	return ADD(ADD(ADD(a, b), c), d);
}
RETf SUB(const __m128 x, const __m128 y) { return _mm_sub_ps(x, y); }
RETf MUL(const __m128 x, const __m128 y) { return _mm_mul_ps(x, y); }
RETf MUL(const __m128 x, const float y) { return MUL(x, SET(y)); }
RETf MUL(const float x, const __m128 y) { return MUL(SET(x), y); }
RETf INC(__m128 &x, const __m128 y) { return x = ADD(x, y); }
RETf INC(float &x, const __m128 y) { __m128 t = ADD(LD(x), y); return STR(x, t); }
RETf DEC(__m128 &x, const __m128 y) { return x = SUB(x, y); }
RETf DEC(float &x, const __m128 y) { __m128 t = SUB(LD(x), y); return STR(x, t); }
RETf MIN(const __m128 x, const __m128 y) { return _mm_min_ps(x, y); }
RETf RCP(const __m128 x) { return _mm_rcp_ps(x); }
RETf RCPSQRT(const __m128 x) { return _mm_rsqrt_ps(x); }

// logical operators
RETf AND(const __m128 x, const __m128 y) { return _mm_and_ps(x, y); }
RETi AND(const __m128i x, const __m128i y) { return _mm_and_si128(x, y); }
RETf ANDNOT(const __m128 x, const __m128 y) { return _mm_andnot_ps(x, y); }
RETf OR(const __m128 x, const __m128 y) { return _mm_or_ps(x, y); }
RETf XOR(const __m128 x, const __m128 y) { return _mm_xor_ps(x, y); }

// comparison operators
RETf CMPGT(const __m128 x, const __m128 y) { return _mm_cmpgt_ps(x, y); }
RETf CMPLT(const __m128 x, const __m128 y) { return _mm_cmplt_ps(x, y); }
RETi CMPGT(const __m128i x, const __m128i y) { return _mm_cmpgt_epi32(x, y); }
RETi CMPLT(const __m128i x, const __m128i y) { return _mm_cmplt_epi32(x, y); }

// conversion operators
RETf CVT(const __m128i x) { return _mm_cvtepi32_ps(x); }
RETi CVT(const __m128 x) { return _mm_cvttps_epi32(x); }

#undef RETf
#undef RETi



inline void* wrMalloc(size_t size) { return malloc(size); }
inline void wrFree(void * ptr) { free(ptr); }


// platform independent aligned memory allocation (see also alFree)
static void* alMalloc(size_t size, int alignment) {
	const size_t pSize = sizeof(void*), a = alignment - 1;
	void *raw = wrMalloc(size + a + pSize);
	void *aligned = (void*)(((size_t)raw + pSize + a) & ~a);
	*(void**)((size_t)aligned - pSize) = raw;
	return aligned;
}


// platform independent alignned memory de-allocation (see also alMalloc)
static void alFree(void* aligned) {
	void* raw = *(void**)((char*)aligned - sizeof(void*));
	wrFree(raw);
}


// compute interpolation values for single column for resapling
static void resampleCoef(int ha, int hb, int &n, int *&yas,
	int *&ybs, float *&wts, int bd[2], int pad = 0)
{
	const float s = float(hb) / float(ha), sInv = 1 / s; float wt, wt0 = float(1e-3)*s;
	bool ds = ha > hb; int nMax; bd[0] = bd[1] = 0;
	if (ds) { n = 0; nMax = ha + (pad > 2 ? pad : 2)*hb; }
	else { n = nMax = hb; }
	// initialize memory
	wts = (float*)alMalloc(nMax*sizeof(float), 16);
	yas = (int*)alMalloc(nMax*sizeof(int), 16);
	ybs = (int*)alMalloc(nMax*sizeof(int), 16);
	if (ds) for (int yb = 0; yb<hb; yb++) {
		// create coefficients for downsampling
		float ya0f = yb*sInv, ya1f = ya0f + sInv, W = 0;
		int ya0 = int(ceil(ya0f)), ya1 = int(ya1f), n1 = 0;
		for (int ya = ya0 - 1; ya<ya1 + 1; ya++) {
			wt = s; if (ya == ya0 - 1) wt = (ya0 - ya0f)*s; else if (ya == ya1) wt = (ya1f - ya1)*s;
			if (wt>wt0 && ya >= 0) { ybs[n] = yb; yas[n] = ya; wts[n] = wt; n++; n1++; W += wt; }
		}
		if (W>1) for (int i = 0; i<n1; i++) wts[n - n1 + i] /= W;
		if (n1>bd[0]) bd[0] = n1;
		while (n1 < pad) { ybs[n] = yb; yas[n] = yas[n - 1]; wts[n] = 0; n++; n1++; }
	}
	else for (int yb = 0; yb < hb; yb++) {
		// create coefficients for upsampling
		float yaf = (float(.5) + yb)*sInv - float(.5); int ya = (int)floor(yaf);
		wt = 1; if (ya >= 0 && ya < ha - 1) wt = 1 - (yaf - ya);
		if (ya < 0) { ya = 0; bd[0]++; } if (ya >= ha - 1) { ya = ha - 1; bd[1]++; }
		ybs[yb] = yb; yas[yb] = ya; wts[yb] = wt;
	}
}

// resample A using bilinear interpolation and and store result in B
//static void _resample_float(float *A, float *B, int ha, int hb, int wa, int wb, int d, float r)
static void _resample_float(float *A, float *B, int ha, int wa, int hb, int wb, int d, float r)
{
	int hn, wn, x, x1, y, z, xa, xb, ya; float *A0, *A1, *A2, *A3, *B0, wt, wt1;
	float *C = (float*)alMalloc((ha + 4)*sizeof(float), 16); for (y = ha; y < ha + 4; y++) C[y] = 0;
	bool sse = 1;
	// get coefficients for resampling along w and h
	int *xas, *xbs, *yas, *ybs; float *xwts, *ywts; int xbd[2], ybd[2];
	resampleCoef(wa, wb, wn, xas, xbs, xwts, xbd, 0);
	resampleCoef(ha, hb, hn, yas, ybs, ywts, ybd, 4);
	if (wa == 2 * wb) r /= 2; if (wa == 3 * wb) r /= 3; if (wa == 4 * wb) r /= 4;
	r /= float(1 + 1e-6); for (y = 0; y < hn; y++) ywts[y] *= r;
	// resample each channel in turn
	for (z = 0; z < d; z++) for (x = 0; x < wb; x++) {
		if (x == 0) x1 = 0; xa = xas[x1]; xb = xbs[x1]; wt = xwts[x1]; wt1 = 1 - wt; y = 0;
		A0 = A + z*ha*wa + xa*ha; A1 = A0 + ha, A2 = A1 + ha, A3 = A2 + ha; B0 = B + z*hb*wb + xb*hb;
		// variables for SSE (simple casts to float)
		float *Af0, *Af1, *Af2, *Af3, *Bf0, *Cf, *ywtsf, wtf, wt1f;
		Af0 = (float*)A0; Af1 = (float*)A1; Af2 = (float*)A2; Af3 = (float*)A3;
		Bf0 = (float*)B0; Cf = (float*)C;
		ywtsf = (float*)ywts; wtf = (float)wt; wt1f = (float)wt1;
		// resample along x direction (A -> C)
#define FORs(X) if(sse) for(; y<ha-4; y+=4) STR(Cf[y],X);
#define FORr(X) for(; y<ha; y++) C[y] = X;
		if (wa == 2 * wb) {
			FORs(ADD(LDu(Af0[y]), LDu(Af1[y])));
			FORr(A0[y] + A1[y]); x1 += 2;
		}
		else if (wa == 3 * wb) {
			FORs(ADD(LDu(Af0[y]), LDu(Af1[y]), LDu(Af2[y])));
			FORr(A0[y] + A1[y] + A2[y]); x1 += 3;
		}
		else if (wa == 4 * wb) {
			FORs(ADD(LDu(Af0[y]), LDu(Af1[y]), LDu(Af2[y]), LDu(Af3[y])));
			FORr(A0[y] + A1[y] + A2[y] + A3[y]); x1 += 4;
		}
		else if (wa > wb) {
			int m = 1; while (x1 + m < wn && xb == xbs[x1 + m]) m++; float wtsf[4];
			for (int x0 = 0; x0 < (m < 4 ? m : 4); x0++) wtsf[x0] = float(xwts[x1 + x0]);
#define U(x) MUL( LDu(*(Af ## x + y)), SET(wtsf[x]) )
#define V(x) *(A ## x + y) * xwts[x1+x]
			if (m == 1) { FORs(U(0));                     FORr(V(0)); }
			if (m == 2) { FORs(ADD(U(0), U(1)));           FORr(V(0) + V(1)); }
			if (m == 3) { FORs(ADD(U(0), U(1), U(2)));      FORr(V(0) + V(1) + V(2)); }
			if (m >= 4) { FORs(ADD(U(0), U(1), U(2), U(3))); FORr(V(0) + V(1) + V(2) + V(3)); }
#undef U
#undef V
			for (int x0 = 4; x0 < m; x0++) {
				A1 = A0 + x0*ha; wt1 = xwts[x1 + x0]; Af1 = (float*)A1; wt1f = float(wt1); y = 0;
				FORs(ADD(LD(Cf[y]), MUL(LDu(Af1[y]), SET(wt1f)))); FORr(C[y] + A1[y] * wt1);
			}
			x1 += m;
		}
		else {
			bool xBd = x < xbd[0] || x >= wb - xbd[1]; x1++;
			if (xBd) memcpy(C, A0, ha*sizeof(float));
			if (!xBd) FORs(ADD(MUL(LDu(Af0[y]), SET(wtf)), MUL(LDu(Af1[y]), SET(wt1f))));
			if (!xBd) FORr(A0[y] * wt + A1[y] * wt1);
		}
#undef FORs
#undef FORr
		// resample along y direction (B -> C)
		if (ha == hb * 2) {
			float r2 = r / 2; int k = ((~((size_t)B0) + 1) & 15) / 4; y = 0;
			for (; y < k; y++)  B0[y] = (C[2 * y] + C[2 * y + 1])*r2;
			if (sse) for (; y < hb - 4; y += 4) STR(Bf0[y], MUL((float)r2, _mm_shuffle_ps(ADD(
				LDu(Cf[2 * y]), LDu(Cf[2 * y + 1])), ADD(LDu(Cf[2 * y + 4]), LDu(Cf[2 * y + 5])), 136)));
			for (; y < hb; y++) B0[y] = (C[2 * y] + C[2 * y + 1])*r2;
		}
		else if (ha == hb * 3) {
			for (y = 0; y < hb; y++) B0[y] = (C[3 * y] + C[3 * y + 1] + C[3 * y + 2])*(r / 3);
		}
		else if (ha == hb * 4) {
			for (y = 0; y<hb; y++) B0[y] = (C[4 * y] + C[4 * y + 1] + C[4 * y + 2] + C[4 * y + 3])*(r / 4);
		}
		else if (ha>hb) {
			y = 0;
			//if( sse && ybd[0]<=4 ) for(; y<hb; y++) // Requires SSE4
			//  STR1(Bf0[y],_mm_dp_ps(LDu(Cf[yas[y*4]]),LDu(ywtsf[y*4]),0xF1));
#define U(o) C[ya+o]*ywts[y*4+o]
			if (ybd[0] == 2) for (; y < hb; y++) { ya = yas[y * 4]; B0[y] = U(0) + U(1); }
			if (ybd[0] == 3) for (; y < hb; y++) { ya = yas[y * 4]; B0[y] = U(0) + U(1) + U(2); }
			if (ybd[0] == 4) for (; y<hb; y++) { ya = yas[y * 4]; B0[y] = U(0) + U(1) + U(2) + U(3); }
			if (ybd[0]>4)  for (; y < hn; y++) { B0[ybs[y]] += C[yas[y]] * ywts[y]; }
#undef U
		}
		else {
			for (y = 0; y < ybd[0]; y++) B0[y] = C[yas[y]] * ywts[y];
			for (; y < hb - ybd[1]; y++) B0[y] = C[yas[y]] * ywts[y] + C[yas[y] + 1] * (r - ywts[y]);
			for (; y < hb; y++)        B0[y] = C[yas[y]] * ywts[y];
		}
	}
	alFree(xas); alFree(xbs); alFree(xwts); alFree(C);
	alFree(yas); alFree(ybs); alFree(ywts);
}


/*
@function    edge_drawing_line_detector
@param       [in]      A:						  input,single channel
@param       [in/out]  B:                         output,single channel
@param       [in]      w1:                        height of A
@param       [in]      h1:                        height of A
@param       [in]      w2:                        height of B
@param       [in]      h2:                        height of B
@param       [in]      d:                         channel, try: d = 1
@param       [in]      r:                         factor: B = B .* r, try r = 1.f
@return£º
@brief£º

*/
void resample_float_to_float(float *A, float *B, int w1, int h1, int w2, int h2, int d, float r)
{
	_resample_float(A, B, w1, h1, w2, h2, d, r);
}
#endif

//****************************************************************

/* Define NULL pointer value */
#ifndef NULL
#ifdef __cplusplus
#define NULL    0
#else  /* __cplusplus */
#define NULL    ((void *)0)
#endif  /* __cplusplus */
#endif  /* NULL */

#ifndef MIN
#define MIN(a,b)  ((a) > (b) ? (b) : (a))
#endif

#ifndef MAX
#define MAX(a,b)  ((a) < (b) ? (b) : (a))
#endif


static void _boxfilter(float *src, float *dst,
	int w, int h, int r)
{
	int hei = h, wid = w;
	int i, j, k, off;
	float *p = NULL, *p1 = NULL, *p2 = NULL, *imCum = NULL;

	imCum = (float*)calloc(w*h*sizeof(float),1);

	//cumulative sum over Y axis
	memcpy(imCum, src, w*sizeof(float));

	p = imCum + w, p1 = imCum; p2 = src + w;
	for (k = 0; k < ((h - 1)*w); k++)
	{
		*p = (*p1) + (*p2);
		p++, p1++, p2++;
	}

	//difference over Y axis
	memcpy(dst, imCum+r*w, (r+1)*w*sizeof(float));

	p = dst + (r+1)*w; p1 = imCum+(2*r+1)*w; p2 = imCum;
	for (k = 0; k < (hei-r-r-1)*w; k++)
	{
		*p = (*p1) - (*p2);
		p++, p1++, p2++;
	}

	p = dst + (hei - r)*w; p2 = imCum + (hei - r - r - 1)*w;
	for (i = hei - r; i < hei;i++)
	{
		p1 = imCum + (h - 1)*w;
		for (j = 0; j < w;j++)
		{
			*p = (*p1) - (*p2);
			p++, p1++, p2++;
		}
	}

	//cumulative sum over X axis
	p = imCum; p1 = dst;
	for (i = 0; i < h;i++)
	{
		*p++ = *p1++;
		for (j = 1; j < w;j++)
		{
			*p = (*p1) + *(p-1);
			p++, p1++;
		}
	}

	//difference over Y axis
	for (i = 0; i < h; i++)
	{
		off = i*w;
		for (j = 0; j < r+1; j++)
		{
			dst[off + j] = imCum[off + j + r];
		}
	}

	for (i = 0; i < h; i++)
	{
		off = i*w;
		for (j = r+1; j < wid-r; j++)
		{
			dst[off + j] = imCum[off+j+r] - imCum[off+j-r-1];
		}
	}

	for (i = 0; i < h; i++)
	{
		off = i*w;
		for (j = wid-r; j < wid; j++)
		{
			dst[off + j] = imCum[off + w - 1] - imCum[off + j - r - 1];
		}
	}

	if (imCum)free(imCum);
}

static int boxfilter(float *src, float *dst,
	int w, int h, int radius)
{
	if (src == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || (radius < 1))
		return 1;

	if ((2*radius+1) > MIN(w,h))
		return 1;

	_boxfilter(src, dst, w, h, radius);

	return 0;
}


static void resizeNN_int8u_to_float(unsigned char *src, float *dst, int w1, int h1, int w2, int h2)
{
	int i, j, offy_src, offy_dst;
	float k;
	int *lutx = NULL, *luty = NULL;
	int *px = NULL, *py = NULL;
	unsigned char *psrc = NULL;
	float *pdst = NULL;

	if (w1 != w2 && h1 != h2)
	{
		lutx = (int*)malloc(w2*sizeof(int));
		luty = (int*)malloc(h2*sizeof(int));
		px = lutx, py = luty;

		for (k = 0; k < w2; k += 1)
			*px++ = (int)(k * (w1 - 1) / (w2 - 1) + 0.5f);

		for (k = 0; k < h2; k += 1)
			*py++ = (int)(k * (h1 - 1) / (h2 - 1) + 0.5f);

		for (i = 0; i < h2; i++)
		{
			offy_dst = i*w2;
			offy_src = luty[i] * w1;
			for (j = 0; j < w2; j++)
			{
				dst[offy_dst + j] = src[offy_src + lutx[j]];
			}
		}
	}
	else
	if (w1 == w2 && h1 != h2)
	{
		luty = (int*)malloc(h2*sizeof(int));
		py = luty;

		for (k = 0; k < h2; k += 1)
			*py++ = (int)(k * (h1 - 1) / (h2 - 1) + 0.5f);

		for (i = 0; i < h2; i++)
		{
			offy_dst = i*w2;
			offy_src = luty[i] * w1;
			for (j = 0; j < w2; j++)
			{
				dst[offy_dst + j] = src[offy_src + j];
			}
		}
	}
	else
	if (w1 != w2 && h1 == h2)
	{
		lutx = (int*)malloc(w2*sizeof(int));
		px = lutx;

		for (k = 0; k < w2; k += 1)
			*px++ = (int)(k * (w1 - 1) / (w2 - 1) + 0.5f);

		for (i = 0; i < h2; i++)
		{
			offy_dst = i * w2;
			offy_src = i * w1;
			for (j = 0; j < w2; j++)
			{
				dst[offy_dst + j] = src[offy_src + lutx[j]];
			}
		}
	}
	else // (w1 == w2 && h1 == h2)
	{
		psrc = src;
		pdst = dst;
		for (i = 0; i < w2*h2; i++)
			*pdst++ = *psrc++;
	}

	if (lutx)free(lutx);
	if (luty)free(luty);

}


//four points takes original value, boundary takes nearest neighbor, other takes bilinear 
static void resizeBilinear_float_to_float(float *src, float *dst, int w1, int h1, int w2, int h2)
{
	int i, j, offy_src, offy_dst;
	float k;

	int val;
	int offx1, offx2;
	int x_w, y_n, y_s;
	float dx, dy;
	float fa, fb;
	float _x, _y;

	float *lutx = NULL, *luty = NULL;
	float *px = NULL, *py = NULL;

	if (w1 != w2 && h1 != h2)
	{
		lutx = (float*)malloc((w2)*sizeof(float));
		luty = (float*)malloc((h2)*sizeof(float));
		px = lutx, py = luty;

		for (k = 0; k < w2; k += 1)
			*px++ = (k * (w1 - 1) / (w2 - 1));

		for (k = 0; k < h2; k += 1)
			*py++ = (k * (h1 - 1) / (h2 - 1));


		for (i = 1; i < h2 - 1; i++)
		{
			offy_dst = i*w2;
			_y = luty[i];
			for (j = 1; j < w2 - 1; j++)
			{
				_x = lutx[j];
				x_w = (int)_x;
				y_n = (int)_y;
				y_s = y_n + 1;
				dx = _x - x_w;
				dy = _y - y_n;
				offx1 = y_n*w1 + x_w;
				offx2 = y_s*w1 + x_w;
				fa = dx*(src[offx1 + 1] - src[offx1]) + src[offx1];
				fb = dx*(src[offx2 + 1] - src[offx2]) + src[offx2];
				dst[offy_dst + j] = dy * (fb - fa) + fa;

			}
		}

		////BOUNDARY
		////top and buttom
		offy_dst = (h2 - 1) * w2;
		offy_src = (h1 - 1) * w1;
		for (j = 1; j < w2 - 1; j++)
		{
			val = (int)(lutx[j] + 0.5f);

			dst[j] = src[val];
			dst[offy_dst + j] = src[offy_src + val];
		}

		////left and right
		for (i = 1; i < h2 - 1; i++)
		{
			val = (int)(luty[i] + 0.5f);

			dst[i*w2] = src[val*w1];
			dst[i*w2 + (w2 - 1)] = src[val*w1 + (w1 - 1)];
		}

		////four points
		dst[0] = src[0];
		dst[w2 - 1] = src[w1 - 1];
		dst[(h2 - 1)*w2] = src[(h1 - 1)*w1];
		dst[h2*w2 - 1] = src[h1*w1 - 1];

	}
	else
	if (w1 == w2 && h1 != h2)
	{
		luty = (float*)malloc((h2)*sizeof(float));
		py = luty;

		for (k = 0; k < h2; k += 1)
			*py++ = (k * (h1 - 1) / (h2 - 1));

		for (i = 1; i < h2 - 1; i++)
		{
			offy_dst = i*w2;
			_y = luty[i];
			for (j = 1; j < w2 - 1; j++)
			{
				y_n = (int)_y;
				y_s = y_n + 1;
				dy = _y - y_n;
				offx1 = y_n*w1 + j;
				offx2 = y_s*w1 + j;
				fa = src[offx1];
				dst[offy_dst + j] = dy * (src[offx2] - fa) + fa;

			}
		}

		////BOUNDARY
		////top and buttom
		memcpy(dst + 1, src + 1, (w2 - 2)*sizeof(float));
		memcpy(dst + (h2 - 1)*w2 + 1, src + (h1 - 1)*w1 + 1, (w2 - 2)*sizeof(float));


		////left and right
		for (i = 1; i < h2 - 1; i++)
		{
			val = (int)(luty[i] + 0.5f);

			dst[i*w2] = src[val*w1];
			dst[i*w2 + (w2 - 1)] = src[val*w1 + (w1 - 1)];
		}

		////four points
		dst[0] = src[0];
		dst[w2 - 1] = src[w1 - 1];
		dst[(h2 - 1)*w2] = src[(h1 - 1)*w1];
		dst[h2*w2 - 1] = src[h1*w1 - 1];

	}
	else
	if (w1 != w2 && h1 == h2)
	{
		lutx = (float*)malloc((w2)*sizeof(float));
		px = lutx;

		for (k = 0; k < w2; k += 1)
			*px++ = (k * (w1 - 1) / (w2 - 1));


		for (i = 1; i < h2 - 1; i++)
		{
			offy_dst = i*w2;
			for (j = 1; j < w2 - 1; j++)
			{
				_x = lutx[j];
				x_w = (int)_x;
				dx = _x - x_w;
				offx1 = i*w1 + x_w;
				offx2 = offx1 + w1;
				dst[offy_dst + j] = dx*(src[offx1 + 1] - src[offx1]) + src[offx1];

			}
		}

		////BOUNDARY
		////top and buttom
		offy_dst = (h2 - 1) * w2;
		offy_src = (h1 - 1) * w1;
		for (j = 1; j < w2 - 1; j++)
		{
			val = (int)(lutx[j] + 0.5f);

			dst[j] = src[val];
			dst[offy_dst + j] = src[offy_src + val];
		}

		////left and right
		for (i = 1; i < h2 - 1; i++)
		{
			dst[i*w2] = src[i*w1];
			dst[i*w2 + (w2 - 1)] = src[i*w1 + (w1 - 1)];
		}

		////four points
		dst[0] = src[0];
		dst[w2 - 1] = src[w1 - 1];
		dst[(h2 - 1)*w2] = src[(h1 - 1)*w1];
		dst[h2*w2 - 1] = src[h1*w1 - 1];
	}
	else // (w1 == w2 && h1 == h2)
	{
		memcpy(dst, src, w1*h1*sizeof(float));
	}

	if (lutx)free(lutx);
	if (luty)free(luty);

}

// a = a.*b
static void dot_mul(float *a, unsigned char *b, int w, int h)
{
	int k, n = w*h;
	float *p;
	unsigned char *p1;
	p = a, p1 = b;

	for (k = 0; k < n; k++)
	{
		*p *= (*p1);
		p++, p1++;
	}

}

static void dot_mul(float *a, float *b, float *c, int w, int h)
{
	int k, n = w*h;
	float *p, *p1, *p2;
	p = c, p1 = a, p2 = b;

	for (k = 0; k < n;k++)
	{
		*p = (*p1) * (*p2);
		p++, p1++, p2++;
	}

}


// c = a-b
static void array_sub(float *a, float *b, float *c, int w, int h)
{
	int k, n = w*h;
	float *p, *p1, *p2;
	p = c, p1 = a, p2 = b;

	for (k = 0; k < n; k++)
	{
		*p = (*p1) - (*p2);
		p++, p1++, p2++;
	}

}

//b = a + num
static void array_add_num_inv(float *a, float *b, float num, int w, int h)
{
	int k, n = w*h;
	float *p, *p1;
	p1 = a, p = b;

	for (k = 0; k < n; k++)
	{
		*p = 1.f / ((*p1) + num);
		p++, p1++;
	}

}

// c = a+b
static void array_add(float *a, float *b, unsigned char *c, int w, int h)
{
	int k, n = w*h;
	float *p1, *p2;
	unsigned char *p;
	p1 = a, p2 = b;
	p = c;

	for (k = 0; k < n; k++)
	{
		*p = (unsigned char)((*p1) + (*p2) + 0.5f);
		p++, p1++, p2++;
	}

}




static void  fast_guided_filter_singlechannel(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int r, float rp, float sr)
{
	float invs = 1.f / sr;
	int r_sub = (int)(r * invs + 0.5f);
	int w_sub = (int)(w * invs + 0.5f);
	int h_sub = (int)(h * invs + 0.5f);
	int nlength = w_sub * h_sub;
	int k;

	float *p = NULL, *p1 = NULL, *p2 = NULL;
	float *invN = NULL;
	float *mtmp = NULL, *mtx = NULL;
	float *I_sub = NULL, *p_sub = NULL;
	float *mean_I = NULL, *mean_p = NULL, *mean_Ip = NULL, *cov_Ip = NULL;
	float *mean_a = NULL, *mean_b = NULL;

	I_sub = (float*)malloc(nlength * sizeof(float));
	p_sub = (float*)malloc(nlength * sizeof(float));
	mtmp = (float*)malloc(nlength * sizeof(float));
	mtx = (float*)malloc(nlength * sizeof(float));
	invN = (float*)malloc(nlength * sizeof(float));
	mean_I = (float*)malloc(nlength * sizeof(float));
	mean_p = (float*)malloc(nlength * sizeof(float));
	mean_Ip = (float*)malloc(nlength * sizeof(float));
	cov_Ip = (float*)malloc(nlength * sizeof(float));

	mean_a = (float*)malloc(w*h * sizeof(float));
	mean_b = (float*)malloc(w*h * sizeof(float));

	//downscale 
	resizeNN_int8u_to_float(guidance,I_sub,w,h,w_sub,h_sub);
	resizeNN_int8u_to_float(src, p_sub, w, h, w_sub, h_sub);

	//set one
	p = mtx;
	for (k = 0; k < nlength; k++)
		*p++ = 1;

	boxfilter(mtx, invN, w_sub, h_sub, r_sub);

	p = invN;
	for (k = 0; k < nlength; k++)
	{
		*p = 1.f / *p;
		p++;
	}

	//mean_I
	boxfilter(I_sub, mtx, w_sub, h_sub, r_sub);
	dot_mul(mtx, invN, mean_I, w_sub, h_sub);

	//mean_p
	boxfilter(p_sub, mtx, w_sub, h_sub, r_sub);
	dot_mul(mtx, invN, mean_p, w_sub, h_sub);

	//mean_Ip
	dot_mul(I_sub, p_sub, mtmp, w_sub, h_sub);
	boxfilter(mtmp, mtx, w_sub, h_sub, r_sub);
	dot_mul(mtx, invN, mean_Ip,w_sub, h_sub);

	//cov_Ip
	dot_mul(mean_I, mean_p, mtx, w_sub, h_sub);
	array_sub(mean_Ip, mtx, cov_Ip, w_sub, h_sub);

	//mean_II(mean_Ip)
	dot_mul(I_sub, I_sub, mtmp, w_sub, h_sub);
	boxfilter(mtmp, mtx, w_sub, h_sub, r_sub);
	dot_mul(mtx, invN, mean_Ip, w_sub, h_sub);

	//var_I(I_sub)
	dot_mul(mean_I, mean_I, mtx, w_sub, h_sub);
	array_sub(mean_Ip, mtx, I_sub, w_sub, h_sub);

	//a(p_sub)
	array_add_num_inv(I_sub, mtx, rp, w_sub, h_sub);
	dot_mul(cov_Ip, mtx, p_sub, w_sub, h_sub);

	//b(I_sub)
	dot_mul(p_sub, mean_I, mtx, w_sub, h_sub);
	array_sub(mean_p, mtx, I_sub, w_sub, h_sub);

	//mean_a(mean_Ip)
	boxfilter(p_sub, mtx, w_sub, h_sub, r_sub);
	dot_mul(mtx, invN, mean_Ip, w_sub, h_sub);

	//mean_b(cov_Ip)
	boxfilter(I_sub, mtx, w_sub, h_sub, r_sub);
	dot_mul(mtx, invN, cov_Ip, w_sub, h_sub);

	//upscale
#if use_sse
	resample_float_to_float(mean_Ip, mean_a, w_sub, h_sub, w, h, 1, 1.f);
	resample_float_to_float(cov_Ip, mean_b, w_sub, h_sub, w, h, 1, 1.f);
#else
	resizeBilinear_float_to_float(mean_Ip, mean_a, w_sub, h_sub, w, h);
	resizeBilinear_float_to_float(cov_Ip, mean_b, w_sub, h_sub, w, h);
#endif

	//q
	dot_mul(mean_a, guidance, w,h);
	array_add(mean_a, mean_b, dst, w, h);
	

	if (mtmp)free(mtmp);
	if (mtx)free(mtx);
	if (invN)free(invN);

	if (mean_a)free(mean_a);
	if (mean_b)free(mean_b);

	if (p_sub)free(p_sub);
	if (I_sub)free(I_sub);

	if (cov_Ip)free(cov_Ip);
	if (mean_Ip)free(mean_Ip);
	if (mean_p)free(mean_p);
	if (mean_I)free(mean_I);

}



/*
@function    fast_guided_filter
@param       [in]      src:						  input image,single channel
@param       [in]      guidance:                  guided image,single channel
@param       [in/out]  dst:                       output image,single channel
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      c:						  channel of image, only c = 1
@param       [in]      r:						  local window radius
@param       [in]      rp:			              regularization parameter:eps
@param       [in]	   sr:                        subsampling ratio, sr>1:downscale, 0<sr<1:upscale
@return£º										  0:ok; 1:error
@brief£º

					eg: r = 4, (try sr = r/4 to sr=r),(try rp=0.1^2, 0.2^2, 0.4^2)
					try:(src,guidance,dst,w,h,1,4,0.01,4)
					condition: (MIN(w, h) / sr) > 1
					condition: (int)(r / sr + 0.5f) >= 1
*/
int fast_guided_filter(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int c, int r, float rp, float sr)
{
	int _r;
	float _scale = 255 * 255;	//if regularization, _scale = 1; if no regularization, _scale =  255*255
	float _size;
	if (src == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || (r < 1) || (!(c == 1 || c == 3)))
		return 1;

	_size = 1.f / sr * MIN(w, h);
	
	if (_size < 1.f)
		return 1;

	_r = (int)(r / sr + 0.5f);
	if (_r < 1)
		return 1;

	switch (c)
	{
	case 1:
		fast_guided_filter_singlechannel(src, guidance, dst,
			w, h, r, _scale*rp, sr);
		break;

	case 3:

		break;
	}

	return 0;
}

