// Image interpolation from scattered samples
// Author: Sergey Chaban <sergey.chaban@gmail.com>

#include "img_sys.hpp"
#include "img_interp.hpp"

template<typename T>
inline
bool lu_decomp(T* pMtx, int N, T* pWk /* [N] */, int* pIdx /* [N] */, T* pDetSgn = nullptr) {
	T dsgn = 1;
	T* pScl = pWk;
	for (int i = 0; i < N; ++i) {
		T scl = 0;
		int ri = i * N;
		for (int j = 0; j < N; ++j) {
			T a = pMtx[ri + j];
			if (a < 0) a = -a;
			if (a > scl) scl = a;
		}
		if (scl == 0) {
			if (pDetSgn) {
				*pDetSgn = 0;
			}
			return false;
		}
		pScl[i] = scl;
	}
	for (int i = 0; i < N; ++i) {
		pScl[i] = 1 / pScl[i];
	}
	int imax = 0;
	for (int k = 0; k < N; ++k) {
		int rk = k * N;
		T amax = 0;
		for (int i = k; i < N; ++i) {
			T a = pMtx[(i * N) + k];
			if (a < 0) a = -a;
			a *= pScl[i];
			if (amax <= a) {
				amax = a;
				imax = i;
			}
		}
		if (k != imax) {
			int rm = imax * N;
			for (int j = 0; j < N; ++j) {
				T t = pMtx[rm + j];
				pMtx[rm + j] = pMtx[rk + j];
				pMtx[rk + j] = t;
			}
			dsgn = -dsgn;
			pScl[imax] = pScl[k];
		}
		if (pIdx) {
			pIdx[k] = imax;
		}
		if (pMtx[rk + k] == 0) {
			pMtx[rk + k] = T(1.0e-16);
		}
		for (int i = k + 1; i < N; ++i) {
			int ri = i * N;
			T s = pMtx[ri + k] / pMtx[rk + k];
			pMtx[ri + k] = s;
			T* pDst = &pMtx[ri + k + 1];
			T* pSrc = &pMtx[rk + k + 1];
			for (int j = k + 1; j < N; ++j) {
				*pDst++ -= *pSrc++ * s;
			}
		}
	}
	if (pDetSgn) {
		*pDetSgn = dsgn;
	}
	return true;
}

template<typename T>
inline
void lu_solve(T* pAns, const T* pLU, const T* pRHS, const int* pIdx, int N) {
	if (pAns != pRHS) {
		for (int i = 0; i < N; ++i) {
			pAns[i] = pRHS[i];
		}
	}
	T s = 0;
	int ii = -1;
	for (int i = 0; i < N; ++i) {
		int ri = i * N;
		int idx = pIdx[i];
		s = pAns[idx];
		pAns[idx] = pAns[i];
		if (ii < 0) {
			if (s != 0) {
				ii = i;
			}
		} else {
			for (int j = ii; j < i; ++j) {
				s -= pLU[ri + j] * pAns[j];
			}
		}
		pAns[i] = s;
	}
	for (int i = N; --i >= 0;) {
		int ri = i * N;
		s = pAns[i];
		for (int j = i + 1; j < N; ++j) {
			s -= pLU[ri + j] * pAns[j];
		}
		pAns[i] = s / pLU[ri + i];
	}
}

template<typename T>
inline
bool symm_ldlt_decomp(T* pMtx, int N, T* pDet /* [N] */, int* pIdx /* [N] */) {
	const T zthr = T(1.0e-16);
	const T aprm = T(1 + ::sqrt(17.0)) / 8;
	int eposi = 0;
	int enega = 0;
	pIdx[N - 1] = N - 1;
	int i = 0;
	while (i < N - 1) {
		int ri = i * N;
		int i1 = i + 1;
		int i2 = i + 2;
		T aii = pMtx[ri + i];
		if (aii < 0) aii = -aii;
		pIdx[i] = i;
		T lambda = pMtx[ri + i1];
		int j = i1;
		for (int k = i2; k < N; ++k) {
			T aik = pMtx[ri + k];
			if (aik < 0) aik = -aik;
			if (aik > lambda) {
				j = k;
				lambda = aik;
			}
		}
		int rj = j * N;
		bool flg = true;
		if (aii < aprm*lambda) {
			T sigma = lambda;
			for (int k = i1; k < j; ++k) {
				T akj = pMtx[k*N + j];
				if (akj < 0) akj = -akj;
				if (akj > sigma) sigma = akj;
			}
			for (int k = j + 1; k < N; ++k) {
				T ajk = pMtx[rj + k];
				if (ajk < 0) ajk = -ajk;
				if (ajk > sigma) sigma = ajk;
			}
			if (aii*sigma < lambda) {
				T ajj = pMtx[rj + j];
				if (ajj < 0) ajj = -ajj;
				if (ajj > aprm*sigma) {
					for (int k = j + 1; k < N; ++k) {
						T t = pMtx[rj + k];
						pMtx[rj + k] = pMtx[ri + k];
						pMtx[ri + k] = t;
					}
					for (int k = i1; k < j; ++k) {
						int rk = k*N;
						T t = pMtx[ri + k];
						pMtx[ri + k] = pMtx[rk + j];
						pMtx[rk + j] = t;
					}
					T t = pMtx[ri + i];
					pMtx[ri + i] = pMtx[rj + j];
					pMtx[rj + j] = t;
					pIdx[i] = j;
				} else {
					int ri1 = i1*N;
					if (j > i1) {
						for (int k = j + 1; k < N; ++k) {
							T t = pMtx[rj + k];
							pMtx[rj + k] = pMtx[ri1 + k];
							pMtx[ri1 + k] = t;
						}
						for (int k = i2; k < j; ++k) {
							int rk = k*N;
							T t = pMtx[ri1 + k];
							pMtx[ri1 + k] = pMtx[rk + j];
							pMtx[rk + j] = t;
						}
						T t = pMtx[ri + i];
						pMtx[ri + i] = pMtx[rj + j];
						pMtx[rj + j] = t;
						t = pMtx[ri + j];
						pMtx[ri + j] = pMtx[ri + i1];
						pMtx[ri + i1] = t;
					}
					T t = pMtx[ri + i1];
					T det = pMtx[ri + i]*pMtx[ri1 + i1] - t*t;
					T idet = T(1) / det;
					aii = pMtx[ri + i] * idet;
					T aii1 = pMtx[ri + i1] * idet;
					T ai1i1 = pMtx[ri1 + i1] * idet;
					pIdx[i] = j;
					pIdx[i1] = -1;
					pDet[i] = 1;
					pDet[i1] = det;
					for (int k = i2; k < N; ++k) {
						T u = aii1*pMtx[ri1 + k] - ai1i1*pMtx[ri + k];
						T v = aii1*pMtx[ri + k] - aii*pMtx[ri1 + k];
						int rk = k*N;
						T* pDst = &pMtx[rk + k];
						T* pSrc = &pMtx[ri + k];
						for (int m = k; m < N; ++m) {
							*pDst++ += *pSrc++ * u;
						}
						pDst = &pMtx[rk + k];
						pSrc = &pMtx[ri1 + k];
						for (int m = k; m < N; ++m) {
							*pDst++ += *pSrc++ * v;
						}
						pMtx[ri + k] = u;
						pMtx[ri1 + k] = v;
					}
					i = i2;
					++eposi;
					++enega;
					flg = false;
				}
			}
		}
		if (flg) {
			aii = pMtx[ri + i];
			if (aii < 0) aii = -aii;
			if (aii > zthr) {
				aii = pMtx[ri + i];
				if (aii > 0) {
					++eposi;
				} else {
					++enega;
				}
				pDet[i] = aii;
			}
			T s = -T(1) / aii;
			for (int k = i1; k < N; ++k) {
				int rk = k*N;
				T t = pMtx[ri + k] * s;
				T* pDst = &pMtx[rk + k];
				T* pSrc = &pMtx[ri + k];
				for (int m = k; m < N; ++m) {
					*pDst++ += *pSrc++ * t;
				}
				pMtx[ri + k] = t;
			}
			i = i1;
		}
	}
	if (i == N - 1) {
		int ri = i*N;
		T aii = pMtx[ri + i];
		if (aii < 0) aii = -aii;
		if (aii > zthr) {
			aii = pMtx[ri + i];
			if (aii > 0) {
				++eposi;
			} else {
				++enega;
			}
		}
		pDet[i] = pMtx[i*N + i];
	}
	return (eposi + enega) == N;
}

template<typename T>
inline
void symm_ldlt_solve(T* pAns, const T* pMtx, const T* pRHS, const T* pDet /* [N] */, const int* pIdx /* [N] */, int N) {
	if (pAns != pRHS) {
		for (int i = 0; i < N; ++i) {
			pAns[i] = pRHS[i];
		}
	}
	int i = 0;
	while (i < N - 1) {
		int ri = i*N;
		int idx = pIdx[i];
		T tp = pAns[idx];
		if (pIdx[i + 1] >= 0) {
			pAns[idx] = pAns[i];
			pAns[i] = tp / pMtx[ri + i];
			for (int j = i + 1; j < N; ++j) {
				pAns[j] += pMtx[ri + j] * tp;
			}
			i = i + 1;
		} else {
			int i1 = i + 1;
			int ri1 = i1*N;
			T t = pAns[i];
			pAns[idx] = pAns[i1];
			T idet = T(1) / pDet[i1];
			pAns[i] = (pMtx[ri1 + i1]*t - pMtx[ri + i1]*tp) * idet;
			pAns[i1] = (pMtx[ri + i]*tp - pMtx[ri + i1]*t) * idet;
			for (int j = i + 2; j < N; ++j) {
				pAns[j] += pMtx[ri + j] * t;
			}
			for (int j = i + 2; j < N; ++j) {
				pAns[j] += pMtx[ri1 + j] * tp;
			}
			i = i + 2;
		}
	}
	if (i == N - 1) {
		pAns[i] /= pMtx[i*N + i];
		i = N - 2;
	} else {
		i = N - 3;
	}
	while (i >= 0) {
		int ii = pIdx[i] >= 0 ? i : i - 1;
		for (int j = ii; j <= i; ++j) {
			int rj = j*N;
			T s = pAns[j];
			for (int k = i + 1; k < N; ++k) {
				s += pMtx[rj + k] * pAns[k];
			}
			pAns[j] = s;
		}
		int idx = pIdx[ii];
		T t = pAns[i];
		pAns[i] = pAns[idx];
		pAns[idx] = t;
		i = ii - 1;
	}
}

static float frand01() {
	VAL32 uf;
	uint32_t r0 = ::rand() & 0xFF;
	uint32_t r1 = ::rand() & 0xFF;
	uint32_t r2 = ::rand() & 0xFF;
	uf.f = 1.0f;
	uf.u |= (r0 | (r1 << 8) | ((r2 & 0x7F) << 16));
	uf.f -= 1.0f;
	return uf.f;
}

static bool ck_pnt(IMG_COORD* pPts, int n, float x, float y) {
	for (int i = 0; i < n; ++i) {
		if (pPts[i].x == x && pPts[i].y == y) return false;
	}
	return true;
}

IMG_COORD* scatter_pts_rand(int w, int h, int n, int seed) {
	IMG_COORD* pPts = reinterpret_cast<IMG_COORD*>(mem_alloc(n * sizeof(IMG_COORD), "ImgPts"));
	::srand(seed);
	if (pPts) {
		float sx = float(w - 1);
		float sy = float(h - 1);
		for (int i = 0; i < n; ++i) {
			float x = ::roundf(frand01() * sx);
			float y = ::roundf(frand01() * sy);
			if (i > 0) {
				while (!ck_pnt(pPts, i, x, y)) {
					x = ::roundf(frand01() * sx);
					y = ::roundf(frand01() * sy);
				}
			}
			pPts[i].set(x, y);
		}
	}
	return pPts;
}

RGBA* get_img_samples(const IMAGE* pImg, const IMG_COORD* pPts, int n) {
	RGBA* pSmps = nullptr;
	if (pImg && pPts && n > 0) {
		pSmps = reinterpret_cast<RGBA*>(mem_alloc(n * sizeof(RGBA), "ImgSmps"));
		if (pSmps) {
			for (int i = 0; i < n; ++i) {
				pSmps[i] = pImg->get_pixel(pPts[i]);
			}
		}
	}
	return pSmps;
}

static float* alloc_solution(int N) {
	/* NxN + wgt[N*3] + tmp[N] + idx[N] */
	int64_t size = N*N * sizeof(float);
	size += N * 3 * sizeof(float);
	size += N * sizeof(float);
	size += N * sizeof(int);
	return reinterpret_cast<float*>(mem_alloc(size, "Solution"));
}

static void calc_dist_mtx(float* pWk, const IMG_COORD* pPts, int N) {
	for (int i = 0; i < N; ++i) {
		int ri = i * N;
		for (int j = 0; j < N; ++j) {
			float d = pt_dist(pPts[i], pPts[j]);
			pWk[ri + j] = d;
		}
	}
}

static bool solve_lu(float* pWk, const IMG_COORD* pPts, const RGBA* pSmps, int N) {
	float* pMtx = pWk;
	float* pWgt = pWk + (N*N);
	float* pWgtR = pWgt;
	float* pWgtG = pWgt + N;
	float* pWgtB = pWgt + N*2;
	float* pTmp = pWgt + N*3;
	int* pIdx = reinterpret_cast<int*>(pWgt + N*4);

	if (!lu_decomp(pMtx, N, pTmp, pIdx)) return false;

	/* Mtx * wgt<RGB> = smp<RGB> */
	float* dsts[] = { pWgtR, pWgtG, pWgtB };
	for (int i = 0; i < 3; ++i) {
		float* pDst = dsts[i];
		for (int j = 0; j < N; ++j) {
			pDst[j] = pSmps[j].ch[i];
		}
		lu_solve(pDst, pMtx, pDst, pIdx, N);
	}
	return true;
}

static bool solve_ldlt(float* pWk, const IMG_COORD* pPts, const RGBA* pSmps, int N) {
	float* pMtx = pWk;
	float* pWgt = pWk + (N*N);
	float* pWgtR = pWgt;
	float* pWgtG = pWgt + N;
	float* pWgtB = pWgt + N*2;
	float* pDet = pWgt + N*3;
	int* pIdx = reinterpret_cast<int*>(pWgt + N*4);

	if (!symm_ldlt_decomp(pMtx, N, pDet, pIdx)) return false;

	float* dsts[] = { pWgtR, pWgtG, pWgtB };
	for (int i = 0; i < 3; ++i) {
		float* pDst = dsts[i];
		for (int j = 0; j < N; ++j) {
			pDst[j] = pSmps[j].ch[i];
		}
		symm_ldlt_solve(pDst, pMtx, pDst, pDet, pIdx, N);
	}
	return true;
}

static bool solve(float* pWk, const IMG_COORD* pPts, const RGBA* pSmps, int N, bool useLDLT = true) {
	calc_dist_mtx(pWk, pPts, N);
	bool res = false;
	if (useLDLT) {
		res = solve_ldlt(pWk, pPts, pSmps, N);
	} else {
		res = solve_lu(pWk, pPts, pSmps, N);
	}
	return res;
}

static void interp(IMAGE* pImg, float* pWk, IMG_COORD* pPts, int N) {
	float* pMtx = pWk;
	float* pWgt = pWk + (N*N);
	float* pWgtR = pWgt;
	float* pWgtG = pWgt + N;
	float* pWgtB = pWgt + N*2;
	float* pTmp = pWgt + N*3;
	int w = pImg->mWidth;
	int h = pImg->mHeight;
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			int idx = y*w + x;
			if (pImg->mPixels[idx].a == 0.0f) {
				IMG_COORD pt;
				float u = float(x);
				float v = float(y);
				pt.set(u, v);
				float sr = 0;
				float sg = 0;
				float sb = 0;
				for (int i = 0; i < N; ++i) {
					float d = pt_dist(pt, pPts[i]);
					sr += pWgtR[i] * d;
					sg += pWgtG[i] * d;
					sb += pWgtB[i] * d;
				}
				pImg->mPixels[idx].set(sr, sg, sb, 1.0f);
			}
		}
	}
}

IMAGE* img_interp(int w, int h, IMG_COORD* pPts, RGBA* pSmps, int nsmp) {
	IMAGE* pImg = img_alloc(w, h);
	if (pImg) {
		pImg->clear(0.0f);
		for (int i = 0; i < nsmp; ++i) {
			pImg->set_pixel(pPts[i], pSmps[i]);
			pImg->set_alpha(pPts[i], 1.0f);
		}
		float* pWk = alloc_solution(nsmp);
		if (pWk) {
			if (solve(pWk, pPts, pSmps, nsmp)) {
				interp(pImg, pWk, pPts, nsmp);
			}
			mem_free(pWk);
			pWk = nullptr;
		}
	}
	return pImg;
}
