#include "img_sys.hpp"
#include "img_interp.hpp"

#define OUT_PATH(_fname) "../../out/" _fname

static bool ck_ext(const char* pStr, const char* pExt) {
	if (pStr && pExt) {
		size_t lenStr = ::strlen(pStr);
		size_t lenPost = ::strlen(pExt);
		if (lenStr < lenPost) return false;
		for (size_t i = 0; i < lenPost; ++i) {
			if (pStr[lenStr - lenPost + i] != pExt[i]) return false;
		}
		return true;
	}
	return false;
}

static void save_samples_png(int w, int h, IMG_COORD* pPts, RGBA* pSmps, int nsmp) {
	IMAGE* pSmpImg = img_alloc(w, h);
	if (pSmpImg) {
		pSmpImg->clear();
		for (int i = 0; i < nsmp; ++i) {
			pSmpImg->set_pixel(pPts[i], pSmps[i]);
		}
		img_save_png(OUT_PATH("img_smp.png"), pSmpImg);
		mem_free(pSmpImg);
	}
}

static void save_samples_geo(int w, int h, IMG_COORD* pPts, RGBA* pSmps, int nsmp) {
	FILE* pOut = file_open(OUT_PATH("img_smp.geo"), "wb");
	if (pOut) {
		::fprintf(pOut, "PGEOMETRY V5\nNPoints %d NPrims 0\nNPointGroups 0 NPrimGroups 0\n", nsmp);
		::fprintf(pOut, "NPointAttrib 1 NVertexAttrib 0 NPrimAttrib 0 NAttrib 0\n");
		::fprintf(pOut, "PointAttrib\n");
		::fprintf(pOut, "Cd 3 float 1 1 1\n");
		float sx = 1.0f / float(w - 1);
		float sy = 1.0f / float(h - 1);
		for (int i = 0; i < nsmp; ++i) {
			::fprintf(pOut, "%.10f %.10f 0 1 (%.10f %.10f %.10f)\n",
				pPts[i].x * sx, 1.0f - (pPts[i].y * sy),
				pSmps[i].r, pSmps[i].g, pSmps[i].b
			);
		}
		::fprintf(pOut, "beginExtra\nendExtra\n");
		file_close(pOut);
	}
}

int main(int argc, char* argv[]) {
	if (argc <= 1) {
		return -1;
	}
	const char* pSrcPath = argv[1];
	IMAGE* pImg = nullptr;
	if (ck_ext(pSrcPath, ".jpg") || ck_ext(pSrcPath, ".jpeg") || ck_ext(pSrcPath, ".JPG")) {
		pImg = img_load_jpg(pSrcPath);
	} else if (ck_ext(pSrcPath, ".png")) {
		pImg = img_load_png(pSrcPath);
	}
	int nsmp = 1000;
	if (pImg) {
		int w = pImg->mWidth;
		int h = pImg->mHeight;
		int npix = w * h;
		::printf("Input image: %dx%d (%d pixels)\n", w, h, npix);
		img_save_png(OUT_PATH("img_org.png"), pImg);
		if (argc > 2) {
			nsmp = ::atoi(argv[2]);
		}
		int smpMin = int(npix / 100.0 * 1.5);
		if (nsmp < 1) {
			nsmp = smpMin;
		}
		if (nsmp > npix) {
			nsmp = int(5.0 * npix / 100.0);
		}
		IMG_COORD* pPts = scatter_pts_rand(w, h, nsmp);
		if (pPts) {
			RGBA* pSmps = get_img_samples(pImg, pPts, nsmp);
			if (pSmps) {
				save_samples_png(w, h, pPts, pSmps, nsmp);
				save_samples_geo(w, h, pPts, pSmps, nsmp);
				::printf("# samples: %d (%f%%)\n", nsmp, nsmp / (npix / 100.0));
				double t0 = time_micros();
				IMAGE* pInterpImg = img_interp(w, h, pPts, pSmps, nsmp);
				double dt = time_micros() - t0;
				::printf("Elapsed time: %f microseconds (%f secs).\n", dt, dt / 1e6);
				if (pInterpImg) {
					img_save_png(OUT_PATH("img_interp.png"), pInterpImg);
				}
			}
		}
		mem_free(pImg);
		pImg = nullptr;
	}
	return 0;
}
