#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

union VAL32 {
	int32_t i;
	uint32_t u;
	float f;
	uint8_t b[4];
};

union VAL64 {
	int64_t i;
	uint64_t u;
	double f;
	uint8_t b[8];
};

union RGBA {
	struct { float r, g, b, a; };
	float ch[4];

	void zero() {
		for (int i = 0; i < 4; ++i) {
			ch[i] = 0.0f;
		}
	}

	void set(float r, float g, float b, float a = 1.0f) {
		this->r = r;
		this->g = g;
		this->b = b;
		this->a = a;
	}

	void scl(float s) {
		for (int i = 0; i < 4; ++i) {
			ch[i] *= s;
		}
	}

	float luma() const {
		return r*0.299f + g*0.587f + b*0.114f;
	}

	float luminance() const {
		return r*0.212671f + g*0.71516f + b*0.072169f;
	}

	float average() const {
		return (r + g + b) / 3.0f;
	}
};

struct IMG_COORD {
	float x, y;

	void set(float x, float y) {
		this->x = x;
		this->y = y;
	}
};

struct IMAGE {
	int32_t mWidth;
	int32_t mHeight;
	int32_t mReserved[2];
	RGBA mPixels[1];

	size_t calc_data_size() const { return mWidth * mHeight * sizeof(RGBA); }

	void clear(float alpha = 1.0f) {
		int n = mWidth * mHeight;
		for (int i = 0; i < n; ++i) {
			mPixels[i].zero();
		}
		if (alpha != 0.0f) {
			for (int i = 0; i < n; ++i) {
				mPixels[i].a = alpha;
			}
		}
	}

	RGBA get_pixel(int x, int y) const {
		return mPixels[y*mWidth + x];
	}

	RGBA get_pixel(const IMG_COORD& pt) const {
		return get_pixel(int(pt.x), int(pt.y));
	}

	void set_pixel(int x, int y, RGBA c) {
		mPixels[y*mWidth + x] = c;
	}

	void set_pixel(const IMG_COORD& pt, RGBA c) {
		set_pixel(int(pt.x), int(pt.y), c);
	}

	void set_alpha(int x, int y, float a) {
		mPixels[y*mWidth + x].a = a;
	}

	void set_alpha(const IMG_COORD& pt, float a) {
		set_alpha(int(pt.x), int(pt.y), a);
	}
};

template<typename T> inline T clamp(T x, T lo, T hi) { return x < lo ? lo : x > hi ? hi : x; }
template<typename T> inline T saturate(T x) { return clamp(x, T(0), T(1)); }

inline float pt_dist(const IMG_COORD& ptA, const IMG_COORD& ptB) {
	float dx = ptA.x - ptB.x;
	float dy = ptA.y - ptB.y;
	return ::sqrtf(dx*dx + dy*dy);
}

double time_micros();

void* mem_alloc(int64_t size, const char* pTag = "Temp");
void mem_free(void* pMem);
template<typename T> inline T* type_alloc(int n, const char* pTag = "Temp") { return reinterpret_cast<T*>(mem_alloc(n * sizeof(T), pTag)); }

FILE* file_open(const char* pPath, const char* mode);
void file_close(FILE* pFile);

#define D_IMG_LIN false

IMAGE* img_alloc(int w, int h);
IMAGE* img_load_png(const char* pPath, bool linearize = D_IMG_LIN);
IMAGE* img_load_jpg(const char* pPath, bool linearize = D_IMG_LIN);
void img_save_png(const char* pPath, const IMAGE* pImg, bool linear = D_IMG_LIN);

