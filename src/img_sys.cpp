// Author: Sergey Chaban <sergey.chaban@gmail.com>

#include "img_sys.hpp"

#include <stdarg.h>
#include <time.h>

#if !defined(_WIN32)
double time_micros() {
	double ms = 0.0f;
	struct timespec t;
	if (clock_gettime(CLOCK_MONOTONIC, &t) != 0) {
		clock_gettime(CLOCK_REALTIME, &t);
	}
	ms = (double)t.tv_nsec*1.0e-3 + (double)t.tv_sec*1.0e6;
	return ms;
}
#endif

struct MEM_INFO {
	MEM_INFO* mpPrev;
	MEM_INFO* mpNext;
	void* mpRaw;
	char* mpTag;
	int64_t mSize;
};

static struct {
	MEM_INFO* mpHead;
	MEM_INFO* mpTail;
	uint32_t mAllocs;
	int64_t mBytes;
	int64_t mPeak;

	void* alloc(int64_t size, const char* pTag) {
		void* pMem = nullptr;
		if (size > 0) {
			int64_t align = 0x10;
			int64_t mask = align - 1;
			size_t tlen = pTag ? ::strlen(pTag) : 0;
			size_t hsize = sizeof(MEM_INFO);
			int64_t bsize = hsize + size + (tlen + 1) + align;
			int64_t asize = (bsize + mask) & (~mask);
			void* pRaw = ::malloc(size_t(asize));
			if (pRaw) {
				MEM_INFO* pInfo;
				pMem = (void*)((intptr_t)((uint8_t*)pRaw + sizeof(MEM_INFO) + mask) & (~mask));
				pInfo = &((MEM_INFO*)pMem)[-1];
				pInfo->mpRaw = pRaw;
				pInfo->mpTag = (char*)pMem + size;
				if (tlen > 0) {
					::memcpy(pInfo->mpTag, pTag, tlen);
				}
				pInfo->mpTag[tlen] = 0;
				pInfo->mSize = size;
				pInfo->mpPrev = nullptr;
				pInfo->mpNext = nullptr;
				if (!mpHead) {
					mpHead = pInfo;
				} else {
					pInfo->mpPrev = mpTail;
				}
				if (mpTail) {
					mpTail->mpNext = pInfo;
				}
				mpTail = pInfo;
				mBytes += size;
				if (mBytes > mPeak) {
					mPeak = mBytes;
				}
				++mAllocs;
			}
		}
		return pMem;
	}

	void free(void* pMem) {
		MEM_INFO* pInfo;
		MEM_INFO* pNext;
		void* pRaw;
		if (!pMem) return;
		pInfo = &((MEM_INFO*)pMem)[-1];
		pNext = pInfo->mpNext;
		if (pInfo->mpPrev) {
			pInfo->mpPrev->mpNext = pNext;
			if (pNext) {
				pNext->mpPrev = pInfo->mpPrev;
			} else {
				mpTail = pInfo->mpPrev;
			}
		} else {
			mpHead = pNext;
			if (pNext) {
				pNext->mpPrev = nullptr;
			} else {
				mpTail = nullptr;
			}
		}
		mBytes -= pInfo->mSize;
		pRaw = pInfo->mpRaw;
		::free(pRaw);
		--mAllocs;
	}
} s_mem = { nullptr, nullptr, 0, 0, 0 };

void* mem_alloc(int64_t size, const char* pTag) {
	return s_mem.alloc(size, pTag);
}

void mem_free(void* pMem) {
	s_mem.free(pMem);
}

FILE* file_open(const char* pPath, const char* mode) {
#if defined(_MSC_VER)
	FILE* f;
	::fopen_s(&f, pPath, mode);
	return f;
#else
	return ::fopen(pPath, mode);
#endif
}

void file_close(FILE* pFile) {
	::fclose(pFile);
}

IMAGE* img_alloc(int w, int h) {
	IMAGE* pImg = nullptr;
	if (w > 0 && h > 0) {
		int size = sizeof(IMAGE) + ((w * h * sizeof(RGBA)) - 1);
		pImg = reinterpret_cast<IMAGE*>(mem_alloc(size, "Image"));
		if (pImg) {
			pImg->mWidth = w;
			pImg->mHeight = h;
		}
	}
	return pImg;
}

