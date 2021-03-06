// Author: Sergey Chaban <sergey.chaban@gmail.com>

#define WIN32_LEAN_AND_MEAN 1
#define NOMINMAX
#define _WIN32_WINNT 0x0500
#include <Windows.h>
#include <tchar.h>
#include <Shlwapi.h>
#include <wincodec.h>

#include "img_sys.hpp"

double time_micros() {
	double ms = 0.0f;
	LARGE_INTEGER frq;
	if (::QueryPerformanceFrequency(&frq)) {
		LARGE_INTEGER ctr;
		::QueryPerformanceCounter(&ctr);
		ms = ((double)ctr.QuadPart / (double)frq.QuadPart) * 1.0e6;
	}
	return ms;
}

void con_locate(int x, int y) {
	HANDLE hstd = ::GetStdHandle(STD_OUTPUT_HANDLE);
	COORD c;
	c.X = (SHORT)x;
	c.Y = (SHORT)y;
	::SetConsoleCursorPosition(hstd, c);
}

void con_text_color(const RGBA& clr) {
	HANDLE hstd = ::GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO info;
	::GetConsoleScreenBufferInfo(hstd, &info);
	WORD attr = 0;
	if (clr.r >= 0.5f) attr |= FOREGROUND_RED;
	if (clr.g >= 0.5f) attr |= FOREGROUND_GREEN;
	if (clr.b >= 0.5f) attr |= FOREGROUND_BLUE;
	if (clr.luma() > 0.5f) attr |= FOREGROUND_INTENSITY;
	::SetConsoleTextAttribute(hstd, attr);
}

static IMAGE* img_load(const char* pPath, const IID& decoderClsId, bool linearize) {
	::CoInitialize(nullptr);
	IMAGE* pImg = nullptr;
	IStream* pStrm = nullptr;
	HRESULT hres = ::SHCreateStreamOnFileA(pPath, STGM_READ, &pStrm);
	if (SUCCEEDED(hres)) {
		IWICBitmapDecoder* pDec = nullptr;
		hres = ::CoCreateInstance(decoderClsId, nullptr, CLSCTX_INPROC_SERVER, __uuidof(IWICBitmapDecoder), (void**)&pDec);
		if (SUCCEEDED(hres)) {
			hres = pDec->Initialize(pStrm, WICDecodeMetadataCacheOnLoad);
			if (SUCCEEDED(hres)) {
				IWICBitmapFrameDecode* pFrm = nullptr;
				hres = pDec->GetFrame(0, &pFrm);
				if (SUCCEEDED(hres)) {
					IWICBitmapSource* pBmp = nullptr;
					hres = ::WICConvertBitmapSource(GUID_WICPixelFormat32bppRGBA, pFrm, &pBmp);
					if (SUCCEEDED(hres)) {
						UINT w = 0;
						UINT h = 0;
						hres = pBmp->GetSize(&w, &h);
						if (SUCCEEDED(hres)) {
							pImg = img_alloc(int(w), int(h));
							if (pImg) {
								VAL32* pTmpPix = type_alloc<VAL32>(w*h);
								if (pTmpPix) {
									WICRect rect;
									rect.X = 0;
									rect.Y = 0;
									rect.Width = w;
									rect.Height = h;
									UINT stride = UINT(w * sizeof(VAL32));
									hres = pBmp->CopyPixels(&rect, stride, (UINT)pImg->calc_data_size(), (BYTE*)pTmpPix);
									if (SUCCEEDED(hres)) {
										RGBA* pClr = pImg->mPixels;
										for (UINT i = 0; i < w*h; ++i) {
											pClr[i].set(float(pTmpPix[i].b[0]), float(pTmpPix[i].b[1]), float(pTmpPix[i].b[2]), float(pTmpPix[i].b[3]));
										}
										for (UINT i = 0; i < w*h; ++i) {
											pClr[i].scl(1.0f / 255);
										}
										if (linearize) {
											for (UINT i = 0; i < w*h; ++i) {
												RGBA* pWk = &pClr[i];
												for (int j = 0; j < 4; ++j) {
													pWk->ch[j] = ::powf(pWk->ch[j], 2.2f);
												}
											}
										}
									}
									mem_free(pTmpPix);
								}
							}
						}
						pBmp->Release();
						pBmp = nullptr;
					}
					pFrm->Release();
					pFrm = nullptr;
				}
			}
			pDec->Release();
			pDec = nullptr;
		}
		pStrm->Release();
		pStrm = nullptr;
	}
	::CoUninitialize();
	return pImg;
}

IMAGE* img_load_png(const char* pPath, bool linearize) {
	return img_load(pPath, CLSID_WICPngDecoder, linearize);
}

IMAGE* img_load_jpg(const char* pPath, bool linearize) {
	return img_load(pPath, CLSID_WICJpegDecoder, linearize);
}

static void img_save(const char* pPath, const IMAGE* pImg, const IID& encoderClsId, bool linear) {
	if (!pPath || !pImg) return;
	::CoInitialize(nullptr);
	IStream* pStrm = nullptr;
	HRESULT hres = ::SHCreateStreamOnFileA(pPath, STGM_CREATE | STGM_WRITE, &pStrm);
	if (SUCCEEDED(hres)) {
		IWICBitmapEncoder* pEnc = nullptr;
		hres = ::CoCreateInstance(encoderClsId, nullptr, CLSCTX_INPROC_SERVER, __uuidof(IWICBitmapEncoder), (void**)&pEnc);
		if (SUCCEEDED(hres)) {
			hres = pEnc->Initialize(pStrm, WICBitmapEncoderNoCache);
			if (SUCCEEDED(hres)) {
				IWICBitmapFrameEncode* pFrm = nullptr;
				hres = pEnc->CreateNewFrame(&pFrm, nullptr);
				if (SUCCEEDED(hres)) {
					pFrm->Initialize(nullptr);
					UINT w = pImg->mWidth;
					UINT h = pImg->mHeight;
					hres = hres = pFrm->SetSize(w, h);
					if (SUCCEEDED(hres)) {
						GUID fmt = GUID_WICPixelFormat32bppRGBA;
						hres = pFrm->SetPixelFormat(&fmt);
						if (SUCCEEDED(hres)) {
							const RGBA* pSrc = pImg->mPixels;
							uint32_t size = w*h * sizeof(uint32_t);
							uint32_t* pDst = (uint32_t*)mem_alloc(size);
							if (pDst) {
								for (UINT i = 0; i < w*h; ++i) {
									RGBA csrc = pSrc[i];
									VAL32 cdst;
									if (linear) {
										for (int i = 0; i < 4; ++i) {
											csrc.ch[i] = ::powf(csrc.ch[i], 1.0f / 2.2f);
										}
									}
									for (int i = 0; i < 4; ++i) {
										cdst.b[i] = uint8_t(saturate(csrc.ch[i]) * 255.0f);
									}
									uint8_t t = cdst.b[0];
									cdst.b[0] = cdst.b[2];
									cdst.b[2] = t;
									pDst[i] = cdst.u;
								}
								UINT stride = w * sizeof(uint32_t);
								hres = pFrm->WritePixels(h, stride, size, (BYTE*)pDst);
								if (SUCCEEDED(hres)) {
									hres = pFrm->Commit();
								}
								mem_free(pDst);
								pDst = nullptr;
							}
						}
					}
					pFrm->Release();
					pFrm = nullptr;
				}
			}
			hres = pEnc->Commit();
			pEnc->Release();
			pEnc = nullptr;
		}
		pStrm->Release();
		pStrm = nullptr;
	}
	::CoUninitialize();
}

void img_save_png(const char* pPath, const IMAGE* pImg, bool linear) {
	img_save(pPath, pImg, CLSID_WICPngEncoder, linear);
}
