
IMG_COORD* scatter_pts_rand(int w, int h, int n, int seed = 7);
RGBA* get_img_samples(const IMAGE* pImg, const IMG_COORD* pPts, int n);
IMAGE* img_interp(int w, int h, IMG_COORD* pPts, RGBA* pSmps, int nsmp);
