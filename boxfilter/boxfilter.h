#ifndef __BOXFILTER_H_
#define  __BOXFILTER_H_

#ifndef use_sse_version
#define use_sse_version (0) // 1: use sse to get faster  0: slower than sse
#endif

/*
@function    BoxfilterFilter
@param       [in]      src:						  input image,single channel
@param       [in/out]  dst:                       output image,single channel
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      c:						  channel of image, only c = 1
@param       [in]      r:						  local window radius
@return£º										  0:ok; 1:error
@brief£º

*/
int BoxfilterFilter(unsigned char *src, unsigned char *dst,
	int w, int h, int c, int r);

#endif