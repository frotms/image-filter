#ifndef __FAST_GUIDED_FILTER_H_
#define __FAST_GUIDED_FILTER_H_

#define use_sse (1)	// 1: use sse to get faster  0: slower than sse


/*******************************************************************
** 函数名:     fast_guided_filter
** 函数描述:   guided filter
** note:       2017.03.16目前只支持单通道
**			   eg: r = 4, (try sr = r/4 to sr = r),(try rp=0.1^2, 0.2^2, 0.4^2)
**
**             (MIN(w, h) / sr) > 1
**			   (int)(r / sr + 0.5f) >= 1
**
** 参数:       [in]      src:                          输入待滤波图像,单通道
**             [in]      guidance:                     输入导向图像,单通道
**             [in/out]  dst:                          输出滤波后图像,单通道
**             [in]      w:                            图像宽度
**             [in]      h:                            图像高度
**             [in]      c:                            图像通道数,目前只支持单通道
**             [in]      r:                            滤波窗口半径
**             [in]      rp:                           规定化参数
**             [in]      sr:                           降采样率
** 返回:                                               错误枚举值
********************************************************************/

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
@return：										  0:ok; 1:error
@brief：

eg: r = 4, (try sr = r/4 to sr=r),(try rp=0.1^2, 0.2^2, 0.4^2)
try:(src,guidance,dst,w,h,1,4,0.01,4)
condition: (MIN(w, h) / sr) > 1
condition: (int)(r / sr + 0.5f) >= 1
*/
int fast_guided_filter(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int c, int r, float rp, float sr);

#endif