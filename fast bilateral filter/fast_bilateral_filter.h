#ifndef __FAST_BILATERAL_FILTER_H_
#define __FAST_BILATERAL_FILTER_H_


/*
@function    FastBilateralFilter
@param       [in]      src:						  input image
@param       [in]      guidance:                  guided image,single channel, only single channel is valid; invalid parameter in three channels 
@param       [in/out]] dst:                       output image
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      c:						  channel of image, only c = 1 or c = 3
@param       [in]      sigma_s:					  filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
@param       [in]      sigma_r:			          filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
@return£º										  0:ok; 1:error
@brief£º
			if guidance is NULL, still could get color filter
*/
int FastBilateralFilter(unsigned char *src, unsigned char *guidance, unsigned char *dst, 
	int w, int h, int c, float sigma_s, float sigma_r);

#endif