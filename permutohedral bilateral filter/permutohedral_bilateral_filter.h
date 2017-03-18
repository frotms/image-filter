#ifndef __PERMUTOHEDRAL_BILATERAL_FILTER_H_
#define __PERMUTOHEDRAL_BILATERAL_FILTER_H_


/*
@function    fast_bilateral_filter
@param       [in]      src:						  input image
@param       [in]      guidance:                  guided image
@param       [in/out]] dst:                       output image
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      c:						  channel of image, only c = 1 or c = 3
@param       [in]      sigma_s:					  filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
@param       [in]      sigma_r:			          filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
@return£º										  0:ok; 1:error
@brief£º
try:(src,guidance,dst,w,h,c,1.6f, 0.6f)
*/
int permutohedral_bilateral_filter(unsigned char *src, unsigned char *guidance, unsigned char *dst,
    int w, int h, int c, float sigma_s, float sigma_r);



#endif