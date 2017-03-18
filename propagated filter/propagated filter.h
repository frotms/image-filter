#ifndef __PROPAGATED_FILTER_H_
#define __PROPAGATED_FILTER_H_


#ifndef use_fixed_type
#define use_fixed_type (0)	// 0: original; 1: use float-type to interger-type; 
#endif


/*
@function    PropagatedFilter
@param       [in]      src:						  input image
@param       [in]      guidance:                  guided image
@param       [in/out]] dst:                       output image
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      c:						  channel of image, only c = 1 or c = 3
@param       [in]      r:						  local window radius
@param       [in]      sigma_s:					  filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
@param       [in]      sigma_r:			          filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
@return£º										  0:ok; 1:error
@brief£º

*/
int PropagatedFilter(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int c, int r, float sigma_s, float sigma_r);

#endif