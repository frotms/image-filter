#include "fast_bilateral_filter.h"


#include <cmath>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#define CHRONO
#include "geom.h"
#include "fast_lbf.h"
#include "fast_color_bf.h"


using namespace std;

typedef Array_2D<double> image_type;
typedef Array_2D<Geometry::Vec3<double> > image_type3;

/*
@function    _fast_bilateral_filter_singlechannel
@param       [in]      src:						  input image,single channel
@param       [in]      guidance:                  guided image,single channel
@param       [in/out]  dst:                       output image,single channel
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      sigma_s:					  filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
@param       [in]      sigma_r:			          filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
@return£º										  0:ok; 1:error
@brief£º

*/
void _fast_bilateral_filter_singlechannel(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, float sigma_s, float sigma_r)
{
	int x, y, offy;
	float _scale = 255 * 255;	//if regularization, _scale = 1; if no regularization, _scale =  255*255
	image_type _src(w, h);
	image_type _guidance(w, h);
	image_type _dst(w, h);

	for (y = 0; y < h; y++)
	{
		offy = y*w;
		for (x = 0; x < w; x++)
		{
			_src(x, y) = src[offy + x];
			_guidance(x, y) = guidance[offy + x];
		}
	}

	Image_filter::fast_LBF(_src, _guidance,
		(double)sigma_s, (double)(_scale*sigma_r),
		0,
		&_dst, &_dst);



	for (y = 0; y < h; y++)
	{
		offy = y*w;
		for (x = 0; x < w; x++)
		{
			dst[offy + x] = (unsigned char)(_dst(x, y) + 0.5f);// *255.0 + 0.5);
		}
	}


}

/*
@function    fast_bilateral_filter_singlechannel
@param       [in]      src:						  input image,single channel
@param       [in]      guidance:                  guided image,single channel
@param       [in/out]  dst:                       output image,single channel
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      c:						  channel of image, only c = 1
@param       [in]      sigma_s:					  filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
@param       [in]      sigma_r:			          filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
@return£º										  0:ok; 1:error
@brief£º

*/
int fast_bilateral_filter_singlechannel(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, int c, float sigma_s, float sigma_r)
{
	if (src == NULL || guidance == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || (c != 1))
		return 1;

	_fast_bilateral_filter_singlechannel(src, guidance, dst, w, h, sigma_s, sigma_r);

	return 0;
}



/*
@function    _fast_bilateral_filter_color
@param       [in]      src:						  input image, bgr...bgr
@param       [in/out]  dst:                       output image, bgr...bgr
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      sigma_s:					  filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
@param       [in]      sigma_r:			          filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
@return£º										  0:ok; 1:error
@brief£º

*/
void _fast_bilateral_filter_color(unsigned char *src, unsigned char *dst, int w, int h, float sigma_s, float sigma_r)
{
	int x, y, offy, offx;
	float _scale = 255 * 255;	//if regularization, _scale = 1; if no regularization, _scale =  255*255
	image_type3 _src(w, h);
	image_type3 _guidance(w, h);
	image_type3 _dst(w, h);

	for (y = 0; y < h; y++)
	{
		offy = y * w * 3;
		for (x = 0; x < w; x++)
		{
			offx = x * 3;
			_src(x, y)[0] = src[offy + offx + 2];
			_src(x, y)[1] = src[offy + offx + 1];
			_src(x, y)[2] = src[offy + offx + 0];

		}
	}

	Image_filter::fast_color_BF(_src,
		(double)sigma_s, (double)(_scale*sigma_r),
		&_dst);

	for (unsigned y = 0; y < h; y++)
	{
		offy = y * w * 3;
		for (unsigned x = 0; x < w; x++)
		{
			offx = x * 3;
			dst[offy + offx + 2] = (unsigned char)(_dst(x, y)[0] + 0.5);
			dst[offy + offx + 1] = (unsigned char)(_dst(x, y)[1] + 0.5);
			dst[offy + offx + 0] = (unsigned char)(_dst(x, y)[2] + 0.5);

		}
	}


}


/*
@function    fast_bilateral_filter_color
@param       [in]      src:						  input image, bgr...bgr
@param       [in/out]  dst:                       output image, bgr...bgr
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      c:						  channel of image, only c = 3
@param       [in]      sigma_s:					  filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
@param       [in]      sigma_r:			          filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
@return£º										  0:ok; 1:error
@brief£º

*/
int fast_bilateral_filter_color(unsigned char *src, unsigned char *dst, int w, int h, int c, float sigma_s, float sigma_r)
{
	if (src == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || (c != 3))
		return 1;

	_fast_bilateral_filter_color(src, dst, w, h, sigma_s, sigma_r);

	return 0;
}


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
int FastBilateralFilter(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, int c, float sigma_s, float sigma_r)
{
	int _code = 0;
	if (src == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || !(c == 1 || c == 3))
		return 1;

	switch (c)
	{
	case 1:
		_code = fast_bilateral_filter_singlechannel(src, guidance, dst, w, h, c, sigma_s, sigma_r);
		break;

	case 3:
		_code = fast_bilateral_filter_color(src, dst, w, h, c, sigma_s, sigma_r);
		break;

	default:
		_code = 1;
		break;
	}
	

	return _code;
}

