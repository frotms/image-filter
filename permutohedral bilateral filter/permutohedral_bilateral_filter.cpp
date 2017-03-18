#include "permutohedral_bilateral_filter.h"

#include <stdio.h>
#include "macros.h"
#include "Image.h"
#include "permutohedral.h"


static int _permutohedral_bilateral_filter_singlechannel(unsigned char *src, unsigned char *guidance, unsigned char *dst,
    int w, int h, float sigma_s, float sigma_r)
{
    int k, x, y, offy;
    float _scale = 255 * 255;	//if regularization, _scale = 1; if no regularization, _scale =  255*255
    float invSpatialStdev = 1.0f / sigma_s;
    float invColorStdev = 1.0f / (sigma_r * _scale);
    float *src_flt = new float[w*h];

    PBFimage input(1, w, h, 1), out;
    PBFimage positions(1, input.width, input.height, 3);//5

    for (k = 0; k < w*h; k++)
        src_flt[k] = src[k];

    input.data = src_flt;

    // Construct the position vectors out of x, y, r, g, and b.
    for (y = 0; y < h; y++)
    {
        offy = y*w;
        for (x = 0; x < w; x++)
        {
            positions(x, y)[0] = invSpatialStdev * x;
            positions(x, y)[1] = invSpatialStdev * y;
            positions(x, y)[2] = invColorStdev * guidance[offy + x];

        }
    }

    // Filter the input with respect to the position vectors. (see permutohedral.h)
    out = PermutohedralLattice::filter(input, positions);

    for (k = 0; k < w*h; k++)
        dst[k] = (out.data[k] + 0.5f);

    return 0;
}

static int permutohedral_bilateral_filter_singlechannel(unsigned char *src, unsigned char *guidance, unsigned char *dst,
    int w, int h, int c, float sigma_s, float sigma_r)
{
    if (src == NULL || guidance == NULL || dst == NULL)
        return 1;

    if (w <= 0 || h <= 0 || (c != 1))
        return 1;

    return _permutohedral_bilateral_filter_singlechannel(src, guidance, dst, w, h, sigma_s, sigma_r);
}

static int _permutohedral_bilateral_filter_color(unsigned char *src, unsigned char *guidance, unsigned char *dst,
    int w, int h, int c, float sigma_s, float sigma_r)
{
    int k, x, y, offy, offx;
    float _scale = 255 * 255;	//if regularization, _scale = 1; if no regularization, _scale =  255*255
    float invSpatialStdev = 1.0f / sigma_s;
    float invColorStdev = 1.0f / (sigma_r * _scale);
    float *src_flt = new float[w*h*c];

    PBFimage input(1, w, h, 3), out;
    PBFimage positions(1, input.width, input.height, 3);//5

    for (k = 0; k < w*h*c; k++)
        src_flt[k] = src[k];

    input.data = src_flt;

    // Construct the position vectors out of x, y, r, g, and b.
    for (y = 0; y < h; y++)
    {
        offy = y*w*c;
        for (x = 0; x < w; x++)
        {
            offx = x*c;
            positions(x, y)[0] = invSpatialStdev * x;
            positions(x, y)[1] = invSpatialStdev * y;
            positions(x, y)[2] = invColorStdev * guidance[offy + offx + 0];
            positions(x, y)[3] = invColorStdev * guidance[offy + offx + 1];
            positions(x, y)[4] = invColorStdev * guidance[offy + offx + 2];
        }
    }

    // Filter the input with respect to the position vectors. (see permutohedral.h)
    out = PermutohedralLattice::filter(input, positions);

    for (k = 0; k < w*h*c; k++)
        dst[k] = (out.data[k] + 0.5f);

    return 0;
}

static int permutohedral_bilateral_filter_color(unsigned char *src, unsigned char *guidance, unsigned char *dst,
    int w, int h, int c, float sigma_s, float sigma_r)
{
    if (src == NULL || guidance == NULL || dst == NULL)
        return 1;

    if (w <= 0 || h <= 0 || (c != 3))
        return 1;

    return _permutohedral_bilateral_filter_color(src, guidance, dst, w, h, c, sigma_s, sigma_r);
}


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
    int w, int h, int c, float sigma_s, float sigma_r)
{
    int _code = 0;
    if (src == NULL || guidance == NULL || dst == NULL)
        return 1;

    if (w <= 0 || h <= 0 || !(c == 1 || c == 3))
        return 1;

    switch (c)
    {
    case 1:
        _code = permutohedral_bilateral_filter_singlechannel(src, guidance, dst, w, h, c, sigma_s, sigma_r);
        break;

    case 3:
        _code = permutohedral_bilateral_filter_color(src, guidance, dst, w, h, c, sigma_s, sigma_r);
        break;

    default:
        _code = 1;
        break;
    }


    return _code;

}