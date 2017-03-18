#include "propagated filter.h"


#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <omp.h>

using namespace std;


#ifdef use_fixed_type
#if use_fixed_type

//resolution factor (float-type to fixed-type)(faster)
#ifndef flt2intFactor
#define flt2intFactor 10	// only flt2intFactor*sigma<65536 could be valid by integer-type
#endif

#endif
#endif

// access matlab 2D matrix
//#define getRef2D(dataPtr,i,j,m,n) dataPtr[(i)+(j)*(m)]
#define getRef2D(dataPtr,i,j,m,n) dataPtr[(j)+(i)*(n)]

// access matlab 3D matrix
//#define getRef3D(dataPtr,i,j,k,m,n,z) dataPtr[(i)+((j)+(k)*(n))*(m)]
#define getRef3D(dataPtr,i,j,k,m,n,z) dataPtr[(k)+((j)+(i)*(n))*(z)]	 // bgr


struct PatternInfo
{
	int fromLocal[2];    
	int toLocal[2];           
};

//float-type to integer-type test
#if use_fixed_type


/**
* Calculate the log relationship of two points in R
* @param RPtr
* @param sigma
* @param fromPoint
* @param toPoint
* @param nRows number of rows in R
* @param nCols number of cols in R
* @return
*/

int calculateLogRelationship(const unsigned char* RPtr, int  sigma, const int(&fromPoint)[2],
	const int(&toPoint)[2], int nRows, int nCols, int zR)
{
	int z;
	int diff;
	int distanceSquare = 0;
	for (z = 0; z < zR; ++z)
	{
		diff = getRef3D(RPtr, fromPoint[0], fromPoint[1], z, nRows, nCols, zR)
			- getRef3D(RPtr, toPoint[0], toPoint[1], z, nRows, nCols, zR);
		distanceSquare += diff*diff;
	}

	return -1 * distanceSquare / (2 * sigma * sigma);
}



/**
* Calculate the log-weight between centerPoint and the toPoint through fromPoint
* @param RPtr
* @param sigma_d
* @param sigma_r
* @param fromPoint
* @param toPoint
* @param nRows number of rows in R
* @param nCols number of cols in R
* @return
*/


int calculateLogWeight(const unsigned char* RPtr, int  sigma_d, int  sigma_r,
	const int* logWeightWindowPtr, int w, const int(&centerPoint)[2], const int(&fromLocal)[2],
	const int(&toLocal)[2], int nRows, int nCols, int zR)
{
	/*计算模板内坐标*/
	int fromPoint[2] = { centerPoint[0] + fromLocal[0], centerPoint[1] + fromLocal[1] };
	int toPoint[2] = { centerPoint[0] + toLocal[0], centerPoint[1] + toLocal[1] };

	/*计算log权重*/
	int pathLogProb = calculateLogRelationship(RPtr, sigma_d, fromPoint, toPoint, nRows, nCols, zR);
	int rangeLogProb = calculateLogRelationship(RPtr, sigma_r, centerPoint, toPoint, nRows, nCols, zR);

	/* w_(s,t-1) + D(t-1,t) + R(s,t) */
	int logWeight = getRef2D(logWeightWindowPtr, fromLocal[0] + w, fromLocal[1] + w, 2 * w + 1, 2 * w + 1)
		+ pathLogProb + rangeLogProb;

	return logWeight;
}

// exp(x/a)
static int exp_int_taylor(float a, int x, int order)
{
	int ret;
	float ia = 1 / a;
	
	int x2, x3, x4;
	float ia2, ia3, ia4;

	switch (order)
	{
	case 1:
		ret = (1 + ia*x)*a;
		break;
	case 2:
		ret = (1 + ia*x + 0.5f*ia*ia*x*x)*a;
		break;
	case 3:
		ia2 = ia*ia;
		x2 = x*x;
		ret = (1 + ia*x + 0.5f*ia2*x2 + 0.1666667f*ia2*ia*x2*x)*a;
		break;
	case 4:
		ia2 = ia*ia;
		x2 = x*x;
		ia3 = ia2*ia;
		x3 = x2*x;
		ret = (1 + ia*x + 0.5f*ia2*x2 + 0.1666667f*ia3*x3 + 0.0416667f*ia3*ia * x3*x)*a;
		break;
	case 5:
		ia2 = ia*ia;
		x2 = x*x;
		ia3 = ia2*ia;
		x3 = x2*x;
		ia4 = ia3*ia;
		x4 = x3*x;
		ret = (1 + ia*x + 0.5f*ia2*x2 + 0.1666667f*ia3*x3 + 0.0416667f*ia4 * x4 + 0.0083333f*ia4*ia*x4*x)*a;
		break;

	default:
		ret = (1 + ia*x + 0.5f*ia*ia*x*x)*a;
		break;
	}


	return ret;
}

/**
* Calculate weighted average of pixels surrounding the point
* @param resultPixelPtr the array (with the size = zA) contains the result pixel values
* @param RPtr
* @param APtr
* @param w
* @param sigma_d
* @param sigma_r
* @param point
* @param m the number of rows of R and A
* @param n the number of columns of R and A
* @return void
*/


void calculateWeightedAverage(unsigned char* resultPixelPtr, const unsigned char* RPtr, const unsigned char* APtr, 
	const PatternInfo *pPatternInfo,
	int w, int  sigma_d, int  sigma_r, const int(&point)[2], int m, int n, int zR, int zA) {

	int k, z;
	int patternSize = 2 * w * (w + 1);
	int wSize = 2 * w + 1;

	int  weight_int;
	int  logWeight_int;

	int * logWeightWindowPtr = new int[wSize * wSize];
	int  totalWeight = 0, itotalWeight = 0;

	int * totalWeightedSum = new int[zA];
	for (int z = 0; z < zA; ++z)
	{
		totalWeightedSum[z] = 0;
	}

	// Calculate the weight of the Center Point
	getRef2D(logWeightWindowPtr, w, w, wSize, wSize) = 0;
	totalWeight += 10;
	for (z = 0; z < zA; ++z)
	{
		totalWeightedSum[z] += getRef3D(APtr, point[0], point[1], z, m, n, zA);
	}

	for (k = 0; k < patternSize; k++)
	{
		if (point[0] + pPatternInfo[k].toLocal[0] < 0 || point[0] + pPatternInfo[k].toLocal[0] > m - 1
			|| point[1] + pPatternInfo[k].toLocal[1] < 0 || point[1] + pPatternInfo[k].toLocal[1] > n - 1)
		{
			continue;
		}


		logWeight_int = calculateLogWeight(RPtr, sigma_d, sigma_r, logWeightWindowPtr,
			w, point, pPatternInfo[k].fromLocal, pPatternInfo[k].toLocal, m, n, zR);

		getRef2D(logWeightWindowPtr, w + pPatternInfo[k].toLocal[0],
			w + pPatternInfo[k].toLocal[1], wSize, wSize) = logWeight_int;

		//weight = exp(logWeight_int * 0.00390625f);// 1 / 256
		weight_int = exp_int_taylor(flt2intFactor, logWeight_int, 2);

		totalWeight += weight_int;

		for (z = 0; z < zA; ++z)
		{
			totalWeightedSum[z] += (weight_int * (int)getRef3D(APtr, point[0] + pPatternInfo[k].toLocal[0],
				point[1] + pPatternInfo[k].toLocal[1], z, m, n, zA)) ;
		}
	}

	//itotalWeight = 1.f / totalWeight;
	// Calculate result pixel value
	for (z = 0; z < zA; ++z)
	{
		resultPixelPtr[z] = (totalWeightedSum[z] / totalWeight);
	}

	delete[] logWeightWindowPtr;
	delete[] totalWeightedSum;

}



/*
@function    FilterPatternInit
@param       [in]      w:                         radius of pattern
@param       [in/out]  pPatternInfo:              pattern
@return：    
@brief：     PropagationFilter pattern generator 
*/
void FilterPatternInit(PatternInfo *pPatternInfo, int w)
{
	int idx = 0;
	int pattern_size = 2 * w * (w + 1);

	// Calculate from distance 1 to window_radius
	for (int r = 1; r < w + 1; ++r)
	{
		for (int dp = 0; dp < r + 1; ++dp)
		{
			for (int pSign = -1; pSign < 2; pSign += 2) // sign = -1, 1
			{
				int p = pSign*dp;


				for (int qSign = -1; qSign < 2; qSign += 2) // sign = -1, 1
				{
					int q = qSign * (r - dp);

					// decide fromLocal (the parent pixel t-1)

					if (p * q == 0) // on the x or y axis
					{
						if (p == 0)
						{
							pPatternInfo[idx].fromLocal[0] = p;
							pPatternInfo[idx].fromLocal[1] = q - qSign;
						}
						else // q == 0
						{
							pPatternInfo[idx].fromLocal[0] = p - pSign;
							pPatternInfo[idx].fromLocal[1] = q;
						}
					}
					else // p*q != 0 (other pixels)
					{
						// if r is odd -> p , else -> q
						if (r % 2 != 0)
						{
							pPatternInfo[idx].fromLocal[0] = p;
							pPatternInfo[idx].fromLocal[1] = q - qSign;
						}
						else
						{
							pPatternInfo[idx].fromLocal[0] = p - pSign;
							pPatternInfo[idx].fromLocal[1] = q;
						}
					}

					pPatternInfo[idx].toLocal[0] = p;
					pPatternInfo[idx].toLocal[1] = q;

					idx++;

					// ensure pixels on the axis is calculated only one time 
					if (q == 0)
						break;

				}

				// ensure pixels on the axis is calculated only one time 
				if (p == 0)
					break;

			}

		}// end dp

	}// end r

}


/**
* Main function of the propagation filter
* @param OutPtr pointer to an array where result will be stored
* @param RPtr pointer to an array containing reference matrix
* @param APtr pointer to an array containing input matrix
* @param nCols the number of columns of R and A
* @param nRows the number of rows of R and A
* @param r window radius
* @param sigma_d
* @param sigma_r
* @param zR the number of channels of R
* @param zA the number of channels of A
*/

void pfilter(unsigned char* OutPtr, const unsigned char* RPtr, const unsigned char* APtr, int nCols, int nRows,
	int r, float sigma_d, float sigma_r, int zR, int zA)
{

# ifdef OPENMP
	// set number of threads used by openmp
	// omp_set_num_threads(floor(omp_get_num_procs()/2));
	omp_set_num_threads(floor(omp_get_num_procs() - 2));
# endif

	PatternInfo *pPatternInfo = NULL;

	pPatternInfo = new PatternInfo[2 * r * (r + 1)]; 

	FilterPatternInit(pPatternInfo, r);
	//#pragma omp parallel default(shared)
	{
		unsigned char* resultPixelPtr = new unsigned char[zA];

		int sigma_d_int = (sigma_d * flt2intFactor);
		int sigma_r_int = (sigma_r * flt2intFactor);

		//#pragma omp for schedule(runtime)
		for (int i = 0; i < nRows; ++i)
		{
			for (int j = 0; j < nCols; ++j)
			{
				int point[2] = { i, j };
				calculateWeightedAverage(resultPixelPtr, RPtr, APtr, pPatternInfo,
					r, sigma_d_int, sigma_r_int, point, nRows, nCols, zR, zA);
				for (int z = 0; z < zA; ++z)
				{
					getRef3D(OutPtr, i, j, z, nRows, nCols, zA) = resultPixelPtr[z];
				}
			}
		}
		delete[] resultPixelPtr;
	}

	if (pPatternInfo)delete[]pPatternInfo;

	return;

}

#else

/**
* Calculate the log relationship of two points in R
* @param RPtr
* @param sigma
* @param fromPoint
* @param toPoint
* @param nRows number of rows in R
* @param nCols number of cols in R
* @return
*/

float calculateLogRelationship(const unsigned char* RPtr, float sigma, const int(&fromPoint)[2], 
	const int(&toPoint)[2], int nRows, int nCols, int zR)
{
	float diff;
	float distanceSquare = 0;
	for (int z = 0; z < zR; ++z)
	{
		diff = (float)getRef3D(RPtr, fromPoint[0], fromPoint[1], z, nRows, nCols, zR)
			- (float)getRef3D(RPtr, toPoint[0], toPoint[1], z, nRows, nCols, zR);
		distanceSquare += diff*diff;
	}

	return -1 * distanceSquare / (2 * sigma * sigma);
}



/**
* Calculate the log-weight between centerPoint and the toPoint through fromPoint
* @param RPtr
* @param sigma_d
* @param sigma_r
* @param fromPoint
* @param toPoint
* @param nRows number of rows in R
* @param nCols number of cols in R
* @return
*/


float calculateLogWeight(const unsigned char* RPtr, float sigma_d, float sigma_r,
	const float* logWeightWindowPtr, int w, const int(&centerPoint)[2], const int(&fromLocal)[2], 
	const int(&toLocal)[2], int nRows, int nCols, int zR)
{
	/*计算模板内坐标*/
	int fromPoint[2] = { centerPoint[0] + fromLocal[0], centerPoint[1] + fromLocal[1] };
	int toPoint[2] = { centerPoint[0] + toLocal[0], centerPoint[1] + toLocal[1] };

	/*计算log权重*/
	float pathLogProb = calculateLogRelationship(RPtr, sigma_d, fromPoint, toPoint, nRows, nCols, zR);
	float rangeLogProb = calculateLogRelationship(RPtr, sigma_r, centerPoint, toPoint, nRows, nCols, zR);

	/* w_(s,t-1) + D(t-1,t) + R(s,t) */
	float logWeight = (float)getRef2D(logWeightWindowPtr, fromLocal[0] + w, fromLocal[1] + w, 2 * w + 1, 2 * w + 1)
		+ pathLogProb + rangeLogProb;

	return logWeight;
}


// exp(x)
static float exp_float_taylor(float x, int order)
{
	float ret;
	float x2, x3, x4;

	switch (order)
	{
	case 1:
		ret = 1 + x;
		break;
	case 2:
		ret = 1 + x + 0.5f*x*x;
		break;
	case 3:
		x2 = x*x;
		ret = 1 + x + 0.5f*x2 + 0.1666667f*x2*x;
		break;
	case 4:
		x2 = x*x;
		x3 = x2*x;
		ret = 1 + x + 0.5f*x2 + 0.1666667f*x3 + 0.0416667f * x3*x;
		break;
	case 5:
		x2 = x*x;
		x3 = x2*x;
		x4 = x3*x;
		ret = 1 + x + 0.5f*x2 + 0.1666667f*x3 + 0.0416667f * x4 + 0.0083333f*x4*x;
		break;
	default:
		ret = 1 + x + 0.5f*x*x;
		break;
	}

	return ret;
}

/**
* Calculate weighted average of pixels surrounding the point
* @param resultPixelPtr the array (with the size = zA) contains the result pixel values
* @param RPtr
* @param APtr
* @param w
* @param sigma_d
* @param sigma_r
* @param point
* @param m the number of rows of R and A
* @param n the number of columns of R and A
* @return void
*/


void calculateWeightedAverage(unsigned char* resultPixelPtr, const unsigned char* RPtr, const unsigned char* APtr, const PatternInfo *pPatternInfo,
	int w, float sigma_d, float sigma_r, const int(&point)[2], int m, int n, int zR, int zA) {

	int k;
	int patternSize = 2 * w * (w + 1);
	int wSize = 2 * w + 1;

	float weight;
	float logWeight;

	float* logWeightWindowPtr = new float[wSize * wSize];
	float totalWeight = 0, itotalWeight = 0;
	float* totalWeightedSum = new float[zA];
	for (int z = 0; z < zA; ++z)
	{
		totalWeightedSum[z] = 0;
	}

	// Calculate the weight of the Center Point
	getRef2D(logWeightWindowPtr, w, w, wSize, wSize) = 0;
	totalWeight += 1.0;
	for (int z = 0; z < zA; ++z)
	{
		totalWeightedSum[z] += (float)getRef3D(APtr, point[0], point[1], z, m, n, zA);
	}

	for (k = 0; k < patternSize; k++)
	{
		if (point[0] + pPatternInfo[k].toLocal[0] < 0 || point[0] + pPatternInfo[k].toLocal[0] > m - 1
			|| point[1] + pPatternInfo[k].toLocal[1] < 0 || point[1] + pPatternInfo[k].toLocal[1] > n - 1)
		{
			continue;
		}


		logWeight = calculateLogWeight(RPtr, sigma_d, sigma_r, logWeightWindowPtr,
			w, point, pPatternInfo[k].fromLocal, pPatternInfo[k].toLocal, m, n, zR);

		getRef2D(logWeightWindowPtr, w + pPatternInfo[k].toLocal[0],
			w + pPatternInfo[k].toLocal[1], wSize, wSize) = logWeight;

		//weight = exp(logWeight);
		weight = exp_float_taylor(logWeight,2);

		totalWeight += weight;

		for (int z = 0; z < zA; ++z)
		{
			totalWeightedSum[z] += weight * getRef3D(APtr, point[0] + pPatternInfo[k].toLocal[0],
				point[1] + pPatternInfo[k].toLocal[1], z, m, n, zA);
		}
	}

	itotalWeight = 1.f / totalWeight;
	// Calculate result pixel value
	for (int z = 0; z < zA; ++z)
	{
		resultPixelPtr[z] = (totalWeightedSum[z] * itotalWeight + 0.5f);
	}

	delete[] logWeightWindowPtr;
	delete[] totalWeightedSum;

}



/*
@function    FilterPatternInit
@param       [in]      w:                         radius of pattern
@param       [in/out]  pPatternInfo:              pattern
@return：    
@brief：     PropagationFilter pattern generator 
*/
void FilterPatternInit(PatternInfo *pPatternInfo, int w)
{
	int idx = 0;
	int pattern_size = 2 * w * (w + 1);

	// Calculate from distance 1 to window_radius
	for (int r = 1; r < w + 1; ++r)
	{
		for (int dp = 0; dp < r + 1; ++dp)
		{
			for (int pSign = -1; pSign < 2; pSign += 2) // sign = -1, 1
			{
				int p = pSign*dp;


				for (int qSign = -1; qSign < 2; qSign += 2) // sign = -1, 1
				{
					int q = qSign * (r - dp);

					// decide fromLocal (the parent pixel t-1)

					if (p * q == 0) // on the x or y axis
					{
						if (p == 0)
						{
							pPatternInfo[idx].fromLocal[0] = p;
							pPatternInfo[idx].fromLocal[1] = q - qSign;
						}
						else // q == 0
						{
							pPatternInfo[idx].fromLocal[0] = p - pSign;
							pPatternInfo[idx].fromLocal[1] = q;
						}
					}
					else // p*q != 0 (other pixels)
					{
						// if r is odd -> p , else -> q
						if (r % 2 != 0)
						{
							pPatternInfo[idx].fromLocal[0] = p;
							pPatternInfo[idx].fromLocal[1] = q - qSign;
						}
						else
						{
							pPatternInfo[idx].fromLocal[0] = p - pSign;
							pPatternInfo[idx].fromLocal[1] = q;
						}
					}

					pPatternInfo[idx].toLocal[0] = p;
					pPatternInfo[idx].toLocal[1] = q;

					idx++;

					// ensure pixels on the axis is calculated only one time 
					if (q == 0)
						break;

				}

				// ensure pixels on the axis is calculated only one time 
				if (p == 0)
					break;

			}

		}// end dp

	}// end r

}


/**
* Main function of the propagation filter
* @param OutPtr pointer to an array where result will be stored
* @param RPtr pointer to an array containing reference matrix
* @param APtr pointer to an array containing input matrix
* @param nCols the number of columns of R and A
* @param nRows the number of rows of R and A
* @param r window radius
* @param sigma_d
* @param sigma_r
* @param zR the number of channels of R
* @param zA the number of channels of A
*/

void pfilter(unsigned char* OutPtr, const unsigned char* RPtr, const unsigned char* APtr, int nCols, int nRows,
	int r, float sigma_d, float sigma_r, int zR, int zA)
{

# ifdef OPENMP
	// set number of threads used by openmp
	// omp_set_num_threads(floor(omp_get_num_procs()/2));
	omp_set_num_threads(floor(omp_get_num_procs() - 2));
# endif

	PatternInfo *pPatternInfo = NULL;

	pPatternInfo = new PatternInfo[2 * r * (r + 1)]; 

	FilterPatternInit(pPatternInfo, r);
//#pragma omp parallel default(shared)
	{
		unsigned char* resultPixelPtr = new unsigned char[zA];



//#pragma omp for schedule(runtime)
		for (int i = 0; i < nRows; ++i)
		{
			for (int j = 0; j < nCols; ++j)
			{
				int point[2] = { i, j };
				calculateWeightedAverage(resultPixelPtr, RPtr, APtr, pPatternInfo,
					r, sigma_d, sigma_r, point, nRows, nCols, zR, zA);
				for (int z = 0; z < zA; ++z)
				{
					getRef3D(OutPtr, i, j, z, nRows, nCols, zA) = resultPixelPtr[z];
				}
			}
		}
		delete[] resultPixelPtr;
	}

	if (pPatternInfo)delete[]pPatternInfo;

	return;

}

#endif


void _propagated_filter_singlechannel(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int r, float sigma_s, float sigma_r)
{
	int c = 1;
	float _scale = 255 * 255;	//if regularization, _scale = 1; if no regularization, _scale =  255*255

	pfilter(dst, src, src, w, h, r, (float)(sigma_s * _scale), (float)(sigma_r* _scale), c, c);

}

int propagated_filter_singlechannel(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int c, int r, float sigma_s, float sigma_r)
{
	if (src == NULL || guidance == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || (c != 1))
		return 1;

	_propagated_filter_singlechannel(src, guidance, dst, w, h, r, sigma_s, sigma_r);

	return 0;
}

void _propagated_filter_color(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int r, float sigma_s, float sigma_r)
{
	int c = 3;
	float _scale = 255 * 255;	//if regularization, _scale = 1; if no regularization, _scale =  255*255

	pfilter(dst, src, src, w, h, r, (float)(sigma_s * _scale), (float)(sigma_r* _scale), c, c);

}


int propagated_filter_color(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int c, int r, float sigma_s, float sigma_r)
{
	if (src == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || (c != 3))
		return 1;

	_propagated_filter_color(src, guidance, dst, w, h, r, sigma_s, sigma_r);

	return 0;
}

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
@return：										  0:ok; 1:error
@brief：

*/
int PropagatedFilter(unsigned char *src, unsigned char *guidance, unsigned char *dst,
	int w, int h, int c, int r, float sigma_s, float sigma_r)
{
	int _code = 0;
	if (src == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || !(c == 1 || c == 3))
		return 1;

	switch (c)
	{
	case 1:
		_code = propagated_filter_singlechannel(src, guidance, dst, w, h, c, r, sigma_s, sigma_r);
		break;

	case 3:
		_code = propagated_filter_color(src, guidance, dst, w, h, c, r, sigma_s, sigma_r);
		break;

	default:
		_code = 1;
		break;
	}


	return _code;
}



