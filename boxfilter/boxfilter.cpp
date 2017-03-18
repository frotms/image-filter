#include "boxfilter.h"

#include <memory.h>
//#include<string.h> 
#include <stdlib.h>
#include <math.h>
//#include "omp.h"


#if use_sse_version
#include <pmmintrin.h>    //sse
#endif

enum RETVAL
{
	RET_OK,
	RET_ERR_OUTOFMEMORY,
	RET_ERR_NULLREFERENCE,
	RET_ERR_ARGUMENTOUTOFRANGE,
	RET_ERR_PARAMISMATCH,
	RET_ERR_DIVIDEBYZERO,
	RET_ERR_NOTSUPPORTED,
	RET_ERR_UNKNOWN
};

enum DEPTH
{
	DEPTH_8U, 				//	unsigned char
	DEPTH_8S, 				//	char
	DEPTH_16S,				//	short
	DEPTH_32S,				//  int
	DEPTH_32F,				//	float
	DEPTH_64F				//	double
};

enum EdgeMode
{
	Tile = 0,				//	duplicate
	Smear = 1				//	mirror
};

struct TMatrix
{
	int Width;
	int Height;
	int Channel;
	int Depth;
	int WidthStep;
	int RowAligned;
	unsigned char *Data;
};

struct Complex
{
	float Real;
	float Imag;
};
//union Cv32suf
//{
//	int i;
//	unsigned u;
//	float f;
//};

#define WidthBytes(bytes) (((bytes * 8) + 31) / 32 * 4)

//core
#if 1


int GetElementSize(int Depth)
{
	int Size;
	switch (Depth)
	{
	case DEPTH_8U:
		Size = sizeof(unsigned char);
		break;
	case DEPTH_8S:
		Size = sizeof(char);
		break;
	case DEPTH_16S:
		Size = sizeof(short);
		break;
	case DEPTH_32S:
		Size = sizeof(int);
		break;
	case DEPTH_32F:
		Size = sizeof(float);
		break;
	case DEPTH_64F:
		Size = sizeof(double);
		break;
	default:
		Size = 0;
		break;
	}
	return Size;
}


void *AllocMemory(unsigned int Size, bool ZeroMemory)	
{
	void *Ptr = _aligned_malloc(Size, 32);					
	if (Ptr != NULL && ZeroMemory == true)
		memset(Ptr, 0, Size);
	return Ptr;
}


void FreeMemory(void *Ptr)
{
	if (Ptr != NULL) _aligned_free(Ptr);		
}


TMatrix *CreateMatrix(int Width, int Height, int Depth, int Channel, int RowAligned)
{
	if (Width < 1 || Height < 1) return NULL;
	if (Depth < DEPTH_8U && Depth > DEPTH_64F) return NULL;

	TMatrix * Matrix = NULL;
	Matrix = (TMatrix *)AllocMemory(sizeof(TMatrix), false);		
	Matrix->Width = Width;
	Matrix->Height = Height;
	Matrix->Depth = Depth;
	Matrix->Channel = Channel;
	Matrix->RowAligned = RowAligned;
	if (RowAligned)
		Matrix->WidthStep = WidthBytes(Width * Channel * GetElementSize(Depth));
	else
		Matrix->WidthStep = Width * Channel * GetElementSize(Depth);
	Matrix->Data = (unsigned char *)AllocMemory(Matrix->Height * Matrix->WidthStep, false);		//	用了True在大循环里换慢很多
	if (Matrix->Data == NULL)
	{
		FreeMemory(Matrix);
		return NULL;
	}
	return Matrix;
}


void FreeMatrix(TMatrix *Matrix)
{
	if (Matrix == NULL) return;
	if (Matrix->Data == NULL)
		FreeMemory(Matrix);
	else
	{
		FreeMemory(Matrix->Data);			
		FreeMemory(Matrix);
	}
}


TMatrix *Clone(TMatrix *Src)
{
	if (Src == NULL || Src->Data == NULL) return NULL;
	TMatrix *Dest = CreateMatrix(Src->Width, Src->Height, Src->Depth, Src->Channel, Src->RowAligned);
	if (Dest != NULL)
	{
		memcpy(Dest->Data, Src->Data, Src->Height * Src->WidthStep);
		return Dest;
	}
	return NULL;
}


#endif

//Utility
#if 1

void SplitRGBA(TMatrix *Src, TMatrix *&Blue, TMatrix *&Green, TMatrix *&Red, TMatrix *&Alpha)
{

	Blue = CreateMatrix(Src->Width, Src->Height, Src->Depth, 1, true);
	Green = CreateMatrix(Src->Width, Src->Height, Src->Depth, 1, true);
	Red = CreateMatrix(Src->Width, Src->Height, Src->Depth, 1, true);
	if (Src->Channel == 4)	Alpha = CreateMatrix(Src->Width, Src->Height, Src->Depth, 1, true);

	int X, Y, Block, Width = Src->Width, Height = Src->Height;
	unsigned char *LinePS = NULL, *LinePB = NULL, *LinePG = NULL, *LinePR = NULL, *LinePA = NULL;
	const int BlockSize = 8;
	Block = Width / BlockSize;

	if (Src->Channel == 3)
	{
		for (Y = 0; Y < Height; Y++)
		{
			LinePS = Src->Data + Y * Src->WidthStep;
			LinePB = Blue->Data + Y * Blue->WidthStep;
			LinePG = Green->Data + Y * Green->WidthStep;
			LinePR = Red->Data + Y * Red->WidthStep;
			for (X = 0; X < Block * BlockSize; X += BlockSize)
			{
				LinePB[0] = LinePS[0];		LinePG[0] = LinePS[1];		LinePR[0] = LinePS[2];
				LinePB[1] = LinePS[3];		LinePG[1] = LinePS[4];		LinePR[1] = LinePS[5];
				LinePB[2] = LinePS[6];		LinePG[2] = LinePS[7];		LinePR[2] = LinePS[8];
				LinePB[3] = LinePS[9];		LinePG[3] = LinePS[10];		LinePR[3] = LinePS[11];
				LinePB[4] = LinePS[12];		LinePG[4] = LinePS[13];		LinePR[4] = LinePS[14];
				LinePB[5] = LinePS[15];		LinePG[5] = LinePS[16];		LinePR[5] = LinePS[17];
				LinePB[6] = LinePS[18];		LinePG[6] = LinePS[19];		LinePR[6] = LinePS[20];
				LinePB[7] = LinePS[21];		LinePG[7] = LinePS[22];		LinePR[7] = LinePS[23];
				LinePB += 8;				LinePG += 8;				LinePR += 8;				LinePS += 24;
			}
			for (; X < Width; X++)
			{
				LinePB[0] = LinePS[0];		LinePG[0] = LinePS[1];		LinePR[0] = LinePS[2];
				LinePB++;					LinePG++;					LinePR++;					LinePS += 3;
			}
		}
	}
	else if (Src->Channel == 4)
	{
		for (Y = 0; Y < Height; Y++)
		{
			LinePS = Src->Data + Y * Src->WidthStep;
			LinePB = Blue->Data + Y * Blue->WidthStep;
			LinePG = Green->Data + Y * Green->WidthStep;
			LinePR = Red->Data + Y * Red->WidthStep;
			LinePA = Alpha->Data + Y * Alpha->WidthStep;
			for (X = 0; X < Block * BlockSize; X += BlockSize)
			{
				LinePB[0] = LinePS[0];		LinePG[0] = LinePS[1];		LinePR[0] = LinePS[2];		LinePA[0] = LinePS[3];
				LinePB[1] = LinePS[4];		LinePG[1] = LinePS[5];		LinePR[1] = LinePS[6];		LinePA[1] = LinePS[7];
				LinePB[2] = LinePS[8];		LinePG[2] = LinePS[9];		LinePR[2] = LinePS[10];		LinePA[2] = LinePS[11];
				LinePB[3] = LinePS[12];		LinePG[3] = LinePS[13];		LinePR[3] = LinePS[14];		LinePA[3] = LinePS[15];
				LinePB[4] = LinePS[16];		LinePG[4] = LinePS[17];		LinePR[4] = LinePS[18];		LinePA[4] = LinePS[19];
				LinePB[5] = LinePS[20];		LinePG[5] = LinePS[21];		LinePR[5] = LinePS[22];		LinePA[5] = LinePS[23];
				LinePB[6] = LinePS[24];		LinePG[6] = LinePS[25];		LinePR[6] = LinePS[26];		LinePA[6] = LinePS[27];
				LinePB[7] = LinePS[28];		LinePG[7] = LinePS[29];		LinePR[7] = LinePS[30];		LinePA[7] = LinePS[31];
				LinePB += 8;				LinePG += 8;				LinePR += 8;				LinePA += 8;				LinePS += 32;
			}
			for (; X < Width; X++)
			{
				LinePB[0] = LinePS[0];		LinePG[0] = LinePS[1];		LinePR[0] = LinePS[2];		LinePA[0] = LinePS[3];
				LinePB++;					LinePG++;					LinePR++;					LinePA++;					LinePS += 4;
			}
		}
	}
}

void CombineRGBA(TMatrix *Dest, TMatrix *Blue, TMatrix *Green, TMatrix *Red, TMatrix *Alpha)
{
	int X, Y, Block, Width = Dest->Width, Height = Dest->Height;
	unsigned char *LinePD = NULL, *LinePB = NULL, *LinePG = NULL, *LinePR = NULL, *LinePA = NULL;
	const int BlockSize = 8;
	Block = Width / BlockSize;

	if (Dest->Channel == 3)
	{
		for (Y = 0; Y < Height; Y++)
		{
			LinePD = Dest->Data + Y * Dest->WidthStep;
			LinePB = Blue->Data + Y * Blue->WidthStep;
			LinePG = Green->Data + Y * Green->WidthStep;
			LinePR = Red->Data + Y * Red->WidthStep;
			for (X = 0; X < Block * BlockSize; X += BlockSize)
			{
				LinePD[0] = LinePB[0];		LinePD[1] = LinePG[0];		LinePD[2] = LinePR[0];
				LinePD[3] = LinePB[1];		LinePD[4] = LinePG[1];		LinePD[5] = LinePR[1];
				LinePD[6] = LinePB[2];		LinePD[7] = LinePG[2];		LinePD[8] = LinePR[2];
				LinePD[9] = LinePB[3];		LinePD[10] = LinePG[3];		LinePD[11] = LinePR[3];
				LinePD[12] = LinePB[4];		LinePD[13] = LinePG[4];		LinePD[14] = LinePR[4];
				LinePD[15] = LinePB[5];		LinePD[16] = LinePG[5];		LinePD[17] = LinePR[5];
				LinePD[18] = LinePB[6];		LinePD[19] = LinePG[6];		LinePD[20] = LinePR[6];
				LinePD[21] = LinePB[7];		LinePD[22] = LinePG[7];		LinePD[23] = LinePR[7];
				LinePB += 8;				LinePG += 8;				LinePR += 8;				LinePD += 24;
			}
			for (; X < Width; X++)
			{
				LinePD[0] = LinePB[0];		LinePD[1] = LinePG[0];		LinePD[2] = LinePR[0];
				LinePB++;					LinePG++;					LinePR++;					LinePD += 3;
			}
		}
	}
	else if (Dest->Channel == 4)
	{
		for (Y = 0; Y < Height; Y++)
		{
			LinePD = Dest->Data + Y * Dest->WidthStep;
			LinePB = Blue->Data + Y * Blue->WidthStep;
			LinePG = Green->Data + Y * Green->WidthStep;
			LinePR = Red->Data + Y * Red->WidthStep;
			LinePA = Alpha->Data + Y * Alpha->WidthStep;
			for (X = 0; X < Block * BlockSize; X += BlockSize)
			{
				LinePD[0] = LinePB[0];		LinePD[1] = LinePG[0];		LinePD[2] = LinePR[0];		LinePD[3] = LinePA[0];
				LinePD[4] = LinePB[1];		LinePD[5] = LinePG[1];		LinePD[6] = LinePR[1];		LinePD[7] = LinePA[1];
				LinePD[8] = LinePB[2];		LinePD[9] = LinePG[2];		LinePD[10] = LinePR[2];		LinePD[11] = LinePA[2];
				LinePD[12] = LinePB[3];		LinePD[13] = LinePG[3];		LinePD[14] = LinePR[3];		LinePD[15] = LinePA[3];
				LinePD[16] = LinePB[4];		LinePD[17] = LinePG[4];		LinePD[18] = LinePR[4];		LinePD[19] = LinePA[4];
				LinePD[20] = LinePB[5];		LinePD[21] = LinePG[5];		LinePD[22] = LinePR[5];		LinePD[23] = LinePA[5];
				LinePD[24] = LinePB[6];		LinePD[25] = LinePG[6];		LinePD[26] = LinePR[6];		LinePD[27] = LinePA[6];
				LinePD[28] = LinePB[7];		LinePD[29] = LinePG[7];		LinePD[30] = LinePR[7];		LinePD[31] = LinePA[7];
				LinePB += 8;				LinePG += 8;				LinePR += 8;				LinePA += 8;				LinePD += 32;
			}
			for (; X < Width; X++)
			{
				LinePD[0] = LinePB[0];		LinePD[1] = LinePG[0];		LinePD[2] = LinePR[0];		LinePD[3] = LinePA[0];
				LinePB++;					LinePG++;					LinePD++;					LinePA++;					LinePD += 4;
				X++;
			}
		}
	}
}

int *GetExpandPos(int Length, int Left, int Right, EdgeMode Edge)
{
	int X, PosX;
	int *Pos = NULL;
	if (Left < 0 || Length < 0 || Right < 0) return NULL;
	Pos = (int *)AllocMemory((Left + Length + Right) * sizeof(int), false);
	if (Pos == NULL) return NULL;

	for (X = -Left; X < Length + Right; X++)
	{
		if (X < 0)
		{
			if (Edge == EdgeMode::Tile)
				Pos[X + Left] = 0;
			else
			{
				PosX = -X;
				while (PosX >= Length) PosX -= Length;
				Pos[X + Left] = PosX;
			}
		}
		else if (X >= Length)
		{
			if (Edge == EdgeMode::Tile)
				Pos[X + Left] = Length - 1;
			else
			{
				PosX = Length - (X - Length + 2);
				while (PosX < 0) PosX += Length;
				Pos[X + Left] = PosX;
			}
		}
		else
		{
			Pos[X + Left] = X;
		}
	}
	return Pos;
}


TMatrix *GetExpandMatrix(TMatrix *Src, int Left, int Top, int Right, int Bottom, EdgeMode Edge)
{
	int X, Y, SrcW, SrcH, DstW, DstH, ElementSize;
	TMatrix *Dest = NULL;

	int *RowPos = NULL;
	int *ColPos = NULL;

	unsigned char *LinePS = NULL;
	unsigned char *LinePD = NULL;

	if (Src == NULL || Src->Data == NULL) return NULL;
	if (Left < 0 || Right < 0 || Top < 0 || Bottom < 0) return NULL;

	SrcW = Src->Width;				SrcH = Src->Height;
	DstW = SrcW + Left + Right;		DstH = SrcH + Top + Bottom;
	ElementSize = Src->Channel * GetElementSize(Src->Depth);

	Dest = CreateMatrix(DstW, DstH, Src->Depth, Src->Channel, Src->RowAligned);
	if (Dest != NULL)
	{
		RowPos = GetExpandPos(SrcW, Left, Right, Edge);
		ColPos = GetExpandPos(SrcH, Top, Bottom, Edge);
		for (Y = 0; Y < SrcH; Y++)
		{
			LinePS = Src->Data + Y * Src->WidthStep;
			LinePD = Dest->Data + (Y + Top) * Dest->WidthStep;
			for (X = 0; X < Left; X++)
				memcpy(LinePD + X * ElementSize, LinePS + RowPos[X] * ElementSize, ElementSize);
			memcpy(LinePD + Left * ElementSize, LinePS, SrcW * ElementSize);
			for (X = Left + SrcW; X < Left + SrcW + Right; X++)
				memcpy(LinePD + X * ElementSize, LinePS + RowPos[X] * ElementSize, ElementSize);
		}
		for (Y = 0; Y < Top; Y++)
			memcpy(Dest->Data + Y * Dest->WidthStep, Dest->Data + (Top + ColPos[Y]) * Dest->WidthStep, Dest->WidthStep);	//	行行直接拷贝

		for (Y = Top + SrcH; Y < Top + SrcH + Bottom; Y++)
			memcpy(Dest->Data + Y * Dest->WidthStep, Dest->Data + (Top + ColPos[Y]) * Dest->WidthStep, Dest->WidthStep);

		FreeMemory(RowPos);
		FreeMemory(ColPos);
	}
	return Dest;
}

unsigned char ClampToByte(int Value)
{
	return ((Value | ((signed int)(255 - Value) >> 31)) & ~((signed int)Value >> 31));
}


#endif

//sse version
#if 1
#if use_sse_version

void __stdcall BoxBlurSSE(TMatrix *Src, TMatrix *Dest, int Radius, EdgeMode Edge)
{

	int X, Y, Z, Width, Height, Channel, Size, Index, X3, Z3;
	int Value, ValueB, ValueG, ValueR;
	int *RowPos, *ColPos, *ColSum, *Diff;

	float Scale;
	int Amount;

	unsigned char *RowData = NULL;
	TMatrix *Sum = NULL;

	unsigned char *LinePS = NULL;
	int *LinePD = NULL;

	unsigned char *AddPos = NULL;
	unsigned char *SubPos = NULL;
	int *AddPos_int = NULL;
	int *SubPos_int = NULL;

	Width = Src->Width, Height = Src->Height, Channel = Src->Channel, Size = 2 * Radius + 1;
	Scale = 1.0 / (Size * Size);

	Amount = Size * Size;

	RowPos = GetExpandPos(Width, Radius, Radius, Edge);
	ColPos = GetExpandPos(Height, Radius, Radius, Edge);
	ColSum = (int *)AllocMemory(Width * Channel * sizeof(int), true);
	Diff = (int *)AllocMemory((Width - 1) * Channel * sizeof(int), true);
	RowData = (unsigned char *)AllocMemory((Width + 2 * Radius) * Channel, true);
	Sum = CreateMatrix(Width, Height, DEPTH_32S, Channel, true);

	for (Y = 0; Y < Height; Y++)					//	水平方向的耗时比垂直方向上的大
	{
		LinePS = Src->Data + Y * Src->WidthStep;
		LinePD = (int *)(Sum->Data + Y * Sum->WidthStep);

		//	拷贝一行数据及边缘部分部分到临时的缓冲区中
		if (Channel == 1)
		{
			for (X = 0; X < Radius; X++)
				RowData[X] = LinePS[RowPos[X]];
			memcpy(RowData + Radius, LinePS, Width);
			for (X = Radius + Width; X < Radius + Width + Radius; X++)
				RowData[X] = LinePS[RowPos[X]];
		}
		else if (Channel == 3)
		{
			for (X = 0; X < Radius; X++)
			{
				X3 = X * 3;
				Index = RowPos[X] * 3;
				RowData[X3] = LinePS[Index];
				RowData[X3 + 1] = LinePS[Index + 1];
				RowData[X3 + 2] = LinePS[Index + 2];
			}
			memcpy(RowData + Radius * 3, LinePS, Width * 3);
			for (X = Radius + Width; X < Radius + Width + Radius; X++)
			{
				X3 = X * 3;
				Index = RowPos[X] * 3;
				RowData[X3 + 0] = LinePS[Index + 0];
				RowData[X3 + 1] = LinePS[Index + 1];
				RowData[X3 + 2] = LinePS[Index + 2];
			}
		}

		AddPos = RowData + Size * Channel;
		SubPos = RowData;
		X = 0;					//	注意这个赋值在下面的循环外部，这可以避免当Width<8时第二个for循环循环变量未初始化			
		__m128i Zero = _mm_setzero_si128();
		for (; X <= (Width - 1) * Channel - 8; X += 8)
		{
			__m128i Add = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i const *)(AddPos + X)), Zero);
			__m128i Sub = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i const *)(SubPos + X)), Zero);
			_mm_store_si128((__m128i *)(Diff + X + 0), _mm_sub_epi32(_mm_unpacklo_epi16(Add, Zero), _mm_unpacklo_epi16(Sub, Zero)));		//	由于采用了_aligned_malloc函数分配内存，可是使用_mm_store_si128
			_mm_store_si128((__m128i *)(Diff + X + 4), _mm_sub_epi32(_mm_unpackhi_epi16(Add, Zero), _mm_unpackhi_epi16(Sub, Zero)));
		}
		for (; X < (Width - 1) * Channel; X++)
			Diff[X] = AddPos[X] - SubPos[X];

		//	第一个点要特殊处理
		if (Channel == 1)
		{
			for (Z = 0, Value = 0; Z < Size; Z++)	Value += RowData[Z];
			LinePD[0] = Value;

			for (X = 1; X < Width; X++)
			{
				Value += Diff[X - 1];
				LinePD[X] = Value;
			}
		}
		else if (Channel == 3)
		{
			for (Z = 0, ValueB = ValueG = ValueR = 0; Z < Size; Z++)
			{
				Z3 = Z * 3;
				ValueB += RowData[Z3 + 0];
				ValueG += RowData[Z3 + 1];
				ValueR += RowData[Z3 + 2];
			}
			LinePD[0] = ValueB;	LinePD[1] = ValueG;	LinePD[2] = ValueR;

			for (X = 1; X < Width; X++)
			{
				Index = X * 3;
				ValueB += Diff[Index - 3];		LinePD[Index + 0] = ValueB;
				ValueG += Diff[Index - 2];		LinePD[Index + 1] = ValueG;
				ValueR += Diff[Index - 1];		LinePD[Index + 2] = ValueR;
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////

	for (Y = 0; Y < Size - 1; Y++)			//	注意没有最后一项哦						//	这里的耗时只占整个的15%左右
	{
		X = 0;
		LinePD = (int *)(Sum->Data + ColPos[Y] * Sum->WidthStep);
		for (; X <= Width * Channel - 4; X += 4)
		{
			__m128i SumP = _mm_load_si128((const __m128i*)(ColSum + X));
			__m128i SrcP = _mm_loadu_si128((const __m128i*)(LinePD + X));
			_mm_store_si128((__m128i *)(ColSum + X), _mm_add_epi32(SumP, SrcP));
		}

		for (; X < Width * Channel; X++)	ColSum[X] += LinePD[X];
	}

	for (Y = 0; Y < Height; Y++)
	{
		LinePS = Dest->Data + Y * Dest->WidthStep;
		AddPos_int = (int*)(Sum->Data + ColPos[Y + Size - 1] * Sum->WidthStep);
		SubPos_int = (int*)(Sum->Data + ColPos[Y] * Sum->WidthStep);

		X = 0;
		const __m128 Inv = _mm_set1_ps(Scale);
		for (; X <= Width * Channel - 8; X += 8)
		{
			__m128i Sub1 = _mm_loadu_si128((const __m128i*)(SubPos_int + X + 0));
			__m128i Sub2 = _mm_loadu_si128((const __m128i*)(SubPos_int + X + 4));

			__m128i Add1 = _mm_loadu_si128((const __m128i*)(AddPos_int + X + 0));
			__m128i Add2 = _mm_loadu_si128((const __m128i*)(AddPos_int + X + 4));
			__m128i Col1 = _mm_load_si128((const __m128i*)(ColSum + X + 0));
			__m128i Col2 = _mm_load_si128((const __m128i*)(ColSum + X + 4));

			__m128i Sum1 = _mm_add_epi32(Col1, Add1);
			__m128i Sum2 = _mm_add_epi32(Col2, Add2);

			__m128i Dest1 = _mm_cvtps_epi32(_mm_mul_ps(Inv, _mm_cvtepi32_ps(Sum1)));
			__m128i Dest2 = _mm_cvtps_epi32(_mm_mul_ps(Inv, _mm_cvtepi32_ps(Sum2)));

			Dest1 = _mm_packs_epi32(Dest1, Dest2);
			_mm_storel_epi64((__m128i *)(LinePS + X), _mm_packus_epi16(Dest1, Dest1));

			_mm_store_si128((__m128i *)(ColSum + X + 0), _mm_sub_epi32(Sum1, Sub1));
			_mm_store_si128((__m128i *)(ColSum + X + 4), _mm_sub_epi32(Sum2, Sub2));
		}
		for (; X < Width * Channel; X++)
		{
			Value = ColSum[X] + AddPos_int[X];
			LinePS[X] = Value * Scale;
			ColSum[X] = Value - SubPos_int[X];
		}
	}

	FreeMemory(RowPos);
	FreeMemory(ColPos);
	FreeMemory(Diff);
	FreeMatrix(Sum);
	FreeMemory(ColSum);
	FreeMemory(RowData);

}
#endif
#endif


//c version
#if 1

void __stdcall BoxBlur(TMatrix *Src, TMatrix *Dest, int Radius, EdgeMode Edge)
{
	int X, Y, Z, Width, Height, Channel, Index, X3, Z3;
	int Value, ValueB, ValueG, ValueR;
	int *RowPos, *ColPos, *ColSum, *Diff;

	int Size, Amount, HalfAmount;

	unsigned char *RowData = NULL;
	TMatrix *Sum = NULL;

	unsigned char *LinePS = NULL;
	int *LinePD = NULL;

	unsigned char *AddPos = NULL;
	unsigned char *SubPos = NULL;
	int *AddPos_int = NULL;
	int *SubPos_int = NULL;

	Width = Src->Width, Height = Src->Height, Channel = Src->Channel;
	Size = 2 * Radius + 1, Amount = Size * Size, HalfAmount = Amount / 2;

	RowPos = GetExpandPos(Width, Radius, Radius, Edge);
	ColPos = GetExpandPos(Height, Radius, Radius, Edge);
	ColSum = (int *)AllocMemory(Width * Channel * sizeof(int), true);
	Diff = (int *)AllocMemory((Width - 1) * Channel * sizeof(int), true);
	RowData = (unsigned char *)AllocMemory((Width + 2 * Radius) * Channel, true);
	Sum = CreateMatrix(Width, Height, DEPTH_32S, Channel, true);

	for (Y = 0; Y < Height; Y++)
	{
		LinePS = Src->Data + Y * Src->WidthStep;
		LinePD = (int *)(Sum->Data + Y * Sum->WidthStep);


		if (Channel == 1)
		{
			for (X = 0; X < Radius; X++)
				RowData[X] = LinePS[RowPos[X]];
			memcpy(RowData + Radius, LinePS, Width);
			for (X = Radius + Width; X < Radius + Width + Radius; X++)
				RowData[X] = LinePS[RowPos[X]];
		}
		else if (Channel == 3)
		{
			for (X = 0; X < Radius; X++)
			{
				X3 = X * 3;
				Index = RowPos[X] * 3;
				RowData[X3] = LinePS[Index];
				RowData[X3 + 1] = LinePS[Index + 1];
				RowData[X3 + 2] = LinePS[Index + 2];
			}
			memcpy(RowData + Radius * 3, LinePS, Width * 3);
			for (X = Radius + Width; X < Radius + Width + Radius; X++)
			{
				X3 = X * 3;
				Index = RowPos[X] * 3;
				RowData[X3 + 0] = LinePS[Index + 0];
				RowData[X3 + 1] = LinePS[Index + 1];
				RowData[X3 + 2] = LinePS[Index + 2];
			}
		}

		AddPos = RowData + Size * Channel;
		SubPos = RowData;

		for (X = 0; X < (Width - 1) * Channel; X++)
			Diff[X] = AddPos[X] - SubPos[X];

		if (Channel == 1)
		{
			for (Z = 0, Value = 0; Z < Size; Z++)	Value += RowData[Z];
			LinePD[0] = Value;

			for (X = 1; X < Width; X++)
			{
				Value += Diff[X - 1];	LinePD[X] = Value;
			}
		}
		else if (Channel == 3)
		{
			for (Z = 0, ValueB = ValueG = ValueR = 0; Z < Size; Z++)
			{
				Z3 = Z * 3;
				ValueB += RowData[Z3 + 0];
				ValueG += RowData[Z3 + 1];
				ValueR += RowData[Z3 + 2];
			}
			LinePD[0] = ValueB;	LinePD[1] = ValueG;	LinePD[2] = ValueR;

			for (X = 1; X < Width; X++)
			{
				Index = X * 3;
				ValueB += Diff[Index - 3];		LinePD[Index + 0] = ValueB;
				ValueG += Diff[Index - 2];		LinePD[Index + 1] = ValueG;
				ValueR += Diff[Index - 1];		LinePD[Index + 2] = ValueR;
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////

	for (Y = 0; Y < Size - 1; Y++)
	{
		LinePD = (int *)(Sum->Data + ColPos[Y] * Sum->WidthStep);
		for (X = 0; X < Width * Channel; X++)	ColSum[X] += LinePD[X];
	}

	for (Y = 0; Y < Height; Y++)
	{
		LinePS = Dest->Data + Y * Dest->WidthStep;
		AddPos_int = (int*)(Sum->Data + ColPos[Y + Size - 1] * Sum->WidthStep);
		SubPos_int = (int*)(Sum->Data + ColPos[Y] * Sum->WidthStep);

		for (X = 0; X < Width * Channel; X++)
		{
			Value = ColSum[X] + AddPos_int[X];
			LinePS[X] = (Value + HalfAmount) / Amount;
			ColSum[X] = Value - SubPos_int[X];
		}
	}

	FreeMemory(RowPos);
	FreeMemory(ColPos);
	FreeMemory(Diff);
	FreeMatrix(Sum);
	FreeMemory(ColSum);
	FreeMemory(RowData);
}
#endif



void _boxfilter_int8u(unsigned char *src, unsigned char *dst,
	int w, int h, int c, int Radius, EdgeMode Edge)
{
	TMatrix *Src = NULL;
	TMatrix *Dst = NULL;

	Src = (TMatrix *)calloc(sizeof(TMatrix), 1);
	Dst = (TMatrix *)calloc(sizeof(TMatrix), 1);

	Src->Channel = c;
	Src->Width = w;
	Src->Height = h;
	Src->RowAligned = 0;
	Src->WidthStep = w*c;
	Src->Depth = DEPTH_8U;
	Src->Data = src;

	Dst->Channel = c;
	Dst->Width = w;
	Dst->Height = h;
	Dst->RowAligned = 0;
	Dst->WidthStep = w*c;
	Dst->Depth = DEPTH_8U;
	Dst->Data = dst;

#if use_sse_version
	BoxBlurSSE(Src, Dst, Radius, Tile);
#else
	BoxBlur(Src, Dst, Radius, Tile);
#endif

	free(Src);
	free(Dst);

}


/*
@function    BoxfilterFilter
@param       [in]      src:						  input image,single channel
@param       [in/out]  dst:                       output image,single channel
@param       [in]      w:                         width of image
@param       [in]      h:                         height of image
@param       [in]      c:						  channel of image, only c = 1
@param       [in]      r:						  local window radius
@return：										  0:ok; 1:error
@brief：

*/
int BoxfilterFilter(unsigned char *src, unsigned char *dst,
	int w, int h, int c, int r)
{
	EdgeMode edge = Tile;
	if (src == NULL || dst == NULL)
		return 1;

	if (w <= 0 || h <= 0 || (r < 1) || (!(c == 1 || c == 3)))
		return 1;

	_boxfilter_int8u(src, dst, w, h, c, r, edge);


	return 0;
}