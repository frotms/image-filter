#ifndef __IMAGE_H
#define __IMAGE_H

#include <string.h>

// Some basic PBFimage classes.
// PBFimage is a reference counting pointer class that refers to an array of floats
// PBFwindow is a subregion of an PBFimage

class PBFwindow {
  public:
    PBFwindow() {
	xstride = ystride = tstride = width = height = frames = channels = 0;
	data = NULL;
    }
    
    PBFwindow(PBFwindow im, int mint_, int minx_, int miny_, int frames_, int width_, int height_) {
	int mint = MAX(0, mint_);
	int maxt = MIN(im.frames, mint_ + frames_);
	int minx = MAX(0, minx_);
	int maxx = MIN(im.width, minx_ + width_);
	int miny = MAX(0, miny_);
	int maxy = MIN(im.height, miny_ + height_);

	xstride = im.xstride;
	ystride = im.ystride;
	tstride = im.tstride;

	width = maxx - minx;
	height = maxy - miny;
	frames = maxt - mint;
	channels = im.channels;

	data = im.data + mint * tstride + miny * ystride + minx * xstride;
    }

    float *operator()(int t, int x, int y) {
	return data + t * tstride + x * xstride + y * ystride;
    }

    float *operator()(int x, int y) {
	return data + x * xstride + y * ystride;
    }

    float *operator()(int x) {
	return data + x * xstride;
    }


    int width, height, frames, channels;
    int xstride, ystride, tstride;
    float *data;    

};

class PBFimage : public PBFwindow {
  public:
    PBFimage() : refCount(NULL) {
	width = frames = height = channels = 0;
	xstride = ystride = tstride = 0;
	data = NULL;
    }
    
    PBFimage(int frames_, int width_, int height_, int channels_, const float *data_ = NULL)
	{
		frames = frames_;
		width = width_;
		height = height_;
		channels = channels_;

		long long memory = ((long long)frames_ *
			(long long)height_ *
			(long long)width_ *
			(long long)channels_);

        data = new float[memory];
        if (!data_) memset(data, 0, memory * sizeof(float));
        else memcpy(data, data_, memory * sizeof(float));
        
        xstride = channels;
        ystride = xstride * width;
        tstride = ystride * height;
        refCount = new int;
        *refCount = 1;
        
        //printf("Making new image "); 
        //debug();
    }
    
    // does not copy data

    PBFimage &operator=(const PBFimage &im) {
	if (refCount) {
	    refCount[0]--;
	    if (*refCount <= 0) {
		delete refCount;
		delete[] data;
	    }
	}

        width = im.width;
        height = im.height;
        channels = im.channels;
        frames = im.frames;
        
        data = im.data;
        
        xstride = channels;
        ystride = xstride * width;
        tstride = ystride * height;	
        
        refCount = im.refCount;       
        if (refCount) refCount[0]++;
        
        return *this;
    }
    
    PBFimage(const PBFimage &im) {
        width = im.width;
        height = im.height;
        channels = im.channels;
        frames = im.frames;
        
        data = im.data;       
        xstride = channels;
        ystride = xstride * width;
        tstride = ystride * height;	
        
        refCount = im.refCount;        
        if (refCount) refCount[0]++;       
    }
    
    // copies data from the window
    PBFimage(PBFwindow im) {
        width = im.width;
        height = im.height;
        channels = im.channels;
        frames = im.frames;
        
        xstride = channels;
        ystride = xstride * width;
        tstride = ystride * height;	
        
        refCount = new int;
        *refCount = 1;
	long long memory = ((long long)width *
			    (long long)height *
			    (long long)channels *
			    (long long)frames);
        data = new float[memory];
        
        for (int t = 0; t < frames; t++) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    for (int c = 0; c < channels; c++) {
                        (*this)(t, x, y)[c] = im(t, x, y)[c];
                    }
                }
            }
        }
        
    }
    
    // makes a new copy of this image
    PBFimage copy() {
        return PBFimage(*((PBFwindow *)this));
    }
    
    ~PBFimage();

    int *refCount;
    
  protected:
    PBFimage &operator=(PBFwindow im) {
        return *this;
    }
};

#endif
