# some mainstream image filter: boxfilter, fast guided filter, fast bilateral filter, permutohedral bilateral filter and propagated filter

- boxfilter 
- fast guided filter
- fast bilateral filter
- permutohedral bilateral filter
- propagated filter

**All original dependencies have been removed. Code could be run  independently.**

## BOXFILTER 

    @param src         				image,single channel.
    
    @param dst      				output image,single channel.
    
    @param w           				width of image.
    
    @param h           				height of image.
    
    @param c      					channel of image, only c = 1.
    
    @param r      					local window radius.
    
    @return            				0:ok; 1:error
    int BoxfilterFilter(unsigned char *src, unsigned char *dst, int w, int h, int c, int r);

## Fast Guided Filter

    @param src         				image,single channel.
    
    @param guidance           		guided image,single channel.
    
    @param dst      				output image,single channel.
    
    @param w           				width of image.
    
    @param h           				height of image.
    
    @param c      					channel of image, only c = 1.
    
    @param r      					local window radius.
    
    @param rp      					regularization parameter:eps.
    
    @param sr      					subsampling ratio, sr>1:downscale; 0<sr<1:upscale.
    
    @return            				0:ok; 1:error
    int FastGuidedFilter(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, int c, int r, float rp, float sr);

## FastBilateralFilter

    @param src         				input image.
    
    @param guidance           		guided image,single channel, only single channel is valid; invalid parameter in three channels.
    
    @param dst      				output image.
    
    @param w           				width of image.
    
    @param h           				height of image.
    
    @param c      					channel of image, only c = 1 or c = 3.
    
    @param sigma_s      			filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
    
    @param sigma_r      			filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
    
    @return            				0:ok; 1:error
    int FastBilateralFilter(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, int c, float sigma_s, float sigma_r);

## PermutohedralBilateralFilter

    @param src         				input image.
    
    @param guidance           		guided image. 
    
    @param dst      				output image.
    
    @param w           				width of image.
    
    @param h           				height of image.
    
    @param c      					channel of image, only c = 1 or c = 3.
    
    @param sigma_s      			filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
    
    @param sigma_r      			filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
    
    @return            				0:ok; 1:error
    int PermutohedralBilateralFilter(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, int c, float sigma_s, float sigma_r);

## PropagatedFilter

    @param src         				input image.
    
    @param guidance           		guided image. 
    
    @param dst      				output image.
    
    @param w           				width of image.
    
    @param h           				height of image.
    
    @param c      					channel of image, only c = 1 or c = 3.
    
    @param r      					local window radius.
    
    @param sigma_s      			filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.
    
    @param sigma_r      			filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.
    
    @return            				0:ok; 1:error
    int PropagatedFilter(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, int c, int r, float sigma_s, float sigma_r);

## Reference

- [OpenCV](https://opencv.org/)
- [Guided Image Filter](http://kaiminghe.com/eccv10/)
- [Propagated Image Filtering](http://mml.citi.sinica.edu.tw/papers/CVPR_2015_Chang.pdf)
- [Propagated Image Filtering code](http://mml.citi.sinica.edu.tw/papers/PFilter_direct_color_MEX.zip)
- [Fast High-Dimensional Filtering Using the Permutohedral Lattice](http://graphics.stanford.edu/papers/permutohedral/)