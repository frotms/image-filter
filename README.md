# some mainstream image filter: boxfilter, fast guided filter, fast bilateral filter, permutohedral bilateral filter and propagated filter

1. __boxfilter__ 
2. __fast guided filter__
3. __fast bilateral filter__
4. __permutohedral bilateral filter__
5. __propagated filter__

__All original dependencies have been removed. Code could be run  independently.__

__BOXFILTER__ Simple Interface

    @param src         				image,single channel.

    @param dst      				output image,single channel.

    @param w           				width of image.

    @param h           				height of image.

    @param c      					channel of image, only c = 1.

    @param r      					local window radius.

    @return            				0:ok; 1:error
                       
int __BoxfilterFilter__(unsigned char *src, unsigned char *dst, int w, int h, int c, int r);



__Fast Guided Filter__ Simple Interface

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
                       
int __FastGuidedFilter__(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, int c, int r, float rp, float sr);

__FastBilateralFilter__ Simple Interface

    @param src         				input image.

    @param guidance           		guided image,single channel, only single channel is valid; invalid parameter in three channels.

    @param dst      				output image.

    @param w           				width of image.

    @param h           				height of image.

    @param c      					channel of image, only c = 1 or c = 3.

    @param sigma_s      			filter sigma in the coordinate space. A larger value of the parameter means that farther pixels will influence each other as long as their colors are close enough (see sigmaColor ). When d>0, it specifies the neighborhood size regardless of sigmaSpace. Otherwise, d is proportional to sigmaSpace.

    @param sigma_r      			filter sigma in the color space. A larger value of the parameter means that farther colors within the pixel neighborhood (see sigmaSpace) will be mixed together, resulting in larger areas of semi-equal color.

    @return            				0:ok; 1:error
                       
int __FastBilateralFilter__(unsigned char *src, unsigned char *guidance, unsigned char *dst, int w, int h, int c, float sigma_s, float sigma_r);

__PermutohedralBilateralFilter__ Simple Interface

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

__PropagatedFilter__ Simple Interface

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


Filter implementations are described in the paper:

  Adams, Andrew, Jongmin Baek, and Myers Abraham Davis. "Fast High‐Dimensional Filtering Using the Permutohedral Lattice." Computer Graphics Forum. Vol. 29. No. 2. Blackwell Publishing Ltd, 2010.

  He, Kaiming, Jian Sun, and Xiaoou Tang. "Guided image filtering." European conference on computer vision. Springer Berlin Heidelberg, 2010.

  He K, Sun J. Fast guided filter[J]. arXiv preprint arXiv:1505.00996, 2015.

  Paris, Sylvain, and Frédo Durand. "A fast approximation of the bilateral filter using a signal processing approach." European conference on computer vision. Springer Berlin Heidelberg, 2006.

  Durand, Frédo, and Julie Dorsey. "Fast bilateral filtering for the display of high-dynamic-range images." ACM transactions on graphics (TOG). Vol. 21. No. 3. ACM, 2002.

  Rick Chang, Jen-Hao, and Yu-Chiang Frank Wang. "Propagated image filtering." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2015.

Copyright (c) 2016-2017 Frotms(frotms@gmail.com)
