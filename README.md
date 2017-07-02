	In this project we implement the algorithm put forward by hirakawa et. al in this paper: "Adaptive homogeneity-directed demosaicing
algorithm" published on IEEE Transactions on Image Processing 2005. For more details, please visit: http://photo-lovers.org/pdf/hirakawa05mndemosaictip.pdf.

The code is devided into five parts:

1. color_array_downsample.m simulate a bayer filter for input image. For more details about bayer filter, please visit: https://en.wikipedia.org/wiki/Color_filter_array

2. rgb2lab.m is a color space transform function, which is used in our implementation and already included by MatlabR2016 but not in the older version.

3. cfa_reconstruct.m offers some basic interpolation methods for comparison with our implementation.

4.homogeneity.m is the main function of our implementation.

5. we also provide some test code in the "test code" folder

If you have any suggestion or find some bugs to report, welcome contact lyxok1@sjtu.edu.cn 
