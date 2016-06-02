/**
    my_img_proc.hpp

    Higher level procedures and functions for image processing, using
    OpenCV library.

    @author Matej Sedlacek
    @version 0.0
*/

#ifndef _MY_IMGPROC_HPP_
#define _MY_IMGPROC_HPP_

/* STANDARD C++ LIBRARIES */
#include <iostream>
#include <string>

/* THIRD PARTY LIBRARIES */
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"


namespace my
{
    /**
	Print message to console. Short version for debugging.

	@param msg - Message to be printed.
	@return void
    */
    template <typename T> void prt(T msg)
    {
	std::cout << msg << std::endl;
    }

    /**
	Gets minimum value of cv::Mat. Short version for debugging.

	@param src - Input image.
	@return Min value of the input image.
    */
    double minMat(const cv::Mat &src)
    {
	double min, max;
	cv::minMaxLoc(src, &min, &max);
	return min;
    }

    /**
	Gets maximum value of cv::Mat. Short version for debugging.

	@param src - Input image.
	@return Max value of the input image.
    */
    double maxMat(const cv::Mat &src)
    {
	double min, max;
	cv::minMaxLoc(src, &min, &max);
	return max;
    }

    /**
	Gets a scaled image to the range [0...255] and casted into
	CV_8U. For details on the transformation see the documentation
	of Mat::convertTo(....).

	@param src - Input image.
	@return Scaled and casted image.
    */
    cv::Mat get_scaled_CV_U8(const cv::Mat &src)
    {
	cv::Mat dst;
	double min_val, max_val;
	cv::minMaxLoc(src, &min_val, &max_val);
	double alpha = 255 / (max_val - min_val);
	double beta = - min_val * alpha;
	src.convertTo(dst, CV_8U, alpha, beta);
	return dst;
    }

    /**
	Displays an unsigned 8bit (CV_8U) image, values [0...255], the
	number of channels allowed is 1 and 3. The image is scaled and
	casted to CV_8U if it is of other type. Short version for
	debugging.

	@param src - Image to be displayed.
	@param window_name - Displayed window name (optional).
	@return void.
    */
    void display(const cv::Mat &src,
		 const std::string &window_name="image")
    {
	CV_Assert(src.channels() == 1 || src.channels() == 3);
	cv::Mat dst;
	// check if src is CV_8U, if not make it
	if (src.type() != CV_8UC1 && src.type() != CV_8UC3) {
	    dst = my::get_scaled_CV_U8(src);
	} else {
	    dst = src;
	}
	cv::namedWindow(window_name, cv::WINDOW_NORMAL);
	cv::imshow(window_name, dst);
	cv::waitKey(0);
    }

    /**
	Gets gradient of an image using Scharr masks.

	@param src - Input image.
	@return A vector of [dx,dy] derivatives.
    */
    std::vector<cv::Mat> get_gradient_scharr(const cv::Mat &src)
    {
	std::vector<cv::Mat> dxdy(2);
	cv::Mat kernel = (cv::Mat_<int>(3,3) << 3, 0, -3,
					       10, 0,-10,
					        3, 0, -3);
	cv::filter2D(src, dxdy[0], CV_32F, kernel);
	cv::transpose(kernel, kernel);
	//cv::flip(kernel,kernel,0);
	/* Flip it if the kernel mask is to be
	applied on the image as it is displayed (in the same way for y as
	for x), otherwise as the y coordinate of cv::Mat image goes from
	top to bottom in an image, it would give results upside down (for
	dy). Counter-clockwise.*/
	cv::filter2D(src, dxdy[1], CV_32F, kernel);
	return dxdy;
    }

  /**
    Gets values that in the absolute value are the maximum of
    the BGR values.

    I(x,y) = B/G/R(x,y)
    if ( abs(B/G/R(x,y)) >= abs(G/R/B(x,y)) &&
         abs(B/G/R(x,y)) >= abs(R/B/G(x,y)) )

    @param src - An input BGR color image.
    @return An image with B/G/R values of max abs B/G/R.
  */
  cv::Mat get_BGR_max_of_abs(const cv::Mat& src)
  {
    CV_Assert(src.channels() == 3);

    cv::Mat src_abs, BgeG, BgeR, GgeR, dst, Bg, Gg, Rg;
    std::vector<cv::Mat> BGR_abs(3), BGR(3);

    src_abs = cv::abs(src);
    cv::split(src_abs, BGR_abs);
    cv::split(src, BGR);

    BgeG = BGR_abs[0] >= BGR_abs[1];
    BgeR = BGR_abs[0] >= BGR_abs[2];
    GgeR = BGR_abs[1] >= BGR_abs[2];
    cv::Mat(((BgeG & BgeR)/255)).convertTo(Bg, CV_32F);
    cv::Mat(((GgeR & ~BgeG)/255)).convertTo(Gg, CV_32F);
    cv::Mat(((~GgeR & ~BgeR)/255)).convertTo(Rg, CV_32F);
    dst = BGR[0].mul(Bg) + BGR[1].mul(Gg) + BGR[2].mul(Rg);
    return  dst;
  }

  /**
    Gets maxima values of BGR values of an image.

    @param src - An input BGR color image.
    @return An image with maxima from input image BGR values.
  */
  cv::Mat get_BGR_maxima(const cv::Mat& src)
  {
    CV_Assert(src.channels() == 3);
    cv::Mat dst;
    std::vector<cv::Mat> BGR(src.channels());
    cv::split(src, BGR);
    return cv::Mat(cv::max(BGR[0],cv::max(BGR[1],BGR[2])));
  }


  /**
    Getting minima values of BGR values of an image.

    @param src - The input BGR color image.
    @return An image with minima from input image BGR values.
  */
  cv::Mat get_BGR_minima(const cv::Mat& src)
  {
    CV_Assert(src.channels() == 3);
    cv::Mat dst;
    std::vector<cv::Mat> BGR(src.channels());
    cv::split(src, BGR);
    return cv::Mat(cv::min(BGR[0],cv::min(BGR[1],BGR[2])));
  }

  /**
    Gets magnitude of a gradient of an image. For a
    multi-channel image, it first computes gradient for all
    channels and then takes the maximum over the channels for
    the computation of gradient.

    @param src - An input image.
    @return A gradient magnitude of the input image.
  */
  cv::Mat get_gradient_magnitude(const cv::Mat& src)
  {
    cv::Mat grad_mag;
    std::vector<cv::Mat> dxdy = my::get_gradient_scharr(src);
    if (src.channels() > 1) {
      dxdy[0] = my::get_BGR_maxima(cv::abs(dxdy[0]));
      dxdy[1] = my::get_BGR_maxima(cv::abs(dxdy[1]));
    }
    cv::magnitude(dxdy[0], dxdy[1], grad_mag);
    return grad_mag;
  }

  /**
    Gets orientation of a gradient of an image. For a
    multi-channel image, it first computes gradient for all
    channels and then takes the maximum over the channels for
    the computation of the orientation.

    @param src - An input image.
    @return A gradient orientation of the input image.
  */
  cv::Mat get_gradient_orientation(const cv::Mat& src)
  {

    cv::Mat grad_orient;
    std::vector<cv::Mat> dxdy = my::get_gradient_scharr(src);
    if (src.channels() > 1) {
      dxdy[0] = my::get_BGR_max_of_abs(dxdy[0]);
      dxdy[1] = my::get_BGR_max_of_abs(dxdy[1]);
    }
    cv::phase(dxdy[0], dxdy[1], grad_orient, true);
    //my::display(grad_orient);
    return grad_orient;
  }

  /**
    Gets a 3-channel image by stacking a single channel image
    along the z-axis (channels).

    @param src - An input single-channel image.
    @return A stacked 3-channel image.
  */
  cv::Mat get_3channel(const cv::Mat& src)
  {
    cv::Mat dst;
    std::vector<cv::Mat> channels;
    channels.push_back(src);
    channels.push_back(src);
    channels.push_back(src);
    cv::merge(channels, dst);
    return dst;
  }


  /**
    Gets image edges based on the gradient magnitude in both
    the intensity and the color of the image.

    @param src - An input image.
    @return An edge image.
  */
  cv::Mat get_edges_color_based(const cv::Mat& src)
  {
    cv::Mat dst, BGR_minima, color, color_grad_mag;
    cv::Mat BGR_minima_grad_mag, mag, mag_edges, edges;
    dst = src;
    //cv::GaussianBlur(src, dst, cv::Size(0,0), 1.0);
    BGR_minima = my::get_BGR_minima(dst);
    BGR_minima = my::get_3channel(BGR_minima);
    color = dst - BGR_minima;

    my::display(BGR_minima);
    cv::Mat minEdges;
    cv::Canny(BGR_minima, minEdges, 100, 150);
    my::display(minEdges);
    //~ my::display(my::get_gradient_scharr(color)[0]);
    //~ my::display(my::get_gradient_scharr(color)[1]);
    //~ my::display(my::get_gradient_scharr(color)[1] + my::get_gradient_scharr(color)[0]);
    //~ my::display(my::get_gradient_scharr(dst)[1] + my::get_gradient_scharr(dst)[0]);
    color_grad_mag = my::get_gradient_magnitude(color);
    //my::display(color_grad_mag);
    BGR_minima_grad_mag = my::get_gradient_magnitude(BGR_minima);
    //my::display(BGR_minima_grad_mag);
    // std::cout << color_grad_mag.type() << " " << BGR_minima_grad_mag.type() << std::endl;
    //mag = BGR_minima_grad_mag;
    mag =  color_grad_mag;//.mul(BGR_minima_grad_mag);
    mag = my::get_scaled_CV_U8(mag);
    //my::display(mag);
    cv::Scalar_<double> mag_mean = cv::mean(mag, mag>0);
    //std::cout << mag_mean << std::endl;
    //cv::Canny(mag, mag_edges, 50, 150);
    //cv::adaptiveThreshold(mag, mag_edges, 255, cv::ADAPTIVE_THRESH_MEAN_C,
	//		cv::THRESH_BINARY, 101, 0);
    cv::threshold(mag, mag_edges, mag_mean[0], 255, cv::THRESH_BINARY);

    mag_edges = my::get_scaled_CV_U8(mag_edges);
    //cv::GaussianBlur(BGR_minima, BGR_minima, cv::Size_<int>(0,0), 1);
    cv::Canny(BGR_minima, edges, 20, 50);

    cv::bitwise_and(edges, mag_edges, dst);
    //my::display(dst);
    dst = edges;
    my::display(edges);
    //LOCAL THRESHOLDING !!!!
    //~ exit(2);
    return dst;
  }

  // COMMENT TODO
  void visualize_points(std::vector<cv::Point_<int> > points,
			cv::Size_<int> size, cv::Point_<int> shift)
  {
    cv::Mat dst(size, CV_8UC1, cv::Scalar_<uchar>(0));
    cv::Rect_<int> area(cv::Point(0,0),size);
    for(auto const& pt : points) {
	if ( (pt+shift).inside(area))
	    dst.at<uchar>(pt+shift) = 255;
    }
    my::display(dst);
  }


//~ #include <X11/Xlib.h>
//~ #include <X11/Xutil.h>
//~ #include <cstdint>
//~ #include <cstring>
//~ #include <vector>
//~ #include <iostream>
//~ #include <opencv2/opencv.hpp>
//~ #include <ctime>
//~
//~ void ImageFromDisplay(std::vector<uint8_t>& Pixels, int& Width, int& Height, int& BitsPerPixel)
//~ {
    //~ Display* display = XOpenDisplay(nullptr);
    //~ Window root = DefaultRootWindow(display);
//~
    //~ XWindowAttributes attributes = {0};
    //~ XGetWindowAttributes(display, root, &attributes);
//~
    //~ Width = attributes.width;
    //~ Height = attributes.height;
//~
    //~ XImage* img = XGetImage(display, root, 0, 0 , Width, Height, AllPlanes, ZPixmap);
    //~ BitsPerPixel = img->bits_per_pixel;
    //~ Pixels.resize(Width * Height * 4);
//~
    //~ memcpy(&Pixels[0], img->data, Pixels.size());
//~
    //~ XFree(img);
    //~ XCloseDisplay(display);
//~ }
//~
//~ int main()
//~ {
        //~ float t0 = clock()*1.0/CLOCKS_PER_SEC;
  //~ for( int ii = 0; ii < 100; ii += 1)
    //~ {
      //~ int Width = 0;
      //~ int Height = 0;
      //~ int Bpp = 0;
      //~ std::vector<std::uint8_t> Pixels;
      //~
      //~ ImageFromDisplay(Pixels, Width, Height, Bpp);
      //~
      //~ if (Width && Height)
	//~ {
	  //~ cv::Mat img = cv::Mat(Height, Width, Bpp > 24 ? CV_8UC4 : CV_8UC3, &Pixels[0]);
	  //~ //Mat(Size(Height, Width), Bpp > 24 ? CV_8UC4 : CV_8UC3, &Pixels[0]);
	  //~ // std::cout << img;
	  //~ cv::namedWindow("WindowTitle", cv::WINDOW_AUTOSIZE);
	  //~ cv::imshow("Display window", img);
//~
	  //~ cv::waitKey(0);
	//~ }
    //~ }
  //~ std::cout << (clock()*1.0/CLOCKS_PER_SEC - t0)/100.0 << std::endl;
   //~ return 0;
//~ }



} /* namespace my */


#endif /* _MY_IMGPROC_HPP_ */
