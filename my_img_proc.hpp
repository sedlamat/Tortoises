/**
    my_imgproc.hpp
    
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
#include <exception>

/* THIRD PARTY LIBRARIES */
#include <opencv2/opencv.hpp>

namespace my 
{
	
	/**
		Scales the values of an input image into the range [0...255] 
		and casts it into CV_8U.

		@param src - The input image.
		@return Scaled and casted image.
	*/
	cv::Mat get_scaled_CV_U8(const cv::Mat& src)
	{
		cv::Mat dst;
		double min_val;
		double max_val;
		cv::minMaxLoc(src, &min_val, &max_val); 
		double alpha = 255/(max_val-min_val);
		double beta = -min_val*alpha;
		src.convertTo(dst,CV_8U,alpha,beta);
		return dst;
	}  	
	
	/**
		Displays an unsigned 8bit (CV_8U) image. The image is scaled
		and casted to CV_8U if it is of other type. 

		@param src - The image to be displayed.
		@param window_name - Displayed window name (optional).
		@return void.
	*/
	void display(const cv::Mat& src, 
	    		 const std::string& window_name="image")
	{
		cv::Mat dst;
		// check if src is CV_8U, if not make it
		if ( src.type() != CV_8UC1 && src.type() != CV_8UC2 && 
		     src.type() != CV_8UC3 && src.type() != CV_8UC4 ) {  
			dst = my::get_scaled_CV_U8(src);
   	  	} else {
			dst = src;
		}
		cv::namedWindow(window_name, cv::WINDOW_NORMAL);
		cv::imshow(window_name, dst);
		cv::waitKey(0);
	}  
 
	/**
		Getting gradient of an image using Scharr masks.

		@param src - The input image.
		@return A vector of [dx,dy] derivatives.
	*/
	std::vector<cv::Mat> get_gradient_scharr(const cv::Mat& src) 
	{	
		std::vector<cv::Mat> dxdy(2);  
		cv::Mat kernel = (cv::Mat_<int>(3,3) << -3, 0, 3, -10, 0, 
												10, -3, 0, 3);
		cv::filter2D(src, dxdy[0], CV_32F, kernel);
		cv::transpose(kernel, kernel);
		cv::filter2D(src, dxdy[1], CV_32F, kernel);
		return dxdy;
	}
 
  	/**
		Getting maxima values of BGR values of an image.

		@param src - The input BGR color image.
		@return An image with maxima from input image BGR values.
	*/ 
	cv::Mat get_BGR_maxima(const cv::Mat& src)
	{
		cv::Mat dst;
		if (src.channels() == 3) {
			std::vector<cv::Mat> BGR(src.channels());
			cv::split(src, BGR);
			return cv::Mat(cv::max(BGR[0],cv::max(BGR[1],BGR[2])));
		} else {
			throw "In get_BGR_maxima( ) invalid num of channels.";
			return src;
		}
	}
	
 	/**
		Getting minima values of BGR values of an image.

		@param src - The input BGR color image.
		@return An image with minima from input image BGR values.
	*/ 
	cv::Mat get_BGR_minima(const cv::Mat& src)
	{
		cv::Mat dst;
		if (src.channels() == 3) {
			std::vector<cv::Mat> BGR(src.channels());
			cv::split(src, BGR);
			return cv::Mat(cv::min(BGR[0],cv::min(BGR[1],BGR[2])));
		} else {
			throw "In get_BGR_minima( ) invalid num of channels.";
			return src;
		}
	}
	
 	/**
		Getting magnitude of a gradient of an image. For a 
		multi-channel image, it first computes gradient for all 
		channels and then takes the maximum over the channels for 
		magnitude computation.

		@param src - The input image.
		@return A gradient magnitude of the input image.
	*/
	cv::Mat get_gradient_magnitude(const cv::Mat& src) 
	{
		std::vector<cv::Mat> dxdy = my::get_gradient_scharr(src);
		my::display(cv::abs(dxdy[0]));
		my::display(dxdy[1]);
		if (src.channels() > 1) {
			dxdy[0] = my::get_BGR_maxima(cv::Mat(cv::abs(dxdy[0])));
			dxdy[1] = my::get_BGR_maxima(cv::abs(dxdy[1]));
		}
		my::display(dxdy[0]);
		my::display(dxdy[1]);
		cv::Mat mag;
		cv::magnitude(dxdy[0],dxdy[1],mag);
		return mag;
	}
	
 	/**
		Getting image edges.

		@param src - The input image.
		@return An edge image.
	*/
	cv::Mat get_edges_color_based(const cv::Mat& src) 
	{
		cv::Mat dst;
		cv::GaussianBlur(src, dst, cv::Size(0,0), 1);
		cv::Mat BGR_minima = my::get_BGR_minima(dst);
		std::vector<cv::Mat> BGR_minima_C3;
		BGR_minima_C3.push_back(BGR_minima);
		BGR_minima_C3.push_back(BGR_minima);
		BGR_minima_C3.push_back(BGR_minima);
		cv::merge(BGR_minima_C3, BGR_minima);
		cv::Mat color = dst - BGR_minima;
		cv::Mat color_grad_mag = my::get_gradient_magnitude(color);
		my::display(color_grad_mag);
		cv::Mat BGR_minima_grad_mag = 
							my::get_gradient_magnitude(BGR_minima);
		cv::Mat mag(color_grad_mag.mul(BGR_minima_grad_mag));
		cv::Mat mag_edges;
		cv::threshold(mag, mag_edges, 25000, 255, cv::THRESH_BINARY);
		mag_edges = my::get_scaled_CV_U8(mag_edges);
		cv::Mat edges;
		cv::Canny(BGR_minima, edges, 0, 50);
		my::display(mag_edges);
		my::display(edges);
		std::cout << edges.type() << std::endl;
		std::cout << mag_edges.type() << std::endl;
		std::cout << edges.type() << std::endl;
		std::cout << edges.type() << std::endl;
		cv::bitwise_and(edges, mag_edges, dst);
		my::display(dst);
		std::cout << edges.type() << std::endl;
		return dst;
	}	
	
}


#endif /* _MY_IMGPROC_HPP_ */
