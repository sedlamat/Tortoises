
#ifndef _IMGPROC_HPP_
#define _IMGPROC_HPP_

/* STANDARD C++ LIBRARIES */
#include <iostream>
#include <string>

/* THIRD PARTY LIBRARIES */
#include <opencv2/opencv.hpp>


using namespace std;
//using namespace cv;

namespace my {

  void display(const cv::Mat& img, const std::string& window_name="image") {
	  
    cv::Mat dst;
    // check whether img is unsigned 8bit, values 0...255
	if ( img.type() != CV_8UC1 && img.type() != CV_8UC2 && 
		 img.type() != CV_8UC3 && img.type() != CV_8UC4 ) {  
		double min_val;
		double max_val;
		cv::minMaxLoc(img, &min_val, &max_val); 
		double alpha = 255/(max_val-min_val);
		double beta = -min_val*alpha;
		img.convertTo(dst,CV_8U,alpha,beta);
	}
	else {
		dst = img;
	}
    cv::namedWindow(window_name, cv::WINDOW_NORMAL);
    cv::imshow(window_name, dst);
    cv::waitKey(0);
    cv::destroyAllWindows();
  }  
  
  void gradient_scharr(const cv::Mat& src, cv::Mat& dst_dx, 
					   cv::Mat& dst_dy)
  {
    cv::Mat kernel = (cv::Mat_<int>(3,3) << -3.f, 0.f, 3.f, -10.f, 
											0.f, 10.f, -3.f, 0.f, 3.f);
	cv::transpose(kernel,kernel);
    std::cout << kernel << std::endl;
    cv::filter2D(src, dst_dx, CV_32F, kernel);
    cv::transpose(kernel,kernel);
    
  }

}


#endif /* _IMGPROC_HPP_ */
