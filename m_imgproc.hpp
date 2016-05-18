
#ifndef _IMGPROC_HPP_
#define _IMGPROC_HPP_

/* STANDARD C++ LIBRARIES */
#include <iostream>
#include <string>

/* THIRD PARTY LIBRARIES */
#include <opencv2/opencv.hpp>


using namespace std;
using namespace cv;

void main() {
	
	Mat img = imread()
}

namespace m_
{

  void display(Mat& img, const string& window_name="image")
  {
    Mat dst;

    
    // check whether img is unsigned 8bit, values 0...255
	if ( img.type() != CV_8UC1 && img.type() != CV_8UC2 && 
		 img.type() != CV_8UC3 && img.type() != CV_8UC4 ) {
			  
		double minVal;
		double maxVal;
		minMaxLoc(img, &minVal, &maxVal); 
		img.convertTo(dst,CV_8U,) 
	}
		
	dst = img;

    cv::namedWindow(window_name, cv::WINDOW_NORMAL);
    cv::imshow(window_name, dst);
    cv::waitKey(0);
    cv::destroyAllWindows();
  }  
  
  void gradient_scharr(const cv::Mat& src, cv::Mat& dst_dx, cv::Mat& dst_dy)
  {
    cv::Mat kernel = (cv::Mat_<int>(3,3) << -3.f, 0.f, 3.f, -10.f, 0.f, 10.f, -3.f, 0.f, 3.f);
    std::cout << kernel << std::endl;
    cv::filter2D(src, dst_dx, CV_32F, kernel);
  }
  
  void scale_values(cv::Mat& dst, const int& norm_min=0, const int& norm_max=255)
  {
    
  }

}

#endif /* _IMGPROC_HPP_ */
