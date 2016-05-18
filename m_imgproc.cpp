/* STANDARD C++ LIBRARIES */
#include <iostream>
#include <string>

/* THIRD PARTY LIBRARIES */
#include <opencv2/opencv.hpp>

/* FIRST PARTY LIBRARIES */
#include "m_imgproc.hpp"

using namespace std;
using namespace cv;

int main()
{
  string a = "HELLo";
  cout << a << endl;
  Mat img = imread("baboon.jpg",1);
  Mat dx, dy, img_GRAY;
  my::gradient_scharr(img, dx, dy);
  //cvtColor(img,img_GRAY, cv::COLOR_BGR2GRAY);
  my::display(dx);
  return 1;
}
