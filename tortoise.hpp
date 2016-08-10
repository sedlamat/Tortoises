/**
    tortoise.hpp

    Tortoise class for tortoise plastron locatization and feature
    extraction.

    @author Matej Sedlacek
    @version 0.0
*/

#ifndef _TORTOISE_HPP_
#define _TORTOISE_HPP_

/* STANDARD C++ LIBRARIES */
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

/* THIRD PARTY LIBRARIES */
#include <opencv2/core/core.hpp>
//#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"
#include "general_hough.hpp"


/**************** class definition *********************************/

class Tortoise {
    // image of the tortoise plastron
    cv::Mat _plastron_img;

    // file name into which info data is saved
    const std::string _tortoise_name, _info_file_name;

    // rotation angle of the plastron
    int _angle;

    // left and right junctions, 0 - 6, from head to tail,
    // Head, GulToHum, HumToPec, PecToAbd, AbdToFem, FemToAna, Tail
    std::vector<cv::Point> _l_juncs, _r_juncs;

    // plastron points
    cv::Point _center, _left_side, _right_side;

    bool _plastron_found, _junctions_found;

    struct BasicColor {
	const cv::Vec3b red = cv::Vec3b(0,0,255);
	const cv::Vec3b green = cv::Vec3b(0,255,0);
	const cv::Vec3b blue = cv::Vec3b(255,0,0);
	const cv::Vec3b white = cv::Vec3b(255,255,255);
	const cv::Vec3b black = cv::Vec3b(0,0,0);
	const cv::Vec3b yellow = cv::Vec3b(0,255,255);
	const cv::Vec3b purple = cv::Vec3b(255,0,255);
	const cv::Vec3b cyan = cv::Vec3b(255,255,0);
    } const _color;

    struct JunctionIndices {
	const int j0Head = 0;
	const int j1GulHum = 1;
	const int j2HumPec = 2;
	const int j3PecAbd = 3;
	const int j4AbdFem = 4;
	const int j5FemAna = 5;
	const int j6Tail = 6;
    } _idx;

public:
    Tortoise(const cv::Mat &plastron_image,
	     const std::string &tortoise_name);
    ~Tortoise(){}

private:
    void locate_plastron();
    void locate_junctions();
    void read_info();
    void write_info() const;
    void print_juncs() const;
    void put_point_on_img(cv::Mat &img,
			  const cv::Rect &img_rect,
			  const cv::Point &pt,
			  const cv::Vec3b &color,
			  const cv::Point &shift);
    void display_points_on_plastron();
    void rotate_img_and_pts(const int angle);
    void rotate_pt(cv::Point &pt, const cv::Point &center,
		   const int angle);
    void locate_central_seam();
    std::vector<cv::Point> get_best_path(const cv::Mat &edges);
    void locate_seams();

};


/*************** Other useful functions *****************************/

void displayImageWithPreprocessing(const cv::Mat &image,
				   const std::string &windowName)
{
    CV_Assert(image.channels() == 1 || image.channels() == 3);

    cv::Mat scaledImage = getScaledImage(image);

    displayImageNoPreprocessing(scaledImage, windowName);
}

// the image is scaled into [0,...,255] and converted into CV_U8
// to control how the image is displayed, all values are in [0,255]
cv::Mat getConvertedScaledImage(const cv::Mat &image)
{
    // finds the image minimum and maximum value
    cv::Mat scaledImage;
    double minValue, maxValue;
    cv::minMaxLoc(image, &minValue, &maxValue);

    // calculates scaling coefficients
    const double alphaCoeff = 255 / (maxValue - minValue);
    const double betaCoeff = -minValue * alpha;

    cv::Mat scaledImage;
    image.convertTo(scaledImage, CV_8U, alphaCoef, betaCoef);

    return scaledImage;
}

void displayImageNoPreprocessing(const cv::Mat &image,
				 const std::string &windowName)
{
    CV_Assert(image.channels() == 1 || image.channels() == 3);

    // creates a resizable(= flag cv::WINDOW_NORMAL) window for the image
    cv::namedWindow(windowName, cv::WINDOW_NORMAL);

    // shows the image until ESC key is pressed

    const int ESC_KEY = 27;
    int keyboardResponse = 0;

    while (keyboardResponse != ESC_KEY) {
	cv::imshow(windowName, image);
	keyboardResponse = cv::waitKey(0);
    }
}


std::istream &operator>>(std::istream &stream, cv::Point &pt);

#endif /* _TORTOISE_HPP_ */
