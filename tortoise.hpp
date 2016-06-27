/**
    tortoise.hpp

    Implementation of the Tortoise class for tortoise plastron
    locatization and feature extraction.

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

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"
#include "general_hough.hpp"

/**************** class definition *********************************/

class Tortoise {
    // image of the tortoise plastron
    const cv::Mat _plastron_img;

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
};

/************* member functions definitions ********************/

Tortoise::Tortoise(const cv::Mat &plastron_image,
		   const std::string &tortoise_name)
	    :
	    _plastron_img(plastron_image),
	    _tortoise_name(tortoise_name),
	    _info_file_name(tortoise_name.substr(0,7) + "_info.txt"),
	    _angle(0),
	    _l_juncs(std::vector<cv::Point>(7)),
	    _r_juncs(std::vector<cv::Point>(7)),
	    _center(),
	    _left_side(),
	    _right_side(),
	    _plastron_found(0),
	    _junctions_found(0)
{
    this->read_info();

    if (!_plastron_found) {
	this->locate_plastron();
	_plastron_found = 1;
	this->write_info();
    }

    if (_plastron_found && !_junctions_found) {
	this->locate_junctions();
	//_junctions_found = 1; enable only if localization is complete
	this->write_info();
    }
}

void Tortoise::read_info()
{
    std::fstream info(_info_file_name, std::fstream::in |
				       std::fstream::out);
    if (info.good()) {
	info >> _plastron_found;
	info >> _junctions_found;
	info >> _angle;
	std::cout << _plastron_found << '\n';
	std::cout << _junctions_found << '\n';
	std::cout << _angle << '\n';
    } else {
	info.close();
	std::cout << "Error when opening/reading Tg#####_info.txt";
	std::cout << " file."<< std::endl;
	std::exit(1);
    }
}

void Tortoise::write_info() const
{
    std::ofstream info(_info_file_name, std::fstream::out);
    if (info.good()) {
	info << _plastron_found << '\n';
	info << _junctions_found << '\n';
	info << _angle << '\n';
	std::cout << _plastron_found << '\n';
	std::cout << _junctions_found << '\n';
	std::cout << _angle << '\n';

	info.close();
    } else {
	info.close();
	std::cout << "Error when opening/reading Tg#####_info.txt";
	std::cout << " file."<< std::endl;
	std::exit(1);
    }
}


void Tortoise::print_juncs() const
{
    std::cout << "Left junctions:" << std::endl;
    for (const auto &junc : _l_juncs) {
	std::cout << junc << std::endl;
    }

    std::cout << "Right junctions:" << std::endl;
    for (const auto &junc : _r_juncs) {
	std::cout << junc << std::endl;
    }
}

void Tortoise::locate_plastron()
{
    const cv::Mat plastron_template
			    = cv::imread("plastron_template.bmp", 0);

    if (plastron_template.empty()) {
	std::cout << "Error: plastron_template.bmp image is not" \
	" inside the directory of the executable file."
	<< std::endl;
	exit(1);
    }

    const std::vector<int> angles = {-5,0,5,-85,-90,-95,180,85,90,95};

    const cv::Point reference_point(65,125);

    std::vector<cv::Point> interest_pts;
    interest_pts.push_back(cv::Point(64,58)); 	// junction 0 - Head
    interest_pts.push_back(cv::Point(64,74));
    interest_pts.push_back(cv::Point(64,97));
    interest_pts.push_back(cv::Point(64,105));
    interest_pts.push_back(cv::Point(64,155));
    interest_pts.push_back(cv::Point(64,168));
    interest_pts.push_back(cv::Point(64,185)); 	// junction 6 - Tail
    interest_pts.push_back(cv::Point(7,131));  	// left side border
    interest_pts.push_back(cv::Point(123,131)); // right side border
    interest_pts.push_back(cv::Point(65,125));  // plastron center

    GeneralHough general_hough(_plastron_img, plastron_template,
			       reference_point,	angles, 20, 200, 1.0,
			       0.3, 50, 100, 0.5, 0,
			       &interest_pts);
    // assigns interest points

    for (int ii = 0; ii < 7; ++ii) {
	_l_juncs[ii] = _r_juncs[ii] = interest_pts[ii];
    }
    _left_side = interest_pts[7];
    _right_side = interest_pts[8];
    _center = interest_pts[9];

    _angle = general_hough.get_best_angle();

    cv::Mat result_img = general_hough.get_result_img();
    sedlamat::display(result_img);
}

void Tortoise::locate_junctions()
{

}

#endif /* _TORTOISE_HPP_ */

/*
    string            m_strLoadDirectory;
    string            m_strSaveDirectory;
    string            m_strTortoiseName;
    string            m_strProcessID;
    string            m_strResultsDataFileName;
//~ //~
    CImg<int>	      m_imgOriginal;
    CImg<int>         m_imgRotatedOriginal;
    CImg<int>         m_imgResultsData;
    CImg<int>         m_imgCentralPathXCoordinates;
//~ //~
    int               m_nWidth, m_nHeight;
    int			          m_nRotationAngleToVerticalInDegrees;
    int		    	      m_imgResultsDataSize;
//~ //~
    plastronJunctions m_sJunctions;
    plastronJunctions m_sLeftJunctions;
    plastronJunctions m_sRightJunctions;
//~ //~
    point2dCoor       m_sPlastronCentreOnRotatedImg;
    point2dCoor       m_sCentralSeamHeadEnd;
    point2dCoor       m_sCentralSeamTailEnd;
//~ //~
   // Methods for tortoise recognition
   void CentralSeamLocalization();
   CImg<int> GetUniformlyResizedImg(CImg<int> img, int resizedMaxSizeInPix);
   inline float GetImgResolutionCoeff(int width, int height);
   CImg<unsigned char> GetEdgeImg(const CImg<int> &img, int lowerThreshold, int upperThreshold, double sigma);
   CImg<unsigned char> GetDirectionOfEdges(const CImg<int> &edgeImg);
   CImg<unsigned char> EliminateOnePixelEdges(CImg<unsigned char> &imgToDoItOn);
   CImg<unsigned char> EliminateTwoPixelsEdges(CImg<unsigned char> &imgToDoItOn);
   void TortoiseContourLocalizationWithGHT();
   void GetAbdJunctionsPreciseLocations(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg);
   void MarkAJunction(CImg<int> &img, const point2dCoor &sJuntionToMark, const bool bSquareFilled, ColourType colour);
   void GetFemToAnaJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions);
   void GetHumToPecJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions);
   void GetGulToHumJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions);
   void GetHeadJunctionLocalizationWithGTH(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg);
   void GetTailJunctionLocalizationWithGTH(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg);
   CImg<unsigned char> SkeletonizationMinDist(CImg<unsigned char> imgEdges);
   CImg<unsigned char> GetCentralSeam(const CImg<unsigned char> &imgEdges);
   void GetCentralPathXCoordinatesForInnerJunctions(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg);
   void GetGulToHumLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
   void GetHumToPecLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
   void GetFemToAnaLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
   void GetAbdLefAndRightJunctions(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions,const  point2dCoor &sPlastronCentreOnRotatedImg);
   void GetCentralPathXCoordinatesForHeadAndTailArea();
   void GetTailJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
   void GetHeadJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
//~ //~
   CImg<double> GetNormalizedJunctionDifferencesAndSeamLengths(const plastronJunctions &sLeftJunctions, const plastronJunctions &sRightJunctions);
 public:
//~ //~
    //Tortoise(string imgDirectory, string tortoiseName, string processID); 	//constructor
    //~Tortoise(){}                                                           // destructor
   void Recognition();
   void Classification();
   * */
 //};
/*
   void AccuracyEvaluation(string strLoadDirectory);
   void kNNclassifier();
   void Classification();
   *
   *
    //~ std::map<std::string, cv::Vec3b> colors;
    //~ std::cout << colors.size() << std::endl;
//~
    //~ colors.insert(std::pair<std::string,
    //~ cv::Vec3b>("a",cv::Vec3b(0,0,255)) );
    //~ colors["red"] = cv::Vec3b(0,0,255);
    //~ color["green"] = cv::Vec3b(0,255,0);
    //~ color["blue"] = cv::Vec3b(255,0,0);
    //~ color["white"] = cv::Vec3b(255,255,255);
    //~ color["black"] = cv::Vec3b(0,0,0);
    //~ color["yellow"] = cv::Vec3b(0,255,255);
    //~ color["purple"] = cv::Vec3b(255,0,255);
    //~ color["cyan"] = cv::Vec3b(255,255,0);

    //~ color["j1Head"] = color["red"];
    //~ color["j2GulToHum"] = color["green"];
    //~ color["j3HumToPec"] = color["blue"];
    //~ color["j4PecToAbd"] = color["white"];
    //~ color["j5AbdToFem"] = color["yellow"];
    //~ color["j6FemToAna"] = color["purple"];
    //~ color["j7Tail"] = color["cyan"];
*/


//~ // using several standard libraries
//~ #include <iostream>
//~ #include <time.h>
//~ #include <cmath>
//~ #include <iomanip>
//~ #include <stdio.h>
//~ #include <fstream>
//~ #include <sstream>
//~ #include <string>
//~ //#include <io.h>
//~ //#include <direct.h>
//~ #include <float.h>
//~ //#include <cassert>
//~
//~ // using libraries for image recognition (CImg and OpenCV)
//~ //#include "CImg.h"
//~ #include <opencv2/opencv.hpp>

//~ #include <opencv2/highgui/highgui.hpp>
//~ #include <opencv2/imgproc/imgproc.hpp>
//~ #include <opencv2/ml/ml.hpp>
//include <opencv2\legacy\legacy.hpp>
