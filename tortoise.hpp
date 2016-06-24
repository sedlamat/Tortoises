
/* STANDARD C++ LIBRARIES */
#include <map>
#include <string>

/* THIRD PARTY LIBRARIES */
#include <opencv2/core/core.hpp>

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

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"
#include "general_hough.hpp"



class Tortoise {

    const cv::Mat _plastron_img; // image of the tortoise plastron

    struct Junctions { // plastron junctions
	cv::Point j1Head, j2GulToHum, j3HumToPec, j4PecToAbd,
		  j5AbdToFem, j6FemToAna, j7Tail;
    };

    Junctions l_juncs, r_juncs; // left and right junctions
//~
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

public:
    Tortoise(const cv::Mat &plastron_image,
	     const cv::Mat &template_image,
	     const cv::Point &template_reference_point);
    virtual ~Tortoise(){}
    void measure(); // locate and measure features

private:
    void locate_plastron(const cv::Mat template_image,
			 const cv::Point reference_point);
};

Tortoise::Tortoise(const cv::Mat &plastron_image,
		   const cv::Mat &template_image,
		   const cv::Point &template_reference_point)
		   : _plastron_img(plastron_image)
{
    this->locate_plastron(template_image, template_reference_point);
}

void Tortoise::locate_plastron(const cv::Mat plastron_template,
			       const cv::Point reference_point)
{
    std::vector<int> angles = {-5,0,5,-85,-90,-95,180,85,90,95};

    std::map<std::string, cv::Point> template_junctions_map;
    template_junctions_map["j1Head"] = cv::Point(64,58);
    template_junctions_map["j2GulToHum"] = cv::Point(64,74);
    template_junctions_map["j3HumToPec"] = cv::Point(64,97);
    template_junctions_map["j4PecToAbd"] = cv::Point(64,105);
    template_junctions_map["j5AbdToFem"] = cv::Point(64,155);
    template_junctions_map["j6FemToAna"] = cv::Point(64,168);
    template_junctions_map["j7Tail"] = cv::Point(64,185);

    GeneralHough general_hough(_plastron_img, plastron_template,
			       reference_point,	angles, 20, 200, 1.0,
			       0.3, 50, 100, 0.5, 0,
			       &template_junctions_map);

    for (auto junc : template_junctions_map) {
	std::cout << junc.second << std::endl;
    }

    cv::Mat result_img = general_hough.get_result_img();
    sedlamat::display(result_img);
    double best_accum_value = general_hough.get_best_accum_val();
    double best_angle = general_hough.get_best_angle();
    double best_scale = general_hough.get_best_scale();
    cv::Point best_reference_point = general_hough.get_best_ref_pt();
}

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
*/
