
/* STANDARD C++ LIBRARIES */
#include <map>
#include <string>

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
#include <opencv2/core/core.hpp>
//~ #include <opencv2/highgui/highgui.hpp>
//~ #include <opencv2/imgproc/imgproc.hpp>
//~ #include <opencv2/ml/ml.hpp>
//include <opencv2\legacy\legacy.hpp>



 // colour that can be used for marking
 // important points in the image
 enum ColourType
 {
   red    = 0,
   green  = 1,
   blue   = 2,
   white  = 3,
   black  = 4,
   yellow = 5, //(red+green)
   purple = 6, //(red+blue)
   cyan   = 7  //(green+blue)
 };

class Tortoise {

    cv::Mat img; // image of the tortoise

    struct Junctions { // plastron junctions
	cv::Point j1Head, j2GulToHum, j3HumToPec, j4PecToAbd,
				      j5AbdToFem, j6FemToAna, j7Tail;
    };

    Junctions l_juncs, r_juncs; // left and right junctions

    std::map<std::string, cv::Vec3b> color;

    color["red"] = cv::Vec3b(0,0,255);
    color["green"] = cv::Vec3b(0,255,0);
    color["blue"] = cv::Vec3b(255,0,0);
    color["white"] = cv::Vec3b(255,255,255);
    color["black"] = cv::Vec3b(0,0,0);
    color["yellow"] = cv::Vec3b(0,255,255);
    color["purple"] = cv::Vec3b(255,0,255);
    color["cyan"] = cv::Vec3b(255,255,0);

    color["j1Head"] = color["red"];
    color["j2GulToHum"] = color["green"];
    color["j3HumToPec"] = color["blue"];
    color["j4PecToAbd"] = color["white"];
    color["j5AbdToFem"] = color["yellow"];
    color["j6FemToAna"] = color["purple"];
    color["j7Tail"] = color["cyan"];

public:
    Tortoise(cv::Mat &img): img(img) {}
    ~Tortoise(){}
    void measure(); // locate and measure features

private:

    string            m_strLoadDirectory;
    string            m_strSaveDirectory;
    string            m_strTortoiseName;
    string            m_strProcessID;
    string            m_strResultsDataFileName;

    CImg<int>	      m_imgOriginal;
    CImg<int>         m_imgRotatedOriginal;
    CImg<int>         m_imgResultsData;
    CImg<int>         m_imgCentralPathXCoordinates;

    int               m_nWidth, m_nHeight;
    int			          m_nRotationAngleToVerticalInDegrees;
    int		    	      m_imgResultsDataSize;

    plastronJunctions m_sJunctions;
    plastronJunctions m_sLeftJunctions;
    plastronJunctions m_sRightJunctions;

    point2dCoor       m_sPlastronCentreOnRotatedImg;
    point2dCoor       m_sCentralSeamHeadEnd;
    point2dCoor       m_sCentralSeamTailEnd;

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

   CImg<double> GetNormalizedJunctionDifferencesAndSeamLengths(const plastronJunctions &sLeftJunctions, const plastronJunctions &sRightJunctions);
 public:

    //Tortoise(string imgDirectory, string tortoiseName, string processID); 	//constructor
    //~Tortoise(){}                                                           // destructor
   void Recognition();
   void Classification();
 };

   void AccuracyEvaluation(string strLoadDirectory);
   void kNNclassifier();
   void Classification();
