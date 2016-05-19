// using several standard libraries
#include <iostream>
#include <time.h>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
//#include <io.h>
//#include <direct.h>
#include <float.h>
//#include <cassert>

// using libraries for image recognition (CImg and OpenCV)
//#include "CImg.h"
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/ml/ml.hpp>
//include <opencv2\legacy\legacy.hpp>

//#include "../precomp.hpp"
// using namespaces for CImg and standard libraries
//using namespace cimg_library;
using namespace std;

// #ifdef PI
// #undef PI
// #endif

// // constant PI in double and float version
// const double PI = 3.141592653589793238462643383;
// const float fPI = 3.141592653589793238462643383f;


// // structure to hold any point on the 2d image
// struct point2dCoor
// {
//   int xCoor;
//   int yCoor;
// };

// // structure to hold the seam junctions
// // on the central seam of the tortoise plastron
// struct plastronJunctions
// {
//   point2dCoor s1Head;
//   point2dCoor s2GulToHum;
//   point2dCoor s3HumToPec;
//   point2dCoor s4PecToAbd;
//   point2dCoor s5AbdToFem;
//   point2dCoor s6FemToAna;
//   point2dCoor s7Tail;
// };

// // colour that can be used for marking 
// // important points in the image
// enum ColourType 
// {
//   red    = 0,
//   green  = 1, 
//   blue   = 2, 
//   white  = 3,
//   black  = 4,
//   yellow = 5, //(red+green)
//   purple = 6, //(red+blue)
//   cyan   = 7  //(green+blue)
// };


// class Tortoise // class containing the methods and the functions used for the tortoise recognition
// {
// private:

//   string            m_strLoadDirectory;
//   string            m_strSaveDirectory;
//   string            m_strTortoiseName;
//   string            m_strProcessID;
//   string            m_strResultsDataFileName;

//   CImg<int>	        m_imgOriginal;	
//   CImg<int>         m_imgRotatedOriginal;
//   CImg<int>         m_imgResultsData;
//   CImg<int>         m_imgCentralPathXCoordinates;

//   int               m_nWidth, m_nHeight;
//   int			          m_nRotationAngleToVerticalInDegrees;
//   int		    	      m_imgResultsDataSize;

//   plastronJunctions m_sJunctions;
//   plastronJunctions m_sLeftJunctions;
//   plastronJunctions m_sRightJunctions;

//   point2dCoor       m_sPlastronCentreOnRotatedImg;
//   point2dCoor       m_sCentralSeamHeadEnd;
//   point2dCoor       m_sCentralSeamTailEnd;

//   // Methods for tortoise recognition
//   void CentralSeamLocalization();
//   CImg<int> GetUniformlyResizedImg(CImg<int> img, int resizedMaxSizeInPix);
//   inline float GetImgResolutionCoeff(int width, int height);
//   CImg<unsigned char> GetEdgeImg(const CImg<int> &img, int lowerThreshold, int upperThreshold, double sigma);
//   CImg<unsigned char> GetDirectionOfEdges(const CImg<int> &edgeImg); 
//   CImg<unsigned char> EliminateOnePixelEdges(CImg<unsigned char> &imgToDoItOn);
//   CImg<unsigned char> EliminateTwoPixelsEdges(CImg<unsigned char> &imgToDoItOn);
//   void TortoiseContourLocalizationWithGHT();
//   void GetAbdJunctionsPreciseLocations(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg);
//   void MarkAJunction(CImg<int> &img, const point2dCoor &sJuntionToMark, const bool bSquareFilled, ColourType colour);
//   void GetFemToAnaJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions);
//   void GetHumToPecJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions);
//   void GetGulToHumJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions);
//   void GetHeadJunctionLocalizationWithGTH(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg);
//   void GetTailJunctionLocalizationWithGTH(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg);
//   CImg<unsigned char> SkeletonizationMinDist(CImg<unsigned char> imgEdges);
//   CImg<unsigned char> GetCentralSeam(const CImg<unsigned char> &imgEdges);
//   void GetCentralPathXCoordinatesForInnerJunctions(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg);
//   void GetGulToHumLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
//   void GetHumToPecLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
//   void GetFemToAnaLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
//   void GetAbdLefAndRightJunctions(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions,const  point2dCoor &sPlastronCentreOnRotatedImg);
//   void GetCentralPathXCoordinatesForHeadAndTailArea();
//   void GetTailJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);
//   void GetHeadJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter);

//   CImg<double> GetNormalizedJunctionDifferencesAndSeamLengths(const plastronJunctions &sLeftJunctions, const plastronJunctions &sRightJunctions);
// public:
                            
// 	Tortoise(string imgDirectory, string tortoiseName, string processID); 	//constructor
//   ~Tortoise(){}                                                           // destructor
//   void Recognition();
//   void Classification();
// };

//   void AccuracyEvaluation(string strLoadDirectory);
//   void kNNclassifier();
//   void Classification();
