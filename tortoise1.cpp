#include "tortoise.hpp"
#include "genHoughTrnf.h"
#include "tortoiseClassification.h"


  //~ //! constructor of Tortoise class
  //~ /**
      //~ \outline
      //~ - Loads file names and the original TG image. Checks if resultsData file exists, which
        //~ will contain data about partial results, and if resultsImg file exists, which will
        //~ show the partial results on an image.
  //~ **/
  //~ Tortoise::Tortoise(string imgDirectory, string tortoiseName, string processID)
  //~ {
    //~ m_strProcessID = processID;
    //~ m_strLoadDirectory.assign(imgDirectory);
    //~ m_strSaveDirectory.assign(imgDirectory);
    //~ m_strTortoiseName.assign(tortoiseName);
    //~ m_strResultsDataFileName.assign(m_strSaveDirectory+"resultsData"+m_strTortoiseName+".cimg");
//~
	  //~ try
	  //~ {
      //~ m_imgOriginal.assign( (m_strLoadDirectory+m_strTortoiseName+".pnm").c_str() ); //	loading original image and its parameters
		  //~ m_nWidth = m_imgOriginal.width();
		  //~ m_nHeight = m_imgOriginal.height();
//~
      //~ m_imgResultsDataSize = 20;
//~
      //~ if (FILE *file = fopen(m_strResultsDataFileName.c_str(), "r")) // if results data exists (data with partial results)
      //~ {
        //~ fclose(file);
        //~ m_imgResultsData.assign(m_strResultsDataFileName.c_str());
      //~ }
      //~ else // if not, create one
      //~ {
        //~ m_imgResultsData.assign(m_imgResultsDataSize,m_imgResultsDataSize,1,1, 0);
      //~ }
	  //~ }
	  //~ catch(CImgIOException&)
	  //~ {
		  //~ cout << "Error: unable to load the image: " << (m_strLoadDirectory+m_strTortoiseName+".pnm").c_str() << endl;
		  //~ throw;
	  //~ }
  //~ }


//------------	 Central Seam Localization	-------------------------------------------------------------------------------/
  CImg<int> Tortoise::GetUniformlyResizedImg(CImg<int> img, int resizedMaxSizeInPix) // resize the original image - with fixed ratio of sides
  {
    int width = img.width();
    int height = img.height();
    int widthResized = 0;
    int heightResized = 0;
    if(width>height)
    {
      double resizeRatio = resizedMaxSizeInPix*1.0/width;
      widthResized  = resizedMaxSizeInPix;
      heightResized = (int) (height*resizeRatio);
    }
    else
    {
      double resizeRatio   = resizedMaxSizeInPix*1.0/height;
      heightResized  = resizedMaxSizeInPix;
      widthResized = (int) (width*resizeRatio);
    }
    return img.get_resize(widthResized,heightResized);
  }

  inline float Tortoise::GetImgResolutionCoeff(int width, int height)
  {
    return (sqrt(1.0f*width*height)/15.55f - 3.67f);
  }

  CImg<unsigned char> Tortoise::GetEdgeImg(const CImg<int> &img, int lowerThreshold, int upperThreshold, double sigma)
{ // canny edge detection, 1-edges, 0-the rest, for any CImg image and the function returns also a CImg image
//	preparing a fileName for saving and loading
	string edgesFileName (m_strSaveDirectory + "edgesWorkingFile" + m_strProcessID + ".pnm");

//	saving CImg image and loading it in OpenCV
	img.save(edgesFileName.c_str());
	cv::Mat imgOpenCV;
	imgOpenCV = cv::imread(edgesFileName.c_str());

//	Pass the image to gray
	cv::Mat src_gray;
    cv::cvtColor( imgOpenCV, src_gray, CV_RGB2GRAY );

//	blur with a gaussian filter
	cv::Size ksize;
	ksize.width = 0;
	ksize.height = 0;
	cv::Mat src_gray_blured;
	if(sigma > 0) cv::GaussianBlur(src_gray, src_gray_blured, ksize, sigma);
	else src_gray_blured = src_gray;

//	Canny edge detector
	cv::Mat edgesOpenCV;
	cv::Canny( src_gray_blured, edgesOpenCV, lowerThreshold, upperThreshold, 3 );

//	saving OpenCV image and loading it in CImg
	cv::imwrite(edgesFileName.c_str(), edgesOpenCV);
	CImg<unsigned char> imgEdges(edgesFileName.c_str());

//	1-edges, 0-the rest
	return imgEdges.threshold(1);
}

  CImg<unsigned char> Tortoise::GetDirectionOfEdges(const CImg<int> &edgeImg)
{
	CImg<unsigned char> imgWithDirectionOfEdges(edgeImg);
  imgWithDirectionOfEdges.fill(0);
	int width = edgeImg.width();
	int height = edgeImg.height();
	for(int x = 2; x < width-2; x++)
	{
		for(int y = 2; y < height-2; y++)
		{
			if(edgeImg(x,y))
			{
				int dirVertical =  edgeImg(x-1,y-1) + edgeImg(x,y-1) + edgeImg(x+1,y-1)
									       + edgeImg(x-1,y+1) + edgeImg(x,y+1) + edgeImg(x+1,y+1);
				int dirHorizontal =  edgeImg(x-1,y-1) + edgeImg(x+1,y-1)
									         + edgeImg(x-1,y  ) +                  edgeImg(x+1,y  )
									         + edgeImg(x-1,y+1) + edgeImg(x+1,y+1);;
				int dirDiagonalRightUp =  edgeImg(x,y-1) + edgeImg(x+1,y-1)
									              + edgeImg(x-1,y  ) +                  edgeImg(x+1,y  )
									              + edgeImg(x-1,y+1) + edgeImg(x,y+1);
				int dirDiagonalLeftUp =  edgeImg(x-1,y-1) + edgeImg(x,y-1)
									             + edgeImg(x-1,y  ) +                  edgeImg(x+1,y  )
									             + edgeImg(x,y+1) + edgeImg(x+1,y+1);
				if (dirHorizontal >= dirVertical && dirHorizontal >= dirDiagonalRightUp && dirHorizontal >= dirDiagonalLeftUp) imgWithDirectionOfEdges(x,y) = 2;
        else if (dirVertical >= dirHorizontal && dirVertical >= dirDiagonalRightUp && dirVertical >= dirDiagonalLeftUp) imgWithDirectionOfEdges(x,y) = 3;
				else if (dirDiagonalRightUp >= dirHorizontal && dirDiagonalRightUp >= dirVertical && dirDiagonalRightUp >= dirDiagonalLeftUp) imgWithDirectionOfEdges(x,y) = 1;
				else if (dirDiagonalLeftUp >= dirHorizontal && dirDiagonalLeftUp >= dirDiagonalRightUp && dirDiagonalLeftUp >= dirVertical) imgWithDirectionOfEdges(x,y) = 1;
			}
		}
	}
	return imgWithDirectionOfEdges;
}

  CImg<unsigned char> Tortoise::EliminateOnePixelEdges(CImg<unsigned char> &imgToDoItOn)
{
	int width = imgToDoItOn.width();
	int height = imgToDoItOn.height();
  imgToDoItOn.threshold(1);
	for(int y = 1; y < height - 1; y++)
	{
		for(int x = 1; x < width-1; x++)
		{
			if(imgToDoItOn(x,y) == 1)
			{ // eliminate one pixel edges
				if(imgToDoItOn(x-1,y)==0 && imgToDoItOn(x-1,y-1)==0 && imgToDoItOn(x-1,y+1)==0 && imgToDoItOn(x,y-1)==0 &&
					imgToDoItOn(x+1,y)==0 && imgToDoItOn(x+1,y-1)==0 && imgToDoItOn(x+1,y+1)==0 && imgToDoItOn(x,y+1)==0) { imgToDoItOn(x,y) = 0; }

			}
		}
	}
	return imgToDoItOn;
}

  CImg<unsigned char> Tortoise::EliminateTwoPixelsEdges(CImg<unsigned char> &imgToDoItOn)
  {
	  int width = imgToDoItOn.width();
	  int height = imgToDoItOn.height();
    imgToDoItOn.threshold(1);
	  for(int y = 2; y < height - 2; y++)
	  {
		  for(int x = 2; x < width-2; x++)
		  {
			  if(imgToDoItOn(x,y) == 1)
			  {
				  CImg<unsigned char> imgCrop(imgToDoItOn.get_crop(x-1,y-1,x+1,y+1));
				  if(imgCrop.sum() <= 2)
				  {
					  imgCrop(1,1)=0;
					  int X = static_cast<int>(imgCrop.get_stats()(8))-1 + x;
					  int Y = static_cast<int>(imgCrop.get_stats()(9))-1 + y;
					  if(imgToDoItOn.get_crop(X-1,Y-1,X+1,Y+1).sum() <=2)
					  {
						  imgToDoItOn(x,y) = 0;
						  imgToDoItOn(X,Y) = 0;
					  }

				  }
			  }
		  }
	  }
	  return imgToDoItOn;
  }

  void Tortoise::CentralSeamLocalization()
  {
    CImg<int> imgResizedOriginal(GetUniformlyResizedImg(m_imgOriginal,500));
	  int nWidth = imgResizedOriginal.width();
	  int nHeight = imgResizedOriginal.height();
    int halfOfTheSizeOfSquareToBeCroped = static_cast<int>( sqrt(2.0)*min(nWidth/2,nHeight/2)/24.0 );
    int heightOfCrop = 12*halfOfTheSizeOfSquareToBeCroped;
    int widthOfCrop = halfOfTheSizeOfSquareToBeCroped/2;
    cout << halfOfTheSizeOfSquareToBeCroped << endl;
    double angleDegrsMax = 180;
    double angleDegrsStep = 5;
    int numOfAngles = static_cast<int>(ceil(angleDegrsMax/angleDegrsStep));
    CImg<int> HoughSpace(numOfAngles,4,1,1, 0);
	  int plastronCenterX = nWidth/2;
	  int plastronCenterY = nHeight/2;
    float koef = GetImgResolutionCoeff(nWidth,nHeight);
    int circleRadius =  max(static_cast<int>(koef*3),1);

    for(int angleNum = 0; angleNum < numOfAngles; ++angleNum)
    {
      int angleDegrs = -90;
      angleDegrs += static_cast<int>(angleNum*angleDegrsStep);
      if((angleDegrs > -75 && angleDegrs < -15) || (angleDegrs > 15 && angleDegrs < 75)) {}
      else
      {
        HoughSpace(angleNum,0) = angleDegrs;
        cout << angleNum <<  " " << angleDegrs << endl;
        CImg<unsigned char> imgRotated(GetDirectionOfEdges(EliminateTwoPixelsEdges(EliminateOnePixelEdges(GetEdgeImg(imgResizedOriginal.get_rotate(static_cast<float>(angleDegrs)),100,0,0)))).threshold(3));
	      for(int x = imgRotated.width()/2-circleRadius; x <= imgRotated.width()/2+circleRadius; ++x)
	      {
          int y = imgRotated.height()/2;
          {
            CImg<unsigned char> imgCropped(imgRotated.get_crop(x-widthOfCrop,y-heightOfCrop,x+widthOfCrop,y+heightOfCrop));
            int mirrorMatch = static_cast<int>(imgCropped.sum());
            if(mirrorMatch > HoughSpace(angleNum,1))
            {
              HoughSpace(angleNum,1) = mirrorMatch;
              HoughSpace(angleNum,2) = x;
              HoughSpace(angleNum,3) = y;
            }
          }
        }
        //imgRotated.display();
      }
	  }
    //HoughSpace.display();
    int maxHough = HoughSpace.get_row(1).max();
    if(maxHough > 0)
    {
      int minHoughXcoor = static_cast<int>(HoughSpace.get_row(1).get_stats()(8));
      int xMin = HoughSpace(minHoughXcoor,2);
      int yMin = HoughSpace(minHoughXcoor,3);
      int minHoughAngle = HoughSpace(minHoughXcoor,0);
      m_nRotationAngleToVerticalInDegrees = minHoughAngle;

      m_sPlastronCentreOnRotatedImg.xCoor = xMin*m_nWidth/nWidth;
      m_sPlastronCentreOnRotatedImg.yCoor = yMin*m_nHeight/nHeight;
    }
    else
    {
      cout << "Error in Hough!" << endl;
    }
  }

  void Tortoise::TortoiseContourLocalizationWithGHT() // GHT = Generalized Hough Transform
  {
    m_imgRotatedOriginal.assign(m_imgOriginal.get_rotate(static_cast<float>(m_nRotationAngleToVerticalInDegrees)));
    int nRotImgWidht = m_imgRotatedOriginal.width();
    int nRotImgHeight = m_imgRotatedOriginal.height();
    int nResizeMaxPixelSize = 400;
    float fSigmaBlurCoef = 0.5f;
    CImg<int> imgResizedRotated(GetUniformlyResizedImg(m_imgRotatedOriginal.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
    int nResRotImgWidth = imgResizedRotated.width();
    int nResRotImgHeight= imgResizedRotated.height();
    const float fImgResCoeff = GetImgResolutionCoeff(m_nWidth,m_nHeight);
    const int nCentralSeamAreaNarrowWidht = static_cast<int>(0.5f*fImgResCoeff);
    float fLeftAccumBndryXCoor = static_cast<float>((m_sPlastronCentreOnRotatedImg.xCoor - nCentralSeamAreaNarrowWidht)*nResRotImgWidth)/nRotImgWidht;
    float fRightAccumBndryXCoor = static_cast<float>((m_sPlastronCentreOnRotatedImg.xCoor + nCentralSeamAreaNarrowWidht)*nResRotImgWidth)/nRotImgWidht;
    CImg<unsigned char> imgEdges(GetEdgeImg(imgResizedRotated,50,100,0));
    imgEdges *= 255;
    string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
    string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
    imgResizedRotated.save(strGHTImgFileName.c_str());
    imgEdges.save(strGHTEdgesFileName.c_str());
    GenHoughTrnf ght("ght\\tortoiseTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
    // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
    ght.Init(150,300,10,4,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
    ght.runGHT(1);
    point2dCoor sAccumMax = {0,0};
    int nRotationAngleMax = 0;
    float fWidthRatioMax  = 0.0f;
    int nPlastronWidht = 0;
    ght.bestCandidate(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht).save(("ght\\GHTRes" + m_strTortoiseName + "WithRestriction.bmp").c_str());
    // the  m_sPlastronCentreOnRotatedImg.xCoor stays the same
    m_sPlastronCentreOnRotatedImg.yCoor = sAccumMax.yCoor*nRotImgHeight/nResRotImgHeight;
    if(nRotationAngleMax>90 && nRotationAngleMax<270)
    {
      m_nRotationAngleToVerticalInDegrees =  (m_nRotationAngleToVerticalInDegrees+180) % 360;
      //m_imgRotatedOriginal.rotate(180);
      m_sPlastronCentreOnRotatedImg.xCoor = nRotImgWidht  - m_sPlastronCentreOnRotatedImg.xCoor;
      m_sPlastronCentreOnRotatedImg.yCoor = nRotImgHeight - m_sPlastronCentreOnRotatedImg.yCoor;
    }
    // plastron junctions position based on template contour
    point2dCoor sCenterTortoiseContour = {64,125};
    plastronJunctions sContourJunctions = {{64,58},{64,74},{64,97},{64,105},{64,155},{64,168},{64,185}};
    // s1Head
    m_sJunctions.s1Head.xCoor = m_sPlastronCentreOnRotatedImg.xCoor;
    m_sJunctions.s1Head.yCoor = m_sPlastronCentreOnRotatedImg.yCoor - static_cast<int>((sCenterTortoiseContour.yCoor-sContourJunctions.s1Head.yCoor)*fWidthRatioMax*nRotImgHeight/nResRotImgHeight);
    m_sJunctions.s1Head.yCoor = max(min(m_sJunctions.s1Head.yCoor,m_imgRotatedOriginal.height()-1),0);
    // s2GulToHum
    m_sJunctions.s2GulToHum.xCoor = m_sPlastronCentreOnRotatedImg.xCoor;
    m_sJunctions.s2GulToHum.yCoor = m_sPlastronCentreOnRotatedImg.yCoor - static_cast<int>((sCenterTortoiseContour.yCoor-sContourJunctions.s2GulToHum.yCoor)*fWidthRatioMax*nRotImgHeight/nResRotImgHeight);
    m_sJunctions.s2GulToHum.yCoor = max(min(m_sJunctions.s2GulToHum.yCoor,m_imgRotatedOriginal.height()-1),0);
    // s3HumToPec
    m_sJunctions.s3HumToPec.xCoor = m_sPlastronCentreOnRotatedImg.xCoor;
    m_sJunctions.s3HumToPec.yCoor = m_sPlastronCentreOnRotatedImg.yCoor - static_cast<int>((sCenterTortoiseContour.yCoor-sContourJunctions.s3HumToPec.yCoor)*fWidthRatioMax*nRotImgHeight/nResRotImgHeight);
    m_sJunctions.s3HumToPec.yCoor = max(min(m_sJunctions.s3HumToPec.yCoor,m_imgRotatedOriginal.height()-1),0);
    // s4PecToAbd
    m_sJunctions.s4PecToAbd.xCoor = m_sPlastronCentreOnRotatedImg.xCoor;
    m_sJunctions.s4PecToAbd.yCoor = m_sPlastronCentreOnRotatedImg.yCoor - static_cast<int>((sCenterTortoiseContour.yCoor-sContourJunctions.s4PecToAbd.yCoor)*fWidthRatioMax*nRotImgHeight/nResRotImgHeight);
    m_sJunctions.s4PecToAbd.yCoor = max(min(m_sJunctions.s4PecToAbd.yCoor,m_imgRotatedOriginal.height()-1),0);
    // s5AbdToFem
    m_sJunctions.s5AbdToFem.xCoor = m_sPlastronCentreOnRotatedImg.xCoor;
    m_sJunctions.s5AbdToFem.yCoor = m_sPlastronCentreOnRotatedImg.yCoor - static_cast<int>((sCenterTortoiseContour.yCoor-sContourJunctions.s5AbdToFem.yCoor)*fWidthRatioMax*nRotImgHeight/nResRotImgHeight);
    m_sJunctions.s5AbdToFem.yCoor = max(min(m_sJunctions.s5AbdToFem.yCoor,m_imgRotatedOriginal.height()-1),0);
    // s6FemToAna
    m_sJunctions.s6FemToAna.xCoor = m_sPlastronCentreOnRotatedImg.xCoor;
    m_sJunctions.s6FemToAna.yCoor = m_sPlastronCentreOnRotatedImg.yCoor - static_cast<int>((sCenterTortoiseContour.yCoor-sContourJunctions.s6FemToAna.yCoor)*fWidthRatioMax*nRotImgHeight/nResRotImgHeight);
    m_sJunctions.s6FemToAna.yCoor = max(min(m_sJunctions.s6FemToAna.yCoor,m_imgRotatedOriginal.height()-1),0);
    // s7Tail
    m_sJunctions.s7Tail.xCoor = m_sPlastronCentreOnRotatedImg.xCoor;
    m_sJunctions.s7Tail.yCoor = m_sPlastronCentreOnRotatedImg.yCoor - static_cast<int>((sCenterTortoiseContour.yCoor-sContourJunctions.s7Tail.yCoor)*fWidthRatioMax*nRotImgHeight/nResRotImgHeight);
    m_sJunctions.s7Tail.yCoor = max(min(m_sJunctions.s7Tail.yCoor,m_imgRotatedOriginal.height()-1),0);
  }

    //------------ Junction localization	-------------------------------------------------------------------------------/
  void Tortoise::GetAbdJunctionsPreciseLocations(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions,
                                                 const point2dCoor &sPlastronCentreOnRotatedImg)
  {
  //using edges on a wider central-seam stripe
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(
      min(m_imgRotatedOriginal.width()-1,max(sPlastronCentreOnRotatedImg.xCoor-nAbdSeamLength,0)),sJunctions.s1Head.yCoor,
      min(m_imgRotatedOriginal.width()-1,max(sPlastronCentreOnRotatedImg.xCoor+nAbdSeamLength,0)),sJunctions.s7Tail.yCoor));
    imgPlastron = GetUniformlyResizedImg(imgPlastron,500);
    int nCentreYCoor = (sPlastronCentreOnRotatedImg.yCoor - sJunctions.s1Head.yCoor)*imgPlastron.height()/(sJunctions.s7Tail.yCoor-sJunctions.s1Head.yCoor);
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,0,0));
    //imgEdges.display();
	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //every edge pix is part of a some path - all edge pixs are assigned an index of the path they belong to
	  CImg<int> sidePathsIndexes(imgEdges.get_fill(0));
  //pixs on the cetral path have 2-column index array, beacause there meet paths from left and right side  - they meet in the middle: x = nImgEdgesWidth/2
	  CImg<int> indexesOnTheCentralPath(2,nImgEdgesHeight,1,1, 0);
  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
  //prefilling - the sides (5 the most left columns) of sidePathsFromLeft and sidePathsIndexes
	  int nPathIndex = 2;
    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
		  int x = 4;
		  {
  //prefilling - for all edge pixs in the area
			  if(imgEdges(x,y))
			  {
  //prefilling - initial value (distance) of sidePathsFromLeft is set to 2
				  sidePathsFromLeft(x,y) = 2;
  //prefilling - edge pixs in the area are all assigned an index in sidePathsIndexes, which they will pass on to pixs that will be in the same path
				  sidePathsIndexes(x,y) = nPathIndex;
			  }
		  }
		  nPathIndex++;
	  }
  //prefilling - the sides (5 the most right columns) of sidePathsFromRight and sidePathsIndexes
	  nPathIndex += 1000;
	  for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      int x = nImgEdgesWidth-5;
		  {
			  if(imgEdges(x,y))
			  {
  //prefilling - initial value (distance) of sidePathsFromRight is set to 2
				  sidePathsFromRight(x,y) = 2;
  //prefilling - edge pixs in the area are all assigned an index in sidePathsIndexes, which they will pass on to pixs that will be in the same path
				  sidePathsIndexes(x,y) = nPathIndex;
			  }
		  }
      nPathIndex++;
	  }
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 3;
      distanceMask(6,y) += 3;
      distanceMask(1,y) += 2;
      distanceMask(5,y) += 2;
      distanceMask(2,y) += 1;
      distanceMask(4,y) += 1;
    }
	  //distanceMask.display();
  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) && x<=nImgEdgesWidth/2)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				    sidePathsIndexes(x,y) = sidePathsIndexes(x-5+static_cast<int>(cropAreaDistances.get_stats()(4)),y-3+static_cast<int>(cropAreaDistances.get_stats()(5)));
  //filling left side paths - if the pix is on the central path, the index is saved in indexesOnTheBestCentralPath
            if(x==nImgEdgesWidth/2) indexesOnTheCentralPath(0,y) = sidePathsIndexes(x,y);
				  }
			  }
		  }
	  }
  //  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
	 distanceMask.rotate(180);
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
        if(imgEdges(x,y) && x >= nImgEdgesWidth/2)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
					  sidePathsIndexes(x,y) = sidePathsIndexes(x+1+static_cast<int>(cropAreaDistances.get_stats()(4)),y-3+static_cast<int>(cropAreaDistances.get_stats()(5)));
            if(x==nImgEdgesWidth/2) indexesOnTheCentralPath(1,y) = sidePathsIndexes(x,y);
				  }
			  }
		  }
	  }
	  //(sidePathsFromLeft,sidePathsFromRight).display();
  //postprocessing - set values 10000 to 0 in sidePathsFromLeft and sidePathsFromRight
    for(int y = 0; y < nImgEdgesHeight; y++)
	  {
      for(int x = 0; x < nImgEdgesWidth; x++)
		  {
        if(sidePathsFromLeft(x,y)  == nMaxPathValue) sidePathsFromLeft(x,y) = 0;
        if(sidePathsFromRight(x,y) == nMaxPathValue) sidePathsFromRight(x,y) = 0;
		  }
	  }
	 //(sidePathsFromLeft+sidePathsFromRight).display("sidePaths");
  //postprocessing - CentralSideCrosses will contain sum of 2x values of distances of pixs on the central seam of the left and right paths
    CImg<int> CentralSideCrosses(1,nImgEdgesHeight,1,1, nMaxPathValue);
  //postprocessing - for all central path pixs near the tortoise center
    for(int y = 6; y < nImgEdgesHeight-6; ++y)
	  {
      int x = nImgEdgesWidth/2;
  //postprocessing - if there are paths from left and right that end in the central-path pix
		  if(sidePathsFromLeft(x,y)>0 && sidePathsFromRight(x,y)>0)
		  {
  //postprocessing - the values are added to CentralSideCrosses
			  CentralSideCrosses(y) = sidePathsFromLeft(x,y) + sidePathsFromRight(x,y);
		  }
  //postprocessing - else, left and right indexes of the central-path pix are set to 0
		  else
		  {
			  indexesOnTheCentralPath(0,y) = 0;
			  indexesOnTheCentralPath(1,y) = 0;
		  }
	  }
   // (CentralSideCrosses,indexesOnTheCentralPath,(sidePathsFromLeft+sidePathsFromRight)).display();
  //finding best side paths - setting number of pixs that are to be discarded, being next (based on index) to a pix where we will find the local minimum (of distances in CentralSideCrosses)
  //                        - with 3000pixs height img time constant 0.025 (based on real images) so that we get ~ 70pixs
  //                        - (it has to be so much that we eliminate all neighbouring paths with similar distances, but so that we do not eliminate other correct side paths)
    int halfIntervalOfIndexes = static_cast<int>(nImgEdgesHeight*0.1);
  //finding best side paths - countOfCrossing counts side paths and is used to index the paths
	  int countOfCrossing = 0;
  //finding best side paths - we want no more than 9 paths (5-side paths, 2-end and beggining of the plastron, 2-ends of tail and head/ends of the shell)
	  const int permittedNumberOfPaths = 9;
	  CImg<int> yCoorOfTheJunctions(1,permittedNumberOfPaths,1,1, 0);
    bool bEnd = 0;
    int nNumOfUpJunctions = 0;
    int nNumOfDownJunctions = 0;
  //finding best side paths - going though CentralSideCrosses until there is no pix with smaller distance than upperSupremum or the is permittedNumberOfPaths reached
   // (CentralSideCrosses,indexesOnTheCentralPath,yCoorOfTheJunctions,valuesOfTheJunctions).display();
	  while(CentralSideCrosses.min()<nMaxPathValue && countOfCrossing < permittedNumberOfPaths && !bEnd)
	  {
  //finding best side paths - finding the current minimum in CentralSideCrosses (distances of left and right paths)
		  int yMin = static_cast<int>(CentralSideCrosses.get_stats()(5));
      int currentMin =  CentralSideCrosses.min();
      if(yMin<nCentreYCoor) nNumOfUpJunctions++;
      else nNumOfDownJunctions++;
      if(nNumOfUpJunctions==2 && nNumOfDownJunctions==2) bEnd=1;
      if((nNumOfUpJunctions<=2&&yMin<nCentreYCoor) || (nNumOfDownJunctions<=2&&yMin>nCentreYCoor)) // spatne neco jinak || &&
      {
		    yCoorOfTheJunctions(countOfCrossing) = yMin;
      }
    //finding best side paths - getting left index of the pix with the minimum and its y coordinate
		    int leftIndexOfMinPix = indexesOnTheCentralPath(0,yMin);
    //finding best side paths - getting right index of the pix with the minimum and its y coordinate
		    int rightIndexOfMinPix = indexesOnTheCentralPath(1,yMin);
    //finding best side paths - deleting indexes close to index of the pix with the current min value
        for(int y = 0; y < nImgEdgesHeight; y++)
		    {
          if( (indexesOnTheCentralPath(0,y) <= leftIndexOfMinPix + halfIntervalOfIndexes) &&
				    (indexesOnTheCentralPath(0,y) >= max(1,leftIndexOfMinPix - halfIntervalOfIndexes)))
			    {
    //finding best side paths - index to be deleted from indexesOnTheBestCentralPath
				    indexesOnTheCentralPath(0,y) = 0;
            indexesOnTheCentralPath(1,y) = 0;
    //finding best side paths - coresponding distances in CentralSideCrosses set to upperSupremum (cannot become minimum in the future)
            CentralSideCrosses(y) = nMaxPathValue;
			    }
    //finding best side paths - the same for the right side indexes
          else if( indexesOnTheCentralPath(1,y) <= rightIndexOfMinPix + halfIntervalOfIndexes &&
				         indexesOnTheCentralPath(1,y) >= max(1,rightIndexOfMinPix - halfIntervalOfIndexes) )
			    {
            indexesOnTheCentralPath(0,y) = 0;
				    indexesOnTheCentralPath(1,y) = 0;
				    CentralSideCrosses(y) = nMaxPathValue;
			    }
		    }
    //finding best side paths - numbering side paths
		    CentralSideCrosses(yMin) = nMaxPathValue + currentMin + 1;

    //finding best side paths - rising the count of side paths
        //(CentralSideCrosses,indexesOnTheCentralPath).display();
       //(CentralSideCrosses,indexesOnTheCentralPath,(sidePathsFromLeft+sidePathsFromRight)).display();
      countOfCrossing++;
	  }
	 //(sidePathsFromLeft+sidePathsFromRight).display("sidePaths");

   //  (CentralSideCrosses,indexesOnTheCentralPath,(sidePathsFromLeft+sidePathsFromRight),yCoorOfTheJunctions).display();
    for(int ii = 0; ii < permittedNumberOfPaths; ++ii)
    {
      if(yCoorOfTheJunctions(ii)>0)
      {
        yCoorOfTheJunctions(ii) = sJunctions.s1Head.yCoor + yCoorOfTheJunctions(ii)
          *(sJunctions.s7Tail.yCoor - sJunctions.s1Head.yCoor)/nImgEdgesHeight;
      }
    }
    const int nPositionOfPecToAbdInTheJunctions = static_cast<int>( (yCoorOfTheJunctions.get_mul((-yCoorOfTheJunctions).get_threshold(-sPlastronCentreOnRotatedImg.yCoor))-sPlastronCentreOnRotatedImg.yCoor).abs().get_stats()(5) );
    const int nPositionOfAbdToFemInTheJunctions = static_cast<int>( (yCoorOfTheJunctions.get_mul((yCoorOfTheJunctions).get_threshold(sPlastronCentreOnRotatedImg.yCoor))-sPlastronCentreOnRotatedImg.yCoor).abs().get_stats()(5) );
    m_sJunctions.s4PecToAbd.yCoor = yCoorOfTheJunctions(nPositionOfPecToAbdInTheJunctions);
    m_sJunctions.s5AbdToFem.yCoor = yCoorOfTheJunctions(nPositionOfAbdToFemInTheJunctions);
    m_sPlastronCentreOnRotatedImg.yCoor = (m_sJunctions.s5AbdToFem.yCoor + m_sJunctions.s4PecToAbd.yCoor)/2;
  }

  void Tortoise::GetFemToAnaJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal,
                                                      const plastronJunctions &sJunctions)
  {
  //using edges on a wider central-seam stripe
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(
      max(0,min(sJunctions.s1Head.xCoor-nAbdSeamLength*3/8,imgRotatedOriginal.width()-1)),max(0,min(sJunctions.s5AbdToFem.yCoor+nAbdSeamLength*4/30,imgRotatedOriginal.height()-1)),
      max(0,min(sJunctions.s1Head.xCoor+nAbdSeamLength*3/8,imgRotatedOriginal.width()-1)),max(0,min(imgRotatedOriginal.height()-1,sJunctions.s5AbdToFem.yCoor+nAbdSeamLength*2/3))));
    //imgPlastron = GetUniformlyResizedImg(imgPlastron,500);
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,0,0));
    //imgEdges.display();
	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 4;
      distanceMask(1,y) += 3;
      distanceMask(2,y) += 2;
      distanceMask(3,y) += 0;
      distanceMask(4,y) += 0;
      distanceMask(5,y) += 1;
      distanceMask(6,y) += 2;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
  //  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
   sidePathsFromLeft += sidePathsFromRight;
   //postprocessing - CentralSideCrosses will contain sum of 2x values of distances of pixs on the central seam of the left and right paths
   CImg<int> CentralSideCrosses(sidePathsFromLeft.get_column(nImgEdgesWidth/2));
   m_sJunctions.s6FemToAna.yCoor = static_cast<int>(CentralSideCrosses.get_stats()(5))+sJunctions.s5AbdToFem.yCoor+nAbdSeamLength*4/30;
  }

  void Tortoise::GetHumToPecJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal,
                                                      const plastronJunctions &sJunctions)
  {
  //using edges on a wider central-seam stripe
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(
      max(0,min(sJunctions.s1Head.xCoor-nAbdSeamLength*2/3,imgRotatedOriginal.width()-1)),max(0,min(sJunctions.s4PecToAbd.yCoor-nAbdSeamLength*2/3,imgRotatedOriginal.height()-1)),
      max(0,min(sJunctions.s1Head.xCoor+nAbdSeamLength*2/3,imgRotatedOriginal.width()-1)),max(0,min(imgRotatedOriginal.height()-1,sJunctions.s4PecToAbd.yCoor))));
    //imgPlastron = GetUniformlyResizedImg(imgPlastron,500);
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,0,0));
    //imgEdges.display();
	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
    //imgEdges.display();
    const int A = nImgEdgesWidth-1;
    const int B = nImgEdgesHeight-1;
	  for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      for(int x = 0; x < nImgEdgesWidth; ++x)
		  {
        if( imgEdges(x,y) && (y>=B*x/(4*A)+B*7/8 || y>=B*(A-x)/(4*A)+B*7/8) )
			  {
  //prefilling - initial value (distance) of sidePathsFromRight is set to 2
				  imgEdges(x,y) = 0;
			  }
		  }
	  }
   // imgEdges.display();
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 3;
      distanceMask(1,y) += 2;
      distanceMask(2,y) += 0;
      distanceMask(3,y) += 0;
      distanceMask(4,y) += 1;
      distanceMask(5,y) += 2;
      distanceMask(6,y) += 3;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
//  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
   sidePathsFromLeft += sidePathsFromRight;
   //postprocessing - CentralSideCrosses will contain sum of 2x values of distances of pixs on the central seam of the left and right paths
   CImg<int> CentralSideCrosses(sidePathsFromLeft.get_column(nImgEdgesWidth/2));
   //(sidePathsFromLeft,CentralSideCrosses).display();
   m_sJunctions.s3HumToPec.yCoor = static_cast<int>(CentralSideCrosses.get_stats()(5))+max(0,min(sJunctions.s4PecToAbd.yCoor-nAbdSeamLength*2/3,imgRotatedOriginal.height()-1));
  }

  void Tortoise::MarkAJunction(CImg<int> &img, const point2dCoor &sJuntionToMark, const bool bSquareFilled, ColourType colour)
  {
    int nDrawnAreaHeight = max(4,4*max(img.height(),img.width())/500);
    int nDrawnAreaContourMargin = max(2,2*max(img.height(),img.width())/500);
    int nDrawnAreaWidth = max(6,6*max(img.height(),img.width())/500);
    int nRed = 0;
    int nGreen = 0;
    int nBlue = 0;
    if(colour==red) nRed = 255;
    else if(colour==green) nGreen = 255;
    else if(colour==blue) nBlue = 255;
    else if(colour==white) {nRed = 255;nGreen = 255;nBlue = 255;}
    else if(colour==black) {}
    else if(colour==yellow) {nRed = 255;nGreen = 255;}
    else if(colour==purple) {nRed = 255;nBlue = 255;}
    else if(colour==cyan) {nGreen = 255;nBlue = 255;}
    for(int ii = -nDrawnAreaWidth; ii <= nDrawnAreaWidth; ++ii)
    {
      for(int jj = -nDrawnAreaHeight; jj <= nDrawnAreaHeight; ++jj)
      {
        if(bSquareFilled)
        {
          int x = sJuntionToMark.xCoor+ii;
          int y = max(0,min(sJuntionToMark.yCoor+jj,img.height()-1));
          img(x,y,0) = nRed;
          img(x,y,1) = nGreen;
          img(x,y,2) = nBlue;
        }
        else if(!(ii>=-nDrawnAreaWidth+nDrawnAreaContourMargin&&ii<=nDrawnAreaWidth-nDrawnAreaContourMargin&&
          jj<=nDrawnAreaHeight-nDrawnAreaContourMargin&&jj>=-nDrawnAreaHeight+nDrawnAreaContourMargin))
        {
          int x = sJuntionToMark.xCoor+ii;
          int y = max(0,min(sJuntionToMark.yCoor+jj,img.height()-1));
          img(x,y,0) = nRed;
          img(x,y,1) = nGreen;
          img(x,y,2) = nBlue;
        }
      }
    }
  }

  void Tortoise::GetGulToHumJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal,
                                                      const plastronJunctions &sJunctions)
  {
  //using edges on a wider central-seam stripe
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(
      max(0,min(sJunctions.s1Head.xCoor-nAbdSeamLength/4,imgRotatedOriginal.width()-1)),max(0,min(sJunctions.s4PecToAbd.yCoor-nAbdSeamLength,imgRotatedOriginal.height()-1)),
      max(0,min(sJunctions.s1Head.xCoor+nAbdSeamLength/4,imgRotatedOriginal.width()-1)),max(0,min(imgRotatedOriginal.height()-1,sJunctions.s3HumToPec.yCoor-nAbdSeamLength/4))));
    //imgPlastron = GetUniformlyResizedImg(imgPlastron,500);
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,0,0));
    //imgEdges.display();
	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
   // imgEdges.display();
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 2;
      distanceMask(1,y) += 1;
      distanceMask(2,y) += 0;
      distanceMask(3,y) += 1;
      distanceMask(4,y) += 2;
      distanceMask(5,y) += 3;
      distanceMask(6,y) += 4;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
//  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
   //for(int x = nImgEdgesWidth-1; x >= 0; --x)
	  //{
   //   for(int y = 0; y < nImgEdgesHeight; ++y)
		 // {
   //     if(sidePathsFromRight(x,y)==nMaxPathValue)
			//  {
   //         sidePathsFromRight(x,y) = 0;
   //
			//  }
   //     if(sidePathsFromLeft(x,y)==nMaxPathValue) sidePathsFromLeft(x,y) = 0;
		 // }
	  //}
   //(sidePathsFromLeft,sidePathsFromRight).display();
   sidePathsFromLeft += sidePathsFromRight;
   //postprocessing - CentralSideCrosses will contain sum of 2x values of distances of pixs on the central seam of the left and right paths
   CImg<int> CentralSideCrosses(sidePathsFromLeft.get_column(nImgEdgesWidth/2));
   //(sidePathsFromLeft,CentralSideCrosses).display();
   m_sJunctions.s2GulToHum.yCoor = static_cast<int>(CentralSideCrosses.get_stats()(5))+max(0,min(sJunctions.s4PecToAbd.yCoor-nAbdSeamLength,imgRotatedOriginal.height()-1));
  }

  void Tortoise::Recognition()  //! Here are the methods for TG recognition called one by one.
  {
    CImg<int> resized(GetUniformlyResizedImg(m_imgOriginal,300)),edgesImg(GetEdgeImg(resized,0,0,0));
    resized.save("resized.bmp");
    edgesImg.save("edges.bmp");
    (resized,edgesImg).display();
    if (m_imgResultsData(0,0) == 1) // central seam localization
    {
      m_nRotationAngleToVerticalInDegrees  =  static_cast<int>(m_imgResultsData(1,0));
      m_sPlastronCentreOnRotatedImg.xCoor  =  static_cast<int>(m_imgResultsData(2,0));
      m_sPlastronCentreOnRotatedImg.yCoor  =  static_cast<int>(m_imgResultsData(3,0));
    }
    //else
    //{
    //  CentralSeamLocalization(); // variables listed below are modified by this method
    //  m_imgResultsData(0,0) = 1;
    //  m_imgResultsData(1,0) = m_nRotationAngleToVerticalInDegrees;
    //  m_imgResultsData(2,0) = m_sPlastronCentreOnRotatedImg.xCoor;
    //  m_imgResultsData(3,0) = m_sPlastronCentreOnRotatedImg.yCoor;

    //  cout << "Central seam localization has finished!" << endl;
    //  m_imgResultsData.save(m_strResultsDataFileName.c_str());
    //}

    if(m_imgResultsData(0,1) == 1)
    {
      m_sPlastronCentreOnRotatedImg.xCoor = static_cast<int>(m_imgResultsData(1,1));
      m_sPlastronCentreOnRotatedImg.yCoor = static_cast<int>(m_imgResultsData(2,1));
      m_nRotationAngleToVerticalInDegrees = static_cast<int>(m_imgResultsData(4,1));
      // s1Head
      m_sJunctions.s1Head.xCoor     = static_cast<int>(m_imgResultsData(5,1));
      m_sJunctions.s1Head.yCoor     = static_cast<int>(m_imgResultsData(6,1));
      // s2GulToHum
      m_sJunctions.s2GulToHum.xCoor = static_cast<int>(m_imgResultsData(7,1));
      m_sJunctions.s2GulToHum.yCoor = static_cast<int>(m_imgResultsData(8,1));
      // s3HumToPec
      m_sJunctions.s3HumToPec.xCoor = static_cast<int>(m_imgResultsData(9,1));
      m_sJunctions.s3HumToPec.yCoor = static_cast<int>(m_imgResultsData(10,1));
      // s4PecToAbd
      m_sJunctions.s4PecToAbd.xCoor = static_cast<int>(m_imgResultsData(11,1));
      m_sJunctions.s4PecToAbd.yCoor = static_cast<int>(m_imgResultsData(12,1));
      // s5AbdToFem
      m_sJunctions.s5AbdToFem.xCoor = static_cast<int>(m_imgResultsData(13,1));
      m_sJunctions.s5AbdToFem.yCoor = static_cast<int>(m_imgResultsData(14,1));
      // s6FemToAna
      m_sJunctions.s6FemToAna.xCoor = static_cast<int>(m_imgResultsData(15,1));
      m_sJunctions.s6FemToAna.yCoor = static_cast<int>(m_imgResultsData(16,1));
      // s7Tail
      m_sJunctions.s7Tail.xCoor     = static_cast<int>(m_imgResultsData(17,1));
      m_sJunctions.s7Tail.yCoor     = static_cast<int>(m_imgResultsData(18,1));
    }
    //else
    //{
    //  TortoiseContourLocalizationWithGHT(); // variables listed below are modified by this method
    //  m_imgResultsData( 0,1) = 1;
    //  m_imgResultsData( 1,1) = m_sPlastronCentreOnRotatedImg.xCoor;
    //  m_imgResultsData( 2,1) = m_sPlastronCentreOnRotatedImg.yCoor;
    //  m_imgResultsData( 4,1) = m_nRotationAngleToVerticalInDegrees;
    //  // s1Head
    //  m_imgResultsData( 5,1) = m_sJunctions.s1Head.xCoor;
    //  m_imgResultsData( 6,1) = m_sJunctions.s1Head.yCoor;
    //  // s2GulToHum
    //  m_imgResultsData( 7,1) = m_sJunctions.s2GulToHum.xCoor;
    //  m_imgResultsData( 8,1) = m_sJunctions.s2GulToHum.yCoor;
    //  // s3HumToPec
    //  m_imgResultsData( 9,1) = m_sJunctions.s3HumToPec.xCoor;
    //  m_imgResultsData(10,1) = m_sJunctions.s3HumToPec.yCoor;
    //  // s4PecToAbd
    //  m_imgResultsData(11,1) = m_sJunctions.s4PecToAbd.xCoor;
    //  m_imgResultsData(12,1) = m_sJunctions.s4PecToAbd.yCoor;
    //  // s5AbdToFem
    //  m_imgResultsData(13,1) = m_sJunctions.s5AbdToFem.xCoor;
    //  m_imgResultsData(14,1) = m_sJunctions.s5AbdToFem.yCoor;
    //  // s6FemToAna
    //  m_imgResultsData(15,1) = m_sJunctions.s6FemToAna.xCoor;
    //  m_imgResultsData(16,1) = m_sJunctions.s6FemToAna.yCoor;
    //  // s7Tail
    //  m_imgResultsData(17,1) = m_sJunctions.s7Tail.xCoor;
    //  m_imgResultsData(18,1) = m_sJunctions.s7Tail.yCoor;

    //  cout << "Tortoise Contour localization has finished!" << endl;
    //  m_imgResultsData.save(m_strResultsDataFileName.c_str());
    //}
    ////cout << 1<<endl;
    //m_imgRotatedOriginal.assign(m_imgOriginal.get_rotate(static_cast<float>(m_nRotationAngleToVerticalInDegrees),1,1));
    ////cout << 2 <<endl;

    if(m_imgResultsData(0,2) == 1) // Abd Seam Junction Localization (PecToAbd and AbdToFem)
    {
      m_sJunctions.s4PecToAbd.yCoor        = m_imgResultsData(1,2);
      m_sJunctions.s5AbdToFem.yCoor        = m_imgResultsData(2,2);
      m_sPlastronCentreOnRotatedImg.yCoor = m_imgResultsData(3,2);
    }
    //else // variables listed below are modified by this method
    //{
    //  GetAbdJunctionsPreciseLocations(m_imgRotatedOriginal, m_sJunctions,m_sPlastronCentreOnRotatedImg);
    //  m_imgResultsData(0,2) = 1;
    //  m_imgResultsData(1,2) = m_sJunctions.s4PecToAbd.yCoor;
    //  m_imgResultsData(2,2) = m_sJunctions.s5AbdToFem.yCoor;
    //  m_imgResultsData(3,2) = m_sPlastronCentreOnRotatedImg.yCoor;
    //  m_imgResultsData.save(m_strResultsDataFileName.c_str());
    //  cout << "Tortoise Abdominal seam junctions precise localization has finished!" << endl;
    //}

    if(m_imgResultsData(0,3) == 1) // FemToAna Localization
    {
      m_sJunctions.s6FemToAna.yCoor = m_imgResultsData(1,3);
    }
    //else
    //{
    //  GetFemToAnaJunctionPreciseLocation(m_imgRotatedOriginal,m_sJunctions); // variables listed below are modified by this method
    //  m_imgResultsData(0,3) = 1;
    //  m_imgResultsData(1,3) = m_sJunctions.s6FemToAna.yCoor;
    //  cout << "Tortoise FemToAna junction precise localization has finished!" << endl;
    //  m_imgResultsData.save(m_strResultsDataFileName.c_str());
    //}

    if(m_imgResultsData(0,4) == 1) // HumToPec Localization
    {
      m_sJunctions.s3HumToPec.yCoor = m_imgResultsData(1,4);
    }
    //else
    //{
    //  GetHumToPecJunctionPreciseLocation(m_imgRotatedOriginal,m_sJunctions); // variables listed below are modified by this method
    //  m_imgResultsData(0,4) = 1;
    //  m_imgResultsData(1,4) = m_sJunctions.s3HumToPec.yCoor;
    //  cout << "Tortoise HumToPec junction precise localization has finished!" << endl;
    //  m_imgResultsData.save(m_strResultsDataFileName.c_str());
    //}

    if(m_imgResultsData(0,5) == 1) // HumToPec Localization
    {
      m_sJunctions.s2GulToHum.yCoor = m_imgResultsData(1,5);
    }
    //else
    //{
    //  GetGulToHumJunctionPreciseLocation(m_imgRotatedOriginal,m_sJunctions); // variables listed below are modified by this method
    //  m_imgResultsData(0,5) = 1;
    //  m_imgResultsData(1,5) = m_sJunctions.s2GulToHum.yCoor;
    //  cout << "Tortoise GulToHum junction precise localization has finished!" << endl;
    //  m_imgResultsData.save(m_strResultsDataFileName.c_str());
    //}
    if(m_imgResultsData(0,6) == 1) // Head and Tail Junction Localization
    {
      m_sJunctions.s1Head.yCoor   = m_imgResultsData(1,6);
      m_sCentralSeamHeadEnd.yCoor = m_imgResultsData(2,6);
      m_sCentralSeamHeadEnd.xCoor = m_imgResultsData(3,6);
      m_sJunctions.s7Tail.yCoor   = m_imgResultsData(4,6);
      m_sCentralSeamTailEnd.yCoor = m_imgResultsData(5,6);
      m_sCentralSeamTailEnd.xCoor = m_imgResultsData(6,6);
    }
    ////else
    //cout << "Hello" << endl;
    //{
    //  // variables listed below are modified by this method
    //  GetHeadJunctionLocalizationWithGTH(m_imgRotatedOriginal,m_sJunctions,m_sPlastronCentreOnRotatedImg);
    //  GetTailJunctionLocalizationWithGTH(m_imgRotatedOriginal,m_sJunctions,m_sPlastronCentreOnRotatedImg);
    //  m_imgResultsData(0,6) = 1;
    //  m_imgResultsData(1,6) = m_sJunctions.s1Head.yCoor;
    //  m_imgResultsData(2,6) = m_sCentralSeamHeadEnd.yCoor;
    //  m_imgResultsData(3,6) = m_sCentralSeamHeadEnd.xCoor;
    //  m_imgResultsData(4,6) = m_sJunctions.s7Tail.yCoor;
    //  m_imgResultsData(5,6) = m_sCentralSeamTailEnd.yCoor;
    //  m_imgResultsData(6,6) = m_sCentralSeamTailEnd.xCoor;

    //  cout << "Tortoise Head and Tail junction localization has finished!" << endl;
    //  m_imgResultsData.save(m_strResultsDataFileName.c_str());
    //}
    //CImg<int> imgShowJunctions(m_imgRotatedOriginal*2/3);
    //MarkAJunction(imgShowJunctions,m_sJunctions.s4PecToAbd,1,red);
    //MarkAJunction(imgShowJunctions,m_sJunctions.s5AbdToFem,1,yellow);
    //MarkAJunction(imgShowJunctions,m_sJunctions.s6FemToAna,1,green);
    //MarkAJunction(imgShowJunctions,m_sJunctions.s3HumToPec,1,purple);
    //MarkAJunction(imgShowJunctions,m_sJunctions.s2GulToHum,1,blue);
    //MarkAJunction(imgShowJunctions,m_sJunctions.s1Head,1,cyan);
    //MarkAJunction(imgShowJunctions,m_sJunctions.s7Tail,1,white);
    //MarkAJunction(imgShowJunctions,m_sCentralSeamHeadEnd,1,green);
    //MarkAJunction(imgShowJunctions,m_sCentralSeamTailEnd,1,green);
    //GetUniformlyResizedImg(imgShowJunctions,500).save(("possibleJunctions\\" + m_strTortoiseName + "HeadAndTailBoundary.bmp").c_str());

    if(m_imgResultsData(0,7) == 1 && m_imgResultsData(0,8) == 1) // All Left And Right Junctions Localization
    {
      // s1Head - JUNCTION
      m_sJunctions.s1Head.xCoor = m_imgResultsData( 1,7);
      m_sJunctions.s1Head.yCoor = m_imgResultsData( 2,7);
      // s2GulToHum
      m_sJunctions.s7Tail.xCoor = m_imgResultsData( 3,7);
      m_sJunctions.s7Tail.yCoor = m_imgResultsData( 4,7);

      // s1Head - LEFT
      m_sLeftJunctions.s1Head.xCoor = m_imgResultsData( 5,7);
      m_sLeftJunctions.s1Head.yCoor = m_imgResultsData( 6,7);
      // s2GulToHum
      m_sLeftJunctions.s2GulToHum.xCoor = m_imgResultsData( 7,7);
      m_sLeftJunctions.s2GulToHum.yCoor = m_imgResultsData( 8,7);
      // s3HumToPec
      m_sLeftJunctions.s3HumToPec.xCoor = m_imgResultsData( 9,7);
      m_sLeftJunctions.s3HumToPec.yCoor = m_imgResultsData(10,7);
      // s4PecToAbd
      m_sLeftJunctions.s4PecToAbd.xCoor = m_imgResultsData(11,7);
      m_sLeftJunctions.s4PecToAbd.yCoor = m_imgResultsData(12,7);
      // s5AbdToFem
      m_sLeftJunctions.s5AbdToFem.xCoor = m_imgResultsData(13,7);
      m_sLeftJunctions.s5AbdToFem.yCoor = m_imgResultsData(14,7);
      // s6FemToAna
      m_sLeftJunctions.s6FemToAna.xCoor = m_imgResultsData(15,7);
      m_sLeftJunctions.s6FemToAna.yCoor = m_imgResultsData(16,7);
      // s7Tail
      m_sLeftJunctions.s7Tail.xCoor = m_imgResultsData(17,7);
      m_sLeftJunctions.s7Tail.yCoor = m_imgResultsData(18,7);

      // s1Head - RIGHT
      m_sRightJunctions.s1Head.xCoor = m_imgResultsData(1,8);
      m_sRightJunctions.s1Head.yCoor = m_imgResultsData(2,8);
      // s2GulToHum
      m_sRightJunctions.s2GulToHum.xCoor = m_imgResultsData(3,8);
      m_sRightJunctions.s2GulToHum.yCoor = m_imgResultsData(4,8);
      // s3HumToPec
      m_sRightJunctions.s3HumToPec.xCoor = m_imgResultsData(5,8);
      m_sRightJunctions.s3HumToPec.yCoor = m_imgResultsData(6,8);
      // s4PecToAbd
      m_sRightJunctions.s4PecToAbd.xCoor = m_imgResultsData(7,8);
      m_sRightJunctions.s4PecToAbd.yCoor = m_imgResultsData(8,8);
      // s5AbdToFem
      m_sRightJunctions.s5AbdToFem.xCoor = m_imgResultsData(9,8);
      m_sRightJunctions.s5AbdToFem.yCoor = m_imgResultsData(10,8);
      // s6FemToAna
      m_sRightJunctions.s6FemToAna.xCoor = m_imgResultsData(11,8);
      m_sRightJunctions.s6FemToAna.yCoor = m_imgResultsData(12,8);
      // s7Tail
      m_sRightJunctions.s7Tail.xCoor = m_imgResultsData(13,8);
      m_sRightJunctions.s7Tail.yCoor = m_imgResultsData(14,8);
    }
      m_sJunctions.s1Head.xCoor = m_sRightJunctions.s1Head.xCoor;
      m_sJunctions.s1Head.yCoor = m_sRightJunctions.s1Head.yCoor;
      // s2GulToHum
      if(abs(m_sJunctions.s2GulToHum.yCoor-m_sRightJunctions.s2GulToHum.yCoor)
        <  abs(m_sJunctions.s2GulToHum.yCoor-m_sLeftJunctions.s2GulToHum.yCoor))  m_sJunctions.s2GulToHum.xCoor = m_sRightJunctions.s2GulToHum.xCoor;
      else m_sJunctions.s2GulToHum.xCoor = m_sLeftJunctions.s2GulToHum.xCoor;
      //m_sJunctions.s2GulToHum.yCoor = m_imgResultsData(4,8);
      // s3HumToPec
      if(abs(m_sJunctions.s3HumToPec.yCoor-m_sRightJunctions.s3HumToPec.yCoor)
        <  abs(m_sJunctions.s3HumToPec.yCoor-m_sLeftJunctions.s3HumToPec.yCoor))  m_sJunctions.s3HumToPec.xCoor = m_sRightJunctions.s3HumToPec.xCoor;
      else m_sJunctions.s3HumToPec.xCoor = m_sLeftJunctions.s3HumToPec.xCoor;
      //m_sJunctions.s3HumToPec.yCoor = m_imgResultsData(6,8);
      // s4PecToAbd
      if(abs(m_sJunctions.s4PecToAbd.yCoor-m_sRightJunctions.s4PecToAbd.yCoor)
        <  abs(m_sJunctions.s4PecToAbd.yCoor-m_sLeftJunctions.s4PecToAbd.yCoor))  m_sJunctions.s4PecToAbd.xCoor = m_sRightJunctions.s4PecToAbd.xCoor;
      else m_sJunctions.s4PecToAbd.xCoor = m_sLeftJunctions.s4PecToAbd.xCoor;
      //m_sJunctions.s4PecToAbd.yCoor = m_imgResultsData(8,8);
      // s5AbdToFem
      if(abs(m_sJunctions.s5AbdToFem.yCoor-m_sRightJunctions.s5AbdToFem.yCoor)
        <  abs(m_sJunctions.s5AbdToFem.yCoor-m_sLeftJunctions.s5AbdToFem.yCoor))  m_sJunctions.s5AbdToFem.xCoor = m_sRightJunctions.s5AbdToFem.xCoor;
      else m_sJunctions.s5AbdToFem.xCoor = m_sLeftJunctions.s5AbdToFem.xCoor;
      //m_sJunctions.s5AbdToFem.yCoor = m_imgResultsData(10,8);
      // s6FemToAna
      if(abs(m_sJunctions.s6FemToAna.yCoor-m_sRightJunctions.s6FemToAna.yCoor)
        <  abs(m_sJunctions.s6FemToAna.yCoor-m_sLeftJunctions.s6FemToAna.yCoor))  m_sJunctions.s6FemToAna.xCoor = m_sRightJunctions.s6FemToAna.xCoor;
      else m_sJunctions.s6FemToAna.xCoor = m_sLeftJunctions.s6FemToAna.xCoor;
      //m_sJunctions.s6FemToAna.yCoor = m_imgResultsData(12,8);
      // s7Tail
      m_sJunctions.s7Tail.xCoor = m_sRightJunctions.s7Tail.xCoor;
      m_sJunctions.s7Tail.yCoor = m_sRightJunctions.s7Tail.yCoor;
    //else
    //{
    //  //cout << "ERRRRRRRROOOOORRRRRRRRRR" << endl;
    //  //// variables listed below are modified by this method
    //  GetCentralPathXCoordinatesForInnerJunctions(m_imgRotatedOriginal, m_sJunctions, m_sPlastronCentreOnRotatedImg);
    //  GetGulToHumLeftAndRightJunctionCoor(m_imgRotatedOriginal, m_sJunctions, m_sPlastronCentreOnRotatedImg);
    //  GetHumToPecLeftAndRightJunctionCoor(m_imgRotatedOriginal, m_sJunctions, m_sPlastronCentreOnRotatedImg);
    //  GetFemToAnaLeftAndRightJunctionCoor(m_imgRotatedOriginal, m_sJunctions, m_sPlastronCentreOnRotatedImg);
    //  GetAbdLefAndRightJunctions(m_imgRotatedOriginal, m_sJunctions, m_sPlastronCentreOnRotatedImg);
    //  GetCentralPathXCoordinatesForHeadAndTailArea();
    //  GetTailJunctionPreciseLocation(m_imgRotatedOriginal,m_sJunctions,m_sPlastronCentreOnRotatedImg);
    //  GetHeadJunctionPreciseLocation(m_imgRotatedOriginal,m_sJunctions,m_sPlastronCentreOnRotatedImg);


      //m_imgResultsData(0,7) = 1;
      //m_imgResultsData(0,8) = 1;
      //// s1Head - JUNCTION
      //m_imgResultsData( 1,7) = m_sJunctions.s1Head.xCoor;
      //m_imgResultsData( 2,7) = m_sJunctions.s1Head.yCoor;
      //// s2GulToHum
      //m_imgResultsData( 3,7) = m_sJunctions.s7Tail.xCoor;
      //m_imgResultsData( 4,7) = m_sJunctions.s7Tail.yCoor;

      //// s1Head - LEFT
      //m_imgResultsData( 5,7) = m_sLeftJunctions.s1Head.xCoor;
      //m_imgResultsData( 6,7) = m_sLeftJunctions.s1Head.yCoor;
      //// s2GulToHum
      //m_imgResultsData( 7,7) = m_sLeftJunctions.s2GulToHum.xCoor;
      //m_imgResultsData( 8,7) = m_sLeftJunctions.s2GulToHum.yCoor;
      //// s3HumToPec
      //m_imgResultsData( 9,7) = m_sLeftJunctions.s3HumToPec.xCoor;
      //m_imgResultsData(10,7) = m_sLeftJunctions.s3HumToPec.yCoor;
      //// s4PecToAbd
      //m_imgResultsData(11,7) = m_sLeftJunctions.s4PecToAbd.xCoor;
      //m_imgResultsData(12,7) = m_sLeftJunctions.s4PecToAbd.yCoor;
      //// s5AbdToFem
      //m_imgResultsData(13,7) = m_sLeftJunctions.s5AbdToFem.xCoor;
      //m_imgResultsData(14,7) = m_sLeftJunctions.s5AbdToFem.yCoor;
      //// s6FemToAna
      //m_imgResultsData(15,7) = m_sLeftJunctions.s6FemToAna.xCoor;
      //m_imgResultsData(16,7) = m_sLeftJunctions.s6FemToAna.yCoor;
      //// s7Tail
      //m_imgResultsData(17,7) = m_sLeftJunctions.s7Tail.xCoor;
      //m_imgResultsData(18,7) = m_sLeftJunctions.s7Tail.yCoor;

      //// s1Head - RIGHT
      //m_imgResultsData(1,8) = m_sRightJunctions.s1Head.xCoor;
      //m_imgResultsData(2,8) = m_sRightJunctions.s1Head.yCoor;
      //// s2GulToHum
      //m_imgResultsData(3,8) = m_sRightJunctions.s2GulToHum.xCoor;
      //m_imgResultsData(4,8) = m_sRightJunctions.s2GulToHum.yCoor;
      //// s3HumToPec
      //m_imgResultsData(5,8) = m_sRightJunctions.s3HumToPec.xCoor;
      //m_imgResultsData(6,8) = m_sRightJunctions.s3HumToPec.yCoor;
      //// s4PecToAbd
      //m_imgResultsData(7,8) = m_sRightJunctions.s4PecToAbd.xCoor;
      //m_imgResultsData(8,8) = m_sRightJunctions.s4PecToAbd.yCoor;
      //// s5AbdToFem
      //m_imgResultsData(9,8) = m_sRightJunctions.s5AbdToFem.xCoor;
      //m_imgResultsData(10,8) = m_sRightJunctions.s5AbdToFem.yCoor;
      //// s6FemToAna
      //m_imgResultsData(11,8) = m_sRightJunctions.s6FemToAna.xCoor;
      //m_imgResultsData(12,8) = m_sRightJunctions.s6FemToAna.yCoor;
      //// s7Tail
      //m_imgResultsData(13,8) = m_sRightJunctions.s7Tail.xCoor;
      //m_imgResultsData(14,8) = m_sRightJunctions.s7Tail.yCoor;

      //cout << "Tortoise Left And Right junction localization has finished!" << endl;
      //m_imgResultsData.save(m_strResultsDataFileName.c_str());
    //}


    //CImg<int> imgShowLeftJunctions(m_imgRotatedOriginal*2/3);
    //MarkAJunction(imgShowLeftJunctions,m_sLeftJunctions.s4PecToAbd,0,red);
    //MarkAJunction(imgShowLeftJunctions,m_sLeftJunctions.s5AbdToFem,0,yellow);
    //MarkAJunction(imgShowLeftJunctions,m_sLeftJunctions.s6FemToAna,0,green);
    //MarkAJunction(imgShowLeftJunctions,m_sLeftJunctions.s3HumToPec,0,purple);
    //MarkAJunction(imgShowLeftJunctions,m_sLeftJunctions.s2GulToHum,0,blue);
    //MarkAJunction(imgShowLeftJunctions,m_sLeftJunctions.s1Head,0,cyan);
    //MarkAJunction(imgShowLeftJunctions,m_sLeftJunctions.s7Tail,0,white);
    //GetUniformlyResizedImg(imgShowLeftJunctions,500).display();//save(("possibleJunctions\\" + m_strTortoiseName + "InnerLeftJunctions.bmp").c_str());

    //CImg<int> imgShowRightJunctions(m_imgRotatedOriginal*2/3);
    //MarkAJunction(imgShowRightJunctions,m_sRightJunctions.s4PecToAbd,0,red);
    //MarkAJunction(imgShowRightJunctions,m_sRightJunctions.s5AbdToFem,0,yellow);
    //MarkAJunction(imgShowRightJunctions,m_sRightJunctions.s6FemToAna,0,green);
    //MarkAJunction(imgShowRightJunctions,m_sRightJunctions.s3HumToPec,0,purple);
    //MarkAJunction(imgShowRightJunctions,m_sRightJunctions.s2GulToHum,0,blue);
    //MarkAJunction(imgShowRightJunctions,m_sRightJunctions.s1Head,0,cyan);
    //MarkAJunction(imgShowRightJunctions,m_sRightJunctions.s7Tail,0,white);
    //GetUniformlyResizedImg(imgShowRightJunctions,500).display();//save(("possibleJunctions\\" + m_strTortoiseName + "InnerRightJunctions.bmp").c_str());
  }


  void Tortoise::GetHeadJunctionLocalizationWithGTH(const CImg<int> &imgRotatedOriginal,
                                                      const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg)
  {
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    cout << "Hello" << endl;
    CImg<int> imgUpperHalfTortoise(imgRotatedOriginal.get_crop(
      max(0,min(sPlastronCentreOnRotatedImg.xCoor-nAbdSeamLength,imgRotatedOriginal.width()-1)),max(0,min(sJunctions.s2GulToHum.yCoor-nAbdSeamLength*65/100,imgRotatedOriginal.height()-1)),
      max(0,min(sPlastronCentreOnRotatedImg.xCoor+nAbdSeamLength,imgRotatedOriginal.width()-1)),max(0,min(imgRotatedOriginal.height()-1,sJunctions.s2GulToHum.yCoor))));
    int nGulToHumYCoor = sJunctions.s2GulToHum.yCoor - max(0,min(sJunctions.s2GulToHum.yCoor-nAbdSeamLength*65/100,imgRotatedOriginal.height()-1));
    //cout << nGulToHumYCoor << endl;
    //imgUpperHalfTortoise.display();
    int nResizeMaxPixelSize = 200;
    float fSigmaBlurCoef = 0.5f;
    CImg<int> imgResizedUpperHalf(GetUniformlyResizedImg(imgUpperHalfTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
    nGulToHumYCoor = nGulToHumYCoor*imgResizedUpperHalf.height()/imgUpperHalfTortoise.height();
    //cout << nGulToHumYCoor << endl;
    //imgResizedUpperHalf.display();
    float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2-8);
    float fRightAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2+8);
    CImg<unsigned char> imgEdgesGHT(GetEdgeImg(imgResizedUpperHalf,50,100,0));
    imgEdgesGHT *= 255;
    string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
    string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
    imgResizedUpperHalf.save(strGHTImgFileName.c_str());
    imgEdgesGHT.save(strGHTEdgesFileName.c_str());
    GenHoughTrnf ght("ght\\1HeadTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
    // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
    ght.Init(90,200,5,5,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
    ght.runGHT(0);
    point2dCoor sAccumMax = {0,0};
    int nRotationAngleMax = 0;
    float fWidthRatioMax  = 0.0f;
    int nPlastronWidht = 0;
    CImg<int> imgPossibleJunctionPreciseLocations(1,1,1,1, nGulToHumYCoor);
    ght.bestCandidateFromGivenPoints(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht,imgPossibleJunctionPreciseLocations);//.display();//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
    cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
    //m_sJunctions.s1Head.yCoor = m_sJunctions.s2GulToHum.yCoor-nLenghtOfAbdSeam/8 - (imgResizedUpperHalf.height()-sAccumMax.yCoor)*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height();
    const int nTemplateHeight = 17; // template contour img
    const int nGularSeamLength = static_cast<int>(nTemplateHeight*fWidthRatioMax)*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height();
    int nUpperHeadJuncBoundary = max(0,min(sJunctions.s2GulToHum.yCoor-nGularSeamLength-nGularSeamLength/3,imgRotatedOriginal.height()-1));
    int nLowerHeadJuncBoundary = max(0,min(sJunctions.s2GulToHum.yCoor-nGularSeamLength+nGularSeamLength/3,imgRotatedOriginal.height()-1));
    m_sCentralSeamHeadEnd.yCoor = max(max(0,min(sJunctions.s2GulToHum.yCoor-nAbdSeamLength*65/100,imgRotatedOriginal.height()-1)),nUpperHeadJuncBoundary);
    m_sJunctions.s1Head.yCoor = min(nLowerHeadJuncBoundary,max(0,min(imgRotatedOriginal.height()-1,sJunctions.s2GulToHum.yCoor-nAbdSeamLength/10)));
    //m_sCentralSeamHeadEnd.xCoor = sPlastronCentreOnRotatedImg.xCoor + (sAccumMax.xCoor-imgResizedUpperHalf.width()/2)*imgUpperHalfTortoise.width()/imgResizedUpperHalf.width();
    //m_sJunctions.s1Head.xCoor = m_sCentralSeamHeadEnd.xCoor;
    m_sCentralSeamHeadEnd.xCoor = m_sJunctions.s1Head.xCoor;
    //cout << m_sCentralSeamHeadEnd.xCoor << " " << m_sCentralSeamHeadEnd.yCoor << endl;
    //imgRotatedOriginal.display();
  }

  void Tortoise::GetTailJunctionLocalizationWithGTH(const CImg<int> &imgRotatedOriginal,
                                                      const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg)
  {
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    CImg<int> imgUpperHalfTortoise(imgRotatedOriginal.get_crop(
      max(0,min(sPlastronCentreOnRotatedImg.xCoor-nAbdSeamLength,imgRotatedOriginal.width()-1)),max(0,min(sJunctions.s6FemToAna.yCoor,imgRotatedOriginal.height()-1)),
      max(0,min(sPlastronCentreOnRotatedImg.xCoor+nAbdSeamLength,imgRotatedOriginal.width()-1)),max(0,min(imgRotatedOriginal.height()-1,sJunctions.s6FemToAna.yCoor+nAbdSeamLength*65/100))));
    int nFemToAnaYCoor = sJunctions.s6FemToAna.yCoor - max(0,min(sJunctions.s6FemToAna.yCoor,imgRotatedOriginal.height()-1));
    //cout << nGulToHumYCoor << endl;
    //imgUpperHalfTortoise.display();
    int nResizeMaxPixelSize = 200;
    float fSigmaBlurCoef = 0.5f;
    CImg<int> imgResizedUpperHalf(GetUniformlyResizedImg(imgUpperHalfTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
    nFemToAnaYCoor = nFemToAnaYCoor*imgResizedUpperHalf.height()/imgUpperHalfTortoise.height();
    //cout << nGulToHumYCoor << endl;
    //imgResizedUpperHalf.display();
    float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2-8);
    float fRightAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2+8);
    CImg<unsigned char> imgEdgesGHT(GetEdgeImg(imgResizedUpperHalf,50,100,0));
    imgEdgesGHT *= 255;
    string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
    string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
    imgResizedUpperHalf.save(strGHTImgFileName.c_str());
    imgEdgesGHT.save(strGHTEdgesFileName.c_str());
    GenHoughTrnf ght("ght\\7TailTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
    // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
    ght.Init(90,200,5,5,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
    ght.runGHT(0);
    point2dCoor sAccumMax = {0,0};
    int nRotationAngleMax = 0;
    float fWidthRatioMax  = 0.0f;
    int nPlastronWidht = 0;
    CImg<int> imgPossibleJunctionPreciseLocations(1,1,1,1, nFemToAnaYCoor);
    ght.bestCandidateFromGivenPoints(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht,imgPossibleJunctionPreciseLocations);//.display();//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
    //cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
    //m_sJunctions.s1Head.yCoor = m_sJunctions.s2GulToHum.yCoor-nLenghtOfAbdSeam/8 - (imgResizedUpperHalf.height()-sAccumMax.yCoor)*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height();
    const int nTemplateHeight = 17; // template contour img
    const int nAnaSeamLength = static_cast<int>(nTemplateHeight*fWidthRatioMax)*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height();
    int nUpperTailJuncBoundary = max(0,min(sJunctions.s6FemToAna.yCoor+nAnaSeamLength-nAnaSeamLength/3,imgRotatedOriginal.height()-1));
    int nLowerTailJuncBoundary = max(0,min(sJunctions.s6FemToAna.yCoor+nAnaSeamLength+nAnaSeamLength/3,imgRotatedOriginal.height()-1));
    //upper
    m_sJunctions.s7Tail.yCoor   = max(max(0,min(sJunctions.s6FemToAna.yCoor,imgRotatedOriginal.height()-1)),nUpperTailJuncBoundary);
    m_sCentralSeamTailEnd.yCoor = min(nLowerTailJuncBoundary,max(0,min(imgRotatedOriginal.height()-1,sJunctions.s6FemToAna.yCoor+nAbdSeamLength*65/100)));
    //m_sCentralSeamTailEnd.xCoor = sPlastronCentreOnRotatedImg.xCoor + (sAccumMax.xCoor-imgResizedUpperHalf.width()/2)*imgUpperHalfTortoise.width()/imgResizedUpperHalf.width();
    //m_sJunctions.s7Tail.xCoor = m_sCentralSeamTailEnd.xCoor;
    m_sCentralSeamTailEnd.xCoor = m_sJunctions.s7Tail.xCoor;
    cout << m_sCentralSeamTailEnd.xCoor << " " << m_sCentralSeamTailEnd.yCoor << endl;

    //imgRotatedOriginal.display();
  }


  CImg<unsigned char> Tortoise::SkeletonizationMinDist(CImg<unsigned char> imgEdges)
{
	imgEdges = EliminateTwoPixelsEdges(EliminateOnePixelEdges(imgEdges));
	CImg<int> res(imgEdges);//,edgesWrapping(getEdgesWrapping(imgToDoItOn)), edgesWrappingInverse;
	int width = res.width();
	int height = res.height();
//sets value 0 for edges in resulting image
	res = (res-1).abs();
//sets value 0 for the area it works
  //edgesWrapping.display();
	//edgesWrappingInverse = (edgesWrapping-1).abs();
  //edgesWrappingInverse.display();
//adds edgesWrapping, so 0s are only in the working area where there is no edge
  //imgToDoItOn.display();
	//imgToDoItOn = imgToDoItOn + edgesWrappingInverse;
  //imgToDoItOn.display();

//
	for(int y = 0; y < height-1; y++)
	{
		for(int x = 0; x < width-1; x++)
		{
			if(imgEdges(x,y) == 0)
			{
				int minDist = 1;
				bool doIncreaseDistance = 1;
				while(doIncreaseDistance && x-minDist>=0 && y-minDist>=0 && x+minDist<width && y+minDist<height)
				{
					if(imgEdges.get_crop(x-minDist,y-minDist,x+minDist,y+minDist).sum() > 0 || minDist >= 3)
					{
						doIncreaseDistance = 0;
					}
					else minDist++;
				}
				res(x,y) = minDist;
			}
		}
	}
//	res = res.get_mul(edgesWrapping);
//eliminates pixs, first those close to edges (cycle 1->3)
	int maxLevel = 3;
	for(int cycle = 1; cycle <= maxLevel; cycle++)
	{
//in each cycle it eliminates also all newly-suitable pixs on lower levels (level 1...cycle)
		for(int level = 1; level <= cycle; level++)
		{
//eliminates until no elimination has occured in previous go-through of the image
			bool change = 1;
			while(change)
			{
				change = 0;
				for(int y = 1; y < height - 1; y++)
				{
					for(int x = 1; x < width-1; x++)
					{
						if(res(x,y) == level)
						{
// eliminate from rigth down corner
							if( res(x-1,y-1)> 0 && res(x,y-1)> 0 && res(x+1,y-1)>=0  &&
								res(x-1,y)  > 0 &&				    res(x+1,y)  ==0  &&
								res(x-1,y+1)>=0 && res(x,y+1)==0 && res(x+1,y+1)>=0     ) { res(x,y) = 0; change = 1; }
// eliminate from left down corner
							else if( res(x-1,y-1)>=0 && res(x,y-1)> 0 && res(x+1,y-1)> 0 &&
									 res(x-1,y)  ==0 &&				    res(x+1,y)  > 0 &&
									 res(x-1,y+1)>=0 && res(x,y+1)==0 && res(x+1,y+1)>=0  ) { res(x,y) = 0; change = 1; }
// eliminate from left up corner
							else if( res(x-1,y-1)>=0 && res(x,y-1)==0 && res(x+1,y-1)>=0 &&
									 res(x-1,y)  ==0 &&				    res(x+1,y)  > 0 &&
									 res(x-1,y+1)>=0 && res(x,y+1)> 0 && res(x+1,y+1)> 0   ) { res(x,y) = 0; change = 1; }
// eliminate from right up corner
							else if( res(x-1,y-1)>=0 && res(x,y-1)==0 && res(x+1,y-1)>=0 &&
									 res(x-1,y)  > 0 &&					res(x+1,y)  ==0 &&
									 res(x-1,y+1)> 0 && res(x,y+1)> 0 && res(x+1,y+1)>=0   ) { res(x,y) = 0; change = 1; }
// eliminate from right side
							else if( res(x-1,y-1)>0 &&  (res(x+1,y-1)==0||res(x+1,y+1)==0) &&
									 res(x-1,y)  >0 &&				 res(x+1,y)  ==0 &&
									 res(x-1,y+1)>0   ) { res(x,y) = 0; change = 1; }
// eliminate from left side
							else if( (res(x-1,y-1)==0||res(x-1,y+1)==0) && res(x+1,y-1)>0 &&
									  res(x-1,y)  ==0 &&					                res(x+1,y) >0 &&
								      res(x+1,y+1)>0  ) { res(x,y) = 0; change = 1; }
// eliminate from up side
							else if( (res(x-1,y-1)==0 || res(x+1,y-1)==0) && res(x,y-1)==0 &&
									 (res(x-1,y)  >0 ||					 res(x+1,y)  >0) &&
									  res(x-1,y+1)>0 && res(x,y+1)>0 && res(x+1,y+1)>0  ) { res(x,y) = 0; change = 1; }
// eliminate from down side
							else if( res(x-1,y-1)>0 && res(x,y-1)>0 && res(x+1,y-1)>0 &&
									(res(x-1,y)  >0 ||					 res(x+1,y)>0) &&
									(res(x-1,y+1)==0 || res(x+1,y+1)==0) && res(x,y+1)==0   ) { res(x,y) = 0; change = 1; }
							else
							{

							}
						}
					}
				}
			}
		}
	}
  //res.display();
	return res.threshold(1);
}

  CImg<unsigned char> Tortoise::GetCentralSeam(const CImg<unsigned char> &imgEdges)
{
  CImg<unsigned char> imgEdgesSkeletonized(SkeletonizationMinDist(imgEdges));
	CImg<double> distanceMask,paths(imgEdgesSkeletonized);
// paths contains skeleton of edge stripe image without large empty spaces
	int width = paths.width();
	int height = paths.height();
//decrLevel ... number of lines that we look backward(forward) to find more (suitable) edge points
	const int decrLevel = 5;
//distanceMask is filled by distances
	distanceMask.assign(5,decrLevel,1,1, 0);
	for(int x = 0; x < 5; x++) for(int y = 0; y < decrLevel; y++) distanceMask(x,y) = abs(x-2)*abs(x-2) + (y+1)*(y+1);
  //distanceMask(1,0) = 1;
  //distanceMask(3,0) = 1;

  //distanceMask.display();
//starting points
  const int nMaxValue = 100000;
//it is minimalization method, so unimportant pixels need to be a lot above other values (0 -> 99999, 1 -> 100000)
  paths += nMaxValue-1;
  //paths.display();
//starting point is set to 2
  for(int x = 0; x < width; ++x) paths(x,decrLevel) = 2;
  for(int x = 0; x < width; ++x) paths(x,height - decrLevel-1) = nMaxValue;
  //paths.display();
	distanceMask.rotate(180);
//starting point should be in y-lower part
	for(int y = decrLevel + 1; y < height; y++)
	{
		for(int x = 2; x < width - 2; x++)
		{
			if(paths(x,y) == nMaxValue)
			{
				CImg<double> cropArea;
				cropArea = paths.get_crop(x-2,y-decrLevel,x+2,y-1);
				if(cropArea.min() < nMaxValue-1)
				{
					cropArea = cropArea + distanceMask;
					paths(x,y) = cropArea.min();
				}
			}
		}
	}
  int nPathXCoor = static_cast<int>(paths.get_row(height - decrLevel-1).get_stats()(4));
  int nPathYCoor = height - decrLevel-1;
  cout << nPathXCoor << " " << nPathYCoor << endl;
  //paths.display();
  //distanceMask.display();
//beacause it is backward mode, the crop area setup stays the same for both scenearios
	while(nPathYCoor > decrLevel)
	{
		paths(nPathXCoor,nPathYCoor) = nMaxValue*2;
		CImg<double> cropArea;
		cropArea = paths.get_crop(nPathXCoor-2,nPathYCoor-decrLevel,nPathXCoor+2,nPathYCoor-1) + distanceMask;
		nPathXCoor = static_cast<int>(cropArea.get_stats()(0,4)) + nPathXCoor-2;
		nPathYCoor = static_cast<int>(cropArea.get_stats()(0,5)) + nPathYCoor-decrLevel;
	}
	return paths.threshold(nMaxValue*2);
}

  void Tortoise::GetCentralPathXCoordinatesForInnerJunctions(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCentreOnRotatedImg)
{
  const int nAbdSeamLenght = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
  CImg<int> imgStripe(imgRotatedOriginal.get_crop(sPlastronCentreOnRotatedImg.xCoor-nAbdSeamLenght/6,m_sJunctions.s1Head.yCoor,
                                                  sPlastronCentreOnRotatedImg.xCoor+nAbdSeamLenght/6,m_sJunctions.s7Tail.yCoor));
  CImg<unsigned char> imgEdges(GetEdgeImg(imgStripe,0,0,0));
  //imgEdges.display();
  CImg<unsigned char> imgCentralPath(GetCentralSeam(imgEdges));
  m_imgCentralPathXCoordinates.assign(1,imgRotatedOriginal.height(),1,1,0);
  for(int jj = 0; jj < imgCentralPath.height(); ++jj)
  {
    int y = jj + m_sJunctions.s1Head.yCoor;
    int x = 0;
    for(int ii = 0; ii < imgCentralPath.width(); ++ii)
    {
      if(imgCentralPath(ii,jj))
      {
        x = ii;
        ii = imgCentralPath.width();
      }
    }
    if(x>0) x += sPlastronCentreOnRotatedImg.xCoor-nAbdSeamLenght/6;
    m_imgCentralPathXCoordinates(y) = x;
  }
  //imgStripe.display();
  //(imgCentralPath,m_imgCentralPathXCoordinates,m_imgRotatedOriginal).display();
  //for(int y = 0; y < m_imgCentralPathXCoordinates.height(); ++y)
  //{
  //  int x = m_imgCentralPathXCoordinates(y);
  //  if(x>0)
  //  {
  //    m_imgRotatedOriginal(x,y,0)=255;
  //    m_imgRotatedOriginal(x,y,1)=0;
  //    m_imgRotatedOriginal(x,y,2)=0;
  //  }
  //}
  //m_imgRotatedOriginal.display();
}

  void Tortoise::GetGulToHumLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter)
  {
  //using edges on a wider central-seam stripe
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    const int nLeftBoundary = max(0,min(sPlastronCenter.xCoor-nAbdSeamLength/4,imgRotatedOriginal.width()-1));
    const int nRightBoundary = max(0,min(sPlastronCenter.xCoor+nAbdSeamLength/4,imgRotatedOriginal.width()-1));
    const int nLowerBoundary = max(0,min(sJunctions.s4PecToAbd.yCoor-nAbdSeamLength,imgRotatedOriginal.height()-1));
    const int nUpperBoundary = max(0,min(imgRotatedOriginal.height()-1,sJunctions.s3HumToPec.yCoor-nAbdSeamLength/4));
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(nLeftBoundary,nLowerBoundary,nRightBoundary,nUpperBoundary));
    //imgPlastron = GetUniformlyResizedImg(imgPlastron,500);
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,30,0));
    //imgEdges.display();
	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      //if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
    //imgEdges.display();
    CImg<int> localCentralPathCoordinates(m_imgCentralPathXCoordinates.get_crop(0,nLowerBoundary,0,nUpperBoundary) - nLeftBoundary);
    for(int jj = 0; jj < localCentralPathCoordinates.height(); ++jj)
    {
      int ii = localCentralPathCoordinates(jj);
      if(ii > 3 && ii < nImgEdgesWidth-4)
      {
        imgEdges(ii,jj) = 1;
        imgEdges(ii-1,jj) = 1;
        imgEdges(ii-2,jj) = 1;
        imgEdges(ii+1,jj) = 1;
        imgEdges(ii+2,jj) = 1;
      }
    }
    //imgEdges.display();
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 2;
      distanceMask(1,y) += 1;
      distanceMask(2,y) += 0;
      distanceMask(3,y) += 1;
      distanceMask(4,y) += 2;
      distanceMask(5,y) += 3;
      distanceMask(6,y) += 4;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
//  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
   for(int x = nImgEdgesWidth-1; x >= 0; --x)
    {
      for(int y = 0; y < nImgEdgesHeight; ++y)
	    {
        if(sidePathsFromRight(x,y)==nMaxPathValue) sidePathsFromRight(x,y) = 0;
        if(sidePathsFromLeft(x,y)==nMaxPathValue) sidePathsFromLeft(x,y) = 0;
        if(x==localCentralPathCoordinates(y))
        {
          sidePathsFromRight(x,y)*=2;
          sidePathsFromLeft(x,y)*=2;
        }
	    }
    }
   int nLocalYCoorOfGulToHumJunction = m_sJunctions.s2GulToHum.yCoor - nLowerBoundary;
   int nOKAreaAroundTheJunction = static_cast<int>(nAbdSeamLength*0.1/2);
   cout<< nOKAreaAroundTheJunction << " " << nLocalYCoorOfGulToHumJunction << endl;
   int nLeftMinYCoor = 0;
   int nLeftMin = nMaxPathValue;
   int nRightMinYCoor = 0;
   int nRightMin = nMaxPathValue;
   for(int jj = max(0,nLocalYCoorOfGulToHumJunction-nOKAreaAroundTheJunction); jj < min(nLocalYCoorOfGulToHumJunction+nOKAreaAroundTheJunction,localCentralPathCoordinates.height()); ++jj)
   {
     int ii = localCentralPathCoordinates(jj);
     if(ii>0 && sidePathsFromLeft(ii,jj)>0 && sidePathsFromLeft(ii,jj)<nLeftMin)
     {
       nLeftMinYCoor = jj;
       nLeftMin = sidePathsFromLeft(ii,jj);
     }
     if(ii>0 && sidePathsFromRight(ii,jj)>0 && sidePathsFromRight(ii,jj)<nRightMin)
     {
       nRightMinYCoor = jj;
       nRightMin = sidePathsFromRight(ii,jj);
     }
   }
//   cout << nRightMinYCoor << " " << nRightMin << " " << nLeftMinYCoor << " " << nLeftMin << endl;
   if(nLeftMin<nMaxPathValue)
   {
    m_sLeftJunctions.s2GulToHum.yCoor = nLeftMinYCoor + nLowerBoundary;
    m_sLeftJunctions.s2GulToHum.xCoor = localCentralPathCoordinates(nLeftMinYCoor) + nLeftBoundary;
   }
   else
   {
    m_sLeftJunctions.s2GulToHum.yCoor = m_sJunctions.s2GulToHum.yCoor;
    m_sLeftJunctions.s2GulToHum.xCoor = m_sJunctions.s2GulToHum.xCoor;
   }
   if(nRightMin<nMaxPathValue)
   {
     m_sRightJunctions.s2GulToHum.yCoor = nRightMinYCoor + nLowerBoundary;
     m_sRightJunctions.s2GulToHum.xCoor = localCentralPathCoordinates(nRightMinYCoor) + nLeftBoundary;
   }
   else
   {
     m_sRightJunctions.s2GulToHum.yCoor = m_sJunctions.s2GulToHum.yCoor;
     m_sRightJunctions.s2GulToHum.xCoor = m_sJunctions.s2GulToHum.xCoor;
   }
 //  (sidePathsFromLeft,sidePathsFromRight).display();
   //cout << m_sLeftJunctions.s2GulToHum.yCoor << " " << m_sLeftJunctions.s2GulToHum.xCoor << " " << m_sRightJunctions.s2GulToHum.yCoor << " " << m_sRightJunctions.s2GulToHum.xCoor << endl;
   //m_imgRotatedOriginal.display();

  }

  void Tortoise::GetHumToPecLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter)
  {
  //using edges on a wider central-seam stripe
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    const int nLeftBoundary = max(0,min(sPlastronCenter.xCoor-nAbdSeamLength*2/3,imgRotatedOriginal.width()-1));
    const int nRightBoundary = max(0,min(sPlastronCenter.xCoor+nAbdSeamLength*2/3,imgRotatedOriginal.width()-1));
    const int nLowerBoundary = max(0,min(sJunctions.s4PecToAbd.yCoor-nAbdSeamLength*2/3,imgRotatedOriginal.height()-1));
    const int nUpperBoundary = max(0,min(imgRotatedOriginal.height()-1,sJunctions.s4PecToAbd.yCoor));
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(nLeftBoundary,nLowerBoundary,nRightBoundary,nUpperBoundary));
    //imgPlastron = GetUniformlyResizedImg(imgPlastron,500);
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,30,0));
    //imgEdges.display();
	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      //if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
    //imgEdges.display();
    CImg<int> localCentralPathCoordinates(m_imgCentralPathXCoordinates.get_crop(0,nLowerBoundary,0,nUpperBoundary) - nLeftBoundary);
    for(int jj = 0; jj < localCentralPathCoordinates.height(); ++jj)
    {
      int ii = localCentralPathCoordinates(jj);
      if(ii > 3 && ii < nImgEdgesWidth-4)
      {
        imgEdges(ii,jj) = 1;
        imgEdges(ii-1,jj) = 1;
        imgEdges(ii-2,jj) = 1;
        imgEdges(ii+1,jj) = 1;
        imgEdges(ii+2,jj) = 1;
      }
    }
    const int A = nImgEdgesWidth-1;
    const int B = nImgEdgesHeight-1;
	  for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      for(int x = 0; x < nImgEdgesWidth; ++x)
		  {
        if( imgEdges(x,y) && (y>=B*x/(4*A)+B*7/8 || y>=B*(A-x)/(4*A)+B*7/8) )
			  {
  //prefilling - initial value (distance) of sidePathsFromRight is set to 2
				  imgEdges(x,y) = 0;
			  }
		  }
	  }
   // imgEdges.display();
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 3;
      distanceMask(1,y) += 2;
      distanceMask(2,y) += 0;
      distanceMask(3,y) += 0;
      distanceMask(4,y) += 1;
      distanceMask(5,y) += 2;
      distanceMask(6,y) += 3;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
//  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }

  for(int x = nImgEdgesWidth-1; x >= 0; --x)
    {
      for(int y = 0; y < nImgEdgesHeight; ++y)
	    {
        if(sidePathsFromRight(x,y)==nMaxPathValue) sidePathsFromRight(x,y) = 0;
        if(sidePathsFromLeft(x,y)==nMaxPathValue) sidePathsFromLeft(x,y) = 0;
        if(x==localCentralPathCoordinates(y))
        {
          sidePathsFromRight(x,y)*=2;
          sidePathsFromLeft(x,y)*=2;
        }
	    }
    }
  //(sidePathsFromRight,sidePathsFromLeft).display();
   int nLocalYCoorJunction = m_sJunctions.s3HumToPec.yCoor - nLowerBoundary;
   int nOKAreaAroundTheJunction = static_cast<int>(nAbdSeamLength*0.1/2);
   cout<< nOKAreaAroundTheJunction << " " << nLocalYCoorJunction << endl;
   int nLeftMinYCoor = 0;
   int nLeftMin = nMaxPathValue;
   int nRightMinYCoor = 0;
   int nRightMin = nMaxPathValue;
   for(int jj = max(0,nLocalYCoorJunction-nOKAreaAroundTheJunction); jj < min(nLocalYCoorJunction+nOKAreaAroundTheJunction,localCentralPathCoordinates.height()); ++jj)
   {
     int ii = localCentralPathCoordinates(jj);
     if(ii>0 && sidePathsFromLeft(ii,jj)>0 && sidePathsFromLeft(ii,jj)<nLeftMin)
     {
       nLeftMinYCoor = jj;
       nLeftMin = sidePathsFromLeft(ii,jj);
     }
     if(ii>0 && sidePathsFromRight(ii,jj)>0 && sidePathsFromRight(ii,jj)<nRightMin)
     {
       nRightMinYCoor = jj;
       nRightMin = sidePathsFromRight(ii,jj);
     }
   }
//   cout << nRightMinYCoor << " " << nRightMin << " " << nLeftMinYCoor << " " << nLeftMin << endl;
   if(nLeftMin<nMaxPathValue)
   {
    m_sLeftJunctions.s3HumToPec.yCoor = nLeftMinYCoor + nLowerBoundary;
    m_sLeftJunctions.s3HumToPec.xCoor = localCentralPathCoordinates(nLeftMinYCoor) + nLeftBoundary;
   }
   else
   {
    m_sLeftJunctions.s3HumToPec.yCoor = m_sJunctions.s3HumToPec.yCoor;
    m_sLeftJunctions.s3HumToPec.xCoor = m_sJunctions.s3HumToPec.xCoor;
   }
   if(nRightMin<nMaxPathValue)
   {
     m_sRightJunctions.s3HumToPec.yCoor = nRightMinYCoor + nLowerBoundary;
     m_sRightJunctions.s3HumToPec.xCoor = localCentralPathCoordinates(nRightMinYCoor) + nLeftBoundary;
   }
   else
   {
     m_sRightJunctions.s3HumToPec.yCoor = m_sJunctions.s3HumToPec.yCoor;
     m_sRightJunctions.s3HumToPec.xCoor = m_sJunctions.s3HumToPec.xCoor;
   }
   //cout << m_sLeftJunctions.s3HumToPec.yCoor << " " << m_sLeftJunctions.s3HumToPec.xCoor << " " << m_sRightJunctions.s3HumToPec.yCoor << " " << m_sRightJunctions.s3HumToPec.xCoor << endl;
   //m_imgRotatedOriginal.display();
  }

  void Tortoise::GetFemToAnaLeftAndRightJunctionCoor(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter)
  {
  //using edges on a wider central-seam stripe
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    const int nLeftBoundary = max(0,min(sJunctions.s1Head.xCoor-nAbdSeamLength*3/8,imgRotatedOriginal.width()-1));
    const int nRightBoundary = max(0,min(sJunctions.s1Head.xCoor+nAbdSeamLength*3/8,imgRotatedOriginal.width()-1));
    const int nLowerBoundary = max(0,min(sJunctions.s5AbdToFem.yCoor+nAbdSeamLength*4/30,imgRotatedOriginal.height()-1));
    const int nUpperBoundary = max(0,min(imgRotatedOriginal.height()-1,sJunctions.s5AbdToFem.yCoor+nAbdSeamLength*2/3));
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(nLeftBoundary,nLowerBoundary,nRightBoundary,nUpperBoundary));
    //imgPlastron = GetUniformlyResizedImg(imgPlastron,500);
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,30,0));
    //imgEdges.display();
	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      //if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
    //imgEdges.display();
    CImg<int> localCentralPathCoordinates(m_imgCentralPathXCoordinates.get_crop(0,nLowerBoundary,0,nUpperBoundary) - nLeftBoundary);
    for(int jj = 0; jj < localCentralPathCoordinates.height(); ++jj)
    {
      int ii = localCentralPathCoordinates(jj);
      if(ii > 3 && ii < nImgEdgesWidth-4)
      {
        imgEdges(ii,jj) = 1;
        imgEdges(ii-1,jj) = 1;
        imgEdges(ii-2,jj) = 1;
        imgEdges(ii+1,jj) = 1;
        imgEdges(ii+2,jj) = 1;
      }
    }
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 4;
      distanceMask(1,y) += 3;
      distanceMask(2,y) += 2;
      distanceMask(3,y) += 0;
      distanceMask(4,y) += 0;
      distanceMask(5,y) += 1;
      distanceMask(6,y) += 2;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
  //  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
   for(int x = nImgEdgesWidth-1; x >= 0; --x)
    {
      for(int y = 0; y < nImgEdgesHeight; ++y)
	    {
        if(sidePathsFromRight(x,y)==nMaxPathValue) sidePathsFromRight(x,y) = 0;
        if(sidePathsFromLeft(x,y)==nMaxPathValue) sidePathsFromLeft(x,y) = 0;
        if(x==localCentralPathCoordinates(y))
        {
          sidePathsFromRight(x,y)*=2;
          sidePathsFromLeft(x,y)*=2;
        }
	    }
    }
   //(sidePathsFromRight,sidePathsFromLeft).display();
   int nLocalYCoorJunction = m_sJunctions.s6FemToAna.yCoor - nLowerBoundary;
   int nOKAreaAroundTheJunction = static_cast<int>(nAbdSeamLength*0.1/2);
   cout<< nOKAreaAroundTheJunction << " " << nLocalYCoorJunction << endl;
   int nLeftMinYCoor = 0;
   int nLeftMin = nMaxPathValue;
   int nRightMinYCoor = 0;
   int nRightMin = nMaxPathValue;
   for(int jj = max(0,nLocalYCoorJunction-nOKAreaAroundTheJunction); jj < min(nLocalYCoorJunction+nOKAreaAroundTheJunction,localCentralPathCoordinates.height()); ++jj)
   {
     int ii = localCentralPathCoordinates(jj);
     if(ii>0 && sidePathsFromLeft(ii,jj)>0 && sidePathsFromLeft(ii,jj)<nLeftMin)
     {
       nLeftMinYCoor = jj;
       nLeftMin = sidePathsFromLeft(ii,jj);
     }
     if(ii>0 && sidePathsFromRight(ii,jj)>0 && sidePathsFromRight(ii,jj)<nRightMin)
     {
       nRightMinYCoor = jj;
       nRightMin = sidePathsFromRight(ii,jj);
     }
   }
//   cout << nRightMinYCoor << " " << nRightMin << " " << nLeftMinYCoor << " " << nLeftMin << endl;
   if(nLeftMin<nMaxPathValue)
   {
    m_sLeftJunctions.s6FemToAna.yCoor = nLeftMinYCoor + nLowerBoundary;
    m_sLeftJunctions.s6FemToAna.xCoor = localCentralPathCoordinates(nLeftMinYCoor) + nLeftBoundary;
   }
   else
   {
    m_sLeftJunctions.s6FemToAna.yCoor = m_sJunctions.s6FemToAna.yCoor;
    m_sLeftJunctions.s6FemToAna.xCoor = m_sJunctions.s6FemToAna.xCoor;
   }
   if(nRightMin<nMaxPathValue)
   {
     m_sRightJunctions.s6FemToAna.yCoor = nRightMinYCoor + nLowerBoundary;
     m_sRightJunctions.s6FemToAna.xCoor = localCentralPathCoordinates(nRightMinYCoor) + nLeftBoundary;
   }
   else
   {
     m_sRightJunctions.s6FemToAna.yCoor = m_sJunctions.s6FemToAna.yCoor;
     m_sRightJunctions.s6FemToAna.xCoor = m_sJunctions.s6FemToAna.xCoor;
   }
   //cout << m_sLeftJunctions.s6FemToAna.yCoor << " " << m_sLeftJunctions.s6FemToAna.xCoor << " " << m_sRightJunctions.s6FemToAna.yCoor << " " << m_sRightJunctions.s6FemToAna.xCoor << endl;
   //m_imgRotatedOriginal.display();
  }

  void Tortoise::GetAbdLefAndRightJunctions(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions,const  point2dCoor &sPlastronCentreOnRotatedImg)
  {
  //using edges on a wider central-seam stripe
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    const int nLeftBoundary = min(m_imgRotatedOriginal.width()-1,max(sPlastronCentreOnRotatedImg.xCoor-nAbdSeamLength,0));
    const int nRightBoundary = min(m_imgRotatedOriginal.width()-1,max(sPlastronCentreOnRotatedImg.xCoor+nAbdSeamLength,0));
    const int nLowerBoundary = sJunctions.s1Head.yCoor;
    const int nUpperBoundary = sJunctions.s7Tail.yCoor;
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(nLeftBoundary,nLowerBoundary,nRightBoundary,nUpperBoundary));
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,30,0));
    //imgEdges.display();
	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      //if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
    //imgEdges.display();
    CImg<int> localCentralPathCoordinates(m_imgCentralPathXCoordinates.get_crop(0,nLowerBoundary,0,nUpperBoundary) - nLeftBoundary);
    for(int jj = 0; jj < localCentralPathCoordinates.height(); ++jj)
    {
      int ii = localCentralPathCoordinates(jj);
      if(ii > 3 && ii < nImgEdgesWidth-4)
      {
        imgEdges(ii,jj) = 1;
        imgEdges(ii-1,jj) = 1;
        imgEdges(ii-2,jj) = 1;
        imgEdges(ii+1,jj) = 1;
        imgEdges(ii+2,jj) = 1;
      }
    }
    //imgEdges.display();
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 3;
      distanceMask(1,y) += 2;
      distanceMask(2,y) += 1;
      distanceMask(3,y) += 0;
      distanceMask(4,y) += 1;
      distanceMask(5,y) += 2;
      distanceMask(6,y) += 3;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
  //  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 5; y < nImgEdgesHeight-5; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
   for(int x = nImgEdgesWidth-1; x >= 0; --x)
    {
      for(int y = 0; y < nImgEdgesHeight; ++y)
	    {
        if(sidePathsFromRight(x,y)==nMaxPathValue) sidePathsFromRight(x,y) = 0;
        if(sidePathsFromLeft(x,y)==nMaxPathValue) sidePathsFromLeft(x,y) = 0;
        if(x==localCentralPathCoordinates(y))
        {
          sidePathsFromRight(x,y)*=2;
          sidePathsFromLeft(x,y)*=2;
        }
	    }
    }
   //(sidePathsFromRight,sidePathsFromLeft).display();
   // first PecToAbd junction left right localization
   int nLocalYCoorJunction = m_sJunctions.s4PecToAbd.yCoor - nLowerBoundary;
   int nOKAreaAroundTheJunction = static_cast<int>(nAbdSeamLength*0.1/2);
   cout<< nOKAreaAroundTheJunction << " " << nLocalYCoorJunction << endl;
   int nLeftMinYCoor = 0;
   int nLeftMin = nMaxPathValue;
   int nRightMinYCoor = 0;
   int nRightMin = nMaxPathValue;
   for(int jj = max(0,nLocalYCoorJunction-nOKAreaAroundTheJunction); jj < min(nLocalYCoorJunction+nOKAreaAroundTheJunction,localCentralPathCoordinates.height()); ++jj)
   {
     int ii = localCentralPathCoordinates(jj);
     if(ii>0 && sidePathsFromLeft(ii,jj)>0 && sidePathsFromLeft(ii,jj)<nLeftMin)
     {
       nLeftMinYCoor = jj;
       nLeftMin = sidePathsFromLeft(ii,jj);
     }
     if(ii>0 && sidePathsFromRight(ii,jj)>0 && sidePathsFromRight(ii,jj)<nRightMin)
     {
       nRightMinYCoor = jj;
       nRightMin = sidePathsFromRight(ii,jj);
     }
   }
//   cout << nRightMinYCoor << " " << nRightMin << " " << nLeftMinYCoor << " " << nLeftMin << endl;
   if(nLeftMin<nMaxPathValue)
   {
    m_sLeftJunctions.s4PecToAbd.yCoor = nLeftMinYCoor + nLowerBoundary;
    m_sLeftJunctions.s4PecToAbd.xCoor = localCentralPathCoordinates(nLeftMinYCoor) + nLeftBoundary;
   }
   else
   {
    m_sLeftJunctions.s4PecToAbd.yCoor = m_sJunctions.s4PecToAbd.yCoor;
    m_sLeftJunctions.s4PecToAbd.xCoor = m_sJunctions.s4PecToAbd.xCoor;
   }
   if(nRightMin<nMaxPathValue)
   {
     m_sRightJunctions.s4PecToAbd.yCoor = nRightMinYCoor + nLowerBoundary;
     m_sRightJunctions.s4PecToAbd.xCoor = localCentralPathCoordinates(nRightMinYCoor) + nLeftBoundary;
   }
   else
   {
     m_sRightJunctions.s4PecToAbd.yCoor = m_sJunctions.s4PecToAbd.yCoor;
     m_sRightJunctions.s4PecToAbd.xCoor = m_sJunctions.s4PecToAbd.xCoor;
   }

   //now localization of left right AbdToFem junction

   nLocalYCoorJunction = m_sJunctions.s5AbdToFem.yCoor - nLowerBoundary;
   nOKAreaAroundTheJunction = static_cast<int>(nAbdSeamLength*0.1/2);
   cout<< nOKAreaAroundTheJunction << " " << nLocalYCoorJunction << endl;
   nLeftMinYCoor = 0;
   nLeftMin = nMaxPathValue;
   nRightMinYCoor = 0;
   nRightMin = nMaxPathValue;
   for(int jj = max(0,nLocalYCoorJunction-nOKAreaAroundTheJunction); jj < min(nLocalYCoorJunction+nOKAreaAroundTheJunction,localCentralPathCoordinates.height()); ++jj)
   {
     int ii = localCentralPathCoordinates(jj);
     if(ii>0 && sidePathsFromLeft(ii,jj)>0 && sidePathsFromLeft(ii,jj)<nLeftMin)
     {
       nLeftMinYCoor = jj;
       nLeftMin = sidePathsFromLeft(ii,jj);
     }
     if(ii>0 && sidePathsFromRight(ii,jj)>0 && sidePathsFromRight(ii,jj)<nRightMin)
     {
       nRightMinYCoor = jj;
       nRightMin = sidePathsFromRight(ii,jj);
     }
   }
//   cout << nRightMinYCoor << " " << nRightMin << " " << nLeftMinYCoor << " " << nLeftMin << endl;
   if(nLeftMin<nMaxPathValue)
   {
    m_sLeftJunctions.s5AbdToFem.yCoor = nLeftMinYCoor + nLowerBoundary;
    m_sLeftJunctions.s5AbdToFem.xCoor = localCentralPathCoordinates(nLeftMinYCoor) + nLeftBoundary;
   }
   else
   {
    m_sLeftJunctions.s5AbdToFem.yCoor = m_sJunctions.s5AbdToFem.yCoor;
    m_sLeftJunctions.s5AbdToFem.xCoor = m_sJunctions.s5AbdToFem.xCoor;
   }
   if(nRightMin<nMaxPathValue)
   {
     m_sRightJunctions.s5AbdToFem.yCoor = nRightMinYCoor + nLowerBoundary;
     m_sRightJunctions.s5AbdToFem.xCoor = localCentralPathCoordinates(nRightMinYCoor) + nLeftBoundary;
   }
   else
   {
     m_sRightJunctions.s5AbdToFem.yCoor = m_sJunctions.s5AbdToFem.yCoor;
     m_sRightJunctions.s5AbdToFem.xCoor = m_sJunctions.s5AbdToFem.xCoor;
   }
   //cout << m_sLeftJunctions.s4PecToAbd.yCoor << " " << m_sLeftJunctions.s4PecToAbd.xCoor << " " << m_sRightJunctions.s4PecToAbd.yCoor << " " << m_sRightJunctions.s4PecToAbd.xCoor << endl;
   //cout << m_sLeftJunctions.s5AbdToFem.yCoor << " " << m_sLeftJunctions.s5AbdToFem.xCoor << " " << m_sRightJunctions.s5AbdToFem.yCoor << " " << m_sRightJunctions.s5AbdToFem.xCoor << endl;

   //m_imgRotatedOriginal.display();
  }

  void Tortoise::GetCentralPathXCoordinatesForHeadAndTailArea()
  {
    const int nLowerAbdToFemJunction = max(m_sLeftJunctions.s6FemToAna.yCoor,m_sRightJunctions.s6FemToAna.yCoor);
    const int nXCoorOfLowerAbdToFemJunction = m_imgCentralPathXCoordinates(nLowerAbdToFemJunction);
    for(int y = nLowerAbdToFemJunction+1; y < m_imgCentralPathXCoordinates.height(); ++y)
    {
      m_imgCentralPathXCoordinates(y) = nXCoorOfLowerAbdToFemJunction;
    }
    const int nUpperGulToHumJunction = min(m_sLeftJunctions.s2GulToHum.yCoor, m_sRightJunctions.s2GulToHum.yCoor);
    const int nXCoorOfUpperGulToHumJunction = m_imgCentralPathXCoordinates(nUpperGulToHumJunction);
    for(int y = nUpperGulToHumJunction-1; y>=0; --y)
    {
      m_imgCentralPathXCoordinates(y) = nXCoorOfUpperGulToHumJunction;
    }
    //for(int y = 0; y < m_imgCentralPathXCoordinates.height(); ++y)
    //{
    //  int x = m_imgCentralPathXCoordinates(y);
    //  m_imgRotatedOriginal(x,y,0) = 255;
    //  m_imgRotatedOriginal(x,y,1) = 0;
    //  m_imgRotatedOriginal(x,y,2) = 0;
    //}

    //(m_imgCentralPathXCoordinates,m_imgRotatedOriginal).display();
  }

  void Tortoise::GetTailJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter)
  {
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    const int nLeftBoundary =  max(0,min(sPlastronCenter.xCoor-nAbdSeamLength/4,imgRotatedOriginal.width()-1));
    const int nRightBoundary = max(0,min(sPlastronCenter.xCoor+nAbdSeamLength/4,imgRotatedOriginal.width()-1));
    const int nLowerBoundary = max(0,min(sJunctions.s6FemToAna.yCoor+nAbdSeamLength/10,imgRotatedOriginal.height()-1));
    const int nUpperBoundary = max(0,min(sJunctions.s6FemToAna.yCoor+nAbdSeamLength*65/100,imgRotatedOriginal.height()-1));
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(nLeftBoundary,nLowerBoundary,nRightBoundary,nUpperBoundary));
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,30,0));

	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);


    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
    //imgEdges.display();
    CImg<int> localCentralPathCoordinates(m_imgCentralPathXCoordinates.get_crop(0,nLowerBoundary,0,nUpperBoundary) - nLeftBoundary);
    for(int jj = 0; jj < localCentralPathCoordinates.height(); ++jj)
    {
      int ii = localCentralPathCoordinates(jj);
      if(ii > 3 && ii < nImgEdgesWidth-4)
      {
        imgEdges(ii,jj) = 1;
        imgEdges(ii-1,jj) = 1;
        imgEdges(ii-2,jj) = 1;
        imgEdges(ii+1,jj) = 1;
        imgEdges(ii+2,jj) = 1;
      }
    }
    //imgEdges.display();
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(7,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),3.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 9;
      distanceMask(1,y) += 4;
      distanceMask(2,y) += 2;
      distanceMask(3,y) += 1;
      distanceMask(4,y) += 0;
      distanceMask(5,y) += 2;
      distanceMask(6,y) += 4;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 3; y < nImgEdgesHeight-3; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
//  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 3; y < nImgEdgesHeight-3; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
   for(int x = nImgEdgesWidth-1; x >= 0; --x)
    {
      for(int y = 0; y < nImgEdgesHeight; ++y)
	    {
        if(sidePathsFromRight(x,y)==nMaxPathValue) sidePathsFromRight(x,y) = 0;
        if(sidePathsFromLeft(x,y)==nMaxPathValue) sidePathsFromLeft(x,y) = 0;
        if(x==localCentralPathCoordinates(y))
        {
          sidePathsFromRight(x,y)*=2;
          sidePathsFromLeft(x,y)*=2;
        }
	    }
    }
   //(sidePathsFromRight,sidePathsFromLeft).display();
   sidePathsFromRight += sidePathsFromLeft;
   //sidePathsFromRight.display();
   const int nLowerLocalYJunctionEstimate = m_sJunctions.s7Tail.yCoor - nLowerBoundary;
   const int nUpperLocalYJunctionEstimate = m_sCentralSeamTailEnd.yCoor - nLowerBoundary;

   int nMinYCoor = 0;
   int nMin = nMaxPathValue;
   for(int jj = max(0,nLowerLocalYJunctionEstimate); jj < min(nUpperLocalYJunctionEstimate,localCentralPathCoordinates.height()); ++jj)
   {
     int ii = localCentralPathCoordinates(jj);
     if(ii>0 && sidePathsFromRight(ii,jj)>0 && sidePathsFromRight(ii,jj)<nMin)
     {
       nMinYCoor = jj;
       nMin = sidePathsFromRight(ii,jj);
     }
   }

   if(nMin<nMaxPathValue)
   {
    m_sJunctions.s7Tail.yCoor = nMinYCoor + nLowerBoundary;
    m_sJunctions.s7Tail.xCoor = localCentralPathCoordinates(nMinYCoor) + nLeftBoundary;
    m_sLeftJunctions.s7Tail.yCoor = m_sJunctions.s7Tail.yCoor;
    m_sLeftJunctions.s7Tail.xCoor = m_sJunctions.s7Tail.xCoor;
    m_sRightJunctions.s7Tail.yCoor = m_sJunctions.s7Tail.yCoor;
    m_sRightJunctions.s7Tail.xCoor = m_sJunctions.s7Tail.xCoor;
   }
   else
   {
    m_sJunctions.s7Tail.yCoor = (m_sJunctions.s7Tail.yCoor + m_sCentralSeamTailEnd.yCoor)/2;
    m_sJunctions.s7Tail.xCoor = m_imgCentralPathXCoordinates(m_sJunctions.s7Tail.yCoor);
    m_sLeftJunctions.s7Tail.yCoor = m_sJunctions.s7Tail.yCoor;
    m_sLeftJunctions.s7Tail.xCoor = m_sJunctions.s7Tail.xCoor;
    m_sRightJunctions.s7Tail.yCoor = m_sJunctions.s7Tail.yCoor;
    m_sRightJunctions.s7Tail.xCoor = m_sJunctions.s7Tail.xCoor;
   }
   //cout << m_sLeftJunctions.s7Tail.xCoor << " " << m_sLeftJunctions.s7Tail.xCoor << endl;
   // cout << nMinYCoor << endl;
   // cout << nLowerLocalYJunctionEstimate << " " << nUpperLocalYJunctionEstimate << endl;
   // (sidePathsFromRight,imgRotatedOriginal).display();
  }

  void Tortoise::GetHeadJunctionPreciseLocation(const CImg<int> &imgRotatedOriginal, const plastronJunctions &sJunctions, const point2dCoor &sPlastronCenter)
  {
    const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
    const int nLeftBoundary =  max(0,min(sJunctions.s1Head.xCoor-nAbdSeamLength/4,imgRotatedOriginal.width()-1));
    const int nRightBoundary = max(0,min(sJunctions.s1Head.xCoor+nAbdSeamLength/4,imgRotatedOriginal.width()-1));
    const int nLowerBoundary = max(0,min(sJunctions.s2GulToHum.yCoor-nAbdSeamLength*65/100,imgRotatedOriginal.height()-1));
    const int nUpperBoundary = max(0,min(sJunctions.s2GulToHum.yCoor-nAbdSeamLength/10,imgRotatedOriginal.height()-1));
    CImg<int> imgPlastron(imgRotatedOriginal.get_crop(nLeftBoundary,nLowerBoundary,nRightBoundary,nUpperBoundary));
    CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,30,0));

	  const int nImgEdgesWidth = imgEdges.width();
	  const int nImgEdgesHeight = imgEdges.height();

  //prefilling - side paths to the value 10000, for future easier processing
    const int nMaxPathValue = 10000;
	  CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);


    for(int y = 0; y < nImgEdgesHeight; ++y)
	  {
      sidePathsFromLeft(4,y) = 2;
      sidePathsFromRight(nImgEdgesWidth-5,y) = 2;
      if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
      if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
	  }
    //imgEdges.display();
    CImg<int> localCentralPathCoordinates(m_imgCentralPathXCoordinates.get_crop(0,nLowerBoundary,0,nUpperBoundary) - nLeftBoundary);
    for(int jj = 0; jj < localCentralPathCoordinates.height(); ++jj)
    {
      int ii = localCentralPathCoordinates(jj);
      if(ii > 3 && ii < nImgEdgesWidth-4)
      {
        imgEdges(ii,jj) = 1;
        imgEdges(ii-1,jj) = 1;
        imgEdges(ii-2,jj) = 1;
        imgEdges(ii+1,jj) = 1;
        imgEdges(ii+2,jj) = 1;
      }
    }
    //imgEdges.display();
    //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
  //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
	  CImg<int> distanceMask(5,5,1,1, 0);
  //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
	  for(int x = 0; x < 5; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
	  for(int y = 0; y < 5; y++)
    {
      distanceMask(0,y) += 2;
      distanceMask(1,y) += 0;
      distanceMask(2,y) += 0;
      distanceMask(3,y) += 0;
      distanceMask(4,y) += 2;
    }

  //filling left side paths - rotating distanceMask
	  distanceMask.rotate(90);
    //distanceMask.display();
  //filling left side paths - going though all img pixs (rows only those where the central path is)
    for(int x = 5; x < nImgEdgesWidth; x++)
	  {
      for(int y = 2; y < nImgEdgesHeight-2; ++y)
		  {
  //filling left side paths - evaluating only edge pixs on the left of the central path
        if(imgEdges(x,y) /*&& x<=nImgEdgesWidth/2*/)
			  {
  //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
				  CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-2,x-5,y+2) + distanceMask);
  //filling left side paths - minimum distance is saved
				  int min = cropAreaDistances.min();
  //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
          if(min<nMaxPathValue)
				  {
  //filling left side paths - minimum distance is saved in the current pix
					  sidePathsFromLeft(x,y) = min;
  //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
//  (sidePathsFromLeft,sidePathsIndexes).display();
  //filling right side paths - rotating distanceMask
  distanceMask.mirror('x');
  //filling right side paths - the same as is in the left side
   for(int x = nImgEdgesWidth-6; x >= 0; --x)
	  {
      for(int y = 2; y < nImgEdgesHeight-2; ++y)
		  {
        if(imgEdges(x,y)/* && x >= nImgEdgesWidth/2*/)
			  {
				  CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-2,x+5,y+2) + distanceMask );
				  int min = cropAreaDistances.min();
          if(min<nMaxPathValue)
				  {
					  sidePathsFromRight(x,y) = min;
				  }
          else
          {
            sidePathsFromRight(x,y) = nMaxPathValue;
            sidePathsFromLeft(x,y) = nMaxPathValue;
            imgEdges(x,y) = 0;
          }
			  }
		  }
	  }
   for(int x = nImgEdgesWidth-1; x >= 0; --x)
    {
      for(int y = 0; y < nImgEdgesHeight; ++y)
	    {
        if(sidePathsFromRight(x,y)==nMaxPathValue) sidePathsFromRight(x,y) = 0;
        if(sidePathsFromLeft(x,y)==nMaxPathValue) sidePathsFromLeft(x,y) = 0;
        if(x==localCentralPathCoordinates(y))
        {
          sidePathsFromRight(x,y)*=2;
          sidePathsFromLeft(x,y)*=2;
        }
	    }
    }
   //(sidePathsFromRight,sidePathsFromLeft).display();
   sidePathsFromRight += sidePathsFromLeft;
   //sidePathsFromRight.display();
   const int nLowerLocalYJunctionEstimate = m_sCentralSeamHeadEnd.yCoor - nLowerBoundary;
   const int nUpperLocalYJunctionEstimate = m_sJunctions.s1Head.yCoor - nLowerBoundary;

   int nMinYCoor = 0;
   int nMin = nMaxPathValue;
   for(int jj = max(0,nLowerLocalYJunctionEstimate); jj < min(nUpperLocalYJunctionEstimate,localCentralPathCoordinates.height()); ++jj)
   {
     int ii = localCentralPathCoordinates(jj);
     if(ii>0 && sidePathsFromRight(ii,jj)>0 && sidePathsFromRight(ii,jj)<nMin)
     {
       nMinYCoor = jj;
       nMin = sidePathsFromRight(ii,jj);
     }
   }

   if(nMin<nMaxPathValue)
   {
    m_sJunctions.s1Head.yCoor = nMinYCoor + nLowerBoundary;
    m_sJunctions.s1Head.xCoor = localCentralPathCoordinates(nMinYCoor) + nLeftBoundary;
    m_sLeftJunctions.s1Head.yCoor = m_sJunctions.s1Head.yCoor;
    m_sLeftJunctions.s1Head.xCoor = m_sJunctions.s1Head.xCoor;
    m_sRightJunctions.s1Head.yCoor = m_sJunctions.s1Head.yCoor;
    m_sRightJunctions.s1Head.xCoor = m_sJunctions.s1Head.xCoor;
   }
   else
   {
    m_sJunctions.s1Head.yCoor = (m_sJunctions.s1Head.yCoor + m_sCentralSeamHeadEnd.yCoor)/2;
    m_sJunctions.s1Head.xCoor = m_imgCentralPathXCoordinates(m_sJunctions.s7Tail.yCoor);
    m_sLeftJunctions.s1Head.yCoor = m_sJunctions.s1Head.yCoor;
    m_sLeftJunctions.s1Head.xCoor = m_sJunctions.s1Head.xCoor;
    m_sRightJunctions.s1Head.yCoor = m_sJunctions.s1Head.yCoor;
    m_sRightJunctions.s1Head.xCoor = m_sJunctions.s1Head.xCoor;
   }
   //cout << m_sLeftJunctions.s1Head.xCoor << " " << m_sLeftJunctions.s1Head.xCoor << endl;
   //cout << nMinYCoor << endl;
   //cout << nLowerLocalYJunctionEstimate << " " << nUpperLocalYJunctionEstimate << endl;
   //(sidePathsFromRight,imgRotatedOriginal).display();
  }

  CImg<double> Tortoise::GetNormalizedJunctionDifferencesAndSeamLengths(const plastronJunctions &sLeftJunctions, const plastronJunctions &sRightJunctions)
  {
    CImg<double> MP(2,7,1,2,0); //measured points
    MP(0,0,0) = sLeftJunctions.s1Head.xCoor;     // left head junction x-coord
    MP(0,0,1) = sLeftJunctions.s1Head.yCoor;     // left head junction y-coord
    MP(0,1,0) = sLeftJunctions.s2GulToHum.xCoor; // left GulToHum junction x-coord
    MP(0,1,1) = sLeftJunctions.s2GulToHum.yCoor; // left GulToHum junction y-coord
    MP(0,2,0) = sLeftJunctions.s3HumToPec.xCoor; // left HumToPec junction x-coord
    MP(0,2,1) = sLeftJunctions.s3HumToPec.yCoor; // left HumToPec junction y-coord
    MP(0,3,0) = sLeftJunctions.s4PecToAbd.xCoor; // left PecToAbd junction x-coord
    MP(0,3,1) = sLeftJunctions.s4PecToAbd.yCoor; // left PecToAbd junction y-coord
    MP(0,4,0) = sLeftJunctions.s5AbdToFem.xCoor; // left AbdToFem junction x-coord
    MP(0,4,1) = sLeftJunctions.s5AbdToFem.yCoor; // left AbdToFem junction y-coord
    MP(0,5,0) = sLeftJunctions.s6FemToAna.xCoor; // left FemToAna junction x-coord
    MP(0,5,1) = sLeftJunctions.s6FemToAna.yCoor; // left FemToAna junction y-coord
    MP(0,6,0) = sLeftJunctions.s7Tail.xCoor;     // left Tail junction x-coord
    MP(0,6,1) = sLeftJunctions.s7Tail.yCoor;     // left Tail junction y-coord

    MP(1,0,0) = sRightJunctions.s1Head.xCoor;     // right head junction x-coord
    MP(1,0,1) = sRightJunctions.s1Head.yCoor;     // right head junction y-coord
    MP(1,1,0) = sRightJunctions.s2GulToHum.xCoor; // right GulToHum junction x-coord
    MP(1,1,1) = sRightJunctions.s2GulToHum.yCoor; // right GulToHum junction y-coord
    MP(1,2,0) = sRightJunctions.s3HumToPec.xCoor; // right HumToPec junction x-coord
    MP(1,2,1) = sRightJunctions.s3HumToPec.yCoor; // right HumToPec junction y-coord
    MP(1,3,0) = sRightJunctions.s4PecToAbd.xCoor; // right PecToAbd junction x-coord
    MP(1,3,1) = sRightJunctions.s4PecToAbd.yCoor; // right PecToAbd junction y-coord
    MP(1,4,0) = sRightJunctions.s5AbdToFem.xCoor; // right AbdToFem junction x-coord
    MP(1,4,1) = sRightJunctions.s5AbdToFem.yCoor; // right AbdToFem junction y-coord
    MP(1,5,0) = sRightJunctions.s6FemToAna.xCoor; // right FemToAna junction x-coord
    MP(1,5,1) = sRightJunctions.s6FemToAna.yCoor; // right FemToAna junction y-coord
    MP(1,6,0) = sRightJunctions.s7Tail.xCoor;     // right Tail junction x-coord
    MP(1,6,1) = sRightJunctions.s7Tail.yCoor;     // right Tail junction y-coord
    //MP.display();
    //calculate segments lenghts and junction differences at the same time
	  CImg<double> seamSegments(2,6,1,1, 0), leftRightJunctionDifferences(1,5,1,1, 0);
    double lastNearestX = MP(0,0,0);
    double lastNearestY = MP(0,0,1);
    for(int i = 0; i < 6; i++)
    {
      double distanceFromLastNearestPointOnLeft = sqrt( (double) (lastNearestX-MP(0,i+1,0))*(lastNearestX-MP(0,i+1,0))
                                                                + (lastNearestY-MP(0,i+1,1))*(lastNearestY-MP(0,i+1,1)) );
      double distanceFromLastNearestPointOnRight = sqrt( (double) (lastNearestX-MP(1,i+1,0))*(lastNearestX-MP(1,i+1,0))
                                                                + (lastNearestY-MP(1,i+1,1))*(lastNearestY-MP(1,i+1,1)) );
      double differenceBetweenLeftRightJunctionPoints = sqrt( (double) (MP(0,i+1,0)-MP(1,i+1,0))*(MP(0,i+1,0)-MP(1,i+1,0))
                                                                      + (MP(0,i+1,1)-MP(1,i+1,1))*(MP(0,i+1,1)-MP(1,i+1,1)) );
      if (distanceFromLastNearestPointOnLeft < distanceFromLastNearestPointOnRight)
      {
        seamSegments(0,i) += distanceFromLastNearestPointOnLeft;
        seamSegments(1,i) += distanceFromLastNearestPointOnLeft + differenceBetweenLeftRightJunctionPoints;
        lastNearestX = MP(1,i+1,0);
        lastNearestY = MP(1,i+1,1);
        if (i!=5)
        {
          seamSegments(0,i+1) += differenceBetweenLeftRightJunctionPoints;
          leftRightJunctionDifferences(0,i) = - differenceBetweenLeftRightJunctionPoints;
        }
      }
      else
      {
        seamSegments(1,i) += distanceFromLastNearestPointOnLeft;
        seamSegments(0,i) += distanceFromLastNearestPointOnLeft + differenceBetweenLeftRightJunctionPoints;
        lastNearestX = MP(0,i+1,0);
        lastNearestY = MP(0,i+1,1);
        if (i!=5)
        {
          seamSegments(1,i+1) += differenceBetweenLeftRightJunctionPoints;
          leftRightJunctionDifferences(0,i) = differenceBetweenLeftRightJunctionPoints;
        }
      }
    }
    //normalizing to the size of the plastron
	  double lengthOfCentralSeam = seamSegments.get_column(0).sum();
	  for(int i = 0; i < 6; i++)
    {
      seamSegments(0,i) /= lengthOfCentralSeam;
      seamSegments(1,i) /= lengthOfCentralSeam;
      if (i!=5) leftRightJunctionDifferences(0,i) /= lengthOfCentralSeam;
    }
    //(seamSegments,leftRightJunctionDifferences).display();
    return seamSegments.get_column(0).append(seamSegments.get_column(1),'y').append(leftRightJunctionDifferences,'y');
  }

  void Tortoise::Classification()
  {
    CImg<int> img("sidePaths.bmp"), res(GetUniformlyResizedImg(img,200)),edges(GetEdgeImg(res,0,0,0));
    //(res,edges).display();
    edges.save("sidePathsEdges.bmp");
    res.save("sidePathsRes.bmp");
    //string strAllTGsFeatures("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGoodNoYoungTo6cm\\FeaturesMeasuredByComputerNoLR.cimg");
    //string strSaveFeatureVectorFileName(m_strLoadDirectory + "resultFeaturesNoLR" + m_strTortoiseName + ".cimg");
    //CImg<double> resultFeatureZeroVec(strSaveFeatureVectorFileName.c_str());
    //CImg<double> allTGsFeatures;

    //if (FILE *file = fopen(strAllTGsFeatures.c_str(), "r")) // if results data exists (data with partial results)
    //{
    //  fclose(file);
    //  allTGsFeatures.assign(strAllTGsFeatures.c_str());
    //  allTGsFeatures.append(resultFeatureZeroVec,'x');
    //}
    //else // if not, create one
    //{
    //  allTGsFeatures.append(resultFeatureZeroVec,'x');
    //}
    //allTGsFeatures.save(strAllTGsFeatures.c_str());



    //CImg<double> resultFeatureZeroVec(1,17,1,1, 999);
    //string strSaveFeatureVectorFileName(m_strLoadDirectory + "resultFeaturesNoLR" + m_strTortoiseName + ".cimg");
    //resultFeatureZeroVec.save(strSaveFeatureVectorFileName.c_str());
    //string strVysledkyFilePath("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\PreciseJunctionLocalization\\\FinalResultsAllJunctionsSorted_DatabaseB_OnlyGoodNoYoungTGs\\GOOD\\");
    //string strResultImgFileName(strVysledkyFilePath+ m_strTortoiseName + "JHoughresultsImg.pnm");
    //CImg<int> imgResult(strResultImgFileName.c_str());
    ////imgResult.display();
    //Recognition();
    //GetNormalizedJunctionDifferencesAndSeamLengths(m_sJunctions,m_sJunctions).save(strSaveFeatureVectorFileName.c_str());//display();

    //GetUniformlyResizedImg(m_imgOriginal,500).save((m_strLoadDirectory+m_strTortoiseName+"Res.bmp").c_str());
    //getTrainingAndTestSetsLabelsAndTGnumbersForAGivenGroupOfFeatures(1,0,0,0);

    //CImg<double> FeaturesMeasuredByHand( (m_strLoadDirectory + "FeaturesMeasuredByHand.cimg").c_str());
    //CImg<int> TgNumbersInINT((m_strLoadDirectory + "numbersOfTGs.cimg").c_str());
    //(FeaturesMeasuredByHand,TgNumbersInINT).display();
  }
  //{
  // for(int grp = 0; grp <=2; ++grp)
  //  {
  //    string grps;
  //    if(grp==0) grps.assign("_group0");
  //    if(grp==1) grps.assign("_group1");
  //    if(grp==2) grps.assign("_group2");

  //    string classificationResultSuffix(grps+"_WithStdz_WithPCA_0.03Thres_11PCs");
  //    string classificationResultsDirectory("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\");
  //    string kNNclassificationResultsFileName (classificationResultsDirectory + "kNNResults" + classificationResultSuffix);
  //    CImgList<double> featureSet(getTrainingAndTestSetsLabelsAndTGnumbersForAGivenGroupOfFeatures(grp,1,1,0.03));
  //    CImg<double> trainingSet(featureSet(0)), trainingLabels(featureSet(1)),
  //                testSet(featureSet(2)), testLabels(featureSet(3)),trainingSetTGnumbers(featureSet(4)),testSetTGnumbers(featureSet(5));
  //    //featureSet.display();
  //    const int K = 100; // max k allowed
  //    int numberOfTrainingData = trainingSet.height();
  //    int numberOfDataDimensions = trainingSet.width();
  //    int numberOfTestData = testSet.height();
  //    cv::Mat trainData ( numberOfTrainingData, numberOfDataDimensions, CV_32FC1 );
  //    cv::Mat trainClasses ( numberOfTrainingData, 1, CV_32FC1 );
  //    for(int i = 0; i < numberOfTrainingData; ++i)
  //    {
  //      for(int j = 0; j < numberOfDataDimensions; ++j)
  //      {
  //        trainData.at<float>(i,j) = trainingSet(j,i);
  //      }
  //    }
  //    for(int i = 0; i < numberOfTrainingData; ++i)
  //    {
  //      trainClasses.at<float>(i,0) = trainingLabels(0,i);
  //    }
  //    CvKNearest knn( trainData, trainClasses, cv::Mat(), false, K ); // learning of the classifier, using training set and traning labels

    //  ofstream kNNClassificationResults (kNNclassificationResultsFileName + ".txt");
    //  ofstream kAndErrors (kNNclassificationResultsFileName + "_Errors" + ".txt");
    //  ofstream outliers (kNNclassificationResultsFileName + "_Outliers" + ".txt");
    //  kAndErrors << setw(20) << "Number k - neighbor" << setw(30) << "totalOptimisticError in %" << setw(30) << "totalPesimisticError in %" << setw(30) << "averageTotalError in %" << setw(30) << "ratting of k" << '\n';
    //  outliers << "Class 0 is the class of different tortoises and Class 1 is the class of same tortoises." << endl;
    //  for(int k = 1; k < K; k+=1)
    //  {
    //    if(k>30 && k < K-10) k+=10;
    //    kAndErrors << setw(20) << k;
    //    cout << k << endl;
    //    //int imgSize = 5000;
    //    //CImg<int> resultsOnImg(imgSize,3,1,2, 10);  // if feature space is 2d then results can be save onto an image
    //    //double sumOf1DfeaturesInClass0 = 0;
    //    //double sumOfSquaresOf1DfeaturesInClass0 = 0;
    //    //double countOfSamplesInClass0 = 0;
    //    //double sumOf1DfeaturesInClass1 = 0;
    //    //double sumOfSquaresOf1DfeaturesInClass1 = 0;
    //    //double countOfSamplesInClass1 = 0;
    //    //for (int i =  0; i < numberOfTrainingData; ++i) // OPTIMISTIC CLASSIFICATION - using training set
    //    //{
    //    //  cv::Mat sampleMat (1,1,CV_32FC1);
    //    //  sampleMat.at<float>(0,0) = trainingSet(0,i);
    //    //  float response = knn.find_nearest(sampleMat,k,0,0,0,0);
    //    //  if (response == 1) //TGs are same
    //    //  {
    //    //    sumOf1DfeaturesInClass1 +=  trainingSet(0,i);
    //    //    sumOfSquaresOf1DfeaturesInClass1 +=  trainingSet(0,i)*trainingSet(0,i);
    //    //    ++countOfSamplesInClass1;
    //    //  }
    //    //  else if (response == 0) //TGs are not same
    //    //  {
    //    //    sumOf1DfeaturesInClass0 +=  trainingSet(0,i);
    //    //    sumOfSquaresOf1DfeaturesInClass0 +=  trainingSet(0,i)*trainingSet(0,i);
    //    //    ++countOfSamplesInClass0;
    //    //  }
    //    //  else
    //    //  {
    //    //    cout << "Error in kNN classification" << endl;
    //    //  }
    //    ////}
    //    //double meanOfClass0 = sumOf1DfeaturesInClass0/countOfSamplesInClass0;
    //    //double stdOfClass0 = sqrt( ( sumOfSquaresOf1DfeaturesInClass0 - (sumOf1DfeaturesInClass0)*(sumOf1DfeaturesInClass0)/countOfSamplesInClass0 ) / (countOfSamplesInClass0-1) );
    //    //double meanOfClass1 = sumOf1DfeaturesInClass1/countOfSamplesInClass1;
    //    //double stdOfClass1 = sqrt( ( sumOfSquaresOf1DfeaturesInClass1 - (sumOf1DfeaturesInClass1)*(sumOf1DfeaturesInClass1)/countOfSamplesInClass1 ) / (countOfSamplesInClass1-1) );
    //    //cout << meanOfClass0 << '\t' << stdOfClass0 << '\t' << meanOfClass1 << '\t' << stdOfClass1 << '\t' << endl;
    //    //{
    //    //  for(int x = meanOfClass0*imgSize/trainingSet.max()+1; x > (meanOfClass0 - stdOfClass0)*imgSize/trainingSet.max()+1; --x)
    //    //  {
    //    //    resultsOnImg(x,0,0) =  5;
    //    //    resultsOnImg(x,0,1) =  5;
    //    //  }
    //    //  for(int x = meanOfClass1*imgSize/trainingSet.max()+1; x < (meanOfClass1 + stdOfClass1)*imgSize/trainingSet.max()+1; ++x)
    //    //  {
    //    //    resultsOnImg(x,0,0) = 15;
    //    //    resultsOnImg(x,0,1) = 15;
    //    //  }
    //    //}
    //    int numberOfCorrectlyClassifiedAsSameTortoises = 0;
    //    int numberOfCorrectlyClassifiedAsDifferentTortoises = 0;
    //    int numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises = 0;
    //    int numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises = 0;
    //    outliers << '\n' << k << " - nearest neigbours" << '\n';
    //    //outliers << setw(20) << "meanClass1" << setw(20) << "stdClass1" << setw(20) << "meanClass0" << setw(20) << "stdClass0"<< '\n';
    //    //outliers << setw(20) << meanOfClass1 << setw(20) << stdOfClass1 << setw(20) << meanOfClass0 << setw(20) << stdOfClass0 << '\n';
    //    outliers << "TrainingSet Classification" << '\n';
    //    outliers << setw(20) <<"TgNumber1" << setw(20) << "TgNumber2" << setw(20) << "sumOfFeatures" << setw(20) << "GroundTruthClass" << setw(20) << "ClassifiedClass" << '\n';
    //
    //    for (int i =  0; i < numberOfTrainingData; ++i) // OPTIMISTIC CLASSIFICATION - using training set
    //    {
    //      cv::Mat sampleMat (1,numberOfDataDimensions,CV_32FC1);
    //      for(int j = 0; j < numberOfDataDimensions; ++j)
    //      {
    //        sampleMat.at<float>(0,j) = trainingSet(j,i);
    //      }
    //      //int x = trainingSet(0,i)*imgSize/trainingSet.max()+1;
    //      //int y = 1;
    //      //int y = trainingSet(1,i)*1000/trainingSet.max()+1;
    //      float response = knn.find_nearest(sampleMat,k,0,0,0,0);
    //      outliers << setw(20) << trainingSetTGnumbers(0,i) << setw(20) << trainingSetTGnumbers(1,i) << setw(20) << trainingSet.get_row(i).sum() << setw(20) << trainingLabels(0,i) << setw(20) << response << '\n';

    //      if (response == 1 && trainingLabels(0,i)==1)
    //      {
    //        ++numberOfCorrectlyClassifiedAsSameTortoises;
    //        //if(x < imgSize && y < imgSize)
    //        //{
    //        //  resultsOnImg(x,y,0) =  15;
    //        //  resultsOnImg(x,y,1) =  15;
    //        //}
    //      }
    //      else if (response == 0 && trainingLabels(0,i)==0)
    //      {
    //        ++numberOfCorrectlyClassifiedAsDifferentTortoises;
    //        // if(x < imgSize && y < imgSize)
    //        // {
    //        //   resultsOnImg(x,y,0) =  5;
    //        //   resultsOnImg(x,y,1) =  5;
    //        //}
    //      }
    //      else if (response == 1 && trainingLabels(0,i)==0) // false positive
    //      {
    //        ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises;
    //        //if(x < imgSize && y < imgSize)
    //        //{
    //        //  resultsOnImg(x,y,0) =  5;
    //        //  resultsOnImg(x,y,1) =  15;
    //        //}
    //      }
    //      else if (response ==0 && trainingLabels(0,i)==1) //false negative
    //      {
    //        ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises;
    //        //outliers << setw(20) << trainingSetTGnumbers(0,i) << setw(20) << trainingSetTGnumbers(1,i) << setw(20) << " " << setw(20) << (trainingSet(0,i)-meanOfClass0)/stdOfClass0 << '\n';
    //        //if(x < imgSize && y < imgSize)
    //        //{
    //        //  resultsOnImg(x,y,0) =  15;
    //        //  resultsOnImg(x,y,1) =  5;
    //        //}
    //      }
    //      else
    //      {
    //        cout << response << " " << trainingLabels(0,i) << " neco se pokazilo" << endl;
    //      }
    //    }
    //    //if(k==51) resultsOnImg.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\classificationImgwithoutPow2.pnm");
    //    kNNClassificationResults << "  " <<  k << " - NN" << '\n';
    //    kNNClassificationResults << "OPTIMISTIC CLASSIFICATION ERROR" << '\n';
    //    kNNClassificationResults << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalOptimisticSuccess in %" << '\n';
    //    kNNClassificationResults << left << setw(30) << numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTrainingData << setw(40) << numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTrainingData;
    //    kNNClassificationResults << left << setw(40) << (numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0 << '\n';
    //    kNNClassificationResults << left << setw(30) << "incorrectlyClassifiedAsSame" << setw(40) << "incorrectlyClassifiedAsDifferent" << setw(40) << "totalOptimisticError in %" << '\n';
    //    kNNClassificationResults << left << setw(30) << numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData << setw(40) << numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData;
    //    kNNClassificationResults << left << setw(40) << (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0 << '\n';
    //    double optimisticTotalError = (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0;
    //    kAndErrors << setw(30) << optimisticTotalError;
    //    kNNClassificationResults << '\n';
    //
    //    numberOfCorrectlyClassifiedAsSameTortoises = 0;
    //    numberOfCorrectlyClassifiedAsDifferentTortoises = 0;
    //    numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises = 0;
    //    numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises = 0;
    //    outliers << '\n' << "TestSet Classification" << '\n';
    //    outliers << setw(20) << "TgNumber1" << setw(20) << "TgNumber2" << setw(20) << "sumOfFeatures" << setw(20) << "GroundTruthClass" << setw(20) << "ClassifiedClass" << '\n';
    //    for (int i =  0; i < numberOfTestData; ++i) //PESSIMISTIC CLASSIFICATION - using test set
    //    {
    //      cv::Mat sampleMat (1,numberOfDataDimensions,CV_32FC1);
    //      for(int j = 0; j < numberOfDataDimensions; ++j)
    //      {
    //        sampleMat.at<float>(0,j) = testSet(j,i);
    //      }
    //      //int x = testSet(0,i)*imgSize/testSet.max()+1;
    //      //int y = 2;
    //      float response = knn.find_nearest(sampleMat,k,0,0,0,0);
    //      outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << testSet.get_row(i).sum() << setw(20) << testLabels(0,i) << setw(20) << response << '\n';
    //      if (response == 1 && testLabels(0,i)==1)
    //      {
    //        ++numberOfCorrectlyClassifiedAsSameTortoises;
    //        //if(x < imgSize && y < imgSize)
    //        //{
    //        //  resultsOnImg(x,y,0) =  15;
    //        //  resultsOnImg(x,y,1) =  15;
    //        //}
    //      }
    //      else if (response ==0 && testLabels(0,i)==0)
    //      {
    //        ++numberOfCorrectlyClassifiedAsDifferentTortoises;
    //        //if(x < imgSize && y < imgSize)
    //        //{
    //        //  resultsOnImg(x,y,0) =  5;
    //        //  resultsOnImg(x,y,1) =  5;
    //        //}
    //      }
    //      else if (response == 1 && testLabels(0,i)==0)
    //      {
    //        ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises;
    //        //outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << 1 << setw(20) << " " << '\n';
    //        //if(x < imgSize && y < imgSize)
    //        //{
    //        //  resultsOnImg(x,y,0) =  5;
    //        //  resultsOnImg(x,y,1) =  15;
    //        //}
    //      }
    //      else if (response == 0 && testLabels(0,i)==1)
    //      {
    //        ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises;
    //        //outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << " " << setw(20) << 1 << '\n';
    //        //if(x < imgSize && y < imgSize)
    //        //{
    //        //  resultsOnImg(x,y,0) =  15;
    //        //  resultsOnImg(x,y,1) =  5;
    //        //}
    //      }
    //      else
    //      {
    //        cout << response << " " << testLabels(0,i) << " neco se pokazilo" << endl;
    //      }
    //    }
    //    //if(k%10==1) resultsOnImg.display();

    //    kNNClassificationResults << "PESSIMISTIC CLASSIFICATION ERROR" << '\n';
    //    kNNClassificationResults << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalPesimisticSuccess in %" << '\n';
    //    kNNClassificationResults << left << setw(30) << numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTestData << setw(40) << numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTestData;
    //    kNNClassificationResults << left << setw(40) << (numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTestData + numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0 << '\n';
    //    kNNClassificationResults << left << setw(30) << "incorrectlyClassifiedAsSame" << setw(40) << "incorrectlyClassifiedAsDifferent" << setw(40) << "totalPesimisticError in %" << '\n';
    //    kNNClassificationResults << left << setw(30) << numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData << setw(40) << numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData;
    //    kNNClassificationResults << left << setw(40) << (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0 << '\n';
    //    double pessimisticTotalError = (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0;
    //    kAndErrors << setw(30) << pessimisticTotalError <<  setw(30) << (optimisticTotalError + pessimisticTotalError)/2 <<  setw(30) << (optimisticTotalError + pessimisticTotalError)*pow(pessimisticTotalError,2.0)/2 << endl;
    //    kNNClassificationResults << '\n' << endl;
    //    outliers << endl;
    //  }
    //}
  //}
