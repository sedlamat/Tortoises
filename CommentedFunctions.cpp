
    //CImg<int> contour("ght\\contour.bmp");
    ////contour.display();
    ////contour.channel(0).threshold(1);
    ////for(int i = 1; i < contour.width()-1; ++i)
    ////{
    ////  for(int j=1; j < contour.height(); ++j )
    ////  {
    ////    if(contour.get_crop(i-1,j-1,i+1,j+1).sum()>2)
    ////    {
    ////      contour(i,j) = 0;
    ////    }
    ////  }
    ////}
    //contour.display();
    //int x = 0;
    //int y = 0;
    //cin >> x;
    //cin >> y;
    ////contour*=255;
    //contour(x,y,0) = 127;
    //contour.channel(0).display();
    //contour.save("ght\\contour.bmp"); 



/*  void drawVerticalStripe();

	CImg<int> getMediumWideStripeFromWideStripe(CImg<int> anyWideStripeImg);
	CImg<int> getNarrowStripeFromWideStripe(CImg<int> anyWideStripeImg);  

  void centerLocalization();

  void junctionLocalization();
  void localizationUsingGHT();

  CImg<int> Tortoise::kMeans( CImg<int> img, int K, int nIterations);

  CImg<int> Tortoise::getBestPathAcrossImgEdges(const CImg<int> &edgeImg);

  CImg<int> getCoordinatesOfCentralPath(CImg<int> imgWithCentralPathAsOnes);
  void drawBestPath();
  void drawAbdominalJunctions();
	CImg<int> skeletonization(CImg<int> imgToDoItOn);
	CImg<int> skeletonizationClassic(CImg<int> imgToDoItOn);
	CImg<int> GetEdgeImgOfSeams(CImg<int> img);

  CImg<unsigned char> Tortoise::skeletonizationMinDist(CImg<unsigned char> imgToDoItOn);
  CImg<unsigned char> Tortoise::getSkeletonizationMinDistWideStripe(CImg<unsigned char> imgToDoItOn);
CImg<unsigned char> Tortoise::getHorizontalSkeletonizationMinDist(CImg<unsigned char> imgToDoItOn);
	CImg<unsigned char> getEdgesWrapping(CImg<unsigned char> imgToDoItOn);
	CImg<int> improveSidePaths(CImg<int> imgToDoItOn);
  bool getYCoordinatesOfAbdominalJunctions(CImg<unsigned char> edgeImgToDoItOn);
	CImg<int> getCoordinatesOfBestCentralPathOnNarrowStripe();
	CImg<int> getCoordinatesOfBestCentralPathOnMediumWideStripe();
	CImg<int> getCoordinatesOfBestCentralPathOnWideStripe();
	CImg<char> getWideStripeWithDirectionOfEdges();
	CImg<char> getWideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges();
	CImg<double> getCentralAndSideSeamPaths();
	CImg<unsigned char> getPossibleEdgesForTheCentralSeamPath(CImg<unsigned char> edgeImgToDoItOn);
	CImg<unsigned char> getBestCentralPath(CImg<unsigned char> edgeImgToDoItOn);
	CImg<int> getJoiningPointsOfsidePathsToTheCentralPath(const double stripeWidthKoef);
	CImg<int> getYCoordinatesOfPectoralFemoralJunctions(CImg<unsigned char> edgeImgToDoItOn);
	CImg<int> correctionOfCentralPath(CImg<unsigned char> edgeImgToDoItOn);
	void getProlongedNarrowStripeBestPathCoordinates();
	CImg<int> drawJunctions();
	CImg<double> getAbdoPectoFemoSidePaths();
	CImg<double> getGroundTruthSegmentsLengths();
	void classifyTortoise();
	CImg<int> visualizeCentralAndSidePaths();*/

 // void HeadJunctionLocalization();
 // void GulToHumJunctionLocalization();
 // void HumToPecJunctionLocalization();
 // void PecToAbdJunctionLocalization();
 // void AbdToFemJunctionLocalizationGHT();
 // void AbdToFemJunctionLocalizationAStar();
 // void FemToAnaJunctionLocalization();
 // void TailJunctionLocalization();
 // CImg<int> DrawJunctions();
 // CImg<int> drawVerticalStripeAndJunctions();
 // void  GetJunctionCoordinates(int nJunctionNumber, int &nXCoor, int &nYCoor);
 // void AbdSeamJunctionPreciseLocalization();
 // void PecToAbdJunctionLocalizationGHTInJunctionPoints(CImg<int> imgPossibleJunctionPreciseLocations);
 // void AbdToFemJunctionLocalizationGHTInJunctionPoints(CImg<int> imgPossibleJunctionPreciseLocations);





 // CImg<int> GetBestSidePathAcrossImgEdges(CImg<int> &imgEdges, bool bMaskIsAsymetric);


 // 
 // CImg<int>     _wideStripe;
 // CImg<int>	    _narrowStripeEdges, _mediumWideStripeEdges, _wideStripeEdges;



 // bool          _isNewResultsImg;






 // //these are scale Variant, they are in original size of the rotated img
	//int			  	   m_nWidthOfRotatedOriginalImg, m_nHeightOfRotatedOriginalImg, _stripeHeight, _wideStripeWidth; 
	//int			  	  _leftXcoordinateOfStripeInTheRotatedImg, _rightXcoordinateOfStripeInTheRotatedImg;



	///** X,Y coordinates and other measures in pixels of the rotated original image rotated, these are used to recalculate pixels to the non-rotated orig. image
	//*/
	//	int				_lowerYEndOfBestCentralPath, _upperYEndOfBestCentralPath;

	//	int				_endOfThePlastronYDownCoordinate, _endOfThePlastronYUpCoordinate;

	//	int				_abdominalJunctionsYDownCoordinate, _abdominalJunctionsYUpCoordinate;

	//	int				_pectoralFemoralJunctionsYDownCoordinate, _pectoralFemoralJunctionsYUpCoordinate;

	//	int				_coorAbdoDownLeft,_coorAbdoDownRight,_coorAbdoUpLeft,_coorAbdoUpRight,
	//					_coorPectoFemoDownLeft,_coorPectoFemoDownRight,_coorPectoFemoUpLeft, _coorPectoFemoUpRight;


	///** all the helping images
	//*/

 // //CImg<int>	    _narrowStripe, _mediumWideStripe, _wideStripe;





	//	CImg<double> _segmentsLengthsAbdPectFemo, _segmentsLengthsAbdPectFemoAllTortoises;

	//	CImg<int>	_joiningPointsOfsidePathsToTheCentralPath, _coordinatesOfBestCentralPathOnNarrowStripe, _bestJoiningPointsOfsidePathsToTheCentralPath,
	//		_coordinatesOfBestCentralPathOnMediumWideStripe, _coordinatesOfBestCentralPathOnWideStripe,
	//		_sidePathsOnLeftOnMediumWideStripe, _sidePathsOnRightOnMediumWideStripe;

	//	CImg<unsigned char>	  _skeletonizedEdges, _possibleEdgesForTheCentralSeamPath, _skeletonizedEdgesHorizontalWideStripe,
	//		_skeletonizedEdgesWideStripe, _skeletonizedEdgesHorizontalMediumWideStripe,
	//		_possibleEdgesForTheCentralSeamPathWideStripe, _edgesWrappingWideStripe;
	//	CImg<char>  _wideStripeWithDirectionOfEdges, _mediumWideStripeWithDirectionOfEdges, _narrowStripeWithDirectionOfEdges;
	//	CImg<char>  _wideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges, _mediumWideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges, _narrowStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges;
	//
	//	CImg<unsigned char>  _bestCentralPath;

	//	CImg<int>	_edgesWrapping;

	//	CImg<int> _stripeWithJunctions,_evaluationOfSidePathsTrajectories, _positionsOfLeftAndRightJunctions, _sidePathsPositions;
   
  //void Tortoise::AbdSeamJunctionPreciseLocalization()
  //{
  //  m_imgRotatedOriginal.assign(m_imgOriginal.get_rotate(static_cast<float>(m_nRotationAngleToVerticalInDegrees)));
  //  PecToAbdJunctionLocalizationGHTInJunctionPoints(m_imgPossibleJunctionPreciseLocations);
  //  AbdToFemJunctionLocalizationGHTInJunctionPoints(m_imgPossibleJunctionPreciseLocations);
  //  CImg<int> imgWithResults(m_imgRotatedOriginal);
  //  MarkAJunction(m_imgRotatedOriginal,m_sJunctions.s4PecToAbd,1);
  //  MarkAJunction(m_imgRotatedOriginal,m_sJunctions.s5AbdToFem,1);
  //  GetUniformlyResizedImg(m_imgRotatedOriginal,500).save(("possibleJunctions\\" + m_strTortoiseName + "AbdJunctions.bmp").c_str());
  //  //m_sPlastronCentreOnRotatedImg.yCoor = (m_sJunctions.s4PecToAbd.yCoor + m_sJunctions.s5AbdToFem.yCoor)/2;
  //}

  //void Tortoise::PecToAbdJunctionLocalizationGHTInJunctionPoints(CImg<int> imgPossibleJunctionPreciseLocations)
  //{
  //  const int nYCoorLowerBound = max(0,m_sJunctions.s4PecToAbd.yCoor-(m_sJunctions.s5AbdToFem.yCoor-m_sJunctions.s4PecToAbd.yCoor)*2/3);
  //  const int nYCoorUpperBound = m_sPlastronCentreOnRotatedImg.yCoor;
  //  CImg<int> imgUpperHalfTortoise(m_imgRotatedOriginal.get_crop(
  //    m_sPlastronCentreOnRotatedImg.xCoor-m_nPlastronWidth/2, nYCoorLowerBound, 
  //    m_sPlastronCentreOnRotatedImg.xCoor+m_nPlastronWidth/2, nYCoorUpperBound));
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations.mul((-imgPossibleJunctionPreciseLocations).get_threshold(-nYCoorUpperBound));
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations.mul((imgPossibleJunctionPreciseLocations).get_threshold(nYCoorLowerBound));
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations -= nYCoorLowerBound;
  //  //imgPossibleJunctionPreciseLocations.display();
  //  //imgUpperHalfTortoise.display();
  //  int nResizeMaxPixelSize = 200;
  //  float fSigmaBlurCoef = 0.5f;
  //  CImg<int> imgResizedUpperHalf(GetUniformlyResizedImg(imgUpperHalfTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
  //  //imgResizedUpperHalf.display();
  //  imgPossibleJunctionPreciseLocations *= imgResizedUpperHalf.height();
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations /= (nYCoorUpperBound-nYCoorLowerBound);
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations.mul(imgPossibleJunctionPreciseLocations.get_threshold(1));
  //  
  //  //(imgResizedUpperHalf,imgPossibleJunctionPreciseLocations).display();
  //  float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2-8);
  //  float fRightAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2+8);
  //  CImg<unsigned char> imgEdges(GetEdgeImg(imgResizedUpperHalf,50,100,0));
  //  imgEdges *= 255;
  //  string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
  //  string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
  //  imgResizedUpperHalf.save(strGHTImgFileName.c_str());
  //  imgEdges.save(strGHTEdgesFileName.c_str());
  //  GenHoughTrnf ght("ght\\4PecToAbdTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
  //  // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
  //  ght.Init(130,200,5,5,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
  //  ght.runGHT(0);
  //  point2dCoor sAccumMax = {0,0};
  //  int nRotationAngleMax = 0;
  //  float fWidthRatioMax  = 0.0f;
  //  int nPlastronWidht = 0;
  //  ght.bestCandidateFromGivenPoints(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht,imgPossibleJunctionPreciseLocations);//.display();//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
  //  //cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
  //  int nYCoorFound = m_sPlastronCentreOnRotatedImg.yCoor - (imgResizedUpperHalf.height()-sAccumMax.yCoor)*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height();
  //  m_sJunctions.s4PecToAbd.yCoor = m_imgPossibleJunctionPreciseLocations((m_imgPossibleJunctionPreciseLocations - nYCoorFound).abs().get_stats()(5));
  //  //m_sJunctions.s3HumToPec.yCoor = m_sJunctions.s4PecToAbd.yCoor - static_cast<int>(9*fWidthRatioMax*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height());
  //}

  //void Tortoise::AbdToFemJunctionLocalizationGHTInJunctionPoints(CImg<int> imgPossibleJunctionPreciseLocations)
  //{
  //  const int nYCoorLowerBound = m_sPlastronCentreOnRotatedImg.yCoor;
  //  const int nYCoorUpperBound = min(m_imgRotatedOriginal.height()-1,m_sJunctions.s5AbdToFem.yCoor+(m_sJunctions.s5AbdToFem.yCoor-m_sJunctions.s4PecToAbd.yCoor)*2/3);
  //  CImg<int> imgBottomHalfTortoise(m_imgRotatedOriginal.get_crop(
  //    m_sPlastronCentreOnRotatedImg.xCoor-m_nPlastronWidth/2,nYCoorLowerBound,
  //    m_sPlastronCentreOnRotatedImg.xCoor+m_nPlastronWidth/2, nYCoorUpperBound));
  //  imgPossibleJunctionPreciseLocations.mul((imgPossibleJunctionPreciseLocations).get_threshold(nYCoorLowerBound));
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations.mul((-imgPossibleJunctionPreciseLocations).get_threshold(-nYCoorUpperBound));
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations -= nYCoorLowerBound;
  //  //imgPossibleJunctionPreciseLocations.display();
  //  //imgBottomHalfTortoise.display();
  //  int nResizeMaxPixelSize = 200;
  //  float fSigmaBlurCoef = 0.5f;
  //  CImg<int> imgResizedBottomHalf(GetUniformlyResizedImg(imgBottomHalfTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
  //  //imgResizedBottomHalf.display();
  //  imgPossibleJunctionPreciseLocations *= imgResizedBottomHalf.height();
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations /= (nYCoorUpperBound-nYCoorLowerBound);
  //  //imgPossibleJunctionPreciseLocations.display();
  //  imgPossibleJunctionPreciseLocations.mul(imgPossibleJunctionPreciseLocations.get_threshold(1));
  //  //imgPossibleJunctionPreciseLocations.display();
  //  float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedBottomHalf.width())/2-8);
  //  float fRightAccumBndryXCoor = (static_cast<float>(imgResizedBottomHalf.width())/2+8);
  //  CImg<unsigned char> imgEdges(GetEdgeImg(imgResizedBottomHalf,50,100,0));
  //  imgEdges *= 255;
  //  string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
  //  string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
  //  imgResizedBottomHalf.save(strGHTImgFileName.c_str());
  //  imgEdges.save(strGHTEdgesFileName.c_str());
  //  GenHoughTrnf ght("ght\\5AbdToFemTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
  //  // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
  //  ght.Init(130,200,5,5,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
  //  ght.runGHT(0);
  //  point2dCoor sAccumMax = {0,0};
  //  int nRotationAngleMax = 0;
  //  float fWidthRatioMax  = 0.0f;
  //  int nPlastronWidht = 0;
  //  ght.bestCandidateFromGivenPoints(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht,imgPossibleJunctionPreciseLocations);//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
  //  //cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
  //  int nYCoorFound = m_sPlastronCentreOnRotatedImg.yCoor + sAccumMax.yCoor*imgBottomHalfTortoise.height()/imgResizedBottomHalf.height();
  //  m_sJunctions.s5AbdToFem.yCoor = m_imgPossibleJunctionPreciseLocations((m_imgPossibleJunctionPreciseLocations - nYCoorFound).abs().get_stats()(5));
  //}




  //void Tortoise::FemToAnaJunctionLocalization()
  //{
  //  int nLenghtOfAbdSeam = m_sJunctions.s5AbdToFem.yCoor - m_sJunctions.s4PecToAbd.yCoor;
  //  CImg<int> imgBottomThirdTortoise(m_imgRotatedOriginal.get_crop(
  //    m_sPlastronCentreOnRotatedImg.xCoor-nLenghtOfAbdSeam*3/4,m_sJunctions.s5AbdToFem.yCoor, 
  //    m_sPlastronCentreOnRotatedImg.xCoor+nLenghtOfAbdSeam*3/4, min(m_sJunctions.s5AbdToFem.yCoor + 2*nLenghtOfAbdSeam/3, m_imgRotatedOriginal.height()-1 )));
  //  //imgBottomThirdTortoise.display();
  //  int nResizeMaxPixelSize = 200;
  //  float fSigmaBlurCoef = 0.5f;
  //  CImg<int> imgResizedBottomThird(GetUniformlyResizedImg(imgBottomThirdTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
  //  //imgResizedBottomThird.display();
  //  float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedBottomThird.width())/2-8);
  //  float fRightAccumBndryXCoor = (static_cast<float>(imgResizedBottomThird.width())/2+8);
  //  CImg<unsigned char> imgEdges(GetEdgeImg(imgResizedBottomThird,50,100,0));
  //  imgEdges *= 255;
  //  string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
  //  string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
  //  imgResizedBottomThird.save(strGHTImgFileName.c_str());
  //  imgEdges.save(strGHTEdgesFileName.c_str());
  //  GenHoughTrnf ght("ght\\6FemToAnaTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
  //  // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
  //  ght.Init(100,200,5,2,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
  //  ght.runGHT(0);
  //  point2dCoor sAccumMax = {0,0};
  //  int nRotationAngleMax = 0;
  //  float fWidthRatioMax  = 0.0f;
  //  int nPlastronWidht = 0;
  //  ght.bestCandidate(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht);//.display();//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
  //  //cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
  //  m_sJunctions.s6FemToAna.yCoor = m_sJunctions.s5AbdToFem.yCoor + sAccumMax.yCoor*imgBottomThirdTortoise.height()/imgResizedBottomThird.height();
  //}

  //void Tortoise::TailJunctionLocalization()
  //{
  //  int nLenghtOfAbdSeam = m_sJunctions.s5AbdToFem.yCoor - m_sJunctions.s4PecToAbd.yCoor;
  //  CImg<int> imgBottomTortoise(m_imgRotatedOriginal.get_crop(
  //    m_sPlastronCentreOnRotatedImg.xCoor-nLenghtOfAbdSeam*2/3,m_sJunctions.s6FemToAna.yCoor, 
  //    m_sPlastronCentreOnRotatedImg.xCoor+nLenghtOfAbdSeam*2/3,
  //    min(m_sJunctions.s5AbdToFem.yCoor + nLenghtOfAbdSeam*2/3 + nLenghtOfAbdSeam*2/12, m_imgRotatedOriginal.height()-1 )));
  //  //imgBottomTortoise.display();
  //  int nResizeMaxPixelSize = 200;
  //  float fSigmaBlurCoef = 0.5f;
  //  CImg<int> imgResizedBottom(GetUniformlyResizedImg(imgBottomTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
  //  //imgResizedBottom.display();
  //  float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedBottom.width())/2-8);
  //  float fRightAccumBndryXCoor = (static_cast<float>(imgResizedBottom.width())/2+8);
  //  CImg<unsigned char> imgEdges(GetEdgeImg(imgResizedBottom,50,100,0));
  //  imgEdges *= 255;
  //  string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
  //  string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
  //  imgResizedBottom.save(strGHTImgFileName.c_str());
  //  imgEdges.save(strGHTEdgesFileName.c_str());
  //  GenHoughTrnf ght("ght\\7TailTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
  //  // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
  //  ght.Init(90,200,5,2,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
  //  ght.runGHT(0);
  //  point2dCoor sAccumMax = {0,0};
  //  int nRotationAngleMax = 0;
  //  float fWidthRatioMax  = 0.0f;
  //  int nPlastronWidht = 0;
  //  ght.bestCandidate(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht);//.display();//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
  //  //cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
  //  m_sJunctions.s7Tail.yCoor = m_sJunctions.s6FemToAna.yCoor + sAccumMax.yCoor*imgBottomTortoise.height()/imgResizedBottom.height();
  //  cout << m_sJunctions.s7Tail.yCoor << endl;
  //}

  //void Tortoise::HumToPecJunctionLocalization()
  //{
  //  int nLenghtOfAbdSeam = m_sJunctions.s5AbdToFem.yCoor - m_sJunctions.s4PecToAbd.yCoor;
  //  CImg<int> imgUpperHalfTortoise(m_imgRotatedOriginal.get_crop(
  //    m_sPlastronCentreOnRotatedImg.xCoor-nLenghtOfAbdSeam*5/4, max(0,m_sJunctions.s4PecToAbd.yCoor - nLenghtOfAbdSeam*2/3),
  //    m_sPlastronCentreOnRotatedImg.xCoor+nLenghtOfAbdSeam*5/4, m_sJunctions.s4PecToAbd.yCoor));
  //  //imgUpperHalfTortoise.display();
  //  int nResizeMaxPixelSize = 200;
  //  float fSigmaBlurCoef = 0.5f;
  //  CImg<int> imgResizedUpperHalf(GetUniformlyResizedImg(imgUpperHalfTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
  //  //imgResizedUpperHalf.display();
  //  float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2-8);
  //  float fRightAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2+8);
  //  CImg<unsigned char> imgEdges(GetEdgeImg(imgResizedUpperHalf,50,100,0));
  //  imgEdges *= 255;
  //  string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
  //  string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
  //  imgResizedUpperHalf.save(strGHTImgFileName.c_str());
  //  imgEdges.save(strGHTEdgesFileName.c_str());
  //  GenHoughTrnf ght("ght\\3HumToPecTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
  //  // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
  //  ght.Init(130,200,5,2,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
  //  ght.runGHT(0);
  //  point2dCoor sAccumMax = {0,0};
  //  int nRotationAngleMax = 0;
  //  float fWidthRatioMax  = 0.0f;
  //  int nPlastronWidht = 0;
  //  ght.bestCandidate(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht);//.display();//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
  //  //cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
  //  m_sJunctions.s3HumToPec.yCoor = m_sJunctions.s4PecToAbd.yCoor - (imgResizedUpperHalf.height()-sAccumMax.yCoor)*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height();
  //}

  //void Tortoise::GulToHumJunctionLocalization()
  //{
  //  int nLenghtOfAbdSeam = m_sJunctions.s5AbdToFem.yCoor - m_sJunctions.s4PecToAbd.yCoor;
  //  CImg<int> imgUpperHalfTortoise(m_imgRotatedOriginal.get_crop(
  //      m_sPlastronCentreOnRotatedImg.xCoor-nLenghtOfAbdSeam/2, max(m_sJunctions.s4PecToAbd.yCoor - nLenghtOfAbdSeam,0),
  //      m_sPlastronCentreOnRotatedImg.xCoor+nLenghtOfAbdSeam/2, m_sJunctions.s3HumToPec.yCoor-nLenghtOfAbdSeam*2/12)
  //    );
  //  //imgUpperHalfTortoise.display();
  //  int nResizeMaxPixelSize = 200;
  //  float fSigmaBlurCoef = 0.5f;
  //  CImg<int> imgResizedUpperHalf(GetUniformlyResizedImg(imgUpperHalfTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
  //  //imgResizedUpperHalf.display();
  //  float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2-8);
  //  float fRightAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2+8);
  //  CImg<unsigned char> imgEdges(GetEdgeImg(imgResizedUpperHalf,50,100,0));
  //  imgEdges *= 255;
  //  string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
  //  string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
  //  imgResizedUpperHalf.save(strGHTImgFileName.c_str());
  //  imgEdges.save(strGHTEdgesFileName.c_str());
  //  GenHoughTrnf ght("ght\\2GulToHumTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
  //  // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
  //  ght.Init(100,200,5,2,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
  //  ght.runGHT(0);
  //  point2dCoor sAccumMax = {0,0};
  //  int nRotationAngleMax = 0;
  //  float fWidthRatioMax  = 0.0f;
  //  int nPlastronWidht = 0;
  //  ght.bestCandidate(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht);//.display();//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
  //  //cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
  //  m_sJunctions.s2GulToHum.yCoor = m_sJunctions.s3HumToPec.yCoor-nLenghtOfAbdSeam*2/12 - (imgResizedUpperHalf.height()-sAccumMax.yCoor)*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height();
  //}

  //void Tortoise::HeadJunctionLocalization()
  //{
  //  int nLenghtOfAbdSeam = m_sJunctions.s5AbdToFem.yCoor - m_sJunctions.s4PecToAbd.yCoor;
  //  CImg<int> imgUpperHalfTortoise(m_imgRotatedOriginal.get_crop(
  //      m_sPlastronCentreOnRotatedImg.xCoor-nLenghtOfAbdSeam/2, max(0,m_sJunctions.s4PecToAbd.yCoor - nLenghtOfAbdSeam*5/4),
  //      m_sPlastronCentreOnRotatedImg.xCoor+nLenghtOfAbdSeam/2, m_sJunctions.s2GulToHum.yCoor-nLenghtOfAbdSeam/8)
  //    );
  //  //imgUpperHalfTortoise.display();
  //  int nResizeMaxPixelSize = 200;
  //  float fSigmaBlurCoef = 0.5f;
  //  CImg<int> imgResizedUpperHalf(GetUniformlyResizedImg(imgUpperHalfTortoise.get_blur(fSigmaBlurCoef),nResizeMaxPixelSize));
  //  //imgResizedUpperHalf.display();
  //  float fLeftAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2-8);
  //  float fRightAccumBndryXCoor = (static_cast<float>(imgResizedUpperHalf.width())/2+8);
  //  CImg<unsigned char> imgEdges(GetEdgeImg(imgResizedUpperHalf,50,100,0));
  //  imgEdges *= 255;
  //  string strGHTImgFileName ("ght\\GHTinputImg" + m_strProcessID + ".bmp");
  //  string strGHTEdgesFileName ("ght\\GHTinputEdges" + m_strProcessID + ".bmp");
  //  imgResizedUpperHalf.save(strGHTImgFileName.c_str());
  //  imgEdges.save(strGHTEdgesFileName.c_str());
  //  GenHoughTrnf ght("ght\\1HeadTemplateContour.bmp","ght\\tortoiseTemplateImg.bmp",strGHTEdgesFileName,strGHTImgFileName);
  //  // Init params : widthMinOfContour, widthMaxOfContour, stepOfContourWidht, windowSizeOfAccumulatorDesity, numOfIntervals, leftAccumBndryXCoor, rightAccumBndryXCoor
  //  ght.Init(150,200,5,2,9,fLeftAccumBndryXCoor,fRightAccumBndryXCoor);
  //  ght.runGHT(0);
  //  point2dCoor sAccumMax = {0,0};
  //  int nRotationAngleMax = 0;
  //  float fWidthRatioMax  = 0.0f;
  //  int nPlastronWidht = 0;
  //  ght.bestCandidate(sAccumMax,nRotationAngleMax,fWidthRatioMax,nPlastronWidht);//.display();//.save((strGHTImgFileName.substr(0,strGHTImgFileName.size()-7) + "AbdToFemRes.bmp").c_str());
  //  //cout << sAccumMax.xCoor << " " << sAccumMax.yCoor << endl;
  //  m_sJunctions.s1Head.yCoor = m_sJunctions.s2GulToHum.yCoor-nLenghtOfAbdSeam/8 - (imgResizedUpperHalf.height()-sAccumMax.yCoor)*imgUpperHalfTortoise.height()/imgResizedUpperHalf.height();
  //}

  //CImg<int> Tortoise::DrawJunctions()
  //{
  //  CImg<int> imgResultsImg(GetUniformlyResizedImg(m_imgOriginal,500));
  //  float fResizeRatio = 1.0f*imgResultsImg.height()/m_imgOriginal.height();
  //  imgResultsImg.rotate(static_cast<float>(m_nRotationAngleToVerticalInDegrees));
  //  int nImgResHeight = imgResultsImg.height();
  //  int nDrawnAreaHeight = 2;
  //  int nDrawnAreaWidht = 4;
  //  for(int nJunctionNum = 4; nJunctionNum <= 5; ++nJunctionNum)
  //  {
  //    point2dCoor sOneJunction = {0,0};
  //    GetJunctionCoordinates(nJunctionNum,sOneJunction.xCoor,sOneJunction.yCoor);
  //    sOneJunction.xCoor = static_cast<int>(sOneJunction.xCoor*fResizeRatio);
  //    sOneJunction.yCoor = static_cast<int>(sOneJunction.yCoor*fResizeRatio);
  //    for(int ii = -nDrawnAreaWidht; ii <= nDrawnAreaWidht; ++ii)
  //    {
  //      for(int jj = -nDrawnAreaHeight; jj <= nDrawnAreaHeight; ++jj)
  //      {
  //        if(!(ii>=-nDrawnAreaWidht+1&&ii<=nDrawnAreaWidht-1&&jj<=nDrawnAreaHeight-1&&jj>=-nDrawnAreaHeight+1))
  //        {
  //          if(nJunctionNum>=5)
  //          {
  //            imgResultsImg(sOneJunction.xCoor+ii,min(sOneJunction.yCoor+jj,nImgResHeight-1),0) = 0;
  //            imgResultsImg(sOneJunction.xCoor+ii,min(sOneJunction.yCoor+jj,nImgResHeight-1),1) = 255;
  //            imgResultsImg(sOneJunction.xCoor+ii,min(sOneJunction.yCoor+jj,nImgResHeight-1),2) = 0;
  //          }
  //          else if(nJunctionNum<=4)
  //          {
  //            imgResultsImg(sOneJunction.xCoor+ii,max(0,sOneJunction.yCoor+jj),0) = 255;
  //            imgResultsImg(sOneJunction.xCoor+ii,max(0,sOneJunction.yCoor+jj),1) = 0;
  //            imgResultsImg(sOneJunction.xCoor+ii,max(0,sOneJunction.yCoor+jj),2) = 0;
  //          }
  //          else
  //          {
  //            imgResultsImg(sOneJunction.xCoor+ii,sOneJunction.yCoor+jj,0) = 0;
  //            imgResultsImg(sOneJunction.xCoor+ii,sOneJunction.yCoor+jj,1) = 255;
  //            imgResultsImg(sOneJunction.xCoor+ii,sOneJunction.yCoor+jj,2) = 0;
  //          }
  //        }
  //      }
  //    }

  //  }
  //  return imgResultsImg;
  //}



  //void Tortoise::GetJunctionCoordinates(int nJunctionNumber, int &nXCoor, int &nYCoor) // junction number 1 to 7
  //{
  //  switch (nJunctionNumber)
  //  {
  //      case 1:
  //        nXCoor = m_sJunctions.s1Head.xCoor;
  //        nYCoor = m_sJunctions.s1Head.yCoor;
  //          break;
  //      case 2:
  //        nXCoor = m_sJunctions.s2GulToHum.xCoor;
  //        nYCoor = m_sJunctions.s2GulToHum.yCoor;
  //          break;
  //      case 3:
  //        nXCoor = m_sJunctions.s3HumToPec.xCoor;
  //        nYCoor = m_sJunctions.s3HumToPec.yCoor;
  //          break;
  //      case 4:
  //        nXCoor = m_sJunctions.s4PecToAbd.xCoor;
  //        nYCoor = m_sJunctions.s4PecToAbd.yCoor;
  //          break;
  //      case 5:
  //        nXCoor = m_sJunctions.s5AbdToFem.xCoor;
  //        nYCoor = m_sJunctions.s5AbdToFem.yCoor;
  //          break;
  //      case 6:
  //        nXCoor = m_sJunctions.s6FemToAna.xCoor;
  //        nYCoor = m_sJunctions.s6FemToAna.yCoor;
  //          break;
  //      case 7:
  //        nXCoor = m_sJunctions.s7Tail.xCoor;
  //        nYCoor = m_sJunctions.s7Tail.yCoor;
  //          break;
  //      default:
  //          cout << "Junction number was not in 1 to 7." << endl;
  //          break;
  //  }
  //}




//
//
//
//  CImg<int> Tortoise::drawVerticalStripeAndJunctions() //m_imgRotatedOriginal has to be already initialized
//  {
//    float fImgResolutionKoef = GetImgResolutionCoeff(m_imgOriginal.width(),m_imgOriginal.height());
//    CImg<int> imgResultsImg(m_imgOriginal.get_rotate(m_nRotationAngleToVerticalInDegrees)); // WILL NOT HAVE IT
//    int nResultsImgWidth = imgResultsImg.width();
//    int nResultsImgHeight = imgResultsImg.height();
//    int nLeftXCoorOfStripe = static_cast<int> (m_sPlastronCentreOnRotatedImg.xCoor - fImgResolutionKoef*m_sCentralSeamAreaWidthCoeff.fWide);
//    int nRightXCoorOfStripe = static_cast<int> (m_sPlastronCentreOnRotatedImg.xCoor + fImgResolutionKoef*m_sCentralSeamAreaWidthCoeff.fWide);
//	  int nHalfOfTheLineWidth = max(static_cast<int>(fImgResolutionKoef*0.01),1);
//    for(int x = nLeftXCoorOfStripe; x <= nRightXCoorOfStripe; ++x)
//    {
//      for(int y = nHalfOfTheLineWidth; y <= 3*nHalfOfTheLineWidth; ++y) 
//      {
//        imgResultsImg(x,y,0) = 255;
//      }
//    }
//    for(int x = nLeftXCoorOfStripe; x <= nRightXCoorOfStripe; ++x)
//    {
//      for(int y = nResultsImgHeight-nHalfOfTheLineWidth; y >= nResultsImgHeight-3*nHalfOfTheLineWidth; --y) 
//      {
//        imgResultsImg(x,y,0) = 255;
//      }
//    }
//    for(int y = 0; y < nResultsImgHeight; ++y) 
//    {
//      for(int x = nLeftXCoorOfStripe-nHalfOfTheLineWidth; x <= nLeftXCoorOfStripe+nHalfOfTheLineWidth; ++x)
//      {
//        imgResultsImg(x,y,0) = 255;
//      }
//      for(int x = nRightXCoorOfStripe-nHalfOfTheLineWidth; x <= nRightXCoorOfStripe+nHalfOfTheLineWidth; ++x)
//      {
//        imgResultsImg(x,y,0) = 255;
//      }
//    }
//
//    return imgResultsImg;
//
//  }
//
// //------------	Stripe Edge Images	-------------------------------------------------------------------------------/
//  CImg<int> Tortoise::getMediumWideStripeFromWideStripe(CImg<int> anyWideStripeImg)
//  {
//	  int wideStripeWidth = anyWideStripeImg.width();
//	  int wideStripeHeight = anyWideStripeImg.height();
//	  int mediumWidthInPixels = wideStripeWidth*m_sCentralSeamAreaWidthCoeff.fMedium/m_sCentralSeamAreaWidthCoeff.fWide;
//	  return anyWideStripeImg.get_columns(wideStripeWidth/2-mediumWidthInPixels/2,wideStripeWidth/2+mediumWidthInPixels/2);
//  }
//
//  CImg<int> Tortoise::getNarrowStripeFromWideStripe(CImg<int> anyWideStripeImg)
//  {
//	  int wideStripeWidth = anyWideStripeImg.width();
//	  int wideStripeHeight = anyWideStripeImg.height();
//	  int narrowWidthInPixels = wideStripeWidth*m_sCentralSeamAreaWidthCoeff.fNarrow/m_sCentralSeamAreaWidthCoeff.fWide;
//	  return anyWideStripeImg.get_columns(wideStripeWidth/2-narrowWidthInPixels/2,wideStripeWidth/2+narrowWidthInPixels/2);
//  }
//
//
//  //void Tortoise::centerLocalization()
//  //{
//  //  CImg<int> stripe(getNarrowStripeFromWideStripe(_wideStripe));
//  //  int width = stripe.width();
//  //  int height = stripe.height();
//  //  //CImg<float> grey(stripe.get_RGBtoHSV().channel(2)), gradx(grey.get_gradient(0,2)(0)), grady(grey.get_gradient(0,2)(1)), dir(gradx.get_atan2(grady)*180/PI), mag((gradx.get_pow(2.0) + grady.get_pow(2.0)).sqrt().normalize(0,255).threshold(1) );
//  //  CImg<float> edges(GetEdgeImg(stripe,0,50,0)),dirVH(edges);
//  //  //for(int x = 0; x < width; ++x)
//  //  //{
//  //  //  for(int y = 0; y < height; ++y)
//  //  //  {
//  //  //    if(dir(x,y) < -125 || (dir(x,y) > -55 && dir(x,y) < 55 && dir(x,y)!=0) || dir(x,y) > 125 )
//  //  //    {
//  //  //      dirVH(x,y) = 1;
//  //  //    }
//  //  //    else if(dir(x,y) !=0) dirVH(x,y) =  2;
//  //  //    else dirVH(x,y) =  0;
//  //  //  }
//  //  //}
//  //  //dirVH = EliminateTwoPixelsEdges(EliminateOnePixelEdges(dirVH.threshold(2).mul(edges)));
//  //  //CImg<int> skeletonizedEdges(skeletonizationMinDist(edges));
//  //  CImg<int> skeletonizedEdges(edges);
//  //  CImg<int> bestPath(getBestPathAcrossImgEdges(skeletonizedEdges));
//  //  //(skeletonizedEdges.mul(bestPath+1), stripe.mul(bestPath+1)).display();
//  //  //int size = width/4;
//  //  //for(int y = 0; y < height; ++y)
//  //  //{
//  //  //  for(int x = 0; x < width; ++x)
//  //  //  {
//  //  //    if(bestPath(x,y))
//  //  //    {
//  //  //      for(int i = -size; i <=size; ++i) if(x+i > 0 && x+i < width) bestPath(x+i,y) = 2;
//  //  //      x=width;
//  //  //    }
//  //  //  }
//  //  //}
//  //  //bestPath.display();
//  //  //bestPath.threshold(2);
//  //  (edges,stripe.mul(bestPath+1)).display();
//  //}
//  
//  CImg<int> Tortoise::GetBestSidePathAcrossImgEdges(CImg<int> &imgEdges, bool bMaskIsAsymetric)
//  {
//    int nImgWidth  = imgEdges.width();
//    int nImgHeight = imgEdges.height();
//    int nSupremum = 10000;
//    
//    CImg<int> evaluatedImgMiddlePart(imgEdges.get_fill(nSupremum)),evaluatedImgMiddlePart2(evaluatedImgMiddlePart);
//	  const int decrLevel = 5;
//	  CImg<int> distanceMask(5,decrLevel,1,1, 0);
//	  for(int x = 0; x < 5; x++) for(int y = 0; y < decrLevel; y++) distanceMask(x,y) = (y+3)*(y+3);
//    distanceMask(1,0) = 1;
//    distanceMask.display();
//    imgEdges.display();
//    int maxRatedValue = distanceMask.max();
//    for(int y = nImgHeight-decrLevel-1; y >= nImgHeight/3; y--)
//		{
//      for(int x = 2; x < nImgWidth - 2; x++)
//			{
//				if(imgEdges(x,y))
//				{
//					CImg<int> cropArea(evaluatedImgMiddlePart.get_crop(x-2,y+1,x+2,y+decrLevel) + distanceMask);
//          if(cropArea.min() < nSupremum)
//					{
//						evaluatedImgMiddlePart(x,y) = cropArea.min();
//            if(evaluatedImgMiddlePart(x,y) > maxRatedValue) maxRatedValue = evaluatedImgMiddlePart(x,y);
//					}
//          else
//          {
//            evaluatedImgMiddlePart(x,y) = maxRatedValue;
//          }
//				}
//			}
//		}
//    for(int y = nImgHeight-1; y >= 0; y--)
//		{
//      for(int x = 0; x < nImgWidth; x++)
//			{
//				if(evaluatedImgMiddlePart(x,y)==nSupremum)
//				{
//          evaluatedImgMiddlePart(x,y) = 0;
//				}
//			}
//		}
//    distanceMask(1,0) = 9;
//    distanceMask(3,0) = 1;
//    //evaluatedImgMiddlePart.display();
//    imgEdges.rotate(180);
//    (evaluatedImgMiddlePart,imgEdges).display();
//    maxRatedValue = distanceMask.max();
//    for(int y = nImgHeight-decrLevel-1; y >= nImgHeight/3; y--)
//		{
//      for(int x = 2; x < nImgWidth - 2; x++)
//			{
//				if(imgEdges(x,y))
//				{
//					CImg<int> cropArea(evaluatedImgMiddlePart2.get_crop(x-2,y+1,x+2,y+decrLevel) + distanceMask);
//          if(cropArea.min() < nSupremum)
//					{
//						evaluatedImgMiddlePart2(x,y) = cropArea.min();
//            if(evaluatedImgMiddlePart2(x,y) > maxRatedValue) maxRatedValue = evaluatedImgMiddlePart2(x,y);
//					}
//          else
//          {
//            evaluatedImgMiddlePart2(x,y) = maxRatedValue;
//          }
//				}
//			}
//		}
//    for(int y = nImgHeight-1; y >= 0; y--)
//		{
//      for(int x = 0; x < nImgWidth; x++)
//			{
//				if(evaluatedImgMiddlePart2(x,y)==nSupremum)
//				{
//          evaluatedImgMiddlePart2(x,y) = 0;
//				}
//			}
//		}
//    (evaluatedImgMiddlePart,evaluatedImgMiddlePart2).display();
//    evaluatedImgMiddlePart2 = evaluatedImgMiddlePart + evaluatedImgMiddlePart2.rotate(180);
//    for(int y = nImgHeight-1; y >= 0; y--)
//		{
//      for(int x = 0; x < nImgWidth; x++)
//			{
//				if(evaluatedImgMiddlePart2(x,y)==2*nSupremum)
//				{
//          evaluatedImgMiddlePart2(x,y) = 0;
//				}
//			}
//		}
//    (evaluatedImgMiddlePart2,imgEdges).display();
//  //  int X = evaluatedImgMiddlePart.get_row(height/5).get_stats()(4);
//  //  int Y = height/5;
//		//while(Y < height*4/5)
//		//{
//		//	evaluatedImgMiddlePart(X,Y) = supremum*2;
//		//	CImg<double> cropArea(evaluatedImgMiddlePart.get_crop(X-2,Y+1,X+2,Y+decrLevel) + distanceMask);
//		//	X = cropArea.get_stats()(0,4) + X-2;
//		//	Y = cropArea.get_stats()(0,5) + Y+1;
//		//}
//  //  evaluatedImgMiddlePart.threshold(supremum*2).display();
//    return evaluatedImgMiddlePart2;
//
//  //  return evaluatedEdgesUpperImgPart;
//
//  //  maxRatedValue = distanceMask.max();
//  //  distanceMask.rotate(180);
//		//for(int y = height*3/8; y < height; y++)
//		//{
//		//	for(int x = 2; x < width - 2; x++) 
//		//	{
//		//		if(edgeImg(x,y))
//		//		{
//		//			CImg<int> cropArea(evaluatedEdgesLowerImgPart.get_crop(x-2,y-decrLevel,x+2,y-1) + distanceMask);
//		//			if(cropArea.min() < supremum)
//		//			{
//		//				evaluatedEdgesLowerImgPart(x,y) = cropArea.min();
//  //          if(evaluatedEdgesLowerImgPart(x,y) > maxRatedValue) maxRatedValue = evaluatedEdgesLowerImgPart(x,y);
//		//			}
//  //        else
//  //        {
//  //          evaluatedEdgesLowerImgPart(x,y) = maxRatedValue;
//  //        }
//		//		}
//		//	}
//		//}
//  //  for(int y = height*7/16+1; y < height*2/3; y++)
//		//{
//		//	for(int x = 2; x < width - 2; x++) 
//  //    {
//  //      evaluatedEdgesUpperImgPart(x,y) = supremum;
//  //    }
//  //  }
//  //  for(int y = height*7/16-1; y > height*2/8; y--)
//		//{
//		//	for(int x = 2; x < width - 2; x++) 
//  //    {
//  //      evaluatedEdgesLowerImgPart(x,y) = supremum;
//  //    }
//  //  }
//
//  //  (evaluatedEdgesLowerImgPart, evaluatedEdgesUpperImgPart).display();
//  //  evaluatedEdges = evaluatedEdgesLowerImgPart.get_min(evaluatedEdgesUpperImgPart);
//		//for(int y = 0; y < height; y++)
//		//{
//  //    if(evaluatedEdges.get_row(y).min() < supremum)
//		//	{
//  //      evaluatedEdges(evaluatedEdges.get_row(y).get_stats()(4),y) = -1;
//		//	}
//		//}
//  //  evaluatedEdges *= -1;
//  //  evaluatedEdges.threshold(1);
//  //  evaluatedEdges.display();
//
//
//  //  return evaluatedEdges;
//  }
//
//  
//CImg<int> Tortoise::kMeans( CImg<int> img, int K, int nIterations)
//{
//	int width = img.width();
//	int height = img.height();
//
//	int pocatecni_rozmezi_K_prumeru = 0;
//
//	if(K>1)
//	{
//		pocatecni_rozmezi_K_prumeru = 255/(K-1);
//	}
//	else
//	{
//		pocatecni_rozmezi_K_prumeru = 127;
//	}
//
//	int *pole_K_prumeru = new int[3*K];
//	for(int i = 0; i<=3*K-3; i+=3)
//	{
//		pole_K_prumeru[i] = (i/3)*pocatecni_rozmezi_K_prumeru;
//		pole_K_prumeru[i+1] = (i/3)*pocatecni_rozmezi_K_prumeru;
//		pole_K_prumeru[i+2] = (i/3)*pocatecni_rozmezi_K_prumeru;
//	}
//
//	int opakovani = 0;
//
//	if(nIterations != 0)
//	{
//		do
//		{
//			////////////////////////// 	Iterace K-prumeru v zelve.    //////////////////////////////////////////
//
//			double *pole_K_suma = new double[3*K];
//			for(int i = 0; i<=3*K-1; i++)pole_K_suma[i] = 0;
//
//			int *pole_poctu = new int[K];
//			for(int i = 0; i<=K-1; i++)pole_poctu[i] = 0;
//
//			for(int y = 0; y < height ; y++)
//			{
//				for(int x = 0; x < width; x++)
//				{
//					int min_vzdalenost = 3*pow(255.0,2); 
//					int index = 0;
//					for(int i = 0; i<=3*K-3; i+=3)
//					{
//						int vzdalenost = pow(img(x,y,0)-pole_K_prumeru[i],2.0) + pow(img(x,y,1)-pole_K_prumeru[i+1],2.0) + pow(img(x,y,2)-pole_K_prumeru[i+2],2.0);
//						if(min_vzdalenost>vzdalenost)
//						{
//							min_vzdalenost = vzdalenost; 
//							index = i;
//						}
//					}
//					pole_K_suma[index] += img(x,y,0);
//					pole_K_suma[index+1] += img(x,y,1);
//					pole_K_suma[index+2] += img(x,y,2);
//					pole_poctu[index/3]++;
//				}
//			}
//		
//
//			for(int i = 0; i<=3*K-3; i+=3)
//			{
//				if(pole_poctu[i/3]!=0)
//				{
//					pole_K_prumeru[i] = pole_K_suma[i]/pole_poctu[i/3];
//					pole_K_prumeru[i+1] = pole_K_suma[i+1]/pole_poctu[i/3];
//					pole_K_prumeru[i+2] = pole_K_suma[i+2]/pole_poctu[i/3];
//				}
//			}
//			opakovani++;
//
//		}while(opakovani < nIterations);
//	}
//  CImg<int> res(width,height,1,1, 0);
//
//  for(int y = 0; y < height ; y++)
//	{
//		for(int x = 0; x < width; x++)
//		{
//			int min_vzdalenost = 3*pow(255.0,2); 
//			int index = 0;
//			for(int i = 0; i<=3*K-3; i+=3)
//			{
//				int vzdalenost = pow(img(x,y,0)-pole_K_prumeru[i],2.0) + pow(img(x,y,1)-pole_K_prumeru[i+1],2.0) + pow(img(x,y,2)-pole_K_prumeru[i+2],2.0);
//				if(min_vzdalenost>vzdalenost)
//				{
//					min_vzdalenost = vzdalenost; 
//					index = i;
//				}
//			}
//			res(x,y) = index;
//		}
//	}
//  return res;
//}
//
//
//  void Tortoise::localizationUsingGHT()
//  {
//
//    //contour.threshold(125).normalize(0,255);
//    //int xMean = 0;
//    //int yMean = 0;
//    //int countOfRefPoints = 0;
//    //for(int x = 0; x < contour.width(); x++)
//    //{
//    //  for(int y = 0; y < contour.height(); y++)
//    //  {
//    //    if(contour(x,y))
//    //    {
//    //      xMean += x;
//    //      yMean += y;
//    //      countOfRefPoints++;
//    //    }
//    //  }
//    //}
//    //xMean /= countOfRefPoints;
//    //yMean /= countOfRefPoints;
//    //contour(xMean,yMean,0) = 127;
//    //contour(xMean,yMean,1) = 127;
//    //contour(xMean,yMean,2) = 127;
//    //contour.display();
//    //contour.save("tortoiseCountour_Ref.bmp");
//
//  }
//
//  void Tortoise::junctionLocalization()
//  {
//
//    CImg<int> stripe(_wideStripe),stripeEdges(GetEdgeImg(stripe, 30,90,0));
//    int width = stripe.width();
//    int height = stripe.height();
//    CImg<float> grey(stripe.get_RGBtoHSV().channel(2)), gradx(grey.get_gradient(0,2)(0)), grady(grey.get_gradient(0,2)(1)), dir(gradx.get_atan2(grady)*180/PI), mag((gradx.get_pow(2.0) + grady.get_pow(2.0)).sqrt().normalize(0,255).threshold(1) );
//    CImg<float> dirVH(dir);
//    for(int x = 0; x < width; ++x)
//    {
//      for(int y = 0; y < height; ++y)
//      {
//        if(dir(x,y) < -105 || (dir(x,y) > -75 && dir(x,y) < 75 && dir(x,y)!=0) || dir(x,y) > 105 )
//        {
//          dirVH(x,y) = 2;
//        }
//        else if(dir(x,y) !=0) dirVH(x,y) =  1;
//        else dirVH(x,y) =  0;
//      }
//    }
//    dirVH.threshold(2);
//    
//    //stripeEdges.display();
//    //(stripeEdges).display();
//    //CImg<int> skeletonizedEdges(skeletonizationMinDist(stripeEdges)), centralPath(getBestCentralPath(skeletonizedEdges));
//    //CImg<int> coordinatesOfCentralPath = getCoordinatesOfCentralPath(centralPath);
//    //(stripeEdges + 2*centralPath).display();
//    //stripeEdges.display();
//    //stripeEdges.mul(dirVH);
//    int absoluteMaxValue = 10000;
//   // CImg<int> ratedWideStripeEdgesUpper(stripeEdges.get_fill(absoluteMaxValue));
//   // CImg<int> ratedWideStripeEdgesLower(stripeEdges.get_fill(absoluteMaxValue)),ratedWideStripeEdges;
//    //decrLevel ... number of lines that we look backward(forward) to find more (suitable) edge points 
//	  const int decrLevel = 5;
//    //distanceMask is filled by distances
//	  CImg<int> distanceMask(5,decrLevel,1,1, 0);
//	  for(int x = 0; x < 5; x++) for(int y = 0; y < decrLevel; y++) distanceMask(x,y) = (x-2)*(x-2) + (y+1)*(y+1)*(y+1);
//    int maxRatedValue = distanceMask.max();
//    distanceMask.display();
//    //distanceMask.display();
//		//for(int y = height/2; y > 0; y--)
//		//{
//		//	for(int x = 2; x < width - 2; x++) 
//		//	{
//		//		if(stripeEdges(x,y))
//		//		{
//		//			CImg<int> cropArea(ratedWideStripeEdgesUpper.get_crop(x-2,y+1,x+2,y+decrLevel) + distanceMask);
//		//			if(cropArea.min() < absoluteMaxValue)
//		//			{
//		//				ratedWideStripeEdgesUpper(x,y) = cropArea.min();
//  //          if(ratedWideStripeEdgesUpper(x,y) > maxRatedValue) maxRatedValue = ratedWideStripeEdgesUpper(x,y);
//		//			}
//  //        else
//  //        {
//  //          ratedWideStripeEdgesUpper(x,y) = maxRatedValue;
//  //        }
//		//		}
//		//	}
//		//}
//  //  maxRatedValue = distanceMask.max();
//  //  distanceMask.rotate(180);//.display();
//		//for(int y = height*3/8; y < height; y++)
//		//{
//		//	for(int x = 2; x < width - 2; x++) 
//		//	{
//		//		if(stripeEdges(x,y))
//		//		{
//		//			CImg<int> cropArea(ratedWideStripeEdgesLower.get_crop(x-2,y-decrLevel,x+2,y-1) + distanceMask);
//		//			if(cropArea.min() < absoluteMaxValue)
//		//			{
//		//				ratedWideStripeEdgesLower(x,y) = cropArea.min();
//  //          if(ratedWideStripeEdgesLower(x,y) > maxRatedValue) maxRatedValue = ratedWideStripeEdgesLower(x,y);
//		//			}
//  //        else
//  //        {
//  //          ratedWideStripeEdgesLower(x,y) = maxRatedValue;
//  //        }
//		//		}
//		//	}
//		//}
//  //  for(int y = height*7/16+1; y < height*2/3; y++)
//		//{
//		//	for(int x = 2; x < width - 2; x++) 
//  //    {
//  //      ratedWideStripeEdgesUpper(x,y) = absoluteMaxValue;
//  //    }
//  //  }
//  //  for(int y = height*7/16-1; y > height*2/8; y--)
//		//{
//		//	for(int x = 2; x < width - 2; x++) 
//  //    {
//  //      ratedWideStripeEdgesLower(x,y) = absoluteMaxValue;
//  //    }
//  //  }
//  //  //(ratedWideStripeEdgesLower, ratedWideStripeEdgesUpper).display();
//  //  ratedWideStripeEdges = ratedWideStripeEdgesLower.get_min(ratedWideStripeEdgesUpper);
//
//		//for(int y = 0; y < height; y++)
//		//{
//  //    int min = absoluteMaxValue;
//		//	for(int x = 0; x < width; x++) 
//		//	{
//  //      if(ratedWideStripeEdges(x,y) < min) 
//  //      {
//  //        min = ratedWideStripeEdges(x,y);
//  //      }
//  //    }
//  //    int countOfMin = 0;
//  //    int sumOfXvalues = 0;
//  //    for(int x = 0; x < width; x++) 
//		//	{
//  //      if(ratedWideStripeEdges(x,y) == min) 
//  //      {
//  //        countOfMin++;
//  //        sumOfXvalues += x;
//  //      }
//  //    }
//  //    ratedWideStripeEdges((int)sumOfXvalues/countOfMin,y) = maxRatedValue*2;
//  //    ratedWideStripeEdges(width-1,y) = (int)sumOfXvalues/countOfMin;
//  //  }
//  //  CImg<int> coordinates(ratedWideStripeEdges.get_columns(width-2,width-1));
//  //  //coordinates.display();
//		//for(int y = 0; y < height; y++)
//		//{
//		//	for(int x = 0; x < width; x++) 
//		//	{
//  //      if(ratedWideStripeEdges(x,y) == absoluteMaxValue) ratedWideStripeEdges(x,y) = 0;
//  //    }
//  //  }
//		//for(int y = 10; y < height-10; y++)
//		//{
//  //    int x = coordinates(1,y);
//  //    int cropSize = GetImgResolutionCoeff(m_nWidth,m_nHeight)/15;
//  //    CImg<int> cropArea1(_wideStripe.get_crop(x-cropSize/2,y-cropSize-1,x+cropSize/2,y-1)), cropArea2(_wideStripe.get_crop(x-cropSize/2,y+1,x+cropSize/2,y+cropSize+1));
//  //    cropArea1 = (cropArea1-cropArea2).abs();
//  //    cropArea1 = cropArea1.get_channel(0) + cropArea1.get_channel(1) + cropArea1.get_channel(2);
//  //    coordinates(0,y) = cropArea1.sum();
//  //  }
//    //(_wideStripe,coordinates).display("VERTIKALNE");
//    //ratedWideStripeEdges.threshold(maxRatedValue*2);
//    //ratedWideStripeEdges.display("VERTIKALNE");
//
//    cout << absoluteMaxValue << " " << maxRatedValue << endl;
//    distanceMask.rotate(180);
//    stripeEdges.rotate(90);
//    CImg<int> ratedWideStripeEdgesUpperRotated(stripeEdges.get_fill(absoluteMaxValue));
//    CImg<int> ratedWideStripeEdgesLowerRotated(stripeEdges.get_fill(absoluteMaxValue)),ratedWideStripeEdgesRotated;
//    int widthRotated= stripeEdges.width();
//	  int heightRotated = stripeEdges.height();
//    maxRatedValue = distanceMask.max();
//		for(int y = decrLevel; y < heightRotated - decrLevel; y++)
//		{
//			for(int x = 2; x < widthRotated - 2; x++) 
//			{
//				if(stripeEdges(x,y))
//				{
//					CImg<int> cropArea(ratedWideStripeEdgesLowerRotated.get_crop(x-2,y-decrLevel,x+2,y-1) + distanceMask);
//					if(cropArea.min() < absoluteMaxValue)
//					{
//						ratedWideStripeEdgesLowerRotated(x,y) = cropArea.min();
//            if(ratedWideStripeEdgesLowerRotated(x,y) > maxRatedValue) maxRatedValue = ratedWideStripeEdgesLowerRotated(x,y);
//					}
//          else
//          {
//            ratedWideStripeEdgesLowerRotated(x,y) = maxRatedValue;
//          }
//				}
//			}
//		}
//    maxRatedValue = distanceMask.max();
//    distanceMask.rotate(180);
//		for(int y = heightRotated - decrLevel -1; y > decrLevel-1; y--)
//		{
//			for(int x = 2; x < widthRotated - 2; x++) 
//			{
//				if(stripeEdges(x,y))
//				{
//					CImg<int> cropArea(ratedWideStripeEdgesUpperRotated.get_crop(x-2,y+1,x+2,y+decrLevel) + distanceMask);
//					if(cropArea.min() < absoluteMaxValue)
//					{
//						ratedWideStripeEdgesUpperRotated(x,y) = cropArea.min();
//            if(ratedWideStripeEdgesUpperRotated(x,y) > maxRatedValue) maxRatedValue = ratedWideStripeEdgesUpperRotated(x,y);
//					}
//          else
//          {
//            ratedWideStripeEdgesUpperRotated(x,y) = maxRatedValue;
//          }
//				}
//			}
//		}
//    ratedWideStripeEdgesRotated = ratedWideStripeEdgesLowerRotated + ratedWideStripeEdgesUpperRotated;
//    //ratedWideStripeEdgesRotated.display();
//    ratedWideStripeEdgesRotated.rotate(-90);
//
//    CImg<int> coordOfHorWdw(1,height,1,1, 1000);
//    //int nWdw = 0;
//    //for(int y = width/4; y < height - width/4-1; y+=1)
//    //{
//    //  CImg<int> crop(ratedWideStripeEdgesRotated.get_crop(0,y-width/16,width-1,y+width/16));
//    //  int sum = crop.sum();
//    //  int nNonZeros = crop.threshold(1).sum();
//    //  coordOfHorWdw(0,y) = sum/nNonZeros;
//    //}
//    //coordOfHorWdw.display_graph();
//    CImg<int> centralStripe(ratedWideStripeEdgesRotated.get_columns(width/2-width/16,width/2+width/16));
//    for(int y = width/4; y < height - width/4-1; y+=1)
//    {
//      coordOfHorWdw(0,y) = centralStripe.get_row(y).min();
//    }
//    for(int x = 0; x < width; ++x)
//    {
//      for(int y = 0; y < height; ++y)
//      {
//        if(ratedWideStripeEdgesRotated(x,y)==2*absoluteMaxValue)ratedWideStripeEdgesRotated(x,y)=0;
//      }
//    }
//    centralStripe.display();
//    coordOfHorWdw(height/2) = coordOfHorWdw.max()*3/2;
//    coordOfHorWdw.display_graph();
//    (stripe,ratedWideStripeEdgesRotated,coordOfHorWdw).display();
//     coordOfHorWdw.display_graph();
//
//
//    //ratedWideStripeEdgesRotated.display();
//    //(ratedWideStripeEdgesLower,ratedWideStripeEdgesUpper,ratedWideStripeEdges).display();
//		//for(int y = 0; y < height; y++)
//		//{
//  //    int xMin = 0;
//  //    int min = absoluteMaxValue;
//		//	for(int x = 0; x < width; x++) 
//		//	{
//  //      if(ratedWideStripeEdges(x,y) < min) 
//  //      {
//  //        xMin = x;
//  //        min = ratedWideStripeEdges(x,y);
//  //      }
//  //    }
//  //    ratedWideStripeEdges(xMin,y) = maxRatedValue*2;
//  //    ratedWideStripeEdges(width-1,y) = min;
//    
///*
//
//    cout << absoluteMaxValue << " " << maxRatedValue << endl;
//    (ratedWideStripeEdges,ratedWideStripeEdgesRotated).display();
//    ratedWideStripeEdges = ratedWideStripeEdges.rotate(90) + ratedWideStripeEdgesRotated;
//    ratedWideStripeEdges.rotate(-90).display();
//    coordinatesOfCentralPath.display();
//    for(int y = height/3; y < height; y++)
//    {
//	    if(coordinatesOfCentralPath(y)==0)
//	    {
//        if(y < height-1)
//        {
//		      if(coordinatesOfCentralPath(y-1)>0 && coordinatesOfCentralPath(y+1)==0) coordinatesOfCentralPath(y) = coordinatesOfCentralPath(y-1);
//		      else if(coordinatesOfCentralPath(y-1)>0 && coordinatesOfCentralPath(y+1)>0) coordinatesOfCentralPath(y) = (coordinatesOfCentralPath(y-1) + coordinatesOfCentralPath(y+1))/2;
//        }
//        else
//        {
//          coordinatesOfCentralPath(y) = coordinatesOfCentralPath(y-1);
//        }
//	    }
//    }
//    for(int y = height*2/3; y >= 0; y--)
//    {
//	    if(coordinatesOfCentralPath(y)==0)
//	    {
//        if(y > 0)
//        {
//		      if(coordinatesOfCentralPath(y-1)==0 && coordinatesOfCentralPath(y+1)>0) coordinatesOfCentralPath(y) = coordinatesOfCentralPath(y+1);
//		      else if(coordinatesOfCentralPath(y-1)>0 && coordinatesOfCentralPath(y+1)>0) coordinatesOfCentralPath(y) = (coordinatesOfCentralPath(y-1) + coordinatesOfCentralPath(y+1))/2;
//        }
//        else
//        {
//          coordinatesOfCentralPath(y) = coordinatesOfCentralPath(y+1);
//        }
//	    }
//    }
//    coordinatesOfCentralPath.display();
//    CImg<int> valuesOnCetralSeam(coordinatesOfCentralPath);
//    for(int y = 0; y < height; y++)
//		{
//      int x = coordinatesOfCentralPath(y);
//      if(x > 0) valuesOnCetralSeam(y) = ratedWideStripeEdgesRotated(x,y);
//    }
//    valuesOnCetralSeam.display();
//    for(int y = height/3; y < height; y++)
//    {
//	    if(valuesOnCetralSeam(y)==0)
//	    {
//        if(y < height-1)
//        {
//		      if(valuesOnCetralSeam(y-1)>0 && valuesOnCetralSeam(y+1)==0) valuesOnCetralSeam(y) = valuesOnCetralSeam(y-1);
//		      else if(valuesOnCetralSeam(y-1)>0 && valuesOnCetralSeam(y+1)>0) valuesOnCetralSeam(y) = (valuesOnCetralSeam(y-1) + valuesOnCetralSeam(y+1))/2;
//        }
//        else
//        {
//          valuesOnCetralSeam(y) = valuesOnCetralSeam(y-1);
//        }
//	    }
//    }
//    for(int y = height*2/3; y >= 0; y--)
//    {
//	    if(valuesOnCetralSeam(y)==0)
//	    {
//        if(y > 0)
//        {
//		      if(valuesOnCetralSeam(y-1)==0 && valuesOnCetralSeam(y+1)>0) valuesOnCetralSeam(y) = valuesOnCetralSeam(y+1);
//		      else if(valuesOnCetralSeam(y-1)>0 && valuesOnCetralSeam(y+1)>0) valuesOnCetralSeam(y) = (valuesOnCetralSeam(y-1) + valuesOnCetralSeam(y+1))/2;
//        }
//        else
//        {
//          valuesOnCetralSeam(y) = valuesOnCetralSeam(y+1);
//        }
//	    }
//    }
//    valuesOnCetralSeam.display();
//
//
//    (ratedWideStripeEdges.rotate(90)).display();
//		for(int y = 0; y < height; y++)
//		{
//			for(int x = 0; x < width; x++) 
//			{
//        if(ratedWideStripeEdges(x,y) == absoluteMaxValue) ratedWideStripeEdges(x,y) = 0;
//      }
//    }
//    (coordinates,_wideStripe).display();*/
//  }
//
//
//
//
//    //  if (m_imgResultsData(0,1) == 1) // coordinates of the the best central path on the wide stripe
//    //{
//    //  _lowerYEndOfBestCentralPath = m_imgResultsData(1,1);
//    //  _upperYEndOfBestCentralPath = m_imgResultsData(2,1);
//    //  //_mediumWideStripeEdges = getMediumWideStripeFromWideStripe(_wideStripeEdges);
//	   // _narrowStripeEdges = getNarrowStripeFromWideStripe(_wideStripeEdges);
//    //  _coordinatesOfBestCentralPathOnNarrowStripe = m_imgResultsData.get_column(m_imgResultsDataSize+_wideStripeWidth);
//	   // _coordinatesOfBestCentralPathOnMediumWideStripe = getCoordinatesOfBestCentralPathOnMediumWideStripe();
//	   // _coordinatesOfBestCentralPathOnWideStripe = getCoordinatesOfBestCentralPathOnWideStripe();
//    //}
//    //else
//    //{
//    //  //_mediumWideStripeEdges = getMediumWideStripeFromWideStripe(_wideStripeEdges);
//	   // _narrowStripeEdges = getNarrowStripeFromWideStripe(_wideStripeEdges);
//	   // _possibleEdgesForTheCentralSeamPath = getPossibleEdgesForTheCentralSeamPath(_narrowStripeEdges);
//	   // _edgesWrapping = edgesWrapping(_possibleEdgesForTheCentralSeamPath);
//	   // _skeletonizedEdges = skeletonizationMinDist(_possibleEdgesForTheCentralSeamPath);
//    //  _bestCentralPath = getBestCentralPath(_skeletonizedEdges);
//    //  _coordinatesOfBestCentralPathOnNarrowStripe = getCoordinatesOfBestCentralPathOnNarrowStripe();
//	   // _coordinatesOfBestCentralPathOnMediumWideStripe = getCoordinatesOfBestCentralPathOnMediumWideStripe();
//	   // _coordinatesOfBestCentralPathOnWideStripe = getCoordinatesOfBestCentralPathOnWideStripe();
//    //  m_imgResultsData(0,1) = 1;
//    //  m_imgResultsData(1,1) = _lowerYEndOfBestCentralPath;
//    //  m_imgResultsData(2,1) = _upperYEndOfBestCentralPath;
//    //  m_imgResultsData.append(_coordinatesOfBestCentralPathOnNarrowStripe,'x');
//
//    //  cout << "PiiP_1" << endl;
//    //}  
//
//  
// //------------	 Best Central Path and its coordinates	-------------------------------------------------------------------------------/
//  CImg<unsigned char> Tortoise::getPossibleEdgesForTheCentralSeamPath(CImg<unsigned char> edgeImgToDoItOn)
//{ //finds edges that are close enough to each other (10 pixs verticaly and 2 pixs horizontaly) 
//	CImg<unsigned char> res(edgeImgToDoItOn); //starting with the stripe edge image
//	unsigned int width = res.width();
//	unsigned int height = res.height();
//// centerY is the approximate tortoise plastron center Y coordinate
//	unsigned int centerY = m_nHeightOfRotatedOriginalImg/2;
////decrLevel is number of pixels we go backward when looking for previous edge pixels
//	const int decrLevel = 10;
//	for(int y = -2; y <= 2; y++) for(int x = 0; x < width; x++) { if(res(x,centerY+y)) res(x,centerY+y) = 2; }
////going through all pixs from the middle upwards
//	for(int y = centerY - 3; y > 0; y--)
//	{
//		for(int x = 2; x < width - 2; x++) 
//		{
////we go though all edge pixels in a row
//			if(res(x,y))
//			{
////edge pixel is marked by number 2, if backwards there is an edge pixel with a number 2 (we look 2 pix on left, 2 on right and 10 pixs back)
//				for(int decr = 1; decr <= decrLevel; decr++)
//				{
//					if(res(x-2,y+decr) == 2 || res(x-1,y+decr) == 2 || res(x,y+decr) == 2 || 
//						res(x+1,y+decr) == 2|| res(x+2,y+decr) == 2) 
//					{ 
//						res(x,y) = 2; 
//						decr = decrLevel+1;
//					}
//				}
//			}
//		}
//	}
////going through all pixs from the middle downwards
//	for(int y = centerY + 3; y < height; y++)
//	{
//		for(int x = 2; x < width - 2; x++) 
//		{
//			if(res(x,y))
//			{
//				for(int decr = 1; decr <= decrLevel; decr++)
//				{
//					if(res(x-2,y-decr) == 2 || res(x-1,y-decr) == 2 || res(x,y-decr) == 2 || 
//						res(x+1,y-decr) == 2 || res(x+2,y-decr) == 2) 
//					{ 
//						res(x,y) = 2; 
//						decr = decrLevel+1;
//
//					}
//				}
//			}
//		}
//	}
//	return 	res.threshold(2);
//}
//
//  CImg<unsigned char> Tortoise::getEdgesWrapping(CImg<unsigned char> imgToDoItOn)
//{
//	int width = imgToDoItOn.width();
//	int height = imgToDoItOn.height();
//	CImg<int> res(imgToDoItOn.get_fill(0));
//	for(int y = 5; y < height - 6; y++)
//	{
//		for(int x = 5; x < width-6; x++)
//		{
//			for(int j = y + 5; j >= y - 5; j--)
//			{
//				for(int i = x+5; i >= x-5; i--)
//				{
//					if(imgToDoItOn(i,j)>0)
//					{
//						res(x,y) = 1;
//						i = x -6;
//						j = y -6;
//					}
//				}
//			}
//		}
//	}
//	return res;
//}
//
//  CImg<unsigned char> Tortoise::skeletonizationMinDist(CImg<unsigned char> imgToDoItOn)
//{	
//	imgToDoItOn = EliminateTwoPixelsEdges(EliminateOnePixelEdges(imgToDoItOn));
//	CImg<int> res(imgToDoItOn),edgesWrapping(getEdgesWrapping(imgToDoItOn)), edgesWrappingInverse;
//	int width = res.width();
//	int height = res.height();
////sets value 0 for edges in resulting image
//	res = (res-1).abs();
////sets value 0 for the area it works
//  //edgesWrapping.display();
//	edgesWrappingInverse = (edgesWrapping-1).abs();
//  //edgesWrappingInverse.display();
////adds edgesWrapping, so 0s are only in the working area where there is no edge
//  //imgToDoItOn.display();
//	imgToDoItOn = imgToDoItOn + edgesWrappingInverse;
//  //imgToDoItOn.display();
//
////
//	for(int y = 0; y < height-1; y++)
//	{
//		for(int x = 0; x < width-1; x++)
//		{
//			if(imgToDoItOn(x,y) == 0)
//			{
//				int minDist = 1;
//				bool doIncreaseDistance = 1;
//				while(doIncreaseDistance && x-minDist>=0 && y-minDist>=0 && x+minDist<width && y+minDist<height)
//				{
//					if(imgToDoItOn.get_crop(x-minDist,y-minDist,x+minDist,y+minDist).sum() > 0 || minDist >= 3)
//					{
//						doIncreaseDistance = 0;
//					}
//					else minDist++;
//				}
//				res(x,y) = minDist;
//			}
//		}
//	}
//	res = res.get_mul(edgesWrapping);
////eliminates pixs, first those close to edges (cycle 1->3)
//	int maxLevel = 3;
//	for(int cycle = 1; cycle <= maxLevel; cycle++)
//	{
////in each cycle it eliminates also all newly-suitable pixs on lower levels (level 1...cycle)
//		for(int level = 1; level <= cycle; level++)
//		{
////eliminates until no elimination has occured in previous go-through of the image
//			bool change = 1;
//			while(change)
//			{
//				change = 0;
//				for(int y = 1; y < height - 1; y++)
//				{
//					for(int x = 1; x < width-1; x++)
//					{
//						if(res(x,y) == level) 
//						{	
//// eliminate from rigth down corner
//							if( res(x-1,y-1)> 0 && res(x,y-1)> 0 && res(x+1,y-1)>=0  &&
//								res(x-1,y)  > 0 &&				    res(x+1,y)  ==0  &&
//								res(x-1,y+1)>=0 && res(x,y+1)==0 && res(x+1,y+1)>=0     ) { res(x,y) = 0; change = 1; }
//// eliminate from left down corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)> 0 && res(x+1,y-1)> 0 &&
//									 res(x-1,y)  ==0 &&				    res(x+1,y)  > 0 &&
//									 res(x-1,y+1)>=0 && res(x,y+1)==0 && res(x+1,y+1)>=0  ) { res(x,y) = 0; change = 1; }
//// eliminate from left up corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)==0 && res(x+1,y-1)>=0 &&
//									 res(x-1,y)  ==0 &&				    res(x+1,y)  > 0 &&
//									 res(x-1,y+1)>=0 && res(x,y+1)> 0 && res(x+1,y+1)> 0   ) { res(x,y) = 0; change = 1; }
//// eliminate from right up corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)==0 && res(x+1,y-1)>=0 &&
//									 res(x-1,y)  > 0 &&					res(x+1,y)  ==0 &&
//									 res(x-1,y+1)> 0 && res(x,y+1)> 0 && res(x+1,y+1)>=0   ) { res(x,y) = 0; change = 1; }
//// eliminate from right side
//							else if( res(x-1,y-1)>0 &&  (res(x+1,y-1)==0||res(x+1,y+1)==0) &&
//									 res(x-1,y)  >0 &&				 res(x+1,y)  ==0 &&
//									 res(x-1,y+1)>0   ) { res(x,y) = 0; change = 1; }
//// eliminate from left side
//							else if( (res(x-1,y-1)==0||res(x-1,y+1)==0) && res(x+1,y-1)>0 &&
//									  res(x-1,y)  ==0 &&					                res(x+1,y) >0 &&
//								      res(x+1,y+1)>0  ) { res(x,y) = 0; change = 1; }
//// eliminate from up side
//							else if( (res(x-1,y-1)==0 || res(x+1,y-1)==0) && res(x,y-1)==0 && 
//									 (res(x-1,y)  >0 ||					 res(x+1,y)  >0) &&
//									  res(x-1,y+1)>0 && res(x,y+1)>0 && res(x+1,y+1)>0  ) { res(x,y) = 0; change = 1; }
//// eliminate from down side
//							else if( res(x-1,y-1)>0 && res(x,y-1)>0 && res(x+1,y-1)>0 &&
//									(res(x-1,y)  >0 ||					 res(x+1,y)>0) &&
//									(res(x-1,y+1)==0 || res(x+1,y+1)==0) && res(x,y+1)==0   ) { res(x,y) = 0; change = 1; }
//							else
//							{		
//
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	return res.threshold(1);
//}
//
//  CImg<unsigned char> Tortoise::getBestCentralPath(CImg<unsigned char> edgeImgToDoItOn)
//{
//	CImg<double> distanceMask,paths(edgeImgToDoItOn),res(paths);
//// paths contains skeleton of edge stripe image without large empty spaces
//	int width = paths.width();
//	int height = paths.height();
//	int centerX = width/2;
//	int centerY = m_nHeightOfRotatedOriginalImg/2;
////decrLevel ... number of lines that we look backward(forward) to find more (suitable) edge points 
//	const int decrLevel = 10;
////distanceMask is filled by distances
//	distanceMask.assign(5,decrLevel,1,1, 0);
//	for(int x = 0; x < 5; x++) for(int y = 0; y < decrLevel; y++) distanceMask(x,y) = (x-2)*(x-2) + (y+1)*(y+1);
////starting point
//	int startX = centerX;
//	int startY = 0;
//	int endX = centerX;
//	int endY = height -1;
////locating starting point
//	while(res(startX,startY)==0)
//	{
//		CImg<double> startApprox;
//		startApprox = res.get_crop(startX-2,startY+1,startX+2,startY+3);
//		if(startApprox.max()==0) { startY = startY+3; }
//		else{ startX = startX + startApprox.get_stats()(0,8) - 2; startY = startY + startApprox.get_stats()(0,9) + 1;}
//	}
////it is minimalization method, so unimportant pixels need to be a lot above other values (0 -> 99999, 1 -> 100000)
//	res = res + 99999;
////starting point is set to 2
//	res(startX,startY) = 2;
//	distanceMask.rotate(180);
////starting point should be in y-lower part
//	if(startY < centerY) // if starting in the y-lower part (CImg images have y-axis inverted)
//	{
//		for(int y = startY + 1; y < height; y++)
//		{
//			for(int x = 2; x < width - 2; x++) 
//			{
//				if(res(x,y) == 100000)
//				{
//					CImg<double> cropArea;
//					cropArea = res.get_crop(x-2,y-decrLevel,x+2,y-1);
//					if(cropArea.min() < 99999)
//					{
//						cropArea = cropArea + distanceMask;
//						res(x,y) = cropArea.min();
//					}
//				}
//			}
//		}
//	}
//	else 
//	{
//			cout << "ERROR: Starting Point in Best path finding is not under _centerRotateY." << endl;
//			return _bestCentralPath.fill(1);
//	}
//	//res.get_rotate(90).display();
////locating ending point
//	bool theBestEnd = 0;
//	while(!theBestEnd)
//	{
//		CImg<double> endApprox;
//		endApprox = res.get_crop(0,endY-2,width-1,endY);
//		if(endApprox.min()>=99999) { endY = endY-3; }
//		else{ endX = endX + endApprox.get_stats()(0,4) - endX; endY = endY + endApprox.get_stats()(0,5) - 2; theBestEnd = 1;}
//	}
////now we want to find the optimal starting point, the same only starting in the upper part so crop area is above (x,y) and distanceMask has to be rotated by 180degrees
//	res = paths;
//	res = res + 99999;
////we begging at the ending point
//	res(endX,endY) = 2;
//	distanceMask.rotate(180);
//	for(int y = endY - 1; y >= startY; y--)
//	{
//		for(int x = 2; x < width - 2; x++) 
//		{
////edge points have 100000
//			if(res(x,y) == 100000)
//			{
////croping 5*decrLevel area under (x,y) point
//				CImg<int> cropArea;
//				cropArea = res.get_crop(x-2,y+1,x+2,y+decrLevel);
////if(there are some edge points in the area)
//				if(cropArea.min() < 99999)
//				{
////area weighted be pixels distance to (x,y)
//					cropArea = cropArea + distanceMask;
////minimum distance is saved in (x,y)
//					res(x,y) = cropArea.min();
//				}
//			}
//		}
//	}
//	//res.get_rotate(90).display();
//	bool theBestStart = 0;
//	while(!theBestStart)
//	{
//		CImg<double> startApprox;
//		startApprox = res.get_crop(0,startY+1,width-1,startY+3);
//		if(startApprox.min()>=99999) { startY = startY+3; }
//		else{ startX = startApprox.get_stats()(0,4); startY = startY + startApprox.get_stats()(0,5) + 1; theBestStart = 1;}
//	}
//	//cout << startX << " " << startY << " " << endX << " " << endY  << endl;
//	int totalDistance = res(startX,startY);
//	CImg<int> copyOfRes(res);
////beacause it is backward mode, the crop area setup stays the same for both scenearios
//	if(startY < centerY)
//	{
//		while(totalDistance > 2)
//		{
//			res(startX,startY) = 10000000;
//			CImg<double> cropArea;
//			cropArea = res.get_crop(startX-2,startY+1,startX+2,startY+decrLevel);
//			cropArea = cropArea + distanceMask;
//			startX = cropArea.get_stats()(0,4) + startX-2;
//			startY = cropArea.get_stats()(0,5) + startY+1;
//			totalDistance = res(startX,startY);
//		}
//	}
//	else {cout << "neco se nepovedlo." << endl;}
//	if(endX == startX && endY == startY) 
//	{		
//		cout << "OK" << endl;
//		cout << startX << " " << startY << " " << endX << " " << endY << " " << totalDistance << endl;
//	}
//	else
//	{
//		cout << "NOT OUKEJ" << endl;
//		cout << startX << " " << startY << " " << endX << " " << endY << " " << totalDistance << endl;
//	}
//	//copyOfRes.mul(res.get_threshold(10000000)).rotate(90).display();
//	//(paths.get_mul(res.get_threshold(10000000)+1) + ((GetEdgeImg(_verticalStripeWithCentralSeamOfTheOriginalImage,0,30,1))-1).abs()).get_rotate(90).display();
//	return res.threshold(10000000);
//}
//  
//  CImg<int> Tortoise::getCoordinatesOfBestCentralPathOnNarrowStripe()
//{
//	int stripeWidth = _bestCentralPath.width();
//	int stripeHeight = _bestCentralPath.height();
//	_lowerYEndOfBestCentralPath = 0;
//	_upperYEndOfBestCentralPath = 0;
//	CImg<int> coordinatesOfBestCentralPath(1,stripeHeight,1,1, 0);
//	for(int y = 0; y < stripeHeight; y++)
//	{
//		for(int x = 0; x < stripeWidth; x++)
//		{
//			if(_bestCentralPath(x,y)) 
//			{
//				if(_lowerYEndOfBestCentralPath == 0) _lowerYEndOfBestCentralPath = y;
//				_upperYEndOfBestCentralPath = y;
//				coordinatesOfBestCentralPath(y) = x;
//				x = stripeWidth;
//			}
//		}
//	}
//	// filling gaps on the best central path
//	for(int y = _lowerYEndOfBestCentralPath; y < stripeHeight/2; y++)
//	{
//		if(coordinatesOfBestCentralPath(y)==0)
//		{
//			if(coordinatesOfBestCentralPath(y-1)>0 && coordinatesOfBestCentralPath(y+1)==0) coordinatesOfBestCentralPath(y) = coordinatesOfBestCentralPath(y-1);
//			else if(coordinatesOfBestCentralPath(y-1)>0 && coordinatesOfBestCentralPath(y+1)>0) coordinatesOfBestCentralPath(y) = (coordinatesOfBestCentralPath(y-1) + coordinatesOfBestCentralPath(y+1))/2;
//		}
//	}
//	for(int y = _upperYEndOfBestCentralPath; y > stripeHeight/2-2; y--)
//	{
//		if(coordinatesOfBestCentralPath(y)==0)
//		{
//			if(coordinatesOfBestCentralPath(y+1)>0 && coordinatesOfBestCentralPath(y-1)==0) coordinatesOfBestCentralPath(y) = coordinatesOfBestCentralPath(y+1);
//			else if(coordinatesOfBestCentralPath(y+1)>0 && coordinatesOfBestCentralPath(y-1)>0) coordinatesOfBestCentralPath(y) = (coordinatesOfBestCentralPath(y+1) + coordinatesOfBestCentralPath(y-1))/2;
//		}
//	}
//	return coordinatesOfBestCentralPath;
//}
//
//  CImg<int> Tortoise::getCoordinatesOfCentralPath(CImg<int> imgWithCentralPathAsOnes)
//{
//	int width = imgWithCentralPathAsOnes.width();
//	int height = imgWithCentralPathAsOnes.height();
//	CImg<int> coordinatesOfCentralPath(1,height,1,1, 0);
//	for(int y = 0; y < height; y++)
//	{
//		for(int x = 0; x < width; x++)
//		{
//			if(imgWithCentralPathAsOnes(x,y)) 
//			{
//				coordinatesOfCentralPath(y) = x;
//				x = width;
//			}
//		}
//	}
//	return coordinatesOfCentralPath;
//}
//
//  CImg<int> Tortoise::getCoordinatesOfBestCentralPathOnMediumWideStripe()
//{
//	CImg<int> coordinatesOfBestCentralPathOnMediumWideStripe(_coordinatesOfBestCentralPathOnNarrowStripe);
//	int shiftFromNarrowToMediumWideStripe = (_mediumWideStripeEdges.width() - _narrowStripeEdges.width())/2;
//	for(int y = _lowerYEndOfBestCentralPath; y <= _upperYEndOfBestCentralPath; y++)
//	{
//		int x = _coordinatesOfBestCentralPathOnNarrowStripe(y);
//		//shifting coordinates of the central seam to fit the wider stripe
//		if(x > 0) coordinatesOfBestCentralPathOnMediumWideStripe(y) = x + shiftFromNarrowToMediumWideStripe;
//	}
//	return coordinatesOfBestCentralPathOnMediumWideStripe;
//}
//
//  CImg<int> Tortoise::getCoordinatesOfBestCentralPathOnWideStripe()
//{
//	CImg<int> coordinatesOfBestCentralPathOnWideStripe(_coordinatesOfBestCentralPathOnNarrowStripe);
//	int shiftFromNarrowToWideStripe = (_wideStripeEdges.width() - _narrowStripeEdges.width())/2;
//	for(int y = _lowerYEndOfBestCentralPath; y <= _upperYEndOfBestCentralPath; y++)
//	{
//		int x = _coordinatesOfBestCentralPathOnNarrowStripe(y);
//		//shifting coordinates of the central seam to fit the wider stripe
//		if(x > 0) coordinatesOfBestCentralPathOnWideStripe(y) = x + shiftFromNarrowToWideStripe;
//	}
//	return coordinatesOfBestCentralPathOnWideStripe;
//}
//
////  void Tortoise::drawBestPath()
////{
////  int x = 0;
////  //float koef = (900.0*(sqrt((float)m_imgResultsImg.width()*m_imgResultsImg.height()))/2799)+ 100 - 368*900.0/2799;
////	//int halfOfTheLineWidth = max((int)(koef*0.006),1);
////  for(int y = 0; y < _coordinatesOfBestCentralPathOnWideStripe.height(); ++y)
////  {
////    x = _coordinatesOfBestCentralPathOnWideStripe(0,y);
////    //cout << x << endl;
////    if (x > 0)
////    {
////      int Y = y*m_imgResultsImg.height()/m_nHeightOfRotatedOriginalImg;
////      x = (x+_leftXcoordinateOfStripeInTheRotatedImg)*m_imgResultsImg.width()/m_nWidthOfRotatedOriginalImg;
////      //for(int i = -halfOfTheLineWidth; i <= halfOfTheLineWidth; ++i)
////      {
////        m_imgResultsImg(x,Y,0) = 0;
////        m_imgResultsImg(x,Y,1) = 255;
////			  m_imgResultsImg(x,Y,2) = 0;
////      }
////    }
////  }
////}
//
////------------	 Abdominal Junctions	-------------------------------------------------------------------------------/
//  bool Tortoise::getYCoordinatesOfAbdominalJunctions(CImg<unsigned char> edgeImgToDoItOn) // returns positions of the best side paths on the best central path of Abdominal Junctions
//{ 
////using edges on a wider central-seam stripe
//	CImg<int> stripeEdges(edgeImgToDoItOn);
//
//	int stripeWidth = stripeEdges.width();
//	int stripeHeight = stripeEdges.height();
//  stripeEdges.display();
//
//	CImg<int> valuesInThirdsOfWidth(stripeWidth,stripeHeight,1,2, 0);
//	CImg<int> coordinatesOfBestCentralPathOnStripe(_coordinatesOfBestCentralPathOnWideStripe);
//	//CImg<int> stripeWithDiscountValuesBasedOnDirectionOfContinuousEdges(_wideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges);
//  coordinatesOfBestCentralPathOnStripe.display();
////every edge pix is part of a some path - all edge pixs are assigned an index of the path they belong to 
//	CImg<int> sidePathsIndexes(stripeEdges.get_fill(0));
////pixs on the cetral path have 2-column index array, beacause there meet paths from left and right side
//	CImg<int> indexesOnTheBestCentralPath(2,coordinatesOfBestCentralPathOnStripe.height(),1,1, 0);
////prefilling - side paths to the value 10000, for future easier processing
//	CImg<int> sidePathsFromLeft(stripeEdges.get_fill(10000)), sidePathsFromRight(sidePathsFromLeft);
//
//	for(int y = 0; y < stripeHeight; y++)
//	{
//		int x = coordinatesOfBestCentralPathOnStripe(y);
//		if(x > 0)
////central path is added to the stripeEdge img and if there are edges (on the original edge img), where the central path is, they are set to 2
//		stripeEdges(x,y)++;
//	}
////prefilling - the sides (5 the most left columns) of sidePathsFromLeft and sidePathsIndexes
//	double index = 1;
//	for(int y = _lowerYEndOfBestCentralPath; y < _upperYEndOfBestCentralPath; y++)
//	{
//		index++;
//		for(int x = 0; x <= 4; x++)
//		{
////prefilling - for all edge pixs in the area
//			if(stripeEdges(x,y)) 
//			{ 
////prefilling - initial value (distance) of sidePathsFromLeft is set to 2
//				sidePathsFromLeft(x,y) = 2; 
////prefilling - edge pixs in the area are all assigned an index in sidePathsIndexes, which they will pass on to pixs that will be in the same path 
//				sidePathsIndexes(x,y) = index; 
//			}
//		}
//	}
////prefilling - the sides (5 the most right columns) of sidePathsFromRight and sidePathsIndexes
//	index += 1000;
//	for(int y = _lowerYEndOfBestCentralPath; y < _upperYEndOfBestCentralPath; y++)
//	{
//		index++;
//		for(int x = stripeWidth-1; x >= stripeWidth-5; x--)
//		{
//			if(stripeEdges(x,y)) 
//			{ 
////prefilling - initial value (distance) of sidePathsFromRight is set to 2
//				sidePathsFromRight(x,y) = 2; 
////prefilling - edge pixs in the area are all assigned an index in sidePathsIndexes, which they will pass on to pixs that will be in the same path 
//				sidePathsIndexes(x,y) = index; 
//			}
//		}
//	}
//
////prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
//	CImg<int> distanceMask(7,5,1,1, 0);
////prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
//	for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = 1*pow((y+1),2.0);
//	for(int y = 0; y < 5; y++) { distanceMask(0,y) += 8; distanceMask(6,y) += 8; distanceMask(1,y) += 4; distanceMask(5,y) += 4; }
//	distanceMask+=3;
//	//distanceMask.display();
////filling left side paths - rotating distanceMask
//	distanceMask.rotate(90);
////filling left side paths - going though all img pixs (rows only those where the central path is)
//	for(int x = 5; x <= stripeWidth-1; x++)
//	{
//		for(int y = max(5,_lowerYEndOfBestCentralPath); y < min(_upperYEndOfBestCentralPath,stripeHeight-5); y++)
//		{
////filling left side paths - evaluating only edge pixs on the left of the central path
//			if(stripeEdges(x,y)>0 && x<=coordinatesOfBestCentralPathOnStripe(y)) 
//			{
////filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs 
//				CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
////filling left side paths - minimum distance is saved
//				int min = cropAreaDistances.min();
////filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
//				if(min<10000) 
//				{
////filling left side paths - minimum distance is saved in the current pix
//					sidePathsFromLeft(x,y) = min;
////filling left side paths - index of the pix with minimum distace is passed on to the current pix
//				    sidePathsIndexes(x,y) = sidePathsIndexes(x-5+cropAreaDistances.get_stats()(4),y-3+cropAreaDistances.get_stats()(5));
////filling left side paths - if the pix is on the central path, the index is saved in indexesOnTheBestCentralPath
//				    if(x==coordinatesOfBestCentralPathOnStripe(y)) indexesOnTheBestCentralPath(0,y) = sidePathsIndexes(x-5+cropAreaDistances.get_stats()(4),y-3+cropAreaDistances.get_stats()(5));
//					if(x>stripeWidth/3-5 && x < stripeWidth/3+5) 
//					{ 
//						valuesInThirdsOfWidth(x,y,0) = sidePathsFromLeft(x,y);
//						valuesInThirdsOfWidth(x,y,1) = sidePathsIndexes(x,y);
//					}
//				}
//			}
//		}
//	}
////filling right side paths - rotating distanceMask 
//	//stripeWithDiscountValuesBasedOnDirectionOfContinuousEdges.display("stripeWithDiscountValuesBasedOnDirectionOfContinuousEdges");
//	distanceMask.rotate(180);
////filling right side paths - the same as is in the left side
//	for(int x = stripeWidth-6; x >= 0; x--)
//	{
//		for(int y = max(5,_lowerYEndOfBestCentralPath); y < min(_upperYEndOfBestCentralPath,stripeHeight-5); y++)
//		{
//			if(stripeEdges(x,y)>0 && x >= coordinatesOfBestCentralPathOnStripe(y)) 
//			{
//				CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
//				int min = cropAreaDistances.min();
//				if(min<10000)
//				{
//					sidePathsFromRight(x,y) = min;
//					sidePathsIndexes(x,y) = sidePathsIndexes(x+1+cropAreaDistances.get_stats()(4),y-3+cropAreaDistances.get_stats()(5));
//					if(x==coordinatesOfBestCentralPathOnStripe(y)) indexesOnTheBestCentralPath(1,y) = sidePathsIndexes(x+1+cropAreaDistances.get_stats()(4),y-3+cropAreaDistances.get_stats()(5));
//					if(x>stripeWidth*2/3-5 && x < stripeWidth*2/3+5) 
//					{ 
//						valuesInThirdsOfWidth(x,y,0) = sidePathsFromRight(x,y);
//						valuesInThirdsOfWidth(x,y,1) = sidePathsIndexes(x,y);
//					}
//				}
//			}
//		}
//	}
//	(sidePathsFromLeft + sidePathsFromRight).display();
////postprocessing - set values 10000 to 0 in sidePathsFromLeft and sidePathsFromRight
//	for(int y = 0; y < stripeHeight; y++) 
//	{
//		for(int x = 0; x < stripeWidth; x++)
//		{
//			if(sidePathsFromLeft(x,y) == 10000) sidePathsFromLeft(x,y) = 0;
//			if(sidePathsFromRight(x,y)== 10000) sidePathsFromRight(x,y) = 0;
//		}
//	}
//	(sidePathsFromLeft+sidePathsFromRight).display("sidePaths");
////postprocessing - setting upperSupremum that will be greater than 4x maximum of maxima in sidePathsFromLeft and sidePathsFromRight
//	int upperSupremum = max(sidePathsFromLeft.max(),sidePathsFromRight.max())*4+1;
////postprocessing - CentralSideCrosses will contain sum of 2x values of distances of pixs on the central seam of the left and right paths
//	CImg<int> CentralSideCrosses(1,coordinatesOfBestCentralPathOnStripe.height(),1,1, upperSupremum);
////postprocessing - for all central path pixs near the tortoise center
//		for(int y = max(max(_lowerYEndOfBestCentralPath,m_sPlastronCentreOnRotatedImg.yCoor - (int)(stripeHeight*0.33)),(int)(stripeHeight*0.20)) ; y < min(min(_upperYEndOfBestCentralPath,m_sPlastronCentreOnRotatedImg.yCoor+(int)(stripeHeight*0.33)), stripeHeight - (int)(stripeHeight*0.20)); y++)
//	{
//		int x = coordinatesOfBestCentralPathOnStripe(y);
////postprocessing - if there are paths from left and right that end in the central-path pix
//		if(sidePathsFromLeft(x,y)>0 && sidePathsFromRight(x,y)>0)
//		{
////postprocessing - multiplying twice the value of central path pixs
//			sidePathsFromLeft(x,y) *= 2;
//			sidePathsFromRight(x,y) *= 2;
////postprocessing - the values are added to CentralSideCrosses
//			CentralSideCrosses(y) = sidePathsFromLeft(x,y) + sidePathsFromRight(x,y);
//		}
////postprocessing - else, left and right indexes of the central-path pix are set to 0
//		else 
//		{
//			indexesOnTheBestCentralPath(0,y) = 0;
//			indexesOnTheBestCentralPath(1,y) = 0;
//		}
//	}
//    CentralSideCrosses.display();
//	CImg<int> CentralSideCrossesCopy(CentralSideCrosses),indexesOnTheBestCentralPathCopy(indexesOnTheBestCentralPath);
////finding best side paths - setting number of pixs that are to be discarded, being next (based on index) to a pix where we will find the local minimum (of distances in CentralSideCrosses) 
////                        - with 3000pixs height img time constant 0.025 (based on real images) so that we get ~ 70pixs 
////                        - (it has to be so much that we eliminate all neighbouring paths with similar distances, but so that we do not eliminate other correct side paths)
//	int halfIntervalOfIndexes = stripeHeight*0.025;
////finding best side paths - countOfCrossing counts side paths and is used to index the paths
//	int countOfCrossing = 0;
////finding best side paths - we want no more than 9 paths (5-side paths, 2-end and beggining of the plastron, 2-ends of tail and head/ends of the shell)
//	int permittedNumberOfPaths = 6;
//	int totalMin = 0;
//	int totalMax = 0;
//	CImg<int> yCoordinateOfTheJunctions(permittedNumberOfPaths,1,1,1, 0);
//	CImg<int> valuesOfTheJunctions(permittedNumberOfPaths,1,1,1, 0);
////finding best side paths - going though CentralSideCrosses until there is no pix with smaller distance than upperSupremum or the is permittedNumberOfPaths reached
//	distanceMask.rotate(180);
//	while(CentralSideCrosses.min()<upperSupremum && countOfCrossing < permittedNumberOfPaths)
//	{
////finding best side paths - finding the current minimum in CentralSideCrosses (distances of left and right paths)
//		int yMin = CentralSideCrosses.get_stats()(5);
//		int currentMin =  CentralSideCrosses.min();
//		yCoordinateOfTheJunctions(countOfCrossing) = yMin;
//		valuesOfTheJunctions(countOfCrossing) = currentMin;
//		//if(countOfCrossing == 1) totalMin = currentMin;
////finding best side paths - getting left index of the pix with the minimum and its y coordinate
//		int leftIndexOfMinPix = indexesOnTheBestCentralPath(0,yMin);
//		int yOfLeftIndexOfMinPix = 0;
//		for(int j = 0; j < stripeHeight; j++) if(sidePathsIndexes(4,j)==leftIndexOfMinPix) yOfLeftIndexOfMinPix = j;
////finding best side paths - getting right index of the pix with the minimum and its y coordinate
//		int rightIndexOfMinPix = indexesOnTheBestCentralPath(1,yMin);
//		int yOfRightIndexOfMinPix = 0;
//		for(int j = 0; j < stripeHeight; j++) if(sidePathsIndexes(stripeWidth-5,j)==rightIndexOfMinPix) yOfRightIndexOfMinPix = j;
////finding best side paths - getting the correct number of indexes that can be deleted on the left of the left side of the index of the minimum pix
//		int localHalfIntervalOfIndexesLeftLeft = 0;
////finding best side paths - counting edge pixs in the interval and dividing it by lenght (in pixs) of that interval, thus getting edge pixs frequency, 
////                        - multiplied by num of pixs to be discarded we get the correct number of indexes to be deleted in CentralSideCrosses while not deleting those belonging to the next correct side path
//		for(int j = yOfLeftIndexOfMinPix; j < min(stripeHeight,yOfLeftIndexOfMinPix+halfIntervalOfIndexes); j++) 
//			if(stripeEdges(0,j)||stripeEdges(1,j)||stripeEdges(2,j)||stripeEdges(3,j)||stripeEdges(4,j)) localHalfIntervalOfIndexesLeftLeft++;
////finding best side paths - the same for the other m_nIntervals
//		int localHalfIntervalOfIndexesLeftRight = 0;
//		for(int j = max(0,yOfLeftIndexOfMinPix-halfIntervalOfIndexes); j < yOfLeftIndexOfMinPix; j++) 
//			if(stripeEdges(0,j)||stripeEdges(1,j)||stripeEdges(2,j)||stripeEdges(3,j)||stripeEdges(4,j)) localHalfIntervalOfIndexesLeftRight++;
//
//		int localHalfIntervalOfIndexesRightLeft = 0;
//		for(int j = max(0,yOfRightIndexOfMinPix-halfIntervalOfIndexes); j < yOfRightIndexOfMinPix; j++) 
//			if(stripeEdges(stripeWidth-1,j)||stripeEdges(stripeWidth-2,j)||stripeEdges(stripeWidth-3,j)||stripeEdges(stripeWidth-4,j)||stripeEdges(stripeWidth-5,j)) localHalfIntervalOfIndexesRightLeft++;
//
//		int localHalfIntervalOfIndexesRightRight = 0;
//		for(int j = yOfRightIndexOfMinPix; j < min(stripeHeight,yOfRightIndexOfMinPix+halfIntervalOfIndexes); j++) 
//			if(stripeEdges(stripeWidth-1,j)||stripeEdges(stripeWidth-2,j)||stripeEdges(stripeWidth-3,j)||stripeEdges(stripeWidth-4,j)||stripeEdges(stripeWidth-5,j)) localHalfIntervalOfIndexesRightRight++;
////finding best side paths - deleting indexes close to index of the pix with the current min value
//		for(int y = _lowerYEndOfBestCentralPath; y < _upperYEndOfBestCentralPath; y++)
//		{
//			if( (indexesOnTheBestCentralPath(0,y) <= leftIndexOfMinPix + localHalfIntervalOfIndexesLeftLeft) &&
//				(indexesOnTheBestCentralPath(0,y) >= leftIndexOfMinPix - localHalfIntervalOfIndexesLeftRight))
//			{
////finding best side paths - index to be deleted from indexesOnTheBestCentralPath
//				indexesOnTheBestCentralPath(0,y) = 0;
////finding best side paths - coresponding distances in CentralSideCrosses set to upperSupremum (cannot become minimum in the future)
//				CentralSideCrosses(y) = upperSupremum;
//			}
////finding best side paths - the same for the right side indexes
//			else if( indexesOnTheBestCentralPath(1,y) <= rightIndexOfMinPix + localHalfIntervalOfIndexesRightRight &&
//				     indexesOnTheBestCentralPath(1,y) >= rightIndexOfMinPix - localHalfIntervalOfIndexesRightLeft )
//			{
//				indexesOnTheBestCentralPath(1,y) = 0;
//				CentralSideCrosses(y) = upperSupremum;
//			}
//		}
////finding best side paths - numbering side paths 
//		CentralSideCrosses(yMin) = upperSupremum + currentMin + 1;
////finding best side paths - rising the count of side paths
//		countOfCrossing++;
//	}
//	int numberOfJunctions = countOfCrossing;
//	CentralSideCrosses -= upperSupremum;
//	yCoordinateOfTheJunctions.crop(0,0,numberOfJunctions-1,0);
//	valuesOfTheJunctions.crop(0,0,numberOfJunctions-1,0);
//	CImg<int> distanceOfTheJunctionsToTheCenterPointOnTheSeam((yCoordinateOfTheJunctions - m_sPlastronCentreOnRotatedImg.yCoor).abs());
//	int middleValueOfTheJunctions =  valuesOfTheJunctions(numberOfJunctions-1,0) + 1;
//	bool isAlreadyUp = 0;
//	bool isAlreadyDown = 0;
//	CImg<int> stripeWithAbdJunc(_mediumWideStripeEdges);
//	cout << numberOfJunctions << endl;
//	for(int i = 0; i < numberOfJunctions; i++) 
//	{
//		cout << i << endl;
//		int currentYCoordinateOfMinJunction = distanceOfTheJunctionsToTheCenterPointOnTheSeam.get_stats()(4);
//		int y = yCoordinateOfTheJunctions(currentYCoordinateOfMinJunction);
//		distanceOfTheJunctionsToTheCenterPointOnTheSeam(currentYCoordinateOfMinJunction,0) = stripeHeight;
//		int valueOfTheJunction = CentralSideCrossesCopy(y);
//		int x = _coordinatesOfBestCentralPathOnWideStripe(y);
//		//if(isAlreadyDown && isAlreadyUp) 
//		//{
//		//	for(int i=-5;i<=5;i++) for(int j=-5;j<=5;j++) 
//		//	{
//		//		stripeShow2(x+i,y+j,0) = 255;
//		//		stripeShow2(x+i,y+j,1) = 0;
//		//		stripeShow2(x+i,y+j,2) = 0;
//		//	}
//		//}
//		if( ( y <= m_sPlastronCentreOnRotatedImg.yCoor && !isAlreadyDown ) || ( y > m_sPlastronCentreOnRotatedImg.yCoor && !isAlreadyUp ) )
//		{
//				//for(int i=-10;i<=10;i++) for(int j=-10;j<=10;j++) 
//				//{
//				//	stripeWithAbdJunc(x+i,y+j,0) = 0;
//				//	stripeWithAbdJunc(x+i,y+j,1) = 255;
//				//	stripeWithAbdJunc(x+i,y+j,2) = 0;
//				//}
//				//stripeWithAbdJunc.display();
//			int numUp = 0;
//			int numDown = 0;
//			int extension = 0;
//			if(y <= m_sPlastronCentreOnRotatedImg.yCoor)
//			{
//				for(int i = 0; i < halfIntervalOfIndexes*2/3 + extension; i++)
//				{
//					if(CentralSideCrossesCopy(y - i) < CentralSideCrossesCopy(y - i - 1)) numUp++;
//					else if(CentralSideCrossesCopy(y - i) > CentralSideCrossesCopy(y - i - 1)) numDown++;
//					else extension++;
//				}
//			}
//			else 
//			{
//				for(int i = 0; i < halfIntervalOfIndexes*2/3 + extension; i++)
//				{
//					if(CentralSideCrossesCopy(y + i) < CentralSideCrossesCopy(y + i + 1)) numUp++;
//					else if(CentralSideCrossesCopy(y + i) > CentralSideCrossesCopy(y + i + 1)) numDown++;
//					else extension++;
//				}
//			}
//			if(valueOfTheJunction < middleValueOfTheJunctions && 2*numUp >= 3*numDown)
//			{
//
//					if(y <= m_sPlastronCentreOnRotatedImg.yCoor) { isAlreadyDown = 1; _abdominalJunctionsYDownCoordinate = y; } 
//					else { isAlreadyUp = 1; _abdominalJunctionsYUpCoordinate = y; }
//			}
//		}
//	}
//	if(isAlreadyUp==0 || isAlreadyDown==0)
//  {
//    cout << "nenaslo abdminal junction" << endl;
//    return 0;
//  }
//  else return 1;
//}
////
////  void Tortoise::drawAbdominalJunctions()
////{
////  float koef = (900.0*(sqrt((float)m_imgResultsImg.width()*m_imgResultsImg.height()))/2799)+ 100 - 368*900.0/2799;
////	int halfOfTheLineWidth = max((int)(koef*0.012),1);
////  int y = _abdominalJunctionsYDownCoordinate;
////  int x = _coordinatesOfBestCentralPathOnWideStripe(0,y);
////  if(x>0)
////  {
////    x = (x+_leftXcoordinateOfStripeInTheRotatedImg)*m_imgResultsImg.width()/m_nWidthOfRotatedOriginalImg;
////    y = y*m_imgResultsImg.height()/m_nHeightOfRotatedOriginalImg;
////    for(int j = -halfOfTheLineWidth; j <= halfOfTheLineWidth; ++j)
////    {
////      for(int i = -halfOfTheLineWidth; i <= halfOfTheLineWidth; ++i)
////      {
////        m_imgResultsImg(x+i,y+j,0) = 0;
////        m_imgResultsImg(x+i,y+j,1) = 0;
////		    m_imgResultsImg(x+i,y+j,2) = 255;
////      }
////    }
////  }
////  y = _abdominalJunctionsYUpCoordinate;
////  x = _coordinatesOfBestCentralPathOnWideStripe(0,y);
////  if(x>0)
////  {
////    x = (x+_leftXcoordinateOfStripeInTheRotatedImg)*m_imgResultsImg.width()/m_nWidthOfRotatedOriginalImg;
////    y = y*m_imgResultsImg.height()/m_nHeightOfRotatedOriginalImg;
////    for(int j = -halfOfTheLineWidth; j <= halfOfTheLineWidth; ++j)
////    {
////      for(int i = -halfOfTheLineWidth; i <= halfOfTheLineWidth; ++i)
////      {
////        m_imgResultsImg(x+i,y+j,0) = 0;
////        m_imgResultsImg(x+i,y+j,1) = 0;
////		    m_imgResultsImg(x+i,y+j,2) = 255;
////      }
////    }
////  }
////}
//
//  CImg<double> Tortoise::getCentralAndSideSeamPaths()
//{
//	correctionOfCentralPath(_skeletonizedEdges);
//
//	_coordinatesOfBestCentralPathOnNarrowStripe = getCoordinatesOfBestCentralPathOnNarrowStripe();
//	getProlongedNarrowStripeBestPathCoordinates();
//	_coordinatesOfBestCentralPathOnMediumWideStripe = getCoordinatesOfBestCentralPathOnMediumWideStripe();
//	_coordinatesOfBestCentralPathOnWideStripe = getCoordinatesOfBestCentralPathOnWideStripe();
//	cout << "getYCoordinatesOfPectoralFemoralJunctions" << endl;
//
//	getYCoordinatesOfPectoralFemoralJunctions(_skeletonizedEdgesHorizontalMediumWideStripe);
//	cout << "getAbdoPectoFemoSidePaths" << endl;
//
//	drawJunctions();
//	//_stripeWithJunctions.display();
//	_segmentsLengthsAbdPectFemo = getAbdoPectoFemoSidePaths();
//	cout << "getGroundTruthSegmentsLengths" << endl;
//
//	_segmentsLengthsAbdPectFemoAllTortoises = getGroundTruthSegmentsLengths();
//	cout << "classifyTortoise" << endl;
//
//	classifyTortoise();
//	cout << "visualizeCentralAndSidePaths" << endl;
//
//	visualizeCentralAndSidePaths();
//	//_stripeWithJunctions.display();
//	//_joiningPointsOfsidePathsToTheCentralPath = getJoiningPointsOfsidePathsToTheCentralPath(m_sCentralSeamAreaWidthCoeff.fWide);
//	
//	return m_imgOriginal;
//}
// 
//
//CImg<unsigned char> Tortoise::getSkeletonizationMinDistWideStripe(CImg<unsigned char> imgToDoItOn)
//{	
//	imgToDoItOn = EliminateTwoPixelsEdges(EliminateOnePixelEdges(imgToDoItOn));
//	CImg<int> res(imgToDoItOn),edgesWrapping(_edgesWrappingWideStripe);
//	int width = res.width();
//	int height = res.height();
////sets value 0 for edges in resulting image
//	res = (res-1).abs();
////adds edgesWrapping, so 0s are only in the working area where there is no edge
//	imgToDoItOn = imgToDoItOn + (edgesWrapping-1).abs();
//	//imgToDoItOn.display("imgToDoItOn");
////
//	for(int y = 0; y < height-1; y++)
//	{
//		for(int x = 0; x < width-1; x++)
//		{
//			if(imgToDoItOn(x,y) == 0)
//			{
//				int minDist = 1;
//				bool doIncreaseDistance = 1;
//				while(doIncreaseDistance && x-minDist>=0 && y-minDist>=0 && x+minDist<width && y+minDist<height)
//				{
//					if(imgToDoItOn.get_crop(x-minDist,y-minDist,x+minDist,y+minDist).sum() > 0 || minDist >= 3)
//					{
//						doIncreaseDistance = 0;
//					}
//					else minDist++;
//				}
//				res(x,y) = minDist;
//			}
//		}
//	}
//	//res.display();
//	res = res.get_mul(edgesWrapping);
//	//res.display();
////eliminates pixs, first those close to edges (cycle 1->3)
//	int maxLevel = 3;
//	for(int cycle = 1; cycle <= maxLevel; cycle++)
//	{
////in each cycle it eliminates also all newly-suitable pixs on lower levels (level 1...cycle)
//		for(int level = 1; level <= cycle; level++)
//		{
////eliminates until no elimination has occured in previous go-through of the image
//			bool change = 1;
//			while(change)
//			{
//				change = 0;
//				for(int y = 1; y < height - 1; y++)
//				{
//					for(int x = 1; x < width-1; x++)
//					{
//						if(res(x,y) == level) 
//						{	
//// eliminate from rigth down corner
//							if( res(x-1,y-1)> 0 && res(x,y-1)> 0 && res(x+1,y-1)>=0  &&
//								res(x-1,y)  > 0 &&				    res(x+1,y)  ==0  &&
//								res(x-1,y+1)>=0 && res(x,y+1)==0 && res(x+1,y+1)>=0     ) { res(x,y) = 0; change = 1; }
//// eliminate from left down corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)> 0 && res(x+1,y-1)> 0 &&
//									 res(x-1,y)  ==0 &&				    res(x+1,y)  > 0 &&
//									 res(x-1,y+1)>=0 && res(x,y+1)==0 && res(x+1,y+1)>=0  ) { res(x,y) = 0; change = 1; }
//// eliminate from left up corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)==0 && res(x+1,y-1)>=0 &&
//									 res(x-1,y)  ==0 &&				    res(x+1,y)  > 0 &&
//									 res(x-1,y+1)>=0 && res(x,y+1)> 0 && res(x+1,y+1)> 0   ) { res(x,y) = 0; change = 1; }
//// eliminate from right up corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)==0 && res(x+1,y-1)>=0 &&
//									 res(x-1,y)  > 0 &&					res(x+1,y)  ==0 &&
//									 res(x-1,y+1)> 0 && res(x,y+1)> 0 && res(x+1,y+1)>=0   ) { res(x,y) = 0; change = 1; }
//// eliminate from right side
//							else if( res(x-1,y-1)>0 &&  (res(x+1,y-1)==0||res(x+1,y+1)==0) &&
//									 res(x-1,y)  >0 &&				       res(x+1,y)  ==0 &&
//									 res(x-1,y+1)>0                                           ) { res(x,y) = 0; change = 1; }
//// eliminate from left side
//							else if( (res(x-1,y-1)==0||res(x-1,y+1)==0) &&  res(x+1,y-1) > 0 &&
//									  res(x-1,y)  ==0  &&                   res(x+1,y)   > 0 &&
//																			res(x+1,y+1) > 0    ) { res(x,y) = 0; change = 1; }
//// eliminate from up side
//							else if( (res(x-1,y-1)==0 || res(x+1,y-1)==0) && res(x,y-1)==0 && 
//									 (res(x-1,y)  >0 ||				   	 res(x+1,y)  >0) &&
//									  res(x-1,y+1)>0 && res(x,y+1)>0 && res(x+1,y+1)>0  ) { res(x,y) = 0; change = 1; }
//// eliminate from down side
//							else if( res(x-1,y-1)>0 && res(x,y-1)>0 && res(x+1,y-1)>0 &&
//									(res(x-1,y)  >0 ||					 res(x+1,y)>0) &&
//									(res(x-1,y+1)==0 || res(x+1,y+1)==0) && res(x,y+1)==0   ) { res(x,y) = 0; change = 1; }
//							else {	}
//						}
//					}
//				}
//			}
//		}
//	}
//	//res.display();
//	return res.threshold(1);
//}
//
//CImg<unsigned char> Tortoise::getHorizontalSkeletonizationMinDist(CImg<unsigned char> imgToDoItOn)
//{	
//	imgToDoItOn = EliminateTwoPixelsEdges(EliminateOnePixelEdges(imgToDoItOn));
//	CImg<int> res(imgToDoItOn),edgesWrapping(_edgesWrappingWideStripe);
//	int width = res.width();
//	int height = res.height();
////sets value 0 for edges in resulting image
//	res = (res-1).abs();
////adds edgesWrapping, so 0s are only in the working area where there is no edge
//	imgToDoItOn = imgToDoItOn + (edgesWrapping-1).abs();
//	//imgToDoItOn.display("imgToDoItOn");
////
//	for(int y = 0; y < height-1; y++)
//	{
//		for(int x = 0; x < width-1; x++)
//		{
//			if(imgToDoItOn(x,y) == 0)
//			{
//				int minDist = 1;
//				bool doIncreaseDistance = 1;
//				while(doIncreaseDistance && x-minDist>=0 && y-minDist>=0 && x+minDist<width && y+minDist<height)
//				{
//					if(imgToDoItOn.get_crop(x-minDist,y-minDist,x+minDist,y+minDist).sum() > 0 || minDist >= 3)
//					{
//						doIncreaseDistance = 0;
//					}
//					else minDist++;
//				}
//				res(x,y) = minDist;
//			}
//		}
//	}
//	//res.display();
//	res = res.get_mul(edgesWrapping);
//	//res.display();
////eliminates pixs, first those close to edges (cycle 1->3)
//	int maxLevel = 3;
//	for(int cycle = 1; cycle <= maxLevel; cycle++)
//	{
////in each cycle it eliminates also all newly-suitable pixs on lower levels (level 1...cycle)
//		for(int level = 1; level <= cycle; level++)
//		{
////eliminates until no elimination has occured in previous go-through of the image
//			bool change = 1;
//			while(change)
//			{
//				change = 0;
//				for(int y = 1; y < height - 1; y++)
//				{
//					for(int x = 1; x < width-1; x++)
//					{
//						if(res(x,y) == level) 
//						{	
//// eliminate from rigth down corner
//							if( res(x-1,y-1)> 0 && res(x,y-1)> 0 && res(x+1,y-1)>=0  &&
//								res(x-1,y)  > 0 &&				    res(x+1,y)  ==0  &&
//								res(x-1,y+1)>=0 && res(x,y+1)==0 && res(x+1,y+1)>=0     ) { res(x,y) = 0; change = 1; }
//// eliminate from left down corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)> 0 && res(x+1,y-1)> 0 &&
//									 res(x-1,y)  ==0 &&				    res(x+1,y)  > 0 &&
//									 res(x-1,y+1)>=0 && res(x,y+1)==0 && res(x+1,y+1)>=0  ) { res(x,y) = 0; change = 1; }
//// eliminate from left up corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)==0 && res(x+1,y-1)>=0 &&
//									 res(x-1,y)  ==0 &&				    res(x+1,y)  > 0 &&
//									 res(x-1,y+1)>=0 && res(x,y+1)> 0 && res(x+1,y+1)> 0   ) { res(x,y) = 0; change = 1; }
//// eliminate from right up corner
//							else if( res(x-1,y-1)>=0 && res(x,y-1)==0 && res(x+1,y-1)>=0 &&
//									 res(x-1,y)  > 0 &&					res(x+1,y)  ==0 &&
//									 res(x-1,y+1)> 0 && res(x,y+1)> 0 && res(x+1,y+1)>=0   ) { res(x,y) = 0; change = 1; }
//// eliminate from right side
//							else if( res(x-1,y-1)>0 &&  (res(x+1,y-1)==0||res(x+1,y+1)==0) &&
//									 res(x-1,y)  >0 &&				       res(x+1,y)  ==0 &&
//									 res(x-1,y+1)>0  && res(x,y-1) > 0 && res(x,y+1) > 0  ) { res(x,y) = 0; change = 1; }
//// eliminate from left side
//							else if( (res(x-1,y-1)==0||res(x-1,y+1)==0) &&  res(x+1,y-1) > 0 &&
//									  res(x-1,y)  ==0  &&                   res(x+1,y)   > 0 &&
//									  res(x,y-1) > 0 && res(x,y+1) > 0 && res(x+1,y+1) > 0    ) { res(x,y) = 0; change = 1; }
//// eliminate from up side
//							else if( (res(x-1,y-1)==0 || res(x+1,y-1)==0) && res(x,y-1)==0 && 
//									
//									  res(x-1,y+1)>0 && res(x,y+1)>0 && res(x+1,y+1)>0  ) { res(x,y) = 0; change = 1; }
//// eliminate from down side
//							else if( res(x-1,y-1)>0 && res(x,y-1)>0 && res(x+1,y-1)>0 &&
//									
//									(res(x-1,y+1)==0 || res(x+1,y+1)==0) && res(x,y+1)==0   ) { res(x,y) = 0; change = 1; }
//							else {	}
//						}
//					}
//				}
//			}
//		}
//	}
//	//res.display();
//	return res.threshold(1);
//}
//
//// direction are 1-vertical, 2-diagonalRightUp, 3-horizontal, 4-diagonalLeftUp
//CImg<char> Tortoise::getWideStripeWithDirectionOfEdges()
//{
//	CImg<char> wideStripeWithDirectionOfEdges(_wideStripeEdges);
//	int wideStripeWidth = _wideStripeEdges.width();
//	for(int x = 2; x < wideStripeWidth-2; x++)
//	{
//		for(int y = 2; y < _stripeHeight-2; y++)
//		{
//			if(_wideStripeEdges(x,y))
//			{
//				int dirVertical =  _wideStripeEdges(x-1,y-2) + _wideStripeEdges(x,y-2) + _wideStripeEdges(x+1,y-2) 
//									   + _wideStripeEdges(x-1,y-1) + _wideStripeEdges(x,y-1) + _wideStripeEdges(x+1,y-1) 
//									   + _wideStripeEdges(x-1,y+1) + _wideStripeEdges(x,y+1) + _wideStripeEdges(x+1,y+1) 
//									   + _wideStripeEdges(x-1,y+2) + _wideStripeEdges(x,y+2) + _wideStripeEdges(x+1,y+2)
//									   + _wideStripeEdges(x+2,y-2) + _wideStripeEdges(x+2,y+2) + _wideStripeEdges(x-2,y-2) 
//									   + _wideStripeEdges(x-2,y+2);
//				int dirHorizontal =   _wideStripeEdges(x-2,y-1) + _wideStripeEdges(x-2,y) + _wideStripeEdges(x-2,y+1) 
//										  + _wideStripeEdges(x-1,y-1) + _wideStripeEdges(x-1,y) + _wideStripeEdges(x-1,y+1)
//										  + _wideStripeEdges(x+1,y-1) + _wideStripeEdges(x+1,y) + _wideStripeEdges(x+1,y+1) 
//									      + _wideStripeEdges(x+2,y-1) + _wideStripeEdges(x+2,y) + _wideStripeEdges(x+2,y+1)
//										  + _wideStripeEdges(x+2,y-2) + _wideStripeEdges(x+2,y+2) + _wideStripeEdges(x-2,y-2) 
//									      + _wideStripeEdges(x-2,y+2);
//				if (dirHorizontal >= dirVertical) wideStripeWithDirectionOfEdges(x,y) = 2;
//				else wideStripeWithDirectionOfEdges(x,y) = 1;
////				int dirVertical =  _wideStripeEdges(x-1,y-2) + _wideStripeEdges(x,y-2) + _wideStripeEdges(x+1,y-2) 
////									   + _wideStripeEdges(x-1,y-1) + _wideStripeEdges(x,y-1) + _wideStripeEdges(x+1,y-1) 
////									   + _wideStripeEdges(x-1,y+1) + _wideStripeEdges(x,y+1) + _wideStripeEdges(x+1,y+1) 
////									   + _wideStripeEdges(x-1,y+2) + _wideStripeEdges(x,y+2) + _wideStripeEdges(x+1,y+2);
////				int dirHorizontal =   _wideStripeEdges(x-2,y-1) + _wideStripeEdges(x-2,y) + _wideStripeEdges(x-2,y+1) 
////										  + _wideStripeEdges(x-1,y-1) + _wideStripeEdges(x-1,y) + _wideStripeEdges(x-1,y+1)
////										  + _wideStripeEdges(x+1,y-1) + _wideStripeEdges(x+1,y) + _wideStripeEdges(x+1,y+1) 
////									      + _wideStripeEdges(x+2,y-1) + _wideStripeEdges(x+2,y) + _wideStripeEdges(x+2,y+1);
////				int dirDiagonalRightUp =  _wideStripeEdges(x-2,y+1) + _wideStripeEdges(x-2,y+2) + _wideStripeEdges(x-1,y) 
////											  + _wideStripeEdges(x-1,y+1) + _wideStripeEdges(x-1,y+2) + _wideStripeEdges(x,y-1)
////											  + _wideStripeEdges(x,y+1) + _wideStripeEdges(x+1,y-2) + _wideStripeEdges(x+1,y-1) 
////											  + _wideStripeEdges(x+1,y) + _wideStripeEdges(x+2,y-2) + _wideStripeEdges(x+2,y-1);
////				int dirDiagonalLeftUp =  _wideStripeEdges(x-2,y-2) + _wideStripeEdges(x-2,y-1) + _wideStripeEdges(x-1,y-2) 
////									         + _wideStripeEdges(x-1,y-1) + _wideStripeEdges(x-1,y) + _wideStripeEdges(x,y-1) 
////											 + _wideStripeEdges(x,y+1) + _wideStripeEdges(x+1,y) + _wideStripeEdges(x+1,y+1) 
////									         + _wideStripeEdges(x+1,y+2) + _wideStripeEdges(x+2,y+1) + _wideStripeEdges(x+2,y+2);
////				if (dirHorizontal >= dirVertical && dirHorizontal >= dirDiagonalRightUp && dirHorizontal >= dirDiagonalLeftUp) wideStripeWithDirectionOfEdges(x,y) = 3;
////				else if (dirDiagonalRightUp >= dirHorizontal && dirDiagonalRightUp >= dirVertical && dirDiagonalRightUp >= dirDiagonalLeftUp) wideStripeWithDirectionOfEdges(x,y) = 2;
////				else if (dirDiagonalLeftUp >= dirHorizontal && dirDiagonalLeftUp >= dirDiagonalRightUp && dirDiagonalLeftUp >= dirVertical) wideStripeWithDirectionOfEdges(x,y) = 4;
////				else if(dirVertical >= dirHorizontal && dirVertical >= dirDiagonalRightUp && dirVertical >= dirDiagonalLeftUp) wideStripeWithDirectionOfEdges(x,y) = 1;			
//			}
//		}
//	}
//	return wideStripeWithDirectionOfEdges;
//}
//
////*/
//
//
//CImg<char> Tortoise::getWideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges()
//{
//	CImg<char> wideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges(_wideStripeEdges.get_fill(0));
//	int wideStripeWidth = _wideStripeEdges.width();
//	for(int x = 2; x < wideStripeWidth-2; x++)
//	{
//		for(int y = 2; y < _stripeHeight-2; y++)
//		{
//			if(_wideStripeWithDirectionOfEdges(x,y)>1)
//			{
//				wideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges(x,y) = -1;
//				//if(_wideStripeWithDirectionOfEdges == 3 || _wideStripeWithDirectionOfEdges == 2 || _wideStripeWithDirectionOfEdges == 4) wideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges(x,y) = -1;
//				int countOfNeighbouringPixelsWithTheSameDirection = 0;
//				int directionInThePixel = _wideStripeWithDirectionOfEdges(x,y);
//				if(_wideStripeWithDirectionOfEdges(x-1,y-1) == directionInThePixel) countOfNeighbouringPixelsWithTheSameDirection++;
//				if(_wideStripeWithDirectionOfEdges(x-1,y) == directionInThePixel) countOfNeighbouringPixelsWithTheSameDirection++;
//				if(_wideStripeWithDirectionOfEdges(x-1,y+1) == directionInThePixel) countOfNeighbouringPixelsWithTheSameDirection++;
//				//if(_wideStripeWithDirectionOfEdges(x,y-1) == directionInThePixel) countOfNeighbouringPixelsWithTheSameDirection++;
//				//if(_wideStripeWithDirectionOfEdges(x,y+1) == directionInThePixel) countOfNeighbouringPixelsWithTheSameDirection++;
//				if(_wideStripeWithDirectionOfEdges(x+1,y-1) == directionInThePixel) countOfNeighbouringPixelsWithTheSameDirection++;
//				if(_wideStripeWithDirectionOfEdges(x+1,y) == directionInThePixel) countOfNeighbouringPixelsWithTheSameDirection++;
//				if(_wideStripeWithDirectionOfEdges(x+1,y+1) == directionInThePixel) countOfNeighbouringPixelsWithTheSameDirection++;
//
//				if(countOfNeighbouringPixelsWithTheSameDirection >= 2) wideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges(x,y) -= 3;
//			}
//			else if(_wideStripeWithDirectionOfEdges(x,y)==1) wideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges(x,y) += 4;
//		}
//	}
//	return wideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges;
//}
//
//CImg<int> Tortoise::correctionOfCentralPath(CImg<unsigned char> edgeImgToDoItOn)
//{
////decrLevel ... number of lines that we look backward(forward) to find more (suitable) edge points 
//	const int decrLevel = 10;
//	CImg<int> distanceMask(5,decrLevel,1,1, 0),paths(edgeImgToDoItOn),centralPath(paths.get_fill(1000000));
//// paths contains skeleton of edge stripe image without large empty spaces
//	int stripeWidth = paths.width();
//	int stripeHeight = paths.height();
////distanceMask is filled by distances
//	for(int x = 0; x < 5; x++) for(int y = 0; y < decrLevel; y++) distanceMask(x,y) = (x-2)*(x-2) + (y+1)*(y+1);
//	distanceMask.rotate(180);
////starting points
//	for(int x = 0; x < stripeWidth; x++)
//	{
//		for(int i = 0; i < 5; i++)
//		{
//			if(paths(x,_pectoralFemoralJunctionsYDownCoordinate-i)==1) centralPath(x,_pectoralFemoralJunctionsYDownCoordinate-i) = 2;
//		}
//	}
//	int Ystart = _pectoralFemoralJunctionsYDownCoordinate;
//	int Yend = _pectoralFemoralJunctionsYUpCoordinate;
//	for(int y = Ystart; y <= Yend; y++)
//	{
//		for(int x = 2; x < stripeWidth-2; x++)
//		{
//			if(paths(x,y)==1)
//			{
//				CImg<int> cropArea;
//				cropArea = centralPath.get_crop(x-2,y-decrLevel,x+2,y-1);
//				if(cropArea.min() < 1000000)
//				{
//					cropArea = cropArea + distanceMask;
//					centralPath(x,y) = cropArea.min();
//				}
//			}
//		}
//	}
//	int X = centralPath.get_row(_pectoralFemoralJunctionsYUpCoordinate).get_stats()(4);
//	int Y = _pectoralFemoralJunctionsYUpCoordinate;
//	_bestCentralPath.fill(0);
//	while(Y > _pectoralFemoralJunctionsYDownCoordinate)
//	{
//		_bestCentralPath(X,Y) = 1;
//		CImg<int> cropArea(centralPath.get_crop(X-2,Y-decrLevel,X+2,Y-1) + distanceMask);
//		X = cropArea.get_stats()(4) + X-2;
//		Y = cropArea.get_stats()(5) + Y-decrLevel;
//	}
//return _bestCentralPath;
//}
//
//void Tortoise::getProlongedNarrowStripeBestPathCoordinates()
//{
//	//CImg<double> linearRegressionX(1,_upperYEndOfBestCentralPath -_lowerYEndOfBestCentralPath+1,1,1, 1);
//	//linearRegressionX.append(_coordinatesOfBestCentralPathOnNarrowStripe.get_crop(0,_lowerYEndOfBestCentralPath,0,_upperYEndOfBestCentralPath),'x');
//	////linearRegressionX.display();
//	//CImg<double> linearRegressionY(1,_upperYEndOfBestCentralPath -_lowerYEndOfBestCentralPath+1,1,1, 0);
//	//for(int i = 0; i < _upperYEndOfBestCentralPath -_lowerYEndOfBestCentralPath+1; i++)
//	//{
//	//	linearRegressionY(i) = _lowerYEndOfBestCentralPath + i;
//	//}
//	////(linearRegressionX,linearRegressionY).display();
//	//CImg<double> regressionCoeficients(((linearRegressionX.get_transpose()*linearRegressionX).invert()*linearRegressionX.get_transpose())*linearRegressionY );
//	////regressionCoeficients.display("regressionCoeficients");
//	//double b = regressionCoeficients(0,0);
//	//double a = regressionCoeficients(0,1);
//	for(int i = _endOfThePlastronYDownCoordinate; i < _lowerYEndOfBestCentralPath; i++)
//	{
//		//_coordinatesOfBestCentralPathOnNarrowStripe(i) = (i-b)*1.0/a;
//		_coordinatesOfBestCentralPathOnNarrowStripe(i) = _coordinatesOfBestCentralPathOnNarrowStripe(_lowerYEndOfBestCentralPath);
//	}
//	for(int i = _upperYEndOfBestCentralPath+1; i <= _endOfThePlastronYUpCoordinate; i++)
//	{
//		_coordinatesOfBestCentralPathOnNarrowStripe(i) = _coordinatesOfBestCentralPathOnNarrowStripe(_upperYEndOfBestCentralPath);
//		//_coordinatesOfBestCentralPathOnNarrowStripe(i) = (i-b)*1.0/a;
//		//cout << _coordinatesOfBestCentralPathOnNarrowStripe(i) << endl;
//	}
//	_lowerYEndOfBestCentralPath = _endOfThePlastronYDownCoordinate;
//	_upperYEndOfBestCentralPath = _endOfThePlastronYUpCoordinate;
//	//for(int y = _lowerYEndOfBestCentralPath; y <= _upperYEndOfBestCentralPath;y++)
//	//{
//	//	int x = _coordinatesOfBestCentralPathOnNarrowStripe(y);
//	//	_narrowStripe(x,y,0) *= 2;
//	//	_narrowStripe(x,y,1) *= 2;
//	//	_narrowStripe(x,y,2) *= 2;
//	//}
//	//_narrowStripe.display();
//}
//
//// returns positions of the best side paths on the best central path of penctoral/femoral junctions
//CImg<int> Tortoise::getYCoordinatesOfPectoralFemoralJunctions(CImg<unsigned char> edgeImgToDoItOn)
//{
////using edges on a medium-wide central-seam stripe
//	CImg<int> stripeEdges(edgeImgToDoItOn);
//	int stripeWidth = stripeEdges.width();
//	int stripeHeight = stripeEdges.height();
//	CImg<int> coordinatesOfBestCentralPathOnStripe(_coordinatesOfBestCentralPathOnMediumWideStripe);
//	CImg<int> stripeWithDiscountValuesBasedOnDirectionOfContinuousEdges(_mediumWideStripeWithDiscountValuesBasedOnDirectionOfContinuousEdges);
////every edge pix is part of a some path - all edge pixs are assigned an index of the path they belong to 
//	CImg<int> sidePathsIndexes(stripeEdges.get_fill(0));
////pixs on the cetral path have 2-column index array, beacause there meet paths from left and right side
//	CImg<int> indexesOnTheBestCentralPath(2,coordinatesOfBestCentralPathOnStripe.height(),1,1, 0);
////prefilling - side paths to the value 10000, for future easier processing
//	CImg<int> sidePathsFromLeft(stripeEdges.get_fill(10000)), sidePathsFromRight(sidePathsFromLeft);
//
//	for(int y = 0; y < stripeHeight; y++)
//	{
//		int x = coordinatesOfBestCentralPathOnStripe(y);
////central path is added to the stripeEdge img and if there are edges (on the original edge img), where the central path is, they are set to 2
//		if(x > 0) stripeEdges(x,y)++;
//	}
////prefilling - the sides (5 the most left columns) of sidePathsFromLeft and sidePathsIndexes
//	double index = 1;
//	for(int y = _lowerYEndOfBestCentralPath; y < _upperYEndOfBestCentralPath; y++)
//	{
//		index++;
//		for(int x = 0; x <= 4; x++)
//		{
////prefilling - for all edge pixs in the area
//			if(stripeEdges(x,y)) 
//			{ 
////prefilling - initial value (distance) of sidePathsFromLeft is set to 2
//				sidePathsFromLeft(x,y) = 2; 
////prefilling - edge pixs in the area are all assigned an index in sidePathsIndexes, which they will pass on to pixs that will be in the same path 
//				sidePathsIndexes(x,y) = index; 
//			}
//		}
//	}
////prefilling - the sides (5 the most right columns) of sidePathsFromRight and sidePathsIndexes
//	index += 1000;
//	for(int y = _lowerYEndOfBestCentralPath; y < _upperYEndOfBestCentralPath; y++)
//	{
//		index++;
//		for(int x = stripeWidth-1; x >= stripeWidth-5; x--)
//		{
//			if(stripeEdges(x,y)) 
//			{ 
////prefilling - initial value (distance) of sidePathsFromRight is set to 2
//				sidePathsFromRight(x,y) = 2; 
////prefilling - edge pixs in the area are all assigned an index in sidePathsIndexes, which they will pass on to pixs that will be in the same path 
//				sidePathsIndexes(x,y) = index; 
//			}
//		}
//	}
//
////prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
//	CImg<int> distanceMask(7,5,1,1, 0);
////prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
//	for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = 1*pow((y+1),2.0);
//	for(int y = 0; y < 5; y++) { distanceMask(0,y) += 8; distanceMask(6,y) += 8; distanceMask(1,y) += 4; distanceMask(5,y) += 4; }
//	distanceMask+=3;
////filling left side paths - rotating distanceMask
//	distanceMask.rotate(90);
////filling left side paths - going though all img pixs (rows only those where the central path is)
//	for(int x = 5; x <= stripeWidth-1; x++)
//	{
//		for(int y = max(5,_lowerYEndOfBestCentralPath); y < min(_upperYEndOfBestCentralPath,stripeHeight-5); y++)
//		{
////filling left side paths - evaluating only edge pixs on the left of the central path
//			if(stripeEdges(x,y)>0 && x<=coordinatesOfBestCentralPathOnStripe(y)) 
//			{
////filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs 
//				CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask + stripeWithDiscountValuesBasedOnDirectionOfContinuousEdges.get_crop(x-1,y-3,x-5,y+3));
////filling left side paths - minimum distance is saved
//				int min = cropAreaDistances.min();
////filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
//				if(min<10000) 
//				{
////filling left side paths - minimum distance is saved in the current pix
//					sidePathsFromLeft(x,y) = min;
////filling left side paths - index of the pix with minimum distace is passed on to the current pix
//				    sidePathsIndexes(x,y) = sidePathsIndexes(x-5+cropAreaDistances.get_stats()(4),y-3+cropAreaDistances.get_stats()(5));
////filling left side paths - if the pix is on the central path, the index is saved in indexesOnTheBestCentralPath
//				    if(x==coordinatesOfBestCentralPathOnStripe(y)) indexesOnTheBestCentralPath(0,y) = sidePathsIndexes(x-5+cropAreaDistances.get_stats()(4),y-3+cropAreaDistances.get_stats()(5));
//				}
//			}
//		}
//	}
////filling right side paths - rotating distanceMask 
//	distanceMask.rotate(180);
////filling right side paths - the same as is in the left side
//	for(int x = stripeWidth-6; x >= 0; x--)
//	{
//		for(int y = max(5,_lowerYEndOfBestCentralPath); y < min(_upperYEndOfBestCentralPath,stripeHeight-5); y++)
//		{
//			if(stripeEdges(x,y)>0 && x >= coordinatesOfBestCentralPathOnStripe(y)) 
//			{
//				CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask + stripeWithDiscountValuesBasedOnDirectionOfContinuousEdges.get_crop(x+1,y-3,x+5,y+3));
//				int min = cropAreaDistances.min();
//				if(min<10000)
//				{
//					sidePathsFromRight(x,y) = min;
//					sidePathsIndexes(x,y) = sidePathsIndexes(x+1+cropAreaDistances.get_stats()(4),y-3+cropAreaDistances.get_stats()(5));
//					if(x==coordinatesOfBestCentralPathOnStripe(y)) indexesOnTheBestCentralPath(1,y) = sidePathsIndexes(x+1+cropAreaDistances.get_stats()(4),y-3+cropAreaDistances.get_stats()(5));
//				}
//			}
//		}
//	}
////postprocessing - set values 10000 to 0 in sidePathsFromLeft and sidePathsFromRight
//	for(int y = 0; y < stripeHeight; y++) 
//	{
//		for(int x = 0; x < stripeWidth; x++)
//		{
//			if(sidePathsFromLeft(x,y) == 10000) sidePathsFromLeft(x,y) = 0;
//			if(sidePathsFromRight(x,y)== 10000) sidePathsFromRight(x,y) = 0;
//		}
//	}
//	_sidePathsOnLeftOnMediumWideStripe = sidePathsFromLeft;
//	_sidePathsOnRightOnMediumWideStripe = sidePathsFromRight;
//	//(sidePathsFromLeft+sidePathsFromRight,stripe).display("sidePaths");
////postprocessing - setting upperSupremum that will be greater than 4x maximum of maxima in sidePathsFromLeft and sidePathsFromRight
//	int upperSupremum = max(sidePathsFromLeft.max(),sidePathsFromRight.max())*4+1;
////postprocessing - CentralSideCrosses will contain sum of 2x values of distances of pixs on the central seam of the left and right paths
//	CImg<int> CentralSideCrosses(1,coordinatesOfBestCentralPathOnStripe.height(),1,1, upperSupremum);
////postprocessing - for all central path pixs near the tortoise center
//	for(int y = _pectoralFemoralJunctionsYDownCoordinate ; y < _pectoralFemoralJunctionsYUpCoordinate; y++)
//	{
//		int x = coordinatesOfBestCentralPathOnStripe(y);
////postprocessing - if there are paths from left and right that end in the central-path pix
//		if(sidePathsFromLeft(x,y)>0 && sidePathsFromRight(x,y)>0)
//		{
////postprocessing - multiplying twice the value of central path pixs
//			sidePathsFromLeft(x,y) *= 2;
//			sidePathsFromRight(x,y) *= 2;
////postprocessing - the values are added to CentralSideCrosses
//			CentralSideCrosses(y) = sidePathsFromLeft(x,y) + sidePathsFromRight(x,y);
//		}
////postprocessing - else, left and right indexes of the central-path pix are set to 0
//		else 
//		{
//			indexesOnTheBestCentralPath(0,y) = 0;
//			indexesOnTheBestCentralPath(1,y) = 0;
//		}
//	}
//	for(int y = _abdominalJunctionsYDownCoordinate - (_abdominalJunctionsYDownCoordinate-_pectoralFemoralJunctionsYDownCoordinate)/5; y < _abdominalJunctionsYUpCoordinate + (_pectoralFemoralJunctionsYUpCoordinate - _abdominalJunctionsYUpCoordinate)/5; y++)
//		CentralSideCrosses(y) = upperSupremum;
//	//(CentralSideCrosses,stripe).display();
//
//	int yMin = CentralSideCrosses.get_stats()(5);
//	if(yMin < m_sPlastronCentreOnRotatedImg.yCoor)
//	{
//		_pectoralFemoralJunctionsYDownCoordinate = yMin;
//		for(int y = m_sPlastronCentreOnRotatedImg.yCoor; y >=0; y--) CentralSideCrosses(y) = upperSupremum;
//		_pectoralFemoralJunctionsYUpCoordinate = CentralSideCrosses.get_stats()(5);
//	}
//	else
//	{
//		_pectoralFemoralJunctionsYUpCoordinate = yMin;
//		for(int y = m_sPlastronCentreOnRotatedImg.yCoor; y < stripeHeight; y++) CentralSideCrosses(y) = upperSupremum;
//		_pectoralFemoralJunctionsYDownCoordinate = CentralSideCrosses.get_stats()(5);
//	}
//	
//	//int y = _pectoralFemoralJunctionsYDownCoordinate;
//	//int x = _coordinatesOfBestCentralPathOnWideStripe(y);	    
//	//for(int i=-10;i<=10;i++) for(int j=-10;j<=10;j++) 
//	//{
//	//	_stripeWithJunctions(x+i,y+j,0) = 0;
//	//	_stripeWithJunctions(x+i,y+j,1) = 255;
//	//	_stripeWithJunctions(x+i,y+j,2) = 0;
//	//}
//	//y = _pectoralFemoralJunctionsYUpCoordinate;
//	//x = _coordinatesOfBestCentralPathOnWideStripe(y);	    
//	//for(int i=-10;i<=10;i++) for(int j=-10;j<=10;j++) 
//	//{
//	//	_stripeWithJunctions(x+i,y+j,0) = 0;
//	//	_stripeWithJunctions(x+i,y+j,1) = 255;
//	//	_stripeWithJunctions(x+i,y+j,2) = 0;
//	//}
//	//(valuesOfTheJunctions.rotate(90),CentralSideCrosses,CentralSideCrossesCopy,stripeWithAbdJunc,valuesInThirdsOfWidth,stripeWithDiscountValuesBasedOnDirectionOfContinuousEdges,stripeEdges).display();
//	//save(_stripeWithJunctions,"AbdoPectFemoJunctions2");
//	return _stripeWithJunctions;
//}
//
//CImg<int> Tortoise::drawJunctions()
//{
//	//_stripeWithJunctions.assign(_wideStripe);
//
//	//int y = _pectoralFemoralJunctionsYDownCoordinate;
//	//int x = _coordinatesOfBestCentralPathOnWideStripe(y);	    
//	//for(int i=-10;i<=10;i++) for(int j=-10;j<=10;j++) 
//	//{
//	//	_stripeWithJunctions(x+i,y+j,0) = 0;
//	//	_stripeWithJunctions(x+i,y+j,1) = 255;
//	//	_stripeWithJunctions(x+i,y+j,2) = 0;
//	//}
//	//y = _pectoralFemoralJunctionsYUpCoordinate;
//	//x = _coordinatesOfBestCentralPathOnWideStripe(y);	    
//	//for(int i=-10;i<=10;i++) for(int j=-10;j<=10;j++) 
//	//{
//	//	_stripeWithJunctions(x+i,y+j,0) = 0;
//	//	_stripeWithJunctions(x+i,y+j,1) = 255;
//	//	_stripeWithJunctions(x+i,y+j,2) = 0;
//	//}
//	//y = _abdominalJunctionsYUpCoordinate;
//	//x = _coordinatesOfBestCentralPathOnWideStripe(y);	    
//	//for(int i=-10;i<=10;i++) for(int j=-10;j<=10;j++) 
//	//{
//	//	_stripeWithJunctions(x+i,y+j,0) = 0;
//	//	_stripeWithJunctions(x+i,y+j,1) = 255;
//	//	_stripeWithJunctions(x+i,y+j,2) = 0;
//	//}
//
//	//y = _abdominalJunctionsYDownCoordinate;
//	//x = _coordinatesOfBestCentralPathOnWideStripe(y);	    
//	//for(int i=-10;i<=10;i++) for(int j=-10;j<=10;j++) 
//	//{
//	//	_stripeWithJunctions(x+i,y+j,0) = 0;
//	//	_stripeWithJunctions(x+i,y+j,1) = 255;
//	//	_stripeWithJunctions(x+i,y+j,2) = 0;
//	//}
//	return _stripeWithJunctions;
//}
//
//CImg<double> Tortoise::getAbdoPectoFemoSidePaths()
//{
//	CImg<int> leftAndRightJunctions(2,4,1,1, 
//		_pectoralFemoralJunctionsYDownCoordinate, _pectoralFemoralJunctionsYDownCoordinate,
//		_abdominalJunctionsYDownCoordinate      , _abdominalJunctionsYDownCoordinate,
//		_abdominalJunctionsYUpCoordinate        , _abdominalJunctionsYUpCoordinate,
//		_pectoralFemoralJunctionsYUpCoordinate  , _pectoralFemoralJunctionsYUpCoordinate );
//
//	//leftAndRightJunctions.display();
//	int areaForSideMinimum = (_pectoralFemoralJunctionsYUpCoordinate-_pectoralFemoralJunctionsYDownCoordinate)/20;
//	//int downIsGreaterIncrease = 0;
//
//
//	for(int i = 0; i < 2; i++)
//	{
//		for(int j = 0; j < 4; j++)
//		{
//			CImg<int> sidePaths;
//			int y = leftAndRightJunctions(i,j);
//			cout << y << endl;
//			int x = _coordinatesOfBestCentralPathOnMediumWideStripe(y);
//			if(i==0) sidePaths.assign(_sidePathsOnLeftOnMediumWideStripe);
//			else sidePaths.assign(_sidePathsOnRightOnMediumWideStripe);
//			int localSideMin = sidePaths(x,y);
//			int yOfLocalSideMin = y;
//			for(int k = -areaForSideMinimum; k <= areaForSideMinimum; k++)
//			{
//				y = leftAndRightJunctions(i,j) + k;
//				x = _coordinatesOfBestCentralPathOnMediumWideStripe(y);
//				if(sidePaths(x,y)<localSideMin && sidePaths(x,y)!=0)
//				{
//					localSideMin = sidePaths(x,y);
//					yOfLocalSideMin = y;
//				}
//			}
//			leftAndRightJunctions(i,j) = yOfLocalSideMin;
//			//if(sidePaths(_coordinatesOfBestCentralPathOnMediumWideStripe(yOfLocalSideMin+areaForSideMinimum),yOfLocalSideMin+areaForSideMinimum)
//			//	< sidePaths(_coordinatesOfBestCentralPathOnMediumWideStripe(yOfLocalSideMin-areaForSideMinimum),yOfLocalSideMin-areaForSideMinimum))
//			//	downIsGreaterIncrease++;
//		}
//	}
//	//cout << downIsGreaterIncrease << " " << areaForSideMinimum << endl;
//	//(leftAndRightJunctions,_mediumWideStripe,_sidePathsOnLeftOnMediumWideStripe,_sidePathsOnRightOnMediumWideStripe).display();
//	_sidePathsPositions = leftAndRightJunctions;
//	CImg<double> AbdPectFemoSegmentsLengths(1,3,1,2, 0);
//	for(int j = 0; j < 2; j++)
//	{
//		for(int i = 0; i < 3 ; i++)
//		{	
//			int y1 = leftAndRightJunctions(j,i);
//			int x1 = _coordinatesOfBestCentralPathOnMediumWideStripe(y1);
//			int y2 = leftAndRightJunctions(j,i+1);
//			int x2 = _coordinatesOfBestCentralPathOnMediumWideStripe(y2);
//			AbdPectFemoSegmentsLengths(0,i,j) = sqrt((double)((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1)));
//		}
//	}
//	double sumOfSegmentsLeft = AbdPectFemoSegmentsLengths.get_channel(0).sum();
//	double sumOfSegmentsRight = AbdPectFemoSegmentsLengths.get_channel(1).sum();
//	for(int i = 0; i < 3 ; i++)
//	{	
//		AbdPectFemoSegmentsLengths(0,i,0) /= sumOfSegmentsLeft;
//		AbdPectFemoSegmentsLengths(0,i,1) /= sumOfSegmentsRight;
//	}
//	//AbdPectFemoSegmentsLengths.display();
//	return AbdPectFemoSegmentsLengths;
//}
//
//CImg<double> Tortoise::getGroundTruthSegmentsLengths()
//{
//	CImg<double> SegmentsLengthsAndJunctionsAngles;
//	for( int i = 0 ; i <= 2 ; i++ )
//	{
//		for( int j = 0 ; j <= 9 ; j++ )
//		{
//			for( int k = 0 ; k <= 9 ; k++ )
//			{
//				if(i==0||(i==1&&j<5)||(i==1&&j==5&&k<=5))
//				{
//					char tortoiseName[10] = {'T' , 'g' , i+48 , j+48 , k+48 };
//					char loadFileDirectoryAndPrefix[1000] = "C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\RucneMereneDelkySvu\\";
//					strncat(loadFileDirectoryAndPrefix, tortoiseName,10);
//					strncat(loadFileDirectoryAndPrefix, "SegmentsLengthsAndJunctionsAngles.cimg",200);
//					//cout << loadFileDirectoryAndPrefix << endl;
//					CImg<double> OneMesurement(loadFileDirectoryAndPrefix);
//					OneMesurement.crop(0,2,0,4);
//					double sumOfSegmentsLeft = OneMesurement.get_channel(0).sum();
//					double sumOfSegmentsRight = OneMesurement.get_channel(1).sum();
//					for(int i = 0; i < 3 ; i++)
//					{	
//						OneMesurement(0,i,0) /= sumOfSegmentsLeft;
//						OneMesurement(0,i,1) /= sumOfSegmentsRight;
//					}
//					SegmentsLengthsAndJunctionsAngles.append(OneMesurement,'x');
//				}
//			}
//		}
//	}
//	CImg<double> segmentsLengthsForAbdPectFemo(SegmentsLengthsAndJunctionsAngles);
//	return segmentsLengthsForAbdPectFemo;
//}
//
//void Tortoise::classifyTortoise()
//{
//	CImg<double> resultArray(156,2,1,1, 0);
//	double minLengthDifference = 10;
//	for(int i = 0; i < 156; i++)
//	{
//		resultArray(i,1) = i;
//		double Pl = _segmentsLengthsAbdPectFemo(0,0,0);
//		double Al = _segmentsLengthsAbdPectFemo(0,1,0);
//		double Fl = _segmentsLengthsAbdPectFemo(0,2,0);
//		double Plall = _segmentsLengthsAbdPectFemoAllTortoises(i,0,0);
//		double Alall = _segmentsLengthsAbdPectFemoAllTortoises(i,1,0);
//		double Flall = _segmentsLengthsAbdPectFemoAllTortoises(i,2,0);
//		double Pr = _segmentsLengthsAbdPectFemo(0,0,1);
//		double Ar = _segmentsLengthsAbdPectFemo(0,1,1);
//		double Fr = _segmentsLengthsAbdPectFemo(0,2,1);
//		double Prall = _segmentsLengthsAbdPectFemoAllTortoises(i,0,1);
//		double Arall = _segmentsLengthsAbdPectFemoAllTortoises(i,1,1);
//		double Frall = _segmentsLengthsAbdPectFemoAllTortoises(i,2,1);
//
//		double lengthDifference1 = sqrt((Pl-Plall)*(Pl-Plall) + (Al-Alall)*(Al-Alall) + (Fl-Flall)*(Fl-Flall) + (Pr-Prall)*(Pr-Prall) + (Ar-Arall)*(Ar-Arall) + (Fr-Frall)*(Fr-Frall));
//		double lengthDifference2 = sqrt((Pl-Flall)*(Pl-Flall) + (Al-Alall)*(Al-Alall) + (Fl-Plall)*(Fl-Plall) + (Pr-Frall)*(Pr-Frall) + (Ar-Arall)*(Ar-Arall) + (Fr-Prall)*(Fr-Prall));
//		double lengthDifference = min(lengthDifference1,lengthDifference2);
//		resultArray(i,0) = lengthDifference;
//	}
//	CImg<double> permutations; /// JE TO PODELENE CELOU DELKOU ... A  TO NENI DOBRE
//	//resultArray.get_row(0).sort(permutations);
//	//_segmentsLengthsAbdPectFemo.append(_segmentsLengthsAbdPectFemoAllTortoises.get_column(81),'x');
//	//for(int i = 0; i < 6; i++)
//	//{
//	//	_segmentsLengthsAbdPectFemo.append(_segmentsLengthsAbdPectFemoAllTortoises.get_column(permutations(i)),'x');
//	//}
//	//(resultArray.get_row(0).sort(permutations), permutations,resultArray,_segmentsLengthsAbdPectFemoAllTortoises,_segmentsLengthsAbdPectFemo).display();
//}
//
//CImg<int> Tortoise::visualizeCentralAndSidePaths()
//{
//	//_sidePathsPositions.display();
//	int stripeWidth = _sidePathsOnLeftOnMediumWideStripe.width();
//	int stripeHeight = _sidePathsOnLeftOnMediumWideStripe.height();
//	_sidePathsOnLeftOnMediumWideStripe.display();
//
//	for(int x = 0; x < stripeWidth; x++) 
//		for(int y = 0; y < stripeHeight; y++) 
//		{
//			if(_sidePathsOnLeftOnMediumWideStripe(x,y) == 0) _sidePathsOnLeftOnMediumWideStripe(x,y) = 10000;
//			if(_sidePathsOnRightOnMediumWideStripe(x,y) == 0) _sidePathsOnRightOnMediumWideStripe(x,y) = 10000;
//		}
//	//_sidePathsOnLeftOnMediumWideStripe.display();
//	CImg<int> stripeEdges(_mediumWideStripeEdges);
//	CImg<int> distanceMask(7,5,1,1, 0);
//	for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = (y+1)*(y+1);
//	for(int y = 0; y < 5; y++) { distanceMask(0,y) += 8; distanceMask(6,y) += 8; distanceMask(1,y) += 4; distanceMask(5,y) += 4; }
//	distanceMask+=2;
//	distanceMask.rotate(-90);
//	//cout << "T" << endl;
//	for(int k = 0; k < 2; k++)
//	{
//		distanceMask.rotate(180);
//		CImg<int> sidePaths;
//		if(k==0) 
//		{
//			sidePaths.assign(_sidePathsOnLeftOnMediumWideStripe);
//			for(int j = 0; j < 4; j++)
//			{
//				int y = _sidePathsPositions(k,j);
//				int x = _coordinatesOfBestCentralPathOnMediumWideStripe(y);
//				//cout << y << endl << x << endl;
//				//_coordinatesOfBestCentralPathOnMediumWideStripe.display();
//
//				while(x > 10)
//				{
//					CImg<int> cropAreaDistances(sidePaths.get_crop(x-5,y-3,x-1,y+3) + distanceMask);
//					//cropAreaDistances.display();
//					int min = cropAreaDistances.min();
//					//cout << min << " " << sidePathsIndexes(X,Y) << " " << leftIndexMin << " " << X << " " << Y << endl;
//					//cropAreaDistances.display();
//					if(min<10000)
//					{
//						x = x + static_cast<int>(cropAreaDistances.get_stats()(4)) - 5;
//						y = y + static_cast<int>(cropAreaDistances.get_stats()(5)) - 3;
//						stripeEdges(x,y)++;
//					}
//				}
//			}
//			//cout << "T" << endl;
//		}
//	
//
//		else 
//		{
//			sidePaths.assign(_sidePathsOnRightOnMediumWideStripe);
//			//_sidePathsOnRightOnMediumWideStripe.display();
//			for(int j = 0; j < 4; j++)
//			{
//				int y = _sidePathsPositions(k,j);
//				int x = _coordinatesOfBestCentralPathOnMediumWideStripe(y);
//				while(x < stripeWidth - 10)
//				{
//					CImg<int> cropAreaDistances(sidePaths.get_crop(x+1,y-3,x+5,y+3) + distanceMask);
//					//cropAreaDistances.display();
//					int min = cropAreaDistances.min();
//					//cout << min << " " << sidePathsIndexes(X,Y) << " " << X << " " << Y << endl;
//					//cropAreaDistances.display();
//					if(min<10000)
//					{
//						//(cropAreaDistances,stripeEdges).display();
//						x = x + cropAreaDistances.get_stats()(4) + 1;
//						y = y + cropAreaDistances.get_stats()(5) - 3;
//						//cout << x << endl << y << endl;
//						stripeEdges(x,y)++;
//
//					}
//				}
//			}
//		}
//	}
//	//cout << "T" << endl;
//
//	//stripeEdges.display();
//	for(int y = _lowerYEndOfBestCentralPath; y <= _upperYEndOfBestCentralPath; y++)
//	{
//		int x = _coordinatesOfBestCentralPathOnMediumWideStripe(y);
//		stripeEdges(x,y) = 2;
//	}
//
//
//	////draw general parabole
//
//	//CImg<double> X(1,5,1,1, 16,30,34,40,56,65,74,80,106),Y(1,5,1,1, 500,550,600,650,700,650,600,550,500),x(X-stripeWidth/2),y(Y-stripeHeight/2),x2(x.get_mul(x)),y2(y.get_mul(y)),xy(x.get_mul(y));
//	//CImg<double> XY(x2.get_append(xy.get_append(y2.get_append(x.get_append(y,'x'),'x'),'x'),'x')), XYsum(5,1,1,1, x2.sum(),xy.sum(),y2.sum(),x.sum(),y.sum());
//	//CImg<double> CoefficientsABCDE(XYsum*((XY.get_transpose()*XY).invert()));
//	//(XY,XYsum,CoefficientsABCDE).display();
//	//for(int x = -stripeWidth/2; x < stripeWidth/2; x++)
//	//	for(int y = -stripeHeight/2; y < stripeHeight/2; y++)
//	//		if( CoefficientsABCDE(0)*1.0*(x)*(x) + CoefficientsABCDE(1)*1.0*(x)*(y) + 
//	//			CoefficientsABCDE(2)*1.0*(y)*(y) + CoefficientsABCDE(3)*1.0*(x) + CoefficientsABCDE(4)*1.0*(y)  - 1 > y )
//	//		{
//	//			stripeEdges(x+stripeWidth/2,y+stripeHeight/2) = 3;
//	//			cout << "YES" << endl;
//	//		}
//	(stripeEdges,_mediumWideStripeEdges).display();
//	return stripeEdges;
//}









//CImg<int> templateImg("ght\\tortoiseTemplateImg.bmp"), contour(GetEdgeImg(templateImg,50,100,0));
    //contour.display();
    //contour.channel(0).display();
    //contour.threshold(1);
    //for(int i = 1; i < contour.width()-1; ++i)
    //{
    //  for(int j=1; j < contour.height(); ++j )
    //  {
    //    if(contour.get_crop(i-1,j-1,i+1,j+1).sum()>2)
    //    {
    //      contour(i,j) = 0;
    //    }
    //  }
    //}
    ////(contour*255).save("ght\\tortoiseTemplateContourNew2.bmp");
    ////for(int i = 1; i < contour.width()-1; ++i)
    ////{
    ////  for(int j=1; j < contour.height(); ++j )
    ////  {
    ////    if(contour.get_crop(i-1,j-1,i+1,j+1).sum()>1)
    ////    {
    ////      contour(i,j) = 0;
    ////    }
    ////  }
    ////}
    //contour.display();
    //int x = 0;
    //int y = 0;
    //cin >> x;
    //cin >> y;
    //contour*=255;
    ////contour.save("ght\\tortoiseTemplateContourNew.bmp");
    ////contour.display();
    ////contour*=255;
    //contour(x,y,0) = 127;
    //contour.display();
    //contour.save("ght\\tortoiseTemplateContour.bmp");



    //CImg<int> contour("ght\\tortoiseTemplateContour.bmp");//,cont(GetEdgeImg(contour,50,100,0));

    //contour.channel(0).threshold(255).display();
    //for(int i = 1; i < contour.width()-1; ++i)
    //{
    //  for(int j=1; j < contour.height(); ++j )
    //  {
    //    if(contour.get_crop(i-1,j-1,i+1,j+1).sum()>2)
    //    {
    //      contour(i,j) = 0;
    //    }
    //  }
    //}
    //for(int i = 1; i < contour.width()-1; ++i)
    //{
    //  for(int j=1; j < contour.height(); ++j )
    //  {
    //    if(contour.get_crop(i-1,j-1,i+1,j+1).sum()>1)
    //    {
    //      contour(i,j) = 0;
    //    }
    //  }
    //}

    //contour.display();
    //contour*=255;
    //contour(65,125,0) = 127;
    //contour.display();
    //contour.save("ght\\tortoiseTemplateContour.bmp");




    //CImg<float> grey(img.get_RGBtoHSV().channel(2)), gradx(grey.get_gradient(0,2)(0)), grady(grey.get_gradient(0,2)(1)), dir(gradx.get_atan2(grady)*180/PI);//, mag((gradx.get_pow(2.0) + grady.get_pow(2.0)).sqrt().normalize(0,255).threshold(1) );







  //void getGoodImgs()
  //{
  //  CImg<double> SomeTgNumbers;
  //  _chdir("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGoodNoYoungTo6cm\\results\\sorted\\OldHough\\");
  //  string directory ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGoodNoYoungTo6cm\\results\\sorted\\OldHough\\");
  //  string goodDir(directory + "GoodOnes\\");
  //  struct _finddata_t c_file;
  //  long hFile;
  //  if ( (hFile = _findfirst("*.pnm", &c_file)) == -1L ) //checks if there is any pnm file in the directory
  //    printf("No *.pnm files in current directory");
  //  else
  //  {
  //    do
  //    {
  //      //int tgNumberINT;
  //      string tgName( c_file.name); //get the file name
  //      cout << tgName.substr(17,4) << endl;
  //      if( tgName.substr(17,4)=="GOOD")// remove( (directory + tgName.substr(0,17) + ".pnm").c_str() );
  //      {
  //       // cout << tgName.substr() << endl;
  //      cout << rename( (directory+tgName.substr(0,17) + ".pnm").c_str(), (goodDir + tgName.substr(0,17) + ".pnm").c_str() ) << endl;
  //      }
  //    } while ( _findnext(hFile, &c_file) == 0 ); //find next file name
  //    _findclose(hFile);
  //  } 
  //}





//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////                                                                                                                   ///////////////////////////////////////////////////
//////////////////////////////  Methods for measuring the seam segments by hand, their evaluation and a tortoise classification based on them.   ///////////////////////////////////////////////////
//////////////////////////////                                                                                                                   ///////////////////////////////////////////////////
//////////////////////////////                   The methods and functions are not part of the Tortoise class.                                   ///////////////////////////////////////////////////
//////////////////////////////                                                                                                                   ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//  //! event Call Back Function used for recording mouse clicks -> for hand measuring 
//  /**
//      \note 
//      - Procedure that retrieves x and y coordinate of a mouse click and saves them into userdata array
//      \outline
//      - When right mouse button is pressed, x and y coordinate of the mouse poiter is saved
//        into userdata array. Last element in userdata stores place of next free elements where
//        next point will be saved.
//        If CTRL key is pressed when pressing right mouse click, then -1 is assign to userdata[0] 
//        (a flag, which will cause skipping the image without saving the measurements).
//        NumberOfPointsToBeMeasured has to be passed here manually.
//  **/
//  void CallBackFunc(int event, int x, int y, int flags, void* userdata)
//  {
//    const int numberOfPointsToBeMeasured = 12;
//	  if ( event == cv::EVENT_RBUTTONDOWN)
//	  {
//		  ((int*)userdata)[((int*)userdata)[2*numberOfPointsToBeMeasured]] = x;
//		  ((int*)userdata)[((int*)userdata)[2*numberOfPointsToBeMeasured]+1] = y;
//		  ((int*)userdata)[2*numberOfPointsToBeMeasured] += 2;
//		  cout << "Right mouse button is released - position (" << x << ", " << y << ")" << endl;
//	  }
//	  if( flags == (cv::EVENT_FLAG_CTRLKEY + cv::EVENT_FLAG_RBUTTON))
//	  {
//		  ((int*)userdata)[0] = -1;
//	  }
//  }
//
//  //! records position of the mouse pointer when clicked on a specified image
//  /**
//      \note 
//      - Function opens an image and by calling CallBackFunction records coordinates of 
//        the right mouse button click.
//      \outline
//      - Loads image.
//        XYcoordinatesOfLeftAndRightJuncitons[24], which contains the current position in the 
//        XYcoordinatesOfLeftAndRightJuncitons array, where recorded coordinates are saved, 
//        is set to 0, so the saving of coordinates starts at the beginning of the array.
//        setMouseCallback -  ensures that while mouse keys are pressed, the CallBackFunc
//        save mouse pointer coordinates into XYcoordinatesOfLeftAndRightJuncitons array.
//        When any other key (i.e. keyboard key) is pressed, the function ends and returns
//        XYcoordinatesOfLeftAndRightJuncitons array.
//  **/
//  void gettingXYcoordinatesOfMouseClick(string loadFileName, int* XYcoordinatesOfLeftAndRightJunctions, const int numberOfPointsToBeMeasured)
//  {
//	  // Read image from file 
//	  cv::Mat img = cv::imread(loadFileName);
//	  //if fail to read the image
//	  if( img.empty() ) 
//	  { 
//		  cout << "Error loading the image" << endl;
//	  }
//	  //Create a window
//	  cv::namedWindow("My window", CV_WINDOW_NORMAL);
//	  XYcoordinatesOfLeftAndRightJunctions[2*numberOfPointsToBeMeasured] = 0;
//	  //set the callback function for any mouse event
//	  cv::setMouseCallback("My window", CallBackFunc, XYcoordinatesOfLeftAndRightJunctions);
//	  //show the image
//	  cv::imshow("My window", img);
//	  // Wait until user press some key
//	  cv::waitKey(0);
//  }
//
//   //! records position of the central seam junctions for all tortoises in loadFileDirectory
//  /**
//      \note 
//      - Calls gettingXYcoordinatesOfMouseClick to record all 12 points.
//      \outline
//      - All tortoises in loadFileDirectory must have their tortoiseNumber, e.g. Tg00301, writen
//        in TortoisesNumbers.txt which has to be in loadFileDirectory
//        If all 12 points are not recorded, it calls gettingXYcoordinatesOfMouseClickto again.
//        Sequence in which points are to be recorded: 
//          1. centralSeamEndHead
//          2.-6. leftJunctions
//          7. centralSeamEndTail
//          8.-12. rightJunctinos
//	      For each point x and y coordinates are recorded, so the array contains 24 numbers in total.
//        If any of the images are not to be measured(recorded), then pres CTRL key and right button key. It 
//        assigns XYcoordinatesOfLeftAndRightJuncitons[0] = -1, a flag that ensures data will not be
//        recorded. To go to anothe image, press any keyboard key.
//        Recorded date is saved in .cimg file measuredPointsByHandOfSeamSegments, in saveFileDirectory.
//  **/
//  void getAndSaveCoordinatesOfCentralSeamJunctionsForAllBenderTortoises(const int numberOfPointsToBeMeasured)
//  {
//    ifstream tortoisesNumbers ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnm\\tortoisesNumbers.txt");
//    string loadFileDirectory ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnm\\");
//    string saveFileDirectory ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizp\\");
//    string tortoiseNumber;
//    for(int i = 0; i < 1304; ++i) // 1304 is the number of tortoises in tortoisesNumbers.txt
//    {
//      getline( tortoisesNumbers, tortoiseNumber);
//      if(i >= 1005) // condition if you want to measure only some of the tortoises
//      {
//        string loadFileName (loadFileDirectory + tortoiseNumber + ".pnm");
//        int* XYcoordinatesOfLeftAndRightJunctions = new int[2*numberOfPointsToBeMeasured+1]; // 2(for x,y)*numberOfPointsToBeMeasured+1(for the count of measured points)
//        gettingXYcoordinatesOfMouseClick( loadFileName, XYcoordinatesOfLeftAndRightJunctions, numberOfPointsToBeMeasured);
//	      while(XYcoordinatesOfLeftAndRightJunctions[2*numberOfPointsToBeMeasured] != 2*numberOfPointsToBeMeasured && XYcoordinatesOfLeftAndRightJunctions[0] != -1)
//	      {
//		      cout << "NEW TRY" << endl;
//		      gettingXYcoordinatesOfMouseClick( loadFileName, XYcoordinatesOfLeftAndRightJunctions, numberOfPointsToBeMeasured);
//	      }
//	      if(XYcoordinatesOfLeftAndRightJunctions[0] != -1) // if == -1 it skips saving for current image
//	      {
//		      CImg<double> measuredPointsByHandOfSeamSegments(1,2*numberOfPointsToBeMeasured,1,1, 0);
//		      for(int i = 0; i < 2*numberOfPointsToBeMeasured; i++)
//		      {
//			      measuredPointsByHandOfSeamSegments(0,i) = XYcoordinatesOfLeftAndRightJunctions[i]; // copying data to CImg that will be saved
//		      }
//          string saveFileName(saveFileDirectory + tortoiseNumber + "measuredPointsByHandOfSeamSegments.cimg");
//          measuredPointsByHandOfSeamSegments.save(saveFileName.c_str());
//          cout << "Measured data for " << tortoiseNumber << " was saved." << endl;
//	      }
//      }
//    }
//  }
//
//  //! calculate seam segments lengths and inner junction differences from measured points
//  /**
//      \note 
//      - Load the meaured points and return seam segments lengths and inner junction differences
//      \outline
//      - First the measured points are loaded from loadFileDirectoryOfMeasuredPoints, the files
//        have form of for example "Tg03301measuredPointsByHandOfSeamSegments.cimg".
//        The directory has to be changed depending on where the measured points are stored.
//        First the loaded points are ordered in a new array MP of size 2(left,right)x7(junctions)x2(x,y-coordinates).
//        The first and last junctions on left and right will have the same cooordinates, for the head and tail 
//        junctions had only one measurement.
//        Seam segments are saved in two columns of 6 elements, where in left column are left segments and in right column 
//        are right segments.
//        Segment is not a mere L2 distance between cosecutive measured point on the left or the right. Because the central
//        seam is not straight always, espectially in inner-junction space.
//        So the segments lenghts are made of the distance between the nearest left and/or right connections of two consecutive 
//        juntions and then the left-right/right-left inner-junction distances are added.
//        In the same time, the leftRight inner junction differences are computed and saved in leftRightJunctionDifferences.
//        Both arrays are then normalized by the central seam lenght and data is saved in the saveFileDirectory directory.
//  **/
//  void calculateSeamLengthsAndJunctionDifferencesFromMeasuredPoints()
//  {
//    ifstream tortoisesNumbers ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnm\\tortoisesNumbers.txt");
//    string loadFileDirectoryOfMeasuredPoints ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizp\\");
//    string saveFileDirectory ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\data\\");
//    string tortoiseNumber;
//    for(int i = 0; i < 1304; ++i)
//    {
//      getline( tortoisesNumbers, tortoiseNumber);
//      if(i>=1005)
//      {
//        string loadNameOfMeasuredPoints(loadFileDirectoryOfMeasuredPoints + tortoiseNumber + "measuredPointsByHandOfSeamSegments.cimg");
//        CImg<double> measuredPoints(loadNameOfMeasuredPoints.c_str()), MP(2,7,1,2,0);
//        for(int i = 0; i <= 1; i++)
//        {
//          MP(i,0,0) = measuredPoints(0, 0); // left/right head junction x-coord
//          MP(i,0,1) = measuredPoints(0, 1); // left/right head junction y-coord
//          MP(i,6,0) = measuredPoints(0,12); // left/right tail junction x-coord
//          MP(i,6,1) = measuredPoints(0,13); // left/right tail junction y-coord   
//        }
//        for(int i = 1; i <= 5; i++)
//        {
//          MP(0,i,0) = measuredPoints(0,i*2   ); // left  i-th junction x-coord
//          MP(0,i,1) = measuredPoints(0,i*2+ 1); // left  i-th junction y-coord
//          MP(1,i,0) = measuredPoints(0,i*2+12); // right i-th junction x-coord
//          MP(1,i,1) = measuredPoints(0,i*2+13); // right i-th junction y-coord       
//        }
//        //calculate segments lenghts and junction differences at the same time
//	      CImg<double> seamSegments(2,6,1,1, 0), leftRightJunctionDifferences(1,5,1,1, 0);
//        double lastNearestX = MP(0,0,0);
//        double lastNearestY = MP(0,0,1);
//        for(int i = 0; i < 6; i++)
//        {
//          double distanceFromLastNearestPointOnLeft = sqrt( (double) (lastNearestX-MP(0,i+1,0))*(lastNearestX-MP(0,i+1,0)) 
//                                                                   + (lastNearestY-MP(0,i+1,1))*(lastNearestY-MP(0,i+1,1)) );
//          double distanceFromLastNearestPointOnRight = sqrt( (double) (lastNearestX-MP(1,i+1,0))*(lastNearestX-MP(1,i+1,0)) 
//                                                                    + (lastNearestY-MP(1,i+1,1))*(lastNearestY-MP(1,i+1,1)) );
//          double differenceBetweenLeftRightJunctionPoints = sqrt( (double) (MP(0,i+1,0)-MP(1,i+1,0))*(MP(0,i+1,0)-MP(1,i+1,0)) 
//                                                                         + (MP(0,i+1,1)-MP(1,i+1,1))*(MP(0,i+1,1)-MP(1,i+1,1)) );
//          if (distanceFromLastNearestPointOnLeft < distanceFromLastNearestPointOnRight)
//          {
//            seamSegments(0,i) += distanceFromLastNearestPointOnLeft;
//            seamSegments(1,i) += distanceFromLastNearestPointOnLeft + differenceBetweenLeftRightJunctionPoints;
//            lastNearestX = MP(1,i+1,0);
//            lastNearestY = MP(1,i+1,1);
//            if (i!=5)
//            {
//             seamSegments(0,i+1) += differenceBetweenLeftRightJunctionPoints;
//             leftRightJunctionDifferences(0,i) = - differenceBetweenLeftRightJunctionPoints;
//            }      
//          }
//          else
//          {
//            seamSegments(1,i) += distanceFromLastNearestPointOnLeft;
//            seamSegments(0,i) += distanceFromLastNearestPointOnLeft + differenceBetweenLeftRightJunctionPoints;
//            lastNearestX = MP(0,i+1,0);
//            lastNearestY = MP(0,i+1,1);
//            if (i!=5)
//            {
//             seamSegments(1,i+1) += differenceBetweenLeftRightJunctionPoints;
//             leftRightJunctionDifferences(0,i) = differenceBetweenLeftRightJunctionPoints;
//            }  
//          }
//        }
//        //normalizing to the size of the plastron
//	      double lengthOfCentralSeam = seamSegments.get_column(0).sum();
//	      for(int i = 0; i < 6; i++) 
//        {
//          seamSegments(0,i) /= lengthOfCentralSeam;
//          seamSegments(1,i) /= lengthOfCentralSeam;
//          if (i!=5) leftRightJunctionDifferences(0,i) /= lengthOfCentralSeam;
//        }
//        //genarate save directory
//        string saveNameForSeamSegments(saveFileDirectory + tortoiseNumber + "seamSegments.cimg");
//        string saveNameForLeftRightJunctionDifferences(saveFileDirectory + tortoiseNumber + "leftRightJunctionDifferences.cimg");
//        //save in C-string directory
//        seamSegments.save(saveNameForSeamSegments.c_str());
//        leftRightJunctionDifferences.save(saveNameForLeftRightJunctionDifferences.c_str());
//      }
//    }
//  }
//
//  //! a procedure that puts all measured seam segments lengths and inner junction differences into one file
//  /**
//      \note 
//      - We want all seam segments lenghts and all inner junction distances in one file 
//        for all 1304 Bender and Cizp tortoises. And we want left and right seam segments lengths and inner
//        junction differences in one column.
//      \outline
//      - With use of numbering of Bender and Cizp Tortoises in tortoisesNumbers we go though all of them (segments lengths 
//        and inner junction distances) and append them to the file MeasuredFeaturesForAllTortoises which are saved at the end.
//        Also seamSegments that are in two columns, are put in one column (column 0 is appended below column 1) and
//        below that the inner junction differences are appended.
//  **/
//  void putTheCalculatedSeamDistancesForAllTortoisesIntoOneFile()
//  {
//    ifstream tortoisesNumbers ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnm\\tortoisesNumbers.txt");
//    CImg<double> MeasuredFeaturesForAllTortoises;
//    string loadFileDirectoryOfMeasuredFeatures ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\data\\");
//    string tortoiseNumber;
//    for(int i = 0; i < 1304; ++i)
//    {
//      getline( tortoisesNumbers, tortoiseNumber);
//
//      string loadSeamSegments(loadFileDirectoryOfMeasuredFeatures + tortoiseNumber + "seamSegments.cimg");
//      string loadJunctionDifferences(loadFileDirectoryOfMeasuredFeatures + tortoiseNumber + "leftRightJunctionDifferences.cimg");
//      CImg<double> seamSegments(loadSeamSegments.c_str()), junctionDifferences(loadJunctionDifferences.c_str()) ;
//      seamSegments.column(0).append(seamSegments.get_column(1),'y').append(junctionDifferences,'y');
//      MeasuredFeaturesForAllTortoises.append(seamSegments,'x');
//    }
//    MeasuredFeaturesForAllTortoises.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\MeasuredFeaturesForAllTortoises.cimg");
//  }
//
//  //! save TG numbers, read from file names in a given directory, into a cimg file 
//  /**
//      \outline 
//      - We want it to read file names, in our case pnm files fo chosen selection of TGs, from a given directory.
//        Then, for each TG file name, extract the TG number and convert it to integer (Tg00301 -> 00301 -> 301) and then save it
//        into CImg variable, which is saved at the end.
//  **/
//  void getAndSaveTGnumbersOfTGfileNamesFromGivenDirectory()
//  {
//    CImg<double> SomeTgNumbers;
//    _chdir("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGoodNoYoungTo6cm\\"); //changes the directory
//    struct _finddata_t c_file;
//    long hFile;
//    if ( (hFile = _findfirst("*.pnm", &c_file)) == -1L ) //checks if there is any pnm file in the directory
//      printf("No *.pnm files in current directory");
//    else
//    {
//      do
//      {
//        int tgNumberINT;
//        string tgNumber( c_file.name); //get the file name
//        istringstream (tgNumber.substr(2,5) ) >> tgNumberINT; //get the number from the name and convert it to integer
//        CImg<int> TGnumber(1,1,1,1, tgNumberINT);
//        SomeTgNumbers.append(TGnumber,'x');
//      } while ( _findnext(hFile, &c_file) == 0 ); //find next file name
//      _findclose(hFile);
//    } 
//    SomeTgNumbers.display();
//    SomeTgNumbers.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\TgNumbersInINTonlyGoodTGsAndNoYoungTGs.cimg");
//  }
//
//  //! given a subset of TG numbers ( < 1304), it gets corresponding subset of  measured features
//  /**
//      \outline 
//      - Based on the set of TG numbers of all TGs and its subset of some TG numbers, it will select
//        the corresponding measured features for the subset of TGs and will save it.
//  **/
//  void getAndSaveTgFeatureSubSetForAnyChosenSubSetOfTGnumbers()
//  {
//    CImg<double> MeasuredFeaturesForAllTortoises("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\MeasuredFeaturesForAllTortoises.cimg");
//    CImg<int> AllTgNumbers("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\AllTgNumbersInINT.cimg");
//    CImg<int> SomeTgNumbers("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\TgNumbersInINTonlyGoodTGs.cimg");
//    int numOfFeatures = MeasuredFeaturesForAllTortoises.height();
//    for(int i = 0; i < AllTgNumbers.width(); ++i)
//    {
//      if( (SomeTgNumbers - AllTgNumbers(i)).abs().min() > 0 ) //checks if a TG number is not in the subset of TGs
//      {
//        for(int j = 0; j < MeasuredFeaturesForAllTortoises.height(); ++j) //for such a TG all its measurements are set to 1000
//        {
//          MeasuredFeaturesForAllTortoises(i,j) = 1000;
//        }
//      }
//    }
//    for(int i = 0; i < MeasuredFeaturesForAllTortoises.width(); ++i) //eliminates measured features of TG not belonging to the subset
//    {
//      if( MeasuredFeaturesForAllTortoises.get_column(i).sum() == numOfFeatures*1000 )
//      {
//        if (i==0)
//        {
//          MeasuredFeaturesForAllTortoises = MeasuredFeaturesForAllTortoises.get_columns(1,MeasuredFeaturesForAllTortoises.width()-1);
//        }
//        else if (i==MeasuredFeaturesForAllTortoises.width()-1)
//        {
//          MeasuredFeaturesForAllTortoises = MeasuredFeaturesForAllTortoises.get_columns(0,MeasuredFeaturesForAllTortoises.width()-2);
//        }
//        else
//        {
//          MeasuredFeaturesForAllTortoises = MeasuredFeaturesForAllTortoises.get_columns(0,i-1).append(MeasuredFeaturesForAllTortoises.get_columns(i+1,MeasuredFeaturesForAllTortoises.width()-1),'x');
//        }
//        --i;
//      }
//    }
//    MeasuredFeaturesForAllTortoises.display();
//    MeasuredFeaturesForAllTortoises.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\MeasuredFeaturesForOnlyGoodTGs.cimg");
//  }
//
//  //! get and save groups of features and groups of numbers of same tortoises and of different tortoises into a CImgList file
//  /**
//      \outline 
//      - For a given selection of TG measured features and coresponding TG numbers, it compose all possible groups of TGs,
//        where in a group there are multiple measurements/names (of different images) of the same individual. In the last 
//        group there are measurements/names for all individual tortoises, for each only one measurement.
//  **/
//  void getAndSaveGroupsOfFeaturesAndGroupsOfNumbersOfSameAndDifferentTGs()
//  {
//    CImg<double> MeasuredFeaturesForSelectedTortoises ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\MeasuredFeaturesForOnlyGoodTGsAndNoYoungTGs.cimg");
//    CImg<int> NumbersOfSelectedTortoises ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\numbersOfTGsInINTonlyGoodTGsAndNoYoungTGs.cimg");
//    CImgList<double> groupsOfFeatures;
//    CImgList<int> groupsOfTGnumbers;
//    int indexOfCurrentTGgroup = 0;  // numbers/measurements index to the TG that already belong to the next group
//    while (indexOfCurrentTGgroup < NumbersOfSelectedTortoises.width()-1)
//    {
//      int currentTGnumber = floor(NumbersOfSelectedTortoises(indexOfCurrentTGgroup)*1.0/100); //gets the number of current TG individual (301 -> 3, 12603 -> 126)
//      int indexOfNextTGgroup = indexOfCurrentTGgroup + 1;
//      int nextTGnumber = floor(NumbersOfSelectedTortoises(indexOfNextTGgroup)*1.0/100); //gets the number of next TG individual
//      if(currentTGnumber==nextTGnumber) //if next measurement belongs to the same individual, than it will forms a new group
//      {
//        while(currentTGnumber==nextTGnumber && indexOfNextTGgroup < NumbersOfSelectedTortoises.width())//finds the last number/measurement that belong to the same individual
//        {
//          ++indexOfNextTGgroup;
//          nextTGnumber = floor(NumbersOfSelectedTortoises(indexOfNextTGgroup)*1.0/100);
//        }
//        groupsOfFeatures.push_back(MeasuredFeaturesForSelectedTortoises.get_columns(indexOfCurrentTGgroup,indexOfNextTGgroup-1)); //stors the new group of measurements
//        groupsOfTGnumbers.push_back(NumbersOfSelectedTortoises.get_columns(indexOfCurrentTGgroup,indexOfNextTGgroup-1)); //stors the new group of TG numbers
//        indexOfCurrentTGgroup = indexOfNextTGgroup; //set index to the next TG individual
//      }
//      else
//      {
//        indexOfCurrentTGgroup = indexOfNextTGgroup; //set index to the next TG individual
//      }
//    }
//    indexOfCurrentTGgroup = 0;
//    CImg<double> groupOfDifferentTortoises; //each individual will be in the group only once
//    CImg<int> numbersOfDifferentTortoises;
//    while (indexOfCurrentTGgroup < NumbersOfSelectedTortoises.width())
//    {
//      if(indexOfCurrentTGgroup==NumbersOfSelectedTortoises.width()-1) //if index is on the last TG it has therefore only single measurement TG so it is stored
//      {
//        groupOfDifferentTortoises.append(MeasuredFeaturesForSelectedTortoises.get_column(indexOfCurrentTGgroup),'x');
//        numbersOfDifferentTortoises.append(NumbersOfSelectedTortoises.get_column(indexOfCurrentTGgroup),'x');
//        ++indexOfCurrentTGgroup;
//      }
//      else //for any other index it checks if the individual has one or more measurements
//      {
//        int currentTGnumber = floor(NumbersOfSelectedTortoises(indexOfCurrentTGgroup)*1.0/100);
//        int indexOfNextTGgroup = indexOfCurrentTGgroup + 1;
//        int nextTGnumber = floor(NumbersOfSelectedTortoises(indexOfNextTGgroup)*1.0/100);
//        if(currentTGnumber==nextTGnumber) //if the individual has more measurements
//        {
//          while(currentTGnumber==nextTGnumber && indexOfNextTGgroup < NumbersOfSelectedTortoises.width())
//          {
//            ++indexOfNextTGgroup;
//            nextTGnumber = floor(NumbersOfSelectedTortoises(indexOfNextTGgroup)*1.0/100);
//          }
//          groupOfDifferentTortoises.append(MeasuredFeaturesForSelectedTortoises.get_column(floor((indexOfCurrentTGgroup+indexOfNextTGgroup-1)/2.0)+1),'x'); //it selects one of the measurements (in the middle)
//          numbersOfDifferentTortoises.append(NumbersOfSelectedTortoises.get_column(floor((indexOfCurrentTGgroup+indexOfNextTGgroup-1)/2.0)+1),'x'); //it selects one of the measurements (in the middle)
//          indexOfCurrentTGgroup = indexOfNextTGgroup;
//        }
//        else//if the individual doesnt have more measurements it is stored 
//        {
//          groupOfDifferentTortoises.append(MeasuredFeaturesForSelectedTortoises.get_column(indexOfCurrentTGgroup),'x');
//          numbersOfDifferentTortoises.append(NumbersOfSelectedTortoises.get_column(indexOfCurrentTGgroup),'x');
//          indexOfCurrentTGgroup = indexOfNextTGgroup;
//        }
//      }
//    }
//    groupsOfFeatures.push_back(groupOfDifferentTortoises);
//    groupsOfTGnumbers.push_back(numbersOfDifferentTortoises);
//    groupsOfTGnumbers.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfTGnumbersForOnlyGoodTGsAndNoYoungTGs.cimg");
//    groupsOfFeatures.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfMeasuredFeaturesForOnlyGoodTGsAndNoYoungTGs.cimg");
//  }
//
//  //! PCA - get eigen values and eigen vectors of correlation matrix of hand measured lengths for all groups (of same and of different TGs) and return transformed feature vectors
//  /**
//      \outline 
//      - For all the measurements (the segments lengths and the junction difference), that are first standardized it will compute the correlation 
//        matrix and that compute and save its eigenvectors and eigenvalues. All is saved in a txt file, plus the correlation
//        matrix is computed and saved there as well.
//        The method return tranformed measurements.
//
//  **/
//  CImg<double> PCA(CImg<double> measurementsStandardized, double PCAthreshold, int whichGroups)
//  {
//    //measurementsStandardized.display();
//	  unsigned int numOfFeatures = measurementsStandardized.width();
//	  unsigned int numOfMeasurements = measurementsStandardized.height();
//	  CImg<double> measurementsCovarianceMatrix ( (measurementsStandardized.get_transpose()*measurementsStandardized)/(numOfMeasurements-1) ); //computing the covariance matrix
//    //measurementsCovarianceMatrix.display();
//	  CImg<double> eigenValues, eigenVectors;
//	  measurementsCovarianceMatrix.symmetric_eigen(eigenValues,eigenVectors); //computing the eigen vals and vecs
//    eigenValues /= eigenValues.sum(); //normalizing the eigen vals
//    //(eigenValues,eigenVectors,measurementsCovarianceMatrix).display();
//    ofstream eigenvaluesLog("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\eigenvalues.txt",ios_base::app);
//    eigenvaluesLog << "Group: " << whichGroups << endl;
//    for(int i = 0; i < numOfFeatures; i++) eigenvaluesLog << right << setw(20) << eigenValues(i) << endl;
//    eigenvaluesLog << endl << endl;
//    //eigenValues.display_graph();
//
//    int numberOfUsefulEigenVectors = 0;
//    for(int i = 0; i < numOfFeatures; ++i)
//    {
//      if(eigenValues(i) > PCAthreshold) // threshold for the eigenvalues
//      {
//        numberOfUsefulEigenVectors = i;
//      }
//    }
//    if(numberOfUsefulEigenVectors>0) eigenVectors.columns(0,numberOfUsefulEigenVectors);
//    else eigenVectors.column(numberOfUsefulEigenVectors);
//    CImg<double> transformedMeasurementsStandardized;
//    transformedMeasurementsStandardized = measurementsStandardized*eigenVectors;
//    //(measurementsStandardized,transformedMeasurementsStandardized).display();
//    return transformedMeasurementsStandardized;
//  }
//
//  //! get the training/test set, the training/test labels and the training/test TG numbers from a given group a selected features
//  /**
//      \outline 
//      - For a given group of feature we want to get a training set and a test set that will contain differences of features of same
//        individuals and of different individuals, in both sets there will be equally differences of pairs of same TGs and different TGs.
//        Then we want to get training labels and test labels that will indicate whether a difference of a pair belong to same TGs or
//        different TGs, will have values 1 or -1. And finally we want to know TG numbers of all pairs that are is training and test sets.
//  **/
//  CImgList<double> getTrainingAndTestSetsLabelsAndTGnumbersForAGivenGroupOfFeatures(int whichGroups, bool doStandardization, bool doPCA, double PCAthreshold)
//  {
//    ofstream log ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\kNNResults_log.txt",ios_base::app);
//    CImgList<double> groupsOfFeatures;
//    CImgList<int> groupsOfTGnumbers;
//    if(whichGroups==0) // choose which group of features will be used
//    {
//      groupsOfFeatures.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfMeasuredFeaturesForAllTortoises.cimg");
//      groupsOfTGnumbers.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfTGnumbersForAllTGs.cimg");
//    }
//    else if (whichGroups==1)
//    {
//      groupsOfFeatures.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfMeasuredFeaturesForOnlyGoodTGs.cimg");
//      groupsOfTGnumbers.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfTGnumbersForOnlyGoodTGs.cimg");
//    }
//    else if (whichGroups==2)
//    {
//      groupsOfFeatures.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfMeasuredFeaturesForOnlyGoodTGsAndNoYoungTGs.cimg");
//      groupsOfTGnumbers.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfTGnumbersForOnlyGoodTGsAndNoYoungTGs.cimg");
//    }
//    log << left <<setw(10) << "Groups:" << setw(20) << whichGroups << setw(15) << "withStandrd:" << setw(10) << doStandardization << setw(10) << "withPCA:" << setw(10) << doPCA << setw(10) << "PCAthres:" << setw(20) << PCAthreshold << endl;
//
//    unsigned int numOfGroups = groupsOfFeatures.size();
//    log << setw(30) << "NumberOfGroups:" << setw(20) << numOfGroups << endl;
//
//    unsigned int numOfPairsInTheSetOfSameTGs = 0;
//
//    //groupsOfFeatures.display();
//    for(unsigned int groupIndex = 0; groupIndex < numOfGroups; ++groupIndex) 
//    {
//      groupsOfFeatures(groupIndex).transpose(); //transpose of the groups
//      if(groupIndex < numOfGroups-1) 
//      {
//        int numOfTGsInGroup = groupsOfFeatures(groupIndex).height();
//        numOfPairsInTheSetOfSameTGs += (numOfTGsInGroup*(numOfTGsInGroup-1)/2);
//      }
//    }
//    //groupsOfFeatures.display();
//    unsigned int numOfFeatures = groupsOfFeatures(0).width();
//    log << setw(30) << "NumOfFeatures:" << setw(20) << numOfFeatures << endl;
//    log << setw(30) << "NumOfSameAndDiffPairs:" << setw(20) << numOfPairsInTheSetOfSameTGs << endl;
//    cout << numOfFeatures << " " << numOfPairsInTheSetOfSameTGs << endl;
//
//    cout << 1 << endl;
//    CImg<double> setOfDifferencesForSameTGs, setOfTGnumbersForSameTGs; // (i,0) - segment differences, (i,1) - junction differences, (i,2) - label
//    for(int groupIndex = 0; groupIndex < numOfGroups -1 ; ++groupIndex) // prepare the set of same TGs and their differences
//    {
//      CImg<double> groupOfTGs(groupsOfFeatures(groupIndex)), groupOfTGnumbers(groupsOfTGnumbers(groupIndex));
//      for(int i = 0; i < groupOfTGs.height(); ++i)
//      {
//        for(int j = i+1; j < groupOfTGs.height(); ++j)
//        {
//          CImg<float> class1 (1,1,1,1, 1.0);
//          setOfDifferencesForSameTGs.append( (groupOfTGs.get_row(i)-groupOfTGs.get_row(j)).abs().append(class1,'x'), 'y');
//          setOfTGnumbersForSameTGs.append(groupOfTGnumbers.get_crop(i,0,i,0).append(groupOfTGnumbers.get_crop(j,0,j,0),'x'),'y');
//        }
//      }
//    }
//
//    cout << 2 << endl;
//    CImg<double> groupOfDifferentTGs(groupsOfFeatures(numOfGroups-1)), groupOfNumbersOfDifferentTGs(groupsOfTGnumbers(numOfGroups-1)), setOfDifferencesForDifferentTGs, setOfTGnumbersForDifferentTGs;
//    unsigned int numOfDifferentTGs = groupOfDifferentTGs.height();
//    unsigned int numOfPairsOfDifferentTGsAllCombinations = numOfDifferentTGs*(numOfDifferentTGs-1)/2;
//    unsigned int step = numOfPairsOfDifferentTGsAllCombinations/numOfPairsInTheSetOfSameTGs;
//    unsigned int counterOfPairs = 0;
//    for(int i = 0; i < numOfDifferentTGs; ++i) // prepare the set of different TGs and their differences
//    {
//      for(int j = i+step; j < numOfDifferentTGs; j+=step)
//      {
//        if(counterOfPairs < numOfPairsInTheSetOfSameTGs)
//        {
//          CImg<float> class2 (1,1,1,1, 0);
//          setOfDifferencesForDifferentTGs.append( (groupOfDifferentTGs.get_row(i)-groupOfDifferentTGs.get_row(j)).abs().append(class2,'x'), 'y');
//          setOfTGnumbersForDifferentTGs.append(groupOfNumbersOfDifferentTGs.get_crop(i,0,i,0).append(groupOfNumbersOfDifferentTGs.get_crop(j,0,j,0),'x'),'y');
//          ++counterOfPairs;
//        }
//      }
//    }
//
//    cout << 3 << endl;
//
//    if(doStandardization)
//    {
//      CImg<double> measurementsStandardized;
//      CImg<double> setOfDifferencesForSameAndDifferentTGs(setOfDifferencesForSameTGs.get_append(setOfDifferencesForDifferentTGs,'y'));
//      //(setOfDifferencesForSameTGs,setOfDifferencesForDifferentTGs).display();
//	    unsigned int numOfMeasurements = setOfDifferencesForSameAndDifferentTGs.height();
//	    CImg<double> meansOfFeatures(numOfFeatures,1);
//      CImg<double> stdsOfFeatures(numOfFeatures,1);
//	    for(int i = 0; i < numOfFeatures; i++) meansOfFeatures(i) = setOfDifferencesForSameAndDifferentTGs.get_column(i).mean();
//      for(int i = 0; i < numOfFeatures; i++) stdsOfFeatures(i) = sqrt(setOfDifferencesForSameAndDifferentTGs.get_column(i).variance());
//	    
//      for(int i = 0; i < numOfFeatures; i++) measurementsStandardized.append((setOfDifferencesForSameAndDifferentTGs.get_column(i)-meansOfFeatures(i))/stdsOfFeatures(i),'x'); //standardizing the features
//      //(setOfDifferencesForSameAndDifferentTGs,measurementsStandardized).display("standardized");
//      if(doPCA)
//      {
//        CImg<double> PCAresult( PCA( measurementsStandardized.get_columns(0,numOfFeatures-1), PCAthreshold, whichGroups ) );
//        setOfDifferencesForSameAndDifferentTGs = PCAresult.append(setOfDifferencesForSameAndDifferentTGs.get_column(numOfFeatures),'x');
//      }
//      else
//      {
//        setOfDifferencesForSameAndDifferentTGs = measurementsStandardized.append(setOfDifferencesForSameAndDifferentTGs.get_column(numOfFeatures),'x');
//      }
//      setOfDifferencesForSameTGs = setOfDifferencesForSameAndDifferentTGs.get_rows(0,numOfPairsInTheSetOfSameTGs-1);
//      setOfDifferencesForDifferentTGs = setOfDifferencesForSameAndDifferentTGs.get_rows(numOfPairsInTheSetOfSameTGs,2*numOfPairsInTheSetOfSameTGs-1);
//
//    }
//    numOfFeatures = setOfDifferencesForDifferentTGs.width()-1;
//    log << setw(30) << "NewNumOfFeatures:" << setw(20) << numOfFeatures << endl;
//    cout << numOfFeatures << " " << setOfDifferencesForSameTGs.height() << " " << setOfDifferencesForDifferentTGs.height() << endl;
//    //(setOfDifferencesForSameTGs,setOfDifferencesForDifferentTGs).display("sets ater trans");
//
//
//    CImg<double> trainingSet, testSet, trainingLabels, testLabels, trainingSetTGnumbers, testSetTGnumbers;
//    for(int i = 0; i < numOfPairsInTheSetOfSameTGs; ++i) // divide sets into training set and test set - each second is taken to test set from sets with differences of same and different TGs 
//    {
//      if(i%2)
//      {
//        trainingSet.append(setOfDifferencesForSameTGs.get_row(i),'y');
//        trainingSet.append(setOfDifferencesForDifferentTGs.get_row(i),'y');
//        trainingSetTGnumbers.append(setOfTGnumbersForSameTGs.get_row(i),'y');
//        trainingSetTGnumbers.append(setOfTGnumbersForDifferentTGs.get_row(i),'y');
//      }
//      else
//      {
//        testSet.append(setOfDifferencesForSameTGs.get_row(i),'y');
//        testSet.append(setOfDifferencesForDifferentTGs.get_row(i),'y');
//        testSetTGnumbers.append(setOfTGnumbersForSameTGs.get_row(i),'y');
//        testSetTGnumbers.append(setOfTGnumbersForDifferentTGs.get_row(i),'y');
//      }
//    }
//    cout << testSet.height() << " " << trainingSetTGnumbers.height() << endl;
//    log << setw(30) << "numOfPairsInTrainingSet:" << setw(20) << trainingSet.height() << endl;
//    log << setw(30) << "numOfPairsInTestSet:" << setw(20) << testSet.height() << endl;
//    trainingLabels = trainingSet.get_column(numOfFeatures);
//    testLabels = testSet.get_column(numOfFeatures);
//    trainingSet.columns(0,numOfFeatures-1);
//    testSet.columns(0,numOfFeatures-1);
//    log << endl << endl;
//
//    //(trainingSet, testSet).display();
//    //for(int i = 0; i < trainingSet.height(); ++i) // here you can choose to sum the features
//    //{
//    //  trainingSet(0,i) = trainingSet.get_row(i).sum(); 
//    //}
//    //trainingSet.column(0);
//    //for(int i = 0; i < testSet.height(); ++i)
//    //{
//    //  testSet(0,i) = testSet.get_row(i).sum(); 
//    //}
//    //testSet.column(0);
//
//    return trainingSet,trainingLabels,testSet,testLabels,trainingSetTGnumbers,testSetTGnumbers;
//  }
//  
//  //! k - nearest neigbour classifier - OpenCV implementation
//  /**
//      \outline 
//      - Classifier designed to classify pairs of TGs based on their differences of features. 
//        As input it takes training set, training labels, test set, test labels, training numbers and test numbers.
//        Successfulness of a classification is save into txt file, where faulse positive and faulse negative errors are.
//        Also total errors for each k are also saved in extra txt file with suffix _Errors.
//        Also TG numbers of outliers are save in extra txt file with suffix _Outliers.
//  **/  
//  void kNNclassifier()
//  {
//
//    for(int grp = 0; grp <=2; ++grp)
//    {
//      string grps;
//      if(grp==0) grps.assign("_group0");
//      if(grp==1) grps.assign("_group1");
//      if(grp==2) grps.assign("_group2");
//
//      string classificationResultSuffix(grps+"_WithStdz_WithPCA_0.03Thres_11PCs");
//      string classificationResultsDirectory("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\");
//      string kNNclassificationResultsFileName (classificationResultsDirectory + "kNNResults" + classificationResultSuffix);
//      CImgList<double> featureSet(getTrainingAndTestSetsLabelsAndTGnumbersForAGivenGroupOfFeatures(grp,1,1,0.03));
//      CImg<double> trainingSet(featureSet(0)), trainingLabels(featureSet(1)),
//                  testSet(featureSet(2)), testLabels(featureSet(3)),trainingSetTGnumbers(featureSet(4)),testSetTGnumbers(featureSet(5));
//      //featureSet.display();
//      const int K = 100; // max k allowed
//      int numberOfTrainingData = trainingSet.height();
//      int numberOfDataDimensions = trainingSet.width();
//      int numberOfTestData = testSet.height();
//      cv::Mat trainData ( numberOfTrainingData, numberOfDataDimensions, CV_32FC1 );
//      cv::Mat trainClasses ( numberOfTrainingData, 1, CV_32FC1 );
//      for(int i = 0; i < numberOfTrainingData; ++i)
//      {
//        for(int j = 0; j < numberOfDataDimensions; ++j)
//        {
//          trainData.at<float>(i,j) = trainingSet(j,i);
//        }
//      }
//      for(int i = 0; i < numberOfTrainingData; ++i)
//      {
//        trainClasses.at<float>(i,0) = trainingLabels(0,i);
//      }
//      CvKNearest knn( trainData, trainClasses, cv::Mat(), false, K ); // learning of the classifier, using training set and traning labels
//    
//      ofstream kNNClassificationResults (kNNclassificationResultsFileName + ".txt");
//      ofstream kAndErrors (kNNclassificationResultsFileName + "_Errors" + ".txt");
//      ofstream outliers (kNNclassificationResultsFileName + "_Outliers" + ".txt");
//      kAndErrors << setw(20) << "Number k - neighbor" << setw(30) << "totalOptimisticError in %" << setw(30) << "totalPesimisticError in %" << setw(30) << "averageTotalError in %" << setw(30) << "ratting of k" << '\n';
//      outliers << "Class 0 is the class of different tortoises and Class 1 is the class of same tortoises." << endl;
//      for(int k = 1; k < K; k+=1)
//      {
//        if(k>30 && k < K-10) k+=10;
//        kAndErrors << setw(20) << k;
//        cout << k << endl;
//        //int imgSize = 5000;
//        //CImg<int> resultsOnImg(imgSize,3,1,2, 10);  // if feature space is 2d then results can be save onto an image
//        //double sumOf1DfeaturesInClass0 = 0;
//        //double sumOfSquaresOf1DfeaturesInClass0 = 0;
//        //double countOfSamplesInClass0 = 0;
//        //double sumOf1DfeaturesInClass1 = 0;
//        //double sumOfSquaresOf1DfeaturesInClass1 = 0;
//        //double countOfSamplesInClass1 = 0;
//        //for (int i =  0; i < numberOfTrainingData; ++i) // OPTIMISTIC CLASSIFICATION - using training set
//        //{
//        //  cv::Mat sampleMat (1,1,CV_32FC1);
//        //  sampleMat.at<float>(0,0) = trainingSet(0,i);
//        //  float response = knn.find_nearest(sampleMat,k,0,0,0,0);
//        //  if (response == 1) //TGs are same
//        //  {
//        //    sumOf1DfeaturesInClass1 +=  trainingSet(0,i);
//        //    sumOfSquaresOf1DfeaturesInClass1 +=  trainingSet(0,i)*trainingSet(0,i);
//        //    ++countOfSamplesInClass1;
//        //  }
//        //  else if (response == 0) //TGs are not same
//        //  {
//        //    sumOf1DfeaturesInClass0 +=  trainingSet(0,i);
//        //    sumOfSquaresOf1DfeaturesInClass0 +=  trainingSet(0,i)*trainingSet(0,i);
//        //    ++countOfSamplesInClass0;
//        //  }
//        //  else
//        //  {
//        //    cout << "Error in kNN classification" << endl;
//        //  }
//        ////}
//        //double meanOfClass0 = sumOf1DfeaturesInClass0/countOfSamplesInClass0;
//        //double stdOfClass0 = sqrt( ( sumOfSquaresOf1DfeaturesInClass0 - (sumOf1DfeaturesInClass0)*(sumOf1DfeaturesInClass0)/countOfSamplesInClass0 ) / (countOfSamplesInClass0-1) );
//        //double meanOfClass1 = sumOf1DfeaturesInClass1/countOfSamplesInClass1;
//        //double stdOfClass1 = sqrt( ( sumOfSquaresOf1DfeaturesInClass1 - (sumOf1DfeaturesInClass1)*(sumOf1DfeaturesInClass1)/countOfSamplesInClass1 ) / (countOfSamplesInClass1-1) );
//        //cout << meanOfClass0 << '\t' << stdOfClass0 << '\t' << meanOfClass1 << '\t' << stdOfClass1 << '\t' << endl;
//        //{
//        //  for(int x = meanOfClass0*imgSize/trainingSet.max()+1; x > (meanOfClass0 - stdOfClass0)*imgSize/trainingSet.max()+1; --x)
//        //  {
//        //    resultsOnImg(x,0,0) =  5;
//        //    resultsOnImg(x,0,1) =  5;
//        //  }
//        //  for(int x = meanOfClass1*imgSize/trainingSet.max()+1; x < (meanOfClass1 + stdOfClass1)*imgSize/trainingSet.max()+1; ++x)
//        //  {
//        //    resultsOnImg(x,0,0) = 15;
//        //    resultsOnImg(x,0,1) = 15;
//        //  }
//        //}
//        int numberOfCorrectlyClassifiedAsSameTortoises = 0;
//        int numberOfCorrectlyClassifiedAsDifferentTortoises = 0;
//        int numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises = 0;
//        int numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises = 0;
//        outliers << '\n' << k << " - nearest neigbours" << '\n';
//        //outliers << setw(20) << "meanClass1" << setw(20) << "stdClass1" << setw(20) << "meanClass0" << setw(20) << "stdClass0"<< '\n';
//        //outliers << setw(20) << meanOfClass1 << setw(20) << stdOfClass1 << setw(20) << meanOfClass0 << setw(20) << stdOfClass0 << '\n';
//        outliers << "TrainingSet Classification" << '\n';
//        outliers << setw(20) <<"TgNumber1" << setw(20) << "TgNumber2" << setw(20) << "sumOfFeatures" << setw(20) << "GroundTruthClass" << setw(20) << "ClassifiedClass" << '\n';
//      
//        for (int i =  0; i < numberOfTrainingData; ++i) // OPTIMISTIC CLASSIFICATION - using training set
//        {
//          cv::Mat sampleMat (1,numberOfDataDimensions,CV_32FC1);
//          for(int j = 0; j < numberOfDataDimensions; ++j)
//          {
//            sampleMat.at<float>(0,j) = trainingSet(j,i);
//          }
//          //int x = trainingSet(0,i)*imgSize/trainingSet.max()+1;
//          //int y = 1;
//          //int y = trainingSet(1,i)*1000/trainingSet.max()+1;
//          float response = knn.find_nearest(sampleMat,k,0,0,0,0);
//          outliers << setw(20) << trainingSetTGnumbers(0,i) << setw(20) << trainingSetTGnumbers(1,i) << setw(20) << trainingSet.get_row(i).sum() << setw(20) << trainingLabels(0,i) << setw(20) << response << '\n';
//
//          if (response == 1 && trainingLabels(0,i)==1)
//          {
//            ++numberOfCorrectlyClassifiedAsSameTortoises;
//            //if(x < imgSize && y < imgSize)
//            //{
//            //  resultsOnImg(x,y,0) =  15;
//            //  resultsOnImg(x,y,1) =  15;
//            //}
//          }
//          else if (response == 0 && trainingLabels(0,i)==0)
//          {
//            ++numberOfCorrectlyClassifiedAsDifferentTortoises;
//            // if(x < imgSize && y < imgSize)
//            // {
//            //   resultsOnImg(x,y,0) =  5;
//            //   resultsOnImg(x,y,1) =  5;
//            //}
//          }
//          else if (response == 1 && trainingLabels(0,i)==0) // false positive
//          {
//            ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises;
//            //if(x < imgSize && y < imgSize)
//            //{
//            //  resultsOnImg(x,y,0) =  5;
//            //  resultsOnImg(x,y,1) =  15;
//            //}
//          }
//          else if (response ==0 && trainingLabels(0,i)==1) //false negative
//          {
//            ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises;
//            //outliers << setw(20) << trainingSetTGnumbers(0,i) << setw(20) << trainingSetTGnumbers(1,i) << setw(20) << " " << setw(20) << (trainingSet(0,i)-meanOfClass0)/stdOfClass0 << '\n';
//            //if(x < imgSize && y < imgSize)
//            //{
//            //  resultsOnImg(x,y,0) =  15;
//            //  resultsOnImg(x,y,1) =  5;
//            //}
//          }
//          else
//          {
//            cout << response << " " << trainingLabels(0,i) << " neco se pokazilo" << endl;
//          }
//        }
//        //if(k==51) resultsOnImg.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\classificationImgwithoutPow2.pnm");
//        kNNClassificationResults << "  " <<  k << " - NN" << '\n';
//        kNNClassificationResults << "OPTIMISTIC CLASSIFICATION ERROR" << '\n';
//        kNNClassificationResults << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalOptimisticSuccess in %" << '\n';
//        kNNClassificationResults << left << setw(30) << numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTrainingData << setw(40) << numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTrainingData;
//        kNNClassificationResults << left << setw(40) << (numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0 << '\n';
//        kNNClassificationResults << left << setw(30) << "incorrectlyClassifiedAsSame" << setw(40) << "incorrectlyClassifiedAsDifferent" << setw(40) << "totalOptimisticError in %" << '\n';
//        kNNClassificationResults << left << setw(30) << numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData << setw(40) << numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData;
//        kNNClassificationResults << left << setw(40) << (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0 << '\n';
//        double optimisticTotalError = (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0;
//        kAndErrors << setw(30) << optimisticTotalError;
//        kNNClassificationResults << '\n';
//      
//        numberOfCorrectlyClassifiedAsSameTortoises = 0;
//        numberOfCorrectlyClassifiedAsDifferentTortoises = 0;
//        numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises = 0;
//        numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises = 0;
//        outliers << '\n' << "TestSet Classification" << '\n';
//        outliers << setw(20) << "TgNumber1" << setw(20) << "TgNumber2" << setw(20) << "sumOfFeatures" << setw(20) << "GroundTruthClass" << setw(20) << "ClassifiedClass" << '\n';
//        for (int i =  0; i < numberOfTestData; ++i) //PESSIMISTIC CLASSIFICATION - using test set
//        {
//          cv::Mat sampleMat (1,numberOfDataDimensions,CV_32FC1);
//          for(int j = 0; j < numberOfDataDimensions; ++j)
//          {
//            sampleMat.at<float>(0,j) = testSet(j,i);
//          }
//          //int x = testSet(0,i)*imgSize/testSet.max()+1;
//          //int y = 2;
//          float response = knn.find_nearest(sampleMat,k,0,0,0,0);
//          outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << testSet.get_row(i).sum() << setw(20) << testLabels(0,i) << setw(20) << response << '\n';
//          if (response == 1 && testLabels(0,i)==1)
//          {
//            ++numberOfCorrectlyClassifiedAsSameTortoises;
//            //if(x < imgSize && y < imgSize)
//            //{
//            //  resultsOnImg(x,y,0) =  15;
//            //  resultsOnImg(x,y,1) =  15;
//            //}
//          }
//          else if (response ==0 && testLabels(0,i)==0)
//          {
//            ++numberOfCorrectlyClassifiedAsDifferentTortoises;
//            //if(x < imgSize && y < imgSize)
//            //{
//            //  resultsOnImg(x,y,0) =  5;
//            //  resultsOnImg(x,y,1) =  5;
//            //}
//          }
//          else if (response == 1 && testLabels(0,i)==0)
//          {
//            ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises;
//            //outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << 1 << setw(20) << " " << '\n';
//            //if(x < imgSize && y < imgSize)
//            //{
//            //  resultsOnImg(x,y,0) =  5;
//            //  resultsOnImg(x,y,1) =  15;
//            //}
//          }
//          else if (response == 0 && testLabels(0,i)==1)
//          {
//            ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises;
//            //outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << " " << setw(20) << 1 << '\n';
//            //if(x < imgSize && y < imgSize)
//            //{
//            //  resultsOnImg(x,y,0) =  15;
//            //  resultsOnImg(x,y,1) =  5;
//            //}
//          }
//          else
//          {
//            cout << response << " " << testLabels(0,i) << " neco se pokazilo" << endl;
//          }
//        }
//        //if(k%10==1) resultsOnImg.display();
//
//        kNNClassificationResults << "PESSIMISTIC CLASSIFICATION ERROR" << '\n';
//        kNNClassificationResults << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalPesimisticSuccess in %" << '\n';
//        kNNClassificationResults << left << setw(30) << numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTestData << setw(40) << numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTestData;
//        kNNClassificationResults << left << setw(40) << (numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTestData + numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0 << '\n';
//        kNNClassificationResults << left << setw(30) << "incorrectlyClassifiedAsSame" << setw(40) << "incorrectlyClassifiedAsDifferent" << setw(40) << "totalPesimisticError in %" << '\n';
//        kNNClassificationResults << left << setw(30) << numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData << setw(40) << numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData;
//        kNNClassificationResults << left << setw(40) << (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0 << '\n';
//        double pessimisticTotalError = (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0;
//        kAndErrors << setw(30) << pessimisticTotalError <<  setw(30) << (optimisticTotalError + pessimisticTotalError)/2 <<  setw(30) << (optimisticTotalError + pessimisticTotalError)*pow(pessimisticTotalError,2.0)/2 << endl;
//        kNNClassificationResults << '\n' << endl;
//        outliers << endl;
//      }
//    }
//  }
//
//