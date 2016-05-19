
using namespace cv;

struct Rpoint
{ 
	int dx;
	int dy;
	float phi;
};

struct Rpoint2
{
	float x;
	float y;
	int phiindex;
};

class GenHoughTrnf
{
private:
  //template file names
  string m_strTmplContourFileName;
  string m_strTmplImgFileName;
  string m_strImgFileName;
  string m_strEdgeImgFileName;
	// accumulator matrix
	CImg<int> m_accum;
  CImg<int> m_accumPara;
	// accumulator matrix
	Mat showimage;
	// contour points:
	std::vector<Rpoint> pts;
	// reference point (inside contour(
	int refPointX;
  int refPointY;
	// R-table of template object:
	std::vector<std::vector<Vec2i>> Rtable;
	// number of m_nIntervals for angles of R-table:
	int m_nIntervals;
  float m_fWidthOfAPhiInterval;
	// width of template contour
	int wtemplate;
  int htemplate;
	// minimum and maximum width of scaled contour
	int m_nTmplWidthMin;
	int m_nTmplWidthMax;
	// dimension in pixels of squares in image
	int m_nRangeXY;
	// interval to increase scale
	int m_nRangeS;	
  // accumulator limiting conditions
  float m_fLeftAccumBndryXCoor;
  float m_fRightAccumBndryXCoor;

public:

  void Init(int nTmplWidthMin, int nTmplWidthMax, int nRangeS, int nRangeXY, int nIntervals, float fLeftAccumBndryXCoor, float fRightAccumBndryXCoor)
  {
    m_nTmplWidthMin = nTmplWidthMin;
    m_nTmplWidthMax = nTmplWidthMax;
    m_nRangeS = nRangeS;
    m_nRangeXY = nRangeXY;
    m_nIntervals = nIntervals;
    m_fLeftAccumBndryXCoor = (fLeftAccumBndryXCoor/(m_nRangeXY))-1;
    m_fRightAccumBndryXCoor = (fRightAccumBndryXCoor/(m_nRangeXY))+1;
    m_fWidthOfAPhiInterval = 180.0f/m_nIntervals;
  }

	GenHoughTrnf(string strTmplContourFileName, string strTmplImgFileName, string strEdgeImgFileName, string strImgFileName)
  {
    m_strTmplContourFileName = strTmplContourFileName;
    m_strTmplImgFileName     = strTmplImgFileName;
    m_strImgFileName         = strImgFileName;
    m_strEdgeImgFileName     = strEdgeImgFileName;
    Init(150,300,5,4,16,-1.0f,INT_MAX*1.0f);
	}

	// show the best candidate detected on image
  CImg<int> bestCandidate(point2dCoor &sAccumMax, int &nRotationAngleMax, float &fWidthRatioMax, int &nTemplateWidth)
  {
    CImg<unsigned char> edges(m_strImgFileName.c_str());
		double maxval;
    int S = static_cast<int>( ceil(static_cast<float>(m_nTmplWidthMax-m_nTmplWidthMin)/m_nRangeS+1.0f) );
    int ch = static_cast<int>( m_accum.get_stats()(11) );
    int nMaxXCoor = static_cast<int>(m_accum.get_stats()(8));
    int nMaxYCoor = static_cast<int>(m_accum.get_stats()(9));
    int r =  m_accumPara(0,ch);
    int s =  m_accumPara(1,ch);
    int id_max[4] = { nMaxXCoor, nMaxYCoor, s, r};
    //cout << s << " "<< r<< endl;
    maxval = m_accum.max();
    //m_accum.get_channel(ch).display();
		int nl= edges.height();
		int nc= edges.width(); 
		Vec2i referenceP = Vec2i(id_max[0]*m_nRangeXY+(m_nRangeXY+1)/2, id_max[1]*m_nRangeXY+(m_nRangeXY+1)/2);
		// rotate and scale points all at once. Then impress them on image
		std::vector<std::vector<Vec2i>> Rtablerotatedscaled(m_nIntervals);
		double cs = cos(r*PI/180.0);
		double sn = sin(r*PI/180.0);
		int w = m_nTmplWidthMin + id_max[2]*m_nRangeS;
		float wratio = (float)w/(wtemplate);
    sAccumMax.xCoor   = referenceP[0];
    sAccumMax.yCoor   = referenceP[1];
    nRotationAngleMax = r;
    fWidthRatioMax    = wratio;
    nTemplateWidth    = w;
    cout << " wratio " <<  wratio <<" angleInDegrees " << r*PI/180.0 << endl;
		for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii){
			for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj){
				int dx = roundToInt(wratio*(cs*Rtable[ii][jj][0] - sn*Rtable[ii][jj][1]));
				int dy = roundToInt(wratio*(sn*Rtable[ii][jj][0] + cs*Rtable[ii][jj][1]));
				int x = referenceP[0] - dx;
				int y = referenceP[1] - dy;
				//Rtablerotatedscaled[ii].push_back(Vec2i( dx, dy));
				if ( (x<nc-1)&&(y<nl-1)&&(x>0)&&(y>0) ){
          for(int i = -1; i < 2; ++i)
          {
            for(int j = -1; j < 2; ++j)
            {
              	edges(x+i,y+j,0) = 255;
					      edges(x+i,y+j,1) = 0;
					      edges(x+i,y+j,2) = 0;
            }
          }

				}
			}
    }
    //edges.display();
    return edges;
	}

  	// show the best candidate detected on image
  CImg<int> bestCandidateFromGivenPoints(point2dCoor &sAccumMax, int &nRotationAngleMax,
    float &fWidthRatioMax, int &nTemplateWidth, const CImg<int> &imgPossibleJunctionPreciseLocations)
  {    
    int nNumOfPossibleJunctions = imgPossibleJunctionPreciseLocations.height();
    CImg<int> selectedPointsAccum(m_accum.width(),nNumOfPossibleJunctions,1,m_accum.spectrum(),0);
    cout << nNumOfPossibleJunctions << endl;
    //accum preprocessing
    for(int ch = 0; ch < m_accum.spectrum(); ++ch)
    {
      for(int jj = 0; jj < nNumOfPossibleJunctions; ++jj)
      {
        int nYJunctionCoor = imgPossibleJunctionPreciseLocations(jj)/(m_nRangeXY);
        for(int ii = 0; ii<m_accum.width();++ii)
        {
          selectedPointsAccum(ii,jj,ch) = m_accum(ii,nYJunctionCoor,ch);
        }
      }
    }
    //selectedPointsAccum.display();
    CImg<unsigned char> edges(m_strImgFileName.c_str());
		double maxval;
    int S = static_cast<int>( ceil(static_cast<float>(m_nTmplWidthMax-m_nTmplWidthMin)/m_nRangeS+1.0f) );
    //get coordinates of the accum minimum 
    int ch = static_cast<int>( m_accum.get_stats()(11) );
    int nMaxXCoor = static_cast<int>(selectedPointsAccum.get_stats()(8));
    int nMaxYCoor = imgPossibleJunctionPreciseLocations(static_cast<int>(selectedPointsAccum.get_stats()(9)));
    int r =  m_accumPara(0,ch);
    int s =  m_accumPara(1,ch);
    int id_max[4] = { nMaxXCoor, nMaxYCoor, s, r};
    //cout << s << " "<< r<< endl;
    maxval = m_accum.max();
    //m_accum.get_channel(ch).display();
		int nl= edges.height();
		int nc= edges.width(); 
		Vec2i referenceP = Vec2i(id_max[0]*m_nRangeXY+(m_nRangeXY+1)/2, id_max[1]);
		// rotate and scale points all at once. Then impress them on image
		std::vector<std::vector<Vec2i>> Rtablerotatedscaled(m_nIntervals);
		double cs = cos(r*PI/180.0);
		double sn = sin(r*PI/180.0);
		int w = m_nTmplWidthMin + id_max[2]*m_nRangeS;
		float wratio = (float)w/(wtemplate);
    sAccumMax.xCoor   = referenceP[0];
    sAccumMax.yCoor   = referenceP[1];
    nRotationAngleMax = r;
    fWidthRatioMax    = wratio;
    nTemplateWidth    = w;
    //cout << " wratio " <<  wratio <<" angleInDegrees " << r*PI/180.0 << endl;
		for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii){
			for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj){
				int dx = roundToInt(wratio*(cs*Rtable[ii][jj][0] - sn*Rtable[ii][jj][1]));
				int dy = roundToInt(wratio*(sn*Rtable[ii][jj][0] + cs*Rtable[ii][jj][1]));
				int x = referenceP[0] - dx;
				int y = referenceP[1] - dy;
				//Rtablerotatedscaled[ii].push_back(Vec2i( dx, dy));
				if ( (x<nc-1)&&(y<nl-1)&&(x>0)&&(y>0) ){
          for(int i = -1; i < 2; ++i)
          {
            for(int j = -1; j < 2; ++j)
            {
              	edges(x+i,y+j,0) = 255;
					      edges(x+i,y+j,1) = 0;
					      edges(x+i,y+j,2) = 0;
            }
          }

				}
			}
    }
    //edges.display();
    return edges;
	}

  void runGHT(bool bMoreAngles)
  {
		createRtable();
		accumulate(bMoreAngles);
  }

private:
  
  void createRtable()
  {
		// code can be improved reading a pre-saved Rtable
		readPoints();
		readRtable();
    //visuRtable();
	}

	// fill accumulator matrix
	void accumulate(bool bMoreAngles)
  {
    //cout <<    m_nTmplWidthMin  << " " <<  m_nTmplWidthMax << " " <<  m_nRangeS << " " <<  m_nRangeXY  << " " <<   m_nIntervals  << " " 
    //  <<   m_fLeftAccumBndryXCoor << " " <<   m_fRightAccumBndryXCoor << " " << m_fWidthOfAPhiInterval << endl;
    //m_nTmplWidthMin = wtemplate;
    Mat input_img = imread(m_strImgFileName.c_str(), 1);
    Mat detected_edges = imread(m_strEdgeImgFileName.c_str(),1);
    CImg<int> edges(m_strEdgeImgFileName.c_str());
    //edges.display();
		// transform image to grayscale:
		Mat src_gray;
		src_gray.create( Size(input_img.cols, input_img.rows), CV_8UC1);
		cvtColor(input_img, src_gray, CV_BGR2GRAY); 
		// get Scharr matrices from image to obtain contour gradients
		Mat dx;
		dx.create( Size(input_img.cols, input_img.rows), CV_16SC1);
		Sobel(src_gray, dx, CV_16S, 1, 0, CV_SCHARR);
		Mat dy;
		dy.create( Size(input_img.cols, input_img.rows), CV_16SC1);
		Sobel(src_gray, dy, CV_16S, 0, 1, CV_SCHARR);
    //cout << dx << endl;
    //waitKey(0);
		// load all points from image all image contours on vector pts2
		int nl= detected_edges.rows;
		int nc= detected_edges.cols ;
		float inv_m_nRangeXY = (float)1/m_nRangeXY;
		std::vector<Rpoint2> pts2;
		for (int j=0; j<nl; ++j) {
			for (int i=0; i<nc; ++i) {  		
        if ( edges(i,j) == 255  ) // consider only white points (contour)
				{
					short vx = dx.at<short>(j,i);
					short vy = dy.at<short>(j,i);
					Rpoint2 rpt;
					rpt.x = i*inv_m_nRangeXY;
					rpt.y = j*inv_m_nRangeXY;
					float a = static_cast<float>(atan2((float)vy, (float)vx)*180.0f/PI);
					float phi = ((a > 0) ? a : a + 180.0f );    
					int angleindex = (int)(phi/m_fWidthOfAPhiInterval); 
					if (angleindex == m_nIntervals) angleindex=m_nIntervals-1;
					rpt.phiindex = angleindex;
					pts2.push_back( rpt );
				}
			}
		}
    //CImg<int> reconstImg(nc,nl,1,1,0);
    //for (vector<Rpoint2>::size_type t = 0; t < pts2.size(); ++t){
    //  reconstImg(pts2[t].x,pts2[t].y) = 255;
    //}
    //cout << pts2.size() << " " << Rtable.size() << endl;
    //reconstImg.display();
    //cout << "accumulate - pts2 loaded." << endl;
		// OpenCv 4-dimensional matrix definition and in general a useful way for defining multidimensional arrays and vectors in c++
		// create accumulator matrix
		int X = static_cast<int>(ceil((float)nc/m_nRangeXY));
		int Y = static_cast<int>(ceil((float)nl/m_nRangeXY));
		int S = static_cast<int>(ceil((float)(m_nTmplWidthMax-m_nTmplWidthMin)/m_nRangeS+1.0f));
    int R = 1;
		if(bMoreAngles) R = 360; // rotate 0->360
    else R = 1; // no rotation
		// icrease accum cells with hits corresponding with slope in Rtable vector rotatated and scaled
		float inv_wtemplate_m_nRangeXY = (float)1/(wtemplate*m_nRangeXY);
		// rotate RTable from minimum to maximum angle
    int sliceNum = 0;
		for (int r=0; r<R; ++r) 
    {  // rotation
      //cout << r << endl;
			std::vector<std::vector<Vec2f>> Rtablerotated(m_nIntervals);
			// cos and sin are computed in the outer loop to reach computational efficiency
      if(isGoodAngle(r))
      {
			  float cs = static_cast<float>(cos(r*PI/180.0));
			  float sn = static_cast<float>(sin(r*PI/180.0));
			  for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii)
        {
				  for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj)
          {
					  int iimod = (int)(ii+r/m_fWidthOfAPhiInterval) % m_nIntervals; // ii is index of phi
					  Rtablerotated[iimod].push_back(Vec2f(cs*Rtable[ii][jj][0] - sn*Rtable[ii][jj][1], sn*Rtable[ii][jj][0] + cs*Rtable[ii][jj][1]));
				  }
			  }
			  // scale the rotated RTable from minimum to maximum scale
			  for (int s=0; s<S; ++s) 
        {  // scale
          //cout << s << endl;
          CImg<int> accumSlice(X,Y,1,1, 0);
				  std::vector<std::vector<Vec2f>> Rtablescaled(m_nIntervals);
				  int w = m_nTmplWidthMin + s*m_nRangeS;
          int h = static_cast<int>(htemplate*(1.0*w/wtemplate)/m_nRangeXY);  
          //cout << h << " " <<  nl << endl;
				  float wratio = (float)w*inv_wtemplate_m_nRangeXY;	
          //cout << "sliceNum " << sliceNum << " wratio " <<  wratio <<" angleInDegrees " <<angleInDegrees<< endl;
				  for (std::vector<std::vector<Vec2f>>::size_type ii = 0; ii < Rtablerotated.size(); ++ii){
					  for (std::vector<Vec2f>::size_type jj= 0; jj < Rtablerotated[ii].size(); ++jj){
						  Rtablescaled[ii].push_back(Vec2f(wratio*Rtablerotated[ii][jj][0], wratio*Rtablerotated[ii][jj][1]));	
					  }
				  }
				  // iterate through each point of edges and hit corresponding cells from rotated and scaled Rtable
				  for (vector<Rpoint2>::size_type t = 0; t < pts2.size(); ++t)
          { // XY plane		
					  int angleindex = pts2[t].phiindex;
					  for (std::vector<Vec2f>::size_type index = 0; index < Rtablescaled[angleindex].size(); ++index){
						  float deltax = Rtablescaled[angleindex][index][0];
						  //float deltax = Rtable[angleindex][index][0];
						  float deltay = Rtablescaled[angleindex][index][1];							
						  //float deltay = Rtable[angleindex][index][1];							
						  int xcell = (int)(pts2[t].x + deltax);
						  int ycell = (int)(pts2[t].y + deltay);
              if(bMoreAngles)
              {
                if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1)&&(xcell<m_fRightAccumBndryXCoor)&&(xcell>m_fLeftAccumBndryXCoor) //)
                  &&(ycell>=h*5/12)&&(ycell<=nl/m_nRangeXY-h*5/12)&&(xcell>=w*5/(2*6*m_nRangeXY))&&(xcell<=nc/m_nRangeXY-w*5/(2*6*m_nRangeXY)) )
                //if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1) )
                {
                  accumSlice(xcell,ycell)++;
						    }
              }
              else
              {
                if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1)&&(xcell<m_fRightAccumBndryXCoor)&&(xcell>m_fLeftAccumBndryXCoor))
                //if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1) )
                {
                  accumSlice(xcell,ycell)++;
						    }
              }
					  }
				  }
          CImg<int> accumSlicePara(2,1,1,1, r,s);
          m_accumPara.append(accumSlicePara,'y');
          m_accum.append(accumSlice,'c');
          ++sliceNum;
          //cout << s << endl;
			  }
        cout << r << endl;
      }
		}
	}
	// load vector pts with all points from the contour
	void readPoints()
  {
		// read original template image and its worked-out contour
    CImg<unsigned char> template_img (m_strTmplContourFileName.c_str());
		// find reference point inside contour image and save it in variable refPoint
		int nl= template_img.height();
		int nc= template_img.width(); 
		for (int j=0; j<nl; ++j) 
    {
			for (int i=0; i<nc; ++i) 
      {  		
				if ( template_img(i,j) == 127 )
        {
					refPointX = i;
          refPointY = j;
				}
			}
		}
    cout << template_img.get_threshold(255).sum() << endl;
		// get Scharr matrices from original template image to obtain contour gradients
    Mat original_img = imread(m_strTmplImgFileName.c_str(), 1);
		Mat input_img_gray;
		input_img_gray.create( Size(original_img.cols, original_img.rows), CV_8UC1);
		cvtColor(original_img, input_img_gray, CV_BGR2GRAY); 
		Mat dx;
		dx.create( Size(original_img.cols, original_img.rows), CV_16SC1);
		Sobel(input_img_gray, dx, CV_16S, 1, 0, CV_SCHARR);
		Mat dy;
		dy.create( Size(original_img.cols, original_img.rows), CV_16SC1);
		Sobel(input_img_gray, dy, CV_16S, 0, 1, CV_SCHARR);
		// load points on vector
		pts.clear();
		int mindx = INT_MAX;
		int maxdx = INT_MIN;
    int mindy = INT_MAX;
    int maxdy = INT_MIN; 
		for (int j=0; j<nl; ++j) 
    {
			for (int i=0; i<nc; ++i)
      {  		
				if ( template_img(i,j) == 255  ) 
				{
					short vx = dx.at<short>(j,i);
					short vy = dy.at<short>(j,i);
					Rpoint rpt;
					rpt.dx = refPointX-i;
					rpt.dy = refPointY-j;
					float a = atan2((float)vy, (float)vx)*180.0f/fPI; // edge direction (with sense) in radians falling in [-180,+180]
          // we need only the direction or the tangent of a given edge point (horizontal,vertival,diagonal,...)
					rpt.phi = ((a > 0) ? a : a + 180.0f ); // 180 is added to all negative directions
					if (rpt.dx < mindx) mindx=rpt.dx;
					if (rpt.dx > maxdx) maxdx=rpt.dx;
					if (rpt.dy < mindy) mindy=rpt.dy;
					if (rpt.dy > maxdy) maxdy=rpt.dy;
					pts.push_back( rpt );
				}
			}
		}
		// maximum width and height of the contour
    htemplate = maxdy-mindy+1;
		wtemplate = maxdx-mindx+1;
	}

	// create Rtable from contour points
	void readRtable()
  {
		Rtable.clear();
		Rtable.resize(m_nIntervals);
		// put points in the right interval, according to discretized angle and range size 
		for (vector<Rpoint>::size_type t = 0; t < pts.size(); ++t)
    {
      // an interval is assigned to the phi angles
			int angleindex = (int)(pts[t].phi/m_fWidthOfAPhiInterval); 
      // if angleindex == m_nIntervals, than phi==180 would form a separate interval, can joint the the m_nIntervals-1
			if (angleindex == m_nIntervals) angleindex=m_nIntervals-1;
			Rtable[angleindex].push_back( Vec2i(pts[t].dx, pts[t].dy) );
		}
	}

	// create Rtable from contour points
	void visuRtable()
  {
    CImg<int> contour(m_strTmplContourFileName.c_str()),RTabelRes(contour.width(),contour.height(),1,1, 0);
    int nl= contour.height();
	  int nc= contour.width();  
		for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii)
    {
			for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj)
      {
				int dx = Rtable[ii][jj][0];
				int dy = Rtable[ii][jj][1];
				int x = refPointX - dx;
				int y = refPointY - dy;
				//Rtablerotatedscaled[ii].push_back(Vec2i( dx, dy));
				if ( (x<nc)&&(y<nl)&&(x>-1)&&(y>-1) ){
					RTabelRes(x, y) = 255;
				}
			}
		}
    (RTabelRes, contour).display("RTable & contour");
	}

	inline int roundToInt(double num) {
		return (num > 0.0) ? (int)(num + 0.5) : (int)(num - 0.5);
	}

  bool isGoodAngle(int angle)
  {
    if( (angle >= 355 || angle < 5) && !(angle%2)) return 1;
    //else if(angle >= -104 && angle <= -74 ) return 1;
    //else if(angle > -5 && angle < 5 && !(angle%2)) return 1;
    //else if(angle >=  74 && angle <= 105 ) return 1;
    else if(angle >  175 && angle <= 185 && !(angle%2)) return 1;
    else return 0;
  }
};