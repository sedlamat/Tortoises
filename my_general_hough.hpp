/**
  my_general_hough.hpp

  Implementation of the generalized Hough transform, using OpenCV 
  library.    

  @author Matej Sedlacek
  @version 0.0
*/

#ifndef _MY_GENERAL_HOUGH_HPP_
#define _MY_GENERAL_HOUGH_HPP_

/* STANDARD C++ LIBRARIES */
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

/* THIRD PARTY LIBRARIES */
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"

namespace my 
{		

   /**
    Gets the hough points, i.e., a vector (indexed by quantized 
    gradient orientations - 0 to 3) of vector of points (the 
    coordinates of edges).    


    @param src - A color image.
    @param src_edges -  An edge image of the color image.
    @return Vector of vector of points generated from src and src_edges.
  */

  std::vector<std::vector<cv::Point_<int> > >  get_hough_points(
					       const cv::Mat& src, 
					       const cv::Mat& src_edges)
  {
    cv::Mat orient = my::get_gradient_orientation(src);
    std::vector<std::vector<cv::Point_<int> > > hough_points(4);

    // orientations in [0, 360) => make [180, 360) to [0, 180), only
    // direction needed (45deg is the same direction as 225deg)
    // this prevents 45deg being 225deg with inverse intensity.
    // inverse intesity changes orientation, do not want that.
    // direction changes only if the edge changes.
    
    cv::Mat orient_adjust = (orient >= 180)/255;
    orient_adjust.convertTo(orient_adjust,orient.depth());
    orient += orient_adjust * -180;
    
    // points are quantized from [0,180) into 0,1,2,3 indices
    for(int yy = 0; yy < src_edges.rows; yy++) {
      const uchar *ptr_src_edges_irow = src_edges.ptr<uchar>(yy);
      for(int xx = 0; xx < src_edges.cols; xx++) {
	if (ptr_src_edges_irow[xx]) {
	  cv::Point_<int> pt(xx, yy);
	  float phi = orient.at<float>(pt);
	  int orient_quant_index;
	  if (phi >= 157.5 || phi < 22.5) orient_quant_index = 0;
	  else if (phi >= 22.5 && phi < 67.5) orient_quant_index = 1;
	  else if (phi >= 67.5 && phi < 112.5) orient_quant_index = 2;
	  else if (phi >= 112.5 && phi < 157.5) orient_quant_index = 3;
	  else CV_Assert( !"Undefined orientation" );
	  hough_points[orient_quant_index].push_back(pt);
	}
      }
    }
    return hough_points;
  }  
  
  /**
    Gets the hough r-table, i.e., the hough points minus reference 
    point of the template.    

    @param src_templ - A template image for the hough tranform.
    @param src_edges - An edge image of the template image.
    @param ref_point - Reference point of the template.
    @return The r-table (vector of vector of [points - ref_point]).
  */

  std::vector<std::vector<cv::Point_<int> > >  get_r_table(
				      const cv::Mat& src_templ, 
				      const cv::Mat& src_templ_edges,
				      const cv::Point_<int>& ref_point)
  {
    std::vector<std::vector<cv::Point_<int> > >  
	  r_table = get_hough_points(src_templ, src_templ_edges);
        
    for(auto & points : r_table) {
      for(auto & point : points) {
	point = ref_point - point;
      }
    }
    return r_table;
  }
 

   /**
    Gets a filled accumulator for one point. 

    @param ref_points - Possible reference points for a given src_pt.
    @param src_pt - One feature point in the src_image.
    @param size - Size of the src_image/accumulator.
    @return One-layer accumulator for ref_points and src_pt.
  */

  cv::Mat get_accum_layer(std::vector<cv::Point_<int> > ref_points,
			  cv::Point_<int> src_pt,
			  cv::Size_<int> size)
  {
    cv::Mat accum_layer(size, CV_32FC1, cv::Scalar_<float>(0.0));
    for(auto const& ref_pt : ref_points) {
      cv::Point_<double> pt1, pt2;
      double a, b;
      if (src_pt.x == ref_pt.x) {
	pt1.x = pt2.x = src_pt.x;
	pt1.y = 0;
	pt2.y = size.height;
      }
      else {
	a = (src_pt.y - ref_pt.y)/(src_pt.x - ref_pt.x);
	b = (ref_pt.y*src_pt.x - src_pt.y*ref_pt.x)/
						(src_pt.x - ref_pt.x);
	if (std::abs(a) < 1) { // then use y = a*x + b
	  pt1.x = 0;
	  pt1.y = a*pt1.x + b;
	  pt2.x = size.width;
	  pt2.y = a*pt2.x + b;
	}
	else { // then use x = (y-b)/a
	  pt1.y = 0;
	  pt1.x = (pt1.y - b) / a;
	  pt2.y = size.height;
	  pt2.x = (pt2.y - b) / a;
	}	 
      } 
      cv::line(accum_layer, pt1, pt2, cv::Scalar_<float>(1.0), 1);
    }
    return accum_layer;
  } 
  
  
  /**
    Gets the hough accumulator source feature points and reference
    points from the hough r-table.    

    @param r_table - Reference point shifts grouped by their gradient
		     orientations.
    @param src_hough_points - Feature points grouped by their gradient
			      orientations.
    @param size - Size of the accumulator/source_image.
    @return The r-table (vector of vector of [points - ref_point]).
  */

  cv::Mat get_accumulator(
	std::vector<std::vector<cv::Point_<int> > > r_table,
	std::vector<std::vector<cv::Point_<int> > > src_hough_points,
	cv::Size_<int> size)
  {
    cv::Mat accum(size, CV_32FC1, cv::Scalar_<float>(0.0));
    for(auto const& quant_idx : {0,1,2,3}) {
      for(auto & src_pt : src_hough_points[quant_idx]) {
	std::vector<cv::Point_<int> > ref_pts;
	for(auto & pt_diff : r_table[quant_idx]) {
	  ref_pts.push_back(pt_diff + src_pt);
	}
	cv::Mat accum_layer = my::get_accum_layer(ref_pts, src_pt,
						  size);
	//std::cout << my::maxMat(accum_layer) << std::endl;
	accum += accum_layer;
      }
      //return accum;
    }
    return accum;
  }
  
  /**
    General Hough transform.    

    @param src_templ - A template image for the hough tranform.
    @param src_templ_edges - Edges for the template image.
    @param ref_point - Reference point of the template.
    @param src - Source image.
    @param src_edges - Edges of the source image.
    @return SO FAR: Displays the hough accumulator.
  */
  void general_hough(const cv::Mat& src_templ, 
		     const cv::Mat& src_templ_edges,
		     const cv::Point_<int>& ref_point,
		     const cv::Mat& src, 
		     const cv::Mat& src_edges)
  {
    
    std::vector<std::vector<cv::Point_<int> > > r_table = 
		my::get_r_table(src_templ, src_templ_edges, ref_point);
	      
    std::vector<std::vector<cv::Point_<int> > > src_hough_points = 
		my::get_hough_points(src, src_edges);
		
    cv::Size_<int> size = src.size();
    
    cv::Mat accumulator = my::get_accumulator(r_table, src_hough_points, 
					      size);
    my::display(accumulator);
    my::display(src);
  }
				      
  
	  
} /* namespace my */


#endif /* _MY_GENERAL_HOUGH_HPP_ */

//~ 
//~ struct Rpoint
//~ { 
	//~ int dx;
	//~ int dy;
	//~ float phi;
//~ };
//~ 
//~ struct Rpoint2
//~ {
	//~ float x;
	//~ float y;
	//~ int phiindex;
//~ };
//~ 
//~ class generalHough
//~ {
//~ private:
  //~ //template file names
  //~ string m_strTmplContourFileName;
  //~ string m_strTmplImgFileName;
  //~ string m_strImgFileName;
  //~ string m_strEdgeImgFileName;
	//~ // accumulator matrix
	//~ CImg<int> m_accum;
  //~ CImg<int> m_accumPara;
	//~ // accumulator matrix
	//~ Mat showimage;
	//~ // contour points:
	//~ std::vector<Rpoint> pts;
	//~ // reference point (inside contour(
	//~ int refPointX;
  //~ int refPointY;
	//~ // R-table of template object:
	//~ std::vector<std::vector<Vec2i>> Rtable;
	//~ // number of m_nIntervals for angles of R-table:
	//~ int m_nIntervals;
  //~ float m_fWidthOfAPhiInterval;
	//~ // width of template contour
	//~ int wtemplate;
  //~ int htemplate;
	//~ // minimum and maximum width of scaled contour
	//~ int m_nTmplWidthMin;
	//~ int m_nTmplWidthMax;
	//~ // dimension in pixels of squares in image
	//~ int m_nRangeXY;
	//~ // interval to increase scale
	//~ int m_nRangeS;	
  //~ // accumulator limiting conditions
  //~ float m_fLeftAccumBndryXCoor;
  //~ float m_fRightAccumBndryXCoor;
//~ 
//~ public:
//~ 
  //~ void Init(int nTmplWidthMin, int nTmplWidthMax, int nRangeS, int nRangeXY, int nIntervals, float fLeftAccumBndryXCoor, float fRightAccumBndryXCoor)
  //~ {
    //~ m_nTmplWidthMin = nTmplWidthMin;
    //~ m_nTmplWidthMax = nTmplWidthMax;
    //~ m_nRangeS = nRangeS;
    //~ m_nRangeXY = nRangeXY;
    //~ m_nIntervals = nIntervals;
    //~ m_fLeftAccumBndryXCoor = (fLeftAccumBndryXCoor/(m_nRangeXY))-1;
    //~ m_fRightAccumBndryXCoor = (fRightAccumBndryXCoor/(m_nRangeXY))+1;
    //~ m_fWidthOfAPhiInterval = 180.0f/m_nIntervals;
  //~ }
//~ 
	//~ GenHoughTrnf(string strTmplContourFileName, string strTmplImgFileName, string strEdgeImgFileName, string strImgFileName)
  //~ {
    //~ m_strTmplContourFileName = strTmplContourFileName;
    //~ m_strTmplImgFileName     = strTmplImgFileName;
    //~ m_strImgFileName         = strImgFileName;
    //~ m_strEdgeImgFileName     = strEdgeImgFileName;
    //~ Init(150,300,5,4,16,-1.0f,INT_MAX*1.0f);
	//~ }
//~ 
	//~ // show the best candidate detected on image
  //~ CImg<int> bestCandidate(point2dCoor &sAccumMax, int &nRotationAngleMax, float &fWidthRatioMax, int &nTemplateWidth)
  //~ {
    //~ CImg<unsigned char> edges(m_strImgFileName.c_str());
		//~ double maxval;
    //~ int S = static_cast<int>( ceil(static_cast<float>(m_nTmplWidthMax-m_nTmplWidthMin)/m_nRangeS+1.0f) );
    //~ int ch = static_cast<int>( m_accum.get_stats()(11) );
    //~ int nMaxXCoor = static_cast<int>(m_accum.get_stats()(8));
    //~ int nMaxYCoor = static_cast<int>(m_accum.get_stats()(9));
    //~ int r =  m_accumPara(0,ch);
    //~ int s =  m_accumPara(1,ch);
    //~ int id_max[4] = { nMaxXCoor, nMaxYCoor, s, r};
    //~ //cout << s << " "<< r<< endl;
    //~ maxval = m_accum.max();
    //~ //m_accum.get_channel(ch).display();
		//~ int nl= edges.height();
		//~ int nc= edges.width(); 
		//~ Vec2i referenceP = Vec2i(id_max[0]*m_nRangeXY+(m_nRangeXY+1)/2, id_max[1]*m_nRangeXY+(m_nRangeXY+1)/2);
		//~ // rotate and scale points all at once. Then impress them on image
		//~ std::vector<std::vector<Vec2i>> Rtablerotatedscaled(m_nIntervals);
		//~ double cs = cos(r*PI/180.0);
		//~ double sn = sin(r*PI/180.0);
		//~ int w = m_nTmplWidthMin + id_max[2]*m_nRangeS;
		//~ float wratio = (float)w/(wtemplate);
    //~ sAccumMax.xCoor   = referenceP[0];
    //~ sAccumMax.yCoor   = referenceP[1];
    //~ nRotationAngleMax = r;
    //~ fWidthRatioMax    = wratio;
    //~ nTemplateWidth    = w;
    //~ cout << " wratio " <<  wratio <<" angleInDegrees " << r*PI/180.0 << endl;
		//~ for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii){
			//~ for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj){
				//~ int dx = roundToInt(wratio*(cs*Rtable[ii][jj][0] - sn*Rtable[ii][jj][1]));
				//~ int dy = roundToInt(wratio*(sn*Rtable[ii][jj][0] + cs*Rtable[ii][jj][1]));
				//~ int x = referenceP[0] - dx;
				//~ int y = referenceP[1] - dy;
				//~ //Rtablerotatedscaled[ii].push_back(Vec2i( dx, dy));
				//~ if ( (x<nc-1)&&(y<nl-1)&&(x>0)&&(y>0) ){
          //~ for(int i = -1; i < 2; ++i)
          //~ {
            //~ for(int j = -1; j < 2; ++j)
            //~ {
              	//~ edges(x+i,y+j,0) = 255;
					      //~ edges(x+i,y+j,1) = 0;
					      //~ edges(x+i,y+j,2) = 0;
            //~ }
          //~ }
//~ 
				//~ }
			//~ }
    //~ }
    //~ //edges.display();
    //~ return edges;
	//~ }
//~ 
  	//~ // show the best candidate detected on image
  //~ CImg<int> bestCandidateFromGivenPoints(point2dCoor &sAccumMax, int &nRotationAngleMax,
    //~ float &fWidthRatioMax, int &nTemplateWidth, const CImg<int> &imgPossibleJunctionPreciseLocations)
  //~ {    
    //~ int nNumOfPossibleJunctions = imgPossibleJunctionPreciseLocations.height();
    //~ CImg<int> selectedPointsAccum(m_accum.width(),nNumOfPossibleJunctions,1,m_accum.spectrum(),0);
    //~ cout << nNumOfPossibleJunctions << endl;
    //~ //accum preprocessing
    //~ for(int ch = 0; ch < m_accum.spectrum(); ++ch)
    //~ {
      //~ for(int jj = 0; jj < nNumOfPossibleJunctions; ++jj)
      //~ {
        //~ int nYJunctionCoor = imgPossibleJunctionPreciseLocations(jj)/(m_nRangeXY);
        //~ for(int ii = 0; ii<m_accum.width();++ii)
        //~ {
          //~ selectedPointsAccum(ii,jj,ch) = m_accum(ii,nYJunctionCoor,ch);
        //~ }
      //~ }
    //~ }
    //~ //selectedPointsAccum.display();
    //~ CImg<unsigned char> edges(m_strImgFileName.c_str());
		//~ double maxval;
    //~ int S = static_cast<int>( ceil(static_cast<float>(m_nTmplWidthMax-m_nTmplWidthMin)/m_nRangeS+1.0f) );
    //~ //get coordinates of the accum minimum 
    //~ int ch = static_cast<int>( m_accum.get_stats()(11) );
    //~ int nMaxXCoor = static_cast<int>(selectedPointsAccum.get_stats()(8));
    //~ int nMaxYCoor = imgPossibleJunctionPreciseLocations(static_cast<int>(selectedPointsAccum.get_stats()(9)));
    //~ int r =  m_accumPara(0,ch);
    //~ int s =  m_accumPara(1,ch);
    //~ int id_max[4] = { nMaxXCoor, nMaxYCoor, s, r};
    //~ //cout << s << " "<< r<< endl;
    //~ maxval = m_accum.max();
    //~ //m_accum.get_channel(ch).display();
		//~ int nl= edges.height();
		//~ int nc= edges.width(); 
		//~ Vec2i referenceP = Vec2i(id_max[0]*m_nRangeXY+(m_nRangeXY+1)/2, id_max[1]);
		//~ // rotate and scale points all at once. Then impress them on image
		//~ std::vector<std::vector<Vec2i>> Rtablerotatedscaled(m_nIntervals);
		//~ double cs = cos(r*PI/180.0);
		//~ double sn = sin(r*PI/180.0);
		//~ int w = m_nTmplWidthMin + id_max[2]*m_nRangeS;
		//~ float wratio = (float)w/(wtemplate);
    //~ sAccumMax.xCoor   = referenceP[0];
    //~ sAccumMax.yCoor   = referenceP[1];
    //~ nRotationAngleMax = r;
    //~ fWidthRatioMax    = wratio;
    //~ nTemplateWidth    = w;
    //~ //cout << " wratio " <<  wratio <<" angleInDegrees " << r*PI/180.0 << endl;
		//~ for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii){
			//~ for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj){
				//~ int dx = roundToInt(wratio*(cs*Rtable[ii][jj][0] - sn*Rtable[ii][jj][1]));
				//~ int dy = roundToInt(wratio*(sn*Rtable[ii][jj][0] + cs*Rtable[ii][jj][1]));
				//~ int x = referenceP[0] - dx;
				//~ int y = referenceP[1] - dy;
				//~ //Rtablerotatedscaled[ii].push_back(Vec2i( dx, dy));
				//~ if ( (x<nc-1)&&(y<nl-1)&&(x>0)&&(y>0) ){
          //~ for(int i = -1; i < 2; ++i)
          //~ {
            //~ for(int j = -1; j < 2; ++j)
            //~ {
              	//~ edges(x+i,y+j,0) = 255;
					      //~ edges(x+i,y+j,1) = 0;
					      //~ edges(x+i,y+j,2) = 0;
            //~ }
          //~ }
//~ 
				//~ }
			//~ }
    //~ }
    //~ //edges.display();
    //~ return edges;
	//~ }
//~ 
  //~ void runGHT(bool bMoreAngles)
  //~ {
		//~ createRtable();
		//~ accumulate(bMoreAngles);
  //~ }
//~ 
//~ private:
  //~ 
  //~ void createRtable()
  //~ {
		//~ // code can be improved reading a pre-saved Rtable
		//~ readPoints();
		//~ readRtable();
    //~ //visuRtable();
	//~ }
//~ 
	//~ // fill accumulator matrix
	//~ void accumulate(bool bMoreAngles)
  //~ {
    //~ //cout <<    m_nTmplWidthMin  << " " <<  m_nTmplWidthMax << " " <<  m_nRangeS << " " <<  m_nRangeXY  << " " <<   m_nIntervals  << " " 
    //~ //  <<   m_fLeftAccumBndryXCoor << " " <<   m_fRightAccumBndryXCoor << " " << m_fWidthOfAPhiInterval << endl;
    //~ //m_nTmplWidthMin = wtemplate;
    //~ Mat input_img = imread(m_strImgFileName.c_str(), 1);
    //~ Mat detected_edges = imread(m_strEdgeImgFileName.c_str(),1);
    //~ CImg<int> edges(m_strEdgeImgFileName.c_str());
    //~ //edges.display();
		//~ // transform image to grayscale:
		//~ Mat src_gray;
		//~ src_gray.create( Size(input_img.cols, input_img.rows), CV_8UC1);
		//~ cvtColor(input_img, src_gray, CV_BGR2GRAY); 
		//~ // get Scharr matrices from image to obtain contour gradients
		//~ Mat dx;
		//~ dx.create( Size(input_img.cols, input_img.rows), CV_16SC1);
		//~ Sobel(src_gray, dx, CV_16S, 1, 0, CV_SCHARR);
		//~ Mat dy;
		//~ dy.create( Size(input_img.cols, input_img.rows), CV_16SC1);
		//~ Sobel(src_gray, dy, CV_16S, 0, 1, CV_SCHARR);
    //~ //cout << dx << endl;
    //~ //waitKey(0);
		//~ // load all points from image all image contours on vector pts2
		//~ int nl= detected_edges.rows;
		//~ int nc= detected_edges.cols ;
		//~ float inv_m_nRangeXY = (float)1/m_nRangeXY;
		//~ std::vector<Rpoint2> pts2;
		//~ for (int j=0; j<nl; ++j) {
			//~ for (int i=0; i<nc; ++i) {  		
        //~ if ( edges(i,j) == 255  ) // consider only white points (contour)
				//~ {
					//~ short vx = dx.at<short>(j,i);
					//~ short vy = dy.at<short>(j,i);
					//~ Rpoint2 rpt;
					//~ rpt.x = i*inv_m_nRangeXY;
					//~ rpt.y = j*inv_m_nRangeXY;
					//~ float a = static_cast<float>(atan2((float)vy, (float)vx)*180.0f/PI);
					//~ float phi = ((a > 0) ? a : a + 180.0f );    
					//~ int angleindex = (int)(phi/m_fWidthOfAPhiInterval); 
					//~ if (angleindex == m_nIntervals) angleindex=m_nIntervals-1;
					//~ rpt.phiindex = angleindex;
					//~ pts2.push_back( rpt );
				//~ }
			//~ }
		//~ }
    //~ //CImg<int> reconstImg(nc,nl,1,1,0);
    //~ //for (vector<Rpoint2>::size_type t = 0; t < pts2.size(); ++t){
    //~ //  reconstImg(pts2[t].x,pts2[t].y) = 255;
    //~ //}
    //~ //cout << pts2.size() << " " << Rtable.size() << endl;
    //~ //reconstImg.display();
    //~ //cout << "accumulate - pts2 loaded." << endl;
		//~ // OpenCv 4-dimensional matrix definition and in general a useful way for defining multidimensional arrays and vectors in c++
		//~ // create accumulator matrix
		//~ int X = static_cast<int>(ceil((float)nc/m_nRangeXY));
		//~ int Y = static_cast<int>(ceil((float)nl/m_nRangeXY));
		//~ int S = static_cast<int>(ceil((float)(m_nTmplWidthMax-m_nTmplWidthMin)/m_nRangeS+1.0f));
    //~ int R = 1;
		//~ if(bMoreAngles) R = 360; // rotate 0->360
    //~ else R = 1; // no rotation
		//~ // icrease accum cells with hits corresponding with slope in Rtable vector rotatated and scaled
		//~ float inv_wtemplate_m_nRangeXY = (float)1/(wtemplate*m_nRangeXY);
		//~ // rotate RTable from minimum to maximum angle
    //~ int sliceNum = 0;
		//~ for (int r=0; r<R; ++r) 
    //~ {  // rotation
      //~ //cout << r << endl;
			//~ std::vector<std::vector<Vec2f>> Rtablerotated(m_nIntervals);
			//~ // cos and sin are computed in the outer loop to reach computational efficiency
      //~ if(isGoodAngle(r))
      //~ {
			  //~ float cs = static_cast<float>(cos(r*PI/180.0));
			  //~ float sn = static_cast<float>(sin(r*PI/180.0));
			  //~ for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii)
        //~ {
				  //~ for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj)
          //~ {
					  //~ int iimod = (int)(ii+r/m_fWidthOfAPhiInterval) % m_nIntervals; // ii is index of phi
					  //~ Rtablerotated[iimod].push_back(Vec2f(cs*Rtable[ii][jj][0] - sn*Rtable[ii][jj][1], sn*Rtable[ii][jj][0] + cs*Rtable[ii][jj][1]));
				  //~ }
			  //~ }
			  //~ // scale the rotated RTable from minimum to maximum scale
			  //~ for (int s=0; s<S; ++s) 
        //~ {  // scale
          //~ //cout << s << endl;
          //~ CImg<int> accumSlice(X,Y,1,1, 0);
				  //~ std::vector<std::vector<Vec2f>> Rtablescaled(m_nIntervals);
				  //~ int w = m_nTmplWidthMin + s*m_nRangeS;
          //~ int h = static_cast<int>(htemplate*(1.0*w/wtemplate)/m_nRangeXY);  
          //~ //cout << h << " " <<  nl << endl;
				  //~ float wratio = (float)w*inv_wtemplate_m_nRangeXY;	
          //~ //cout << "sliceNum " << sliceNum << " wratio " <<  wratio <<" angleInDegrees " <<angleInDegrees<< endl;
				  //~ for (std::vector<std::vector<Vec2f>>::size_type ii = 0; ii < Rtablerotated.size(); ++ii){
					  //~ for (std::vector<Vec2f>::size_type jj= 0; jj < Rtablerotated[ii].size(); ++jj){
						  //~ Rtablescaled[ii].push_back(Vec2f(wratio*Rtablerotated[ii][jj][0], wratio*Rtablerotated[ii][jj][1]));	
					  //~ }
				  //~ }
				  //~ // iterate through each point of edges and hit corresponding cells from rotated and scaled Rtable
				  //~ for (vector<Rpoint2>::size_type t = 0; t < pts2.size(); ++t)
          //~ { // XY plane		
					  //~ int angleindex = pts2[t].phiindex;
					  //~ for (std::vector<Vec2f>::size_type index = 0; index < Rtablescaled[angleindex].size(); ++index){
						  //~ float deltax = Rtablescaled[angleindex][index][0];
						  //~ //float deltax = Rtable[angleindex][index][0];
						  //~ float deltay = Rtablescaled[angleindex][index][1];							
						  //~ //float deltay = Rtable[angleindex][index][1];							
						  //~ int xcell = (int)(pts2[t].x + deltax);
						  //~ int ycell = (int)(pts2[t].y + deltay);
              //~ if(bMoreAngles)
              //~ {
                //~ if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1)&&(xcell<m_fRightAccumBndryXCoor)&&(xcell>m_fLeftAccumBndryXCoor) //)
                  //~ &&(ycell>=h*5/12)&&(ycell<=nl/m_nRangeXY-h*5/12)&&(xcell>=w*5/(2*6*m_nRangeXY))&&(xcell<=nc/m_nRangeXY-w*5/(2*6*m_nRangeXY)) )
                //~ //if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1) )
                //~ {
                  //~ accumSlice(xcell,ycell)++;
						    //~ }
              //~ }
              //~ else
              //~ {
                //~ if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1)&&(xcell<m_fRightAccumBndryXCoor)&&(xcell>m_fLeftAccumBndryXCoor))
                //~ //if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1) )
                //~ {
                  //~ accumSlice(xcell,ycell)++;
						    //~ }
              //~ }
					  //~ }
				  //~ }
          //~ CImg<int> accumSlicePara(2,1,1,1, r,s);
          //~ m_accumPara.append(accumSlicePara,'y');
          //~ m_accum.append(accumSlice,'c');
          //~ ++sliceNum;
          //~ //cout << s << endl;
			  //~ }
        //~ cout << r << endl;
      //~ }
		//~ }
	//~ }
	//~ // load vector pts with all points from the contour
	//~ void readPoints()
  //~ {
		//~ // read original template image and its worked-out contour
    //~ CImg<unsigned char> template_img (m_strTmplContourFileName.c_str());
		//~ // find reference point inside contour image and save it in variable refPoint
		//~ int nl= template_img.height();
		//~ int nc= template_img.width(); 
		//~ for (int j=0; j<nl; ++j) 
    //~ {
			//~ for (int i=0; i<nc; ++i) 
      //~ {  		
				//~ if ( template_img(i,j) == 127 )
        //~ {
					//~ refPointX = i;
          //~ refPointY = j;
				//~ }
			//~ }
		//~ }
    //~ cout << template_img.get_threshold(255).sum() << endl;
		//~ // get Scharr matrices from original template image to obtain contour gradients
    //~ Mat original_img = imread(m_strTmplImgFileName.c_str(), 1);
		//~ Mat input_img_gray;
		//~ input_img_gray.create( Size(original_img.cols, original_img.rows), CV_8UC1);
		//~ cvtColor(original_img, input_img_gray, CV_BGR2GRAY); 
		//~ Mat dx;
		//~ dx.create( Size(original_img.cols, original_img.rows), CV_16SC1);
		//~ Sobel(input_img_gray, dx, CV_16S, 1, 0, CV_SCHARR);
		//~ Mat dy;
		//~ dy.create( Size(original_img.cols, original_img.rows), CV_16SC1);
		//~ Sobel(input_img_gray, dy, CV_16S, 0, 1, CV_SCHARR);
		//~ // load points on vector
		//~ pts.clear();
		//~ int mindx = INT_MAX;
		//~ int maxdx = INT_MIN;
    //~ int mindy = INT_MAX;
    //~ int maxdy = INT_MIN; 
		//~ for (int j=0; j<nl; ++j) 
    //~ {
			//~ for (int i=0; i<nc; ++i)
      //~ {  		
				//~ if ( template_img(i,j) == 255  ) 
				//~ {
					//~ short vx = dx.at<short>(j,i);
					//~ short vy = dy.at<short>(j,i);
					//~ Rpoint rpt;
					//~ rpt.dx = refPointX-i;
					//~ rpt.dy = refPointY-j;
					//~ float a = atan2((float)vy, (float)vx)*180.0f/fPI; // edge direction (with sense) in radians falling in [-180,+180]
          //~ // we need only the direction or the tangent of a given edge point (horizontal,vertival,diagonal,...)
					//~ rpt.phi = ((a > 0) ? a : a + 180.0f ); // 180 is added to all negative directions
					//~ if (rpt.dx < mindx) mindx=rpt.dx;
					//~ if (rpt.dx > maxdx) maxdx=rpt.dx;
					//~ if (rpt.dy < mindy) mindy=rpt.dy;
					//~ if (rpt.dy > maxdy) maxdy=rpt.dy;
					//~ pts.push_back( rpt );
				//~ }
			//~ }
		//~ }
		//~ // maximum width and height of the contour
    //~ htemplate = maxdy-mindy+1;
		//~ wtemplate = maxdx-mindx+1;
	//~ }
//~ 
	//~ // create Rtable from contour points
	//~ void readRtable()
  //~ {
		//~ Rtable.clear();
		//~ Rtable.resize(m_nIntervals);
		//~ // put points in the right interval, according to discretized angle and range size 
		//~ for (vector<Rpoint>::size_type t = 0; t < pts.size(); ++t)
    //~ {
      //~ // an interval is assigned to the phi angles
			//~ int angleindex = (int)(pts[t].phi/m_fWidthOfAPhiInterval); 
      //~ // if angleindex == m_nIntervals, than phi==180 would form a separate interval, can joint the the m_nIntervals-1
			//~ if (angleindex == m_nIntervals) angleindex=m_nIntervals-1;
			//~ Rtable[angleindex].push_back( Vec2i(pts[t].dx, pts[t].dy) );
		//~ }
	//~ }
//~ 
	//~ // create Rtable from contour points
	//~ void visuRtable()
  //~ {
    //~ CImg<int> contour(m_strTmplContourFileName.c_str()),RTabelRes(contour.width(),contour.height(),1,1, 0);
    //~ int nl= contour.height();
	  //~ int nc= contour.width();  
		//~ for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii)
    //~ {
			//~ for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj)
      //~ {
				//~ int dx = Rtable[ii][jj][0];
				//~ int dy = Rtable[ii][jj][1];
				//~ int x = refPointX - dx;
				//~ int y = refPointY - dy;
				//~ //Rtablerotatedscaled[ii].push_back(Vec2i( dx, dy));
				//~ if ( (x<nc)&&(y<nl)&&(x>-1)&&(y>-1) ){
					//~ RTabelRes(x, y) = 255;
				//~ }
			//~ }
		//~ }
    //~ (RTabelRes, contour).display("RTable & contour");
	//~ }
//~ 
	//~ inline int roundToInt(double num) {
		//~ return (num > 0.0) ? (int)(num + 0.5) : (int)(num - 0.5);
	//~ }
//~ 
  //~ bool isGoodAngle(int angle)
  //~ {
    //~ if( (angle >= 355 || angle < 5) && !(angle%2)) return 1;
    //~ //else if(angle >= -104 && angle <= -74 ) return 1;
    //~ //else if(angle > -5 && angle < 5 && !(angle%2)) return 1;
    //~ //else if(angle >=  74 && angle <= 105 ) return 1;
    //~ else if(angle >  175 && angle <= 185 && !(angle%2)) return 1;
    //~ else return 0;
  //~ }
//~ };
