#include "tortoise.h"



struct Rpoint{ 
	int dx;
	int dy;
	float phi;
};

struct Rpoint2{
	int x;
	int y;
	int phiindex;
};

class GenHoughTrnf{

private:
	// accumulator matrix
	CImg<int> accum;
	// accumulator matrix
	Mat showimage;
	// contour points:
	std::vector<Rpoint> pts;
	// reference point (inside contour(
	Vec2i refPoint;
	// R-table of template object:
	std::vector<std::vector<Vec2i>> Rtable;
	// number of intervals for angles of R-table:
	int intervals;
	// threasholds of canny edge detector
	int thr1;
	int thr2;
	// width of template contour
	int wtemplate;
  int htemplate;
	// minimum and maximum width of scaled contour
	int wmin;
	int wmax;
	// minimum and maximum rotation allowed for template
	float phimin;
	float phimax;
	// dimension in pixels of squares in image
	int rangeXY;
	// interval to increase scale
	int rangeS;	

public:

	GenHoughTrnf(){
		// default values

		// canny threasholds
		thr1 = 50;
		thr2 = 150;

		// minimun and maximum width of the searched template
		wmin = 50;
		wmax = 20;
		// increasing step in pixels of the width
		rangeS = 5;
		// side of the squares in which the image is divided
		rangeXY = 6;

		// min value allowed is -pi
		phimin = -PI;
		// max value allowed is +pi
		phimax = +PI;
		// number of slices (angles) in R-table
		intervals = 24;
	}

	void setTresholds(int t1, int t2){
		thr1 = t1;
		thr2 = t2;
	}

	void setLinearPars(int w1, int w2, int rS, int rXY){
		wmin = w1;
		wmax = w2;
		rangeS = rS;
		rangeXY = rXY;
	}

	void setAngularPars(float p1, float p2, int ints){
		if (p1<=p2){
			if (p1>-PI){
				phimin = p1;	
			}
			if (p2<+PI){
				phimax = p2;
			}
		}
		intervals = ints;
	}

	// save file with canny edge of the original image
	void createTemplate(){
		Mat input_img = imread("files\\template_original.jpg", 1);
		Mat src_gray;
		Mat detected_edges;
		src_gray.create( Size(input_img.cols, input_img.rows), CV_8UC1);
		cvtColor(input_img, src_gray, CV_BGR2GRAY); 
		blur( src_gray, detected_edges, Size(3,3) );
		Canny( detected_edges, detected_edges, 1, 100, 3 );
		imwrite("files\\contour_rough.bmp", detected_edges);
	}

	void createRtable(){
		// code can be improved reading a pre-saved Rtable
		readPoints();
		readRtable();
	}

	// fill accumulator matrix
	void accumulate(cv::Mat& input_img,cv::Mat& detected_edges, int leftX, int rightX){
    wmin = wtemplate;
		showimage = detected_edges;
		// transform image to grayscale:
		Mat src_gray;
		src_gray.create( Size(input_img.cols, input_img.rows), CV_8UC1);
		cvtColor(input_img, src_gray, CV_BGR2GRAY); 
		// reduce noise with a kernel 3x3 and get cannyedge image:
		//Mat detected_edges;
		//blur( src_gray, detected_edges, Size(3,3) );
		//Canny( detected_edges, detected_edges, thr1, thr2, 3 );
		//imshow("detected_edges", detected_edges);
    //waitKey(0);
		// get Scharr matrices from image to obtain contour gradients
		Mat dx;
		dx.create( Size(input_img.cols, input_img.rows), CV_16SC1);
		Sobel(src_gray, dx, CV_16S, 1, 0, CV_SCHARR);
		Mat dy;
		dy.create( Size(input_img.cols, input_img.rows), CV_16SC1);
		Sobel(src_gray, dy, CV_16S, 0, 1, CV_SCHARR);
		// load all points from image all image contours on vector pts2
    //imshow("detected_edges",detected_edges);
		int nl= detected_edges.rows;
		int nc= detected_edges.cols ;
		float deltaphi = PI/intervals;
		float inv_deltaphi = (float)intervals/PI;
		float inv_rangeXY = (float)1/rangeXY;
		float pi_half = PI*0.5f;
		std::vector<Rpoint2> pts2;
		for (int j=0; j<nl; ++j) {
			uchar* data= (uchar*)(detected_edges.data + detected_edges.step.p[0]*j);
			for (int i=0; i<nc; ++i) {  		
				if ( data[i]==255  ) // consider only white points (contour)
				{
					short vx = dx.at<short>(j,i);
					short vy = dy.at<short>(j,i);
					Rpoint2 rpt;
					rpt.x = i*inv_rangeXY;
					rpt.y = j*inv_rangeXY;
					float a = atan2((float)vy, (float)vx);              //	gradient angle in radians
					float phi = ((a > 0) ? a-pi_half : a+pi_half);      // contour angle with respect to x axis
					int angleindex = (int)((phi+PI*0.5f)*inv_deltaphi); // index associated with angle (0 index = -90 degrees)
					if (angleindex == intervals) angleindex=0;// -90°angle and +90° has same effect
					rpt.phiindex = angleindex;
					pts2.push_back( rpt );
				}
			}
		}
    cout << "accumulate - pts2 loaded." << endl;
		// OpenCv 4-dimensional matrix definition and in general a useful way for defining multidimensional arrays and vectors in c++
		// create accumulator matrix
		int X = ceil((float)nc/rangeXY);
		int Y = ceil((float)nl/rangeXY);
		//int S = ceil((float)(wmax-wmin)/rangeS+1.0f);
		int S = 1;
		//int R = ceil(phimax/deltaphi)-floor(phimin/deltaphi);
		int R = 1;
    cout << R << '\t' << deltaphi << '\t' << phimax << '\t' << phimin << '\t' << ceil(phimax/deltaphi) << '\t' << floor(phimin/deltaphi) << endl;
		if (phimax==PI && phimin==-PI) R--;	
		//int r0 = -floor(phimin/deltaphi);
		int r0 = 0;
		int matSizep_S[] = {X, Y, S, R};
    accum.assign(X,Y,1,S*R, 0);
		//accum.create(4, matSizep_S, CV_16S);
		//accum = Scalar::all(0);
		// icrease accum cells with hits corresponding with slope in Rtable vector rotatated and scaled
		float inv_wtemplate_rangeXY = (float)1/(wtemplate*rangeXY);
		// rotate RTable from minimum to maximum angle
		for (int r=0; r<R; ++r) {  // rotation
			int reff = r-r0;
			std::vector<std::vector<Vec2f>> Rtablerotated(intervals);
			// cos and sin are computed in the outer loop to reach computational efficiency
      double angleInDegrees = 1.0*reff*deltaphi*180.0/PI;
      cout << angleInDegrees << endl;
      cout << r << " " << r0 << " " << reff << " " << deltaphi << " " << 1.0*reff*deltaphi*180.0/PI << " "<< angleInDegrees << endl;
      if(isGoodAngle(angleInDegrees))
      {
			  float cs = cos(reff*deltaphi);
			  float sn = sin(reff*deltaphi);
			  for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii){
				  for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj){
					  int iimod = (ii+reff) % intervals; // ii is index of phi
					  Rtablerotated[iimod].push_back(Vec2f(cs*Rtable[ii][jj][0] - sn*Rtable[ii][jj][1], sn*Rtable[ii][jj][0] + cs*Rtable[ii][jj][1]));
				  }
			  }

			  // scale the rotated RTable from minimum to maximum scale
			  for (int s=0; s<S; ++s) {  // scale
				  std::vector<std::vector<Vec2f>> Rtablescaled(intervals);
				  //int w = wmin + s*rangeS;
				  int w = wtemplate;
          int h = htemplate*(1.0*w/wtemplate)/rangeXY;  
          cout << h << " " <<  nl << endl;
				  float wratio = (float)w*inv_wtemplate_rangeXY;	
				  for (std::vector<std::vector<Vec2f>>::size_type ii = 0; ii < Rtablerotated.size(); ++ii){
					  for (std::vector<Vec2f>::size_type jj= 0; jj < Rtablerotated[ii].size(); ++jj){
						  Rtablescaled[ii].push_back(Vec2f(wratio*Rtablerotated[ii][jj][0], wratio*Rtablerotated[ii][jj][1]));	
					  }
				  }
          cout << Rtable[0][0] << endl;
          int cnt = 0;
				  // iterate through each point of edges and hit corresponding cells from rotated and scaled Rtable
				  for (vector<Rpoint2>::size_type t = 0; t < pts2.size(); ++t){ // XY plane				
					  int angleindex = pts2[t].phiindex;
					  for (std::vector<Vec2f>::size_type index = 0; index < Rtable[angleindex].size(); ++index){
						  //float deltax = Rtablescaled[angleindex][index][0];
						  float deltax = Rtable[angleindex][index][0];
						  //float deltay = Rtablescaled[angleindex][index][1];							
						  float deltay = Rtable[angleindex][index][1];							
						  int xcell = (int)(pts2[t].x + deltax);
						  int ycell = (int)(pts2[t].y + deltay);
              //if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1)&&(xcell<rightX)&&(xcell>leftX)&&(ycell>=h/2)&&(ycell<=nl/rangeXY-h/2)&&(xcell>=w/(2*rangeXY))&&(xcell<=nc/rangeXY-w/(2*rangeXY)) ){
              if ( (xcell<X)&&(ycell<Y)&&(xcell>-1)&&(ycell>-1)){
							  //(*( (short*)(accum.data + xcell*accum.step.p[0] + ycell*accum.step.p[1] + s*accum.step.p[2]+ r*accum.step.p[3])))++;
                if(xcell>150&&xcell<160&&ycell>212&&ycell<223) cnt++;
                  //cout << xcell << " " << ycell << endl;
                accum(xcell,ycell,0,r*S+s)++;
							  //(*ptrat4D(accum, xcell, ycell, s, r))++;
						  }
              //imshow("accum",accum);
              //waitKey(0);
					  }
                   
				  }
           cout << cnt << endl;
          cout << s << endl;
			  }

      }
      cout << r << " " << R << endl;


		}
    accum.display();
    //imshow("accum",accum);
    //waitKey(0);
	}


	// show the best candidate detected on image
	Mat bestCandidate(){
		double minval;
		double maxval;
		int id_min[4] = { 0, 0, 0, 0};
    int S = ceil((float)(wmax-wmin)/rangeS+1.0f);
    int r =  ((int)accum.get_stats()(11))/S;
    int s = ((int)accum.get_stats()(11)) % S;
    int id_max[4] = { accum.get_stats()(8), accum.get_stats()(9), s, r};
    maxval = accum.max();
    cout << maxval << " " << id_max[0] << " " << id_max[1] << " " << id_max[2] << " " << id_max[3] << endl;
		//minMaxIdx(accum, &minval, &maxval, id_min, id_max);
     accum.get_channel(accum.get_stats()(11)).display();
    for(int i = 0; i < accum.spectrum(); ++i)
    {
      if(accum.get_channel(i)!=0)
      {
        //accum.get_channel(i).display();
      }

    }
    accum.get_channel(accum.get_stats()(11)).display();
		int nl= showimage.rows;
		int nc= showimage.cols; 
		Mat	input_img2 = showimage.clone();

		Vec2i referenceP = Vec2i(id_max[0]*rangeXY+(rangeXY+1)/2, id_max[1]*rangeXY+(rangeXY+1)/2);
		
		// rotate and scale points all at once. Then impress them on image
		std::vector<std::vector<Vec2i>> Rtablerotatedscaled(intervals);
		float deltaphi = PI/intervals;
		//int r0 = -floor(phimin/deltaphi);
		int r0 = 0;
		int reff = id_max[3]-r0;
		float cs = cos(reff*deltaphi);
		float sn = sin(reff*deltaphi);
		int w = wmin + id_max[2]*rangeS;
		float wratio = (float)w/(wtemplate);	
		for (std::vector<std::vector<Vec2i>>::size_type ii = 0; ii < Rtable.size(); ++ii){
			for (std::vector<Vec2i>::size_type jj= 0; jj < Rtable[ii].size(); ++jj){
				int iimod = (ii+reff) % intervals;
				int dx = roundToInt(wratio*(cs*Rtable[ii][jj][0] - sn*Rtable[ii][jj][1]));
				int dy = roundToInt(wratio*(sn*Rtable[ii][jj][0] + cs*Rtable[ii][jj][1]));
				int x = referenceP[0] - dx;
				int y = referenceP[1] - dy;
				//Rtablerotatedscaled[ii].push_back(Vec2i( dx, dy));
				if ( (x<nc)&&(y<nl)&&(x>-1)&&(y>-1) ){
					input_img2.at<Vec3b>(y, x) = Vec3b(0, 255, 255);
				}
			}
		}
    cout << id_min << endl;
    imshow("input_img", input_img2);

    return input_img2;
    //waitKey(0);
    //imshow("input_img", showimage);
    //waitKey(0);
		// show result
		//bool alt = false; 
		//for(;;)
		//{
		//	char c = (char)waitKey(0);
		//	if( c == 27 )
		//		break;
		//	if (alt){
		//		imshow("input_img", input_img2);
		//	}
		//	else{
		//		imshow("input_img", showimage);
		//	}
		//	alt = !alt;
		//}
	}


	// finds local minima above a certain threashold
	void localMaxima(){
		// to bve implemented ...
	}


private:

	// load vector pts with all points from the contour
	void readPoints(){
		// read original template image and its worked-out contour
		Mat original_img = imread("ght\\tortoiseTemplateImg.bmp", 1);
		Mat input_img_gray;
		input_img_gray.create( Size(original_img.cols, original_img.rows), CV_8UC1);
		cvtColor(original_img, input_img_gray, CV_BGR2GRAY); 
	    //Mat template_img = imread("files\\contour_def.bmp", 1);
	    Mat template_img = imread("ght\\tortoiseTemplateContour.bmp", 1);
		// find reference point inside contour image and save it in variable refPoint
		int nl= template_img.rows;
		int nc= template_img.cols; 
		for (int j=0; j<nl; ++j) {
			Vec3b* data= (Vec3b*)(template_img.data + template_img.step.p[0]*j);
			for (int i=0; i<nc; ++i) {  		
				if ( data[i]==Vec3b(127,127,127)  ){
					refPoint = Vec2i(i,j);
				}
			}
		}
		// get Scharr matrices from original template image to obtain contour gradients
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
		for (int j=0; j<nl; ++j) {
			Vec3b* data= (Vec3b*)(template_img.data + template_img.step.p[0]*j);
			for (int i=0; i<nc; ++i) {  		
				if ( data[i]==Vec3b(255,255,255)  ) 
				{

					//short vx = dx.at<short>(j,i);
					//short vy = dy.at<short>(j,i);
					//Rpoint2 rpt;
					//rpt.x = i*inv_rangeXY;
					//rpt.y = j*inv_rangeXY;
					//float a = atan2((float)vy, (float)vx);              //	gradient angle in radians
					//float phi = ((a > 0) ? a-pi_half : a+pi_half);      // contour angle with respect to x axis
					//int angleindex = (int)((phi+PI*0.5f)*inv_deltaphi); // index associated with angle (0 index = -90 degrees)
					//if (angleindex == intervals) angleindex=intervals-1;// -90°angle and +90° has same effect

					short vx = dx.at<short>(j,i);
					short vy = dy.at<short>(j,i);
					Rpoint rpt;
					//float mag = std::sqrt( float(vx*vx+vy*vy) );
					rpt.dx = refPoint(0)-i;
					rpt.dy = refPoint(1)-j;
					float a = atan2((float)vy, (float)vx); //radians
					rpt.phi = ((a > 0) ? a-PI/2.0 : a+PI/2.0);
					//float a = atan2((float)vy, (float)vx) * 180/3.14159265358979f; //degrees
					//rpt.phi = ((a > 0) ? a-90 : a+90);
					// update further right and left dx
					if (rpt.dx < mindx) mindx=rpt.dx;
					if (rpt.dx > maxdx) maxdx=rpt.dx;
					if (rpt.dy < mindy) mindy=rpt.dy;
					if (rpt.dy > maxdy) maxdy=rpt.dy;
					pts.push_back( rpt );
				}
			}
		}
		// maximum width of the contour
    htemplate = maxdy-mindy+1;
		wtemplate = maxdx-mindx+1;
	}

	// create Rtable from contour points
	void readRtable(){
		Rtable.clear();
		Rtable.resize(intervals);
		// put points in the right interval, according to discretized angle and range size 
		float range = PI/intervals;
    cout <<  pts.size() << endl;
		for (vector<Rpoint>::size_type t = 0; t < pts.size(); ++t){
			int angleindex = (int)((pts[t].phi+PI/2.0)/range);
			if (angleindex == intervals) angleindex=0;
			Rtable[angleindex].push_back( Vec2i(pts[t].dx, pts[t].dy) );
		}
	}

	inline int roundToInt(float num) {
		return (num > 0.0) ? (int)(num + 0.5f) : (int)(num - 0.5f);
	}

	inline short at4D(Mat &mt, int i0, int i1, int i2, int i3){
		//(short)(mt.data + i0*mt.step.p[0] + i1*mt.step.p[1] + i2*mt.step.p[2] + i3*mt.step.p[3]);	
		return *( (short*)(mt.data + i0*mt.step.p[0] + i1*mt.step.p[1] + i2*mt.step.p[2] + i3*mt.step.p[3]));
	}

	inline short* ptrat4D(Mat &mt, int i0, int i1, int i2, int i3){
		return (short*)(mt.data + i0*mt.step.p[0] + i1*mt.step.p[1] + i2*mt.step.p[2] + i3*mt.step.p[3]);
	}

  bool isGoodAngle(double angle)
  {
    if(angle >= -180 && angle < -175 ) return 1;
    //else if(angle >= -104 && angle <= -74 ) return 1;
    else if(angle > -5 && angle < 5 ) return 1;
    //else if(angle >=  74 && angle <= 105 ) return 1;
    else if(angle >  175 && angle <= 180 ) return 1;
    else return 0;
  }
};



void runGHT(char c, string ImgFileName, string edgesFileName, int leftX, int rightX){
	//c = 't'; // create template
	//c = 'r'; // run algorithm
	if (c == 't'){   
		GenHoughTrnf ght;
		ght.createTemplate();
	}
	else if (c == 'r'){
		GenHoughTrnf ght;

		int thr1 = 50; // canny thresholds
		int thr2 = 150;
		ght.setTresholds(thr1, thr2); 
		int wmin = 150; // minimun and maximum width of the searched template
		int wmax = 300;
		int rangeS = 200; // increasing step in pixels of the width
		int rangeXY = 1; // side of the squares in which the image is divided 
    ght.setLinearPars(wmin, wmax, rangeS, rangeXY);
		float phimin = -PI; // min value allowed is -pi
		float phimax = +PI; // max value allowed is +pi
		int intervals = 180; // number of slices (angles) in R-table
    ght.setAngularPars(phimin, phimax, intervals);

		ght.createRtable();
    cout << "R-table created." << endl;
		//Mat detect_img = imread("files\\Img_01.png", 1);
		Mat detect_img = imread(ImgFileName.c_str(), 1);
    Mat detected_edges = imread(edgesFileName.c_str(),1);

		//Mat detect_img = imread("files\\Img_03.png", 1);
		ght.accumulate(detect_img, detected_edges, (leftX/(rangeXY))-1, (rightX/(rangeXY))+2 );
    cout << "Accumulation done." << endl;
    Mat res;
		res = ght.bestCandidate();
    imshow("detected_edges",res);
    waitKey(0);
    imwrite(ImgFileName.substr(0,ImgFileName.size()-7) + "Res.bmp", res );

	}
}