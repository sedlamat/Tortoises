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



std::istream &operator>>(std::istream &stream, cv::Point &pt);

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
    void locate_Abd_juncs();
};

/************* member functions definitions ************************/

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
	    _junctions_found(0),
	    _color()
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
	//this->write_info();
    }
}

std::istream &operator>>(std::istream &stream, cv::Point &pt)
{
    // cv::Point is printed in format [pt.x, pt.y]
    int a = 0;
    // moves through the stream until an integer is ahead
    while (!std::isdigit(a = stream.peek()) && a != '-') stream.get();
    // loads pt.x
    stream >> pt.x;
    // again moves through the stream until an integer is ahead
    while (!std::isdigit(a = stream.peek()) && a != '-') stream.get();
    // loads pt.y
    stream >> pt.y;
    // gets beyond the cv::Point expression
    while (stream.get() != ' ') ;
    return stream;
}

void Tortoise::read_info()
{
    std::ifstream info(_info_file_name);
    if (info.good()) {
	info >> _plastron_found;
	info >> _junctions_found;
	info >> _angle;
	info >> _center;
	info >> _left_side;
	info >> _right_side;
	for (auto &junc : _l_juncs) {
	    info >> junc;
	}
	for (auto &junc : _r_juncs) {
	    info >> junc;
	}
	std::cout << _plastron_found << '\n';
	std::cout << _junctions_found << '\n';
	std::cout << _angle << '\n';
	std::cout << _center << '\n';
	std::cout << _left_side << '\n';
	std::cout << _right_side << '\n';

	this->print_juncs();

	if (!info.good()) { // error when reading -> recompute data
	    _plastron_found = 0;
	    _junctions_found = 0;
	}
    } else {
	std::cout << "File Tg#####_info.txt does not exists. It will";
	std::cout << " be created." << std::endl;
    }

    info.close();
}

void Tortoise::write_info() const
{
    std::ofstream info(_info_file_name, std::fstream::out);
    if (info.good()) {
	info << _plastron_found << " ";
	info << _junctions_found << " ";
	info << _angle << " ";
	info << _center << " ";
	info << _left_side << " ";
	info << _right_side << " ";
	for (const auto &junc : _l_juncs) {
	    info << junc  << " ";
	}
	for (const auto &junc : _r_juncs) {
	    info << junc << " ";
	}

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

void Tortoise::put_point_on_img(cv::Mat &img,
				const cv::Rect &img_rect,
				const cv::Point &pt,
				const cv::Vec3b &color,
				const cv::Point &shift)
{
    if (img_rect.contains(pt)) {
	cv::rectangle(img,
		      pt-shift,
		      pt+shift,
		      cv::Scalar(color), CV_FILLED);
    }
}

void Tortoise::display_points_on_plastron()
{
    cv::Mat img;
    _plastron_img.copyTo(img);
    cv::Rect img_rect(cv::Point(0,0), img.size());

    const double inv_resize_coeff = std::max(img.cols, img.rows) / 200;
    cv::Point shift(inv_resize_coeff,inv_resize_coeff);

    for (const auto &junc : _l_juncs) {
	this->put_point_on_img(img, img_rect, junc, _color.red, shift);
    }
    for (const auto &junc : _r_juncs) {
	this->put_point_on_img(img, img_rect, junc, _color.green,
								shift);
    }

    this->put_point_on_img(img, img_rect, _center, _color.blue,
								shift);
    this->put_point_on_img(img, img_rect, _left_side, _color.purple,
								shift);
    this->put_point_on_img(img, img_rect, _right_side, _color.cyan,
								shift);
    this->put_point_on_img(img, img_rect, _r_juncs[0], _color.red,
								shift);

    sedlamat::display(img);
}

void Tortoise::rotate_pt(cv::Point &pt, const cv::Point &center,
			 const int angle)
{
    const double cs = std::cos(angle * M_PI / 180.0);
    const double sn = std::sin(angle * M_PI / 180.0);

    pt -= center;

    const int x = pt.x;
    const int y = pt.y;
    pt.x = cs*x - sn*y;
    pt.y = sn*x + cs*y;

    pt += center;
}

void Tortoise::rotate_img_and_pts(const int angle)
{

    // makes _plastron_img square matrix with wide borders
    const int w = _plastron_img.size().width;
    const int h = _plastron_img.size().height;
    const int max_dist = std::max(std::max(_center.x, _center.y),
			          std::max(w-_center.x, h-_center.y));
    const int top = max_dist - _center.y;
    const int bottom = max_dist - (h-_center.y);
    const int left = max_dist - _center.x;
    const int right = max_dist - (w-_center.x);
    cv::copyMakeBorder(_plastron_img, _plastron_img, top, bottom, left,
		       right, cv::BORDER_CONSTANT);

    // rotates points and _plastron_img around _center
    const cv::Point shift = cv::Point(max_dist, max_dist) - _center;

    for (auto &junc : _l_juncs) {
	this->rotate_pt(junc, _center, -angle);
	junc += shift;
    }
    for (auto &junc : _r_juncs) {
	this->rotate_pt(junc, _center, -angle);
	junc += shift;
    }
    this->rotate_pt(_left_side, _center, -angle);
    _left_side += shift;
    this->rotate_pt(_right_side, _center, -angle);
    _right_side += shift;

    _center += shift;

    cv::Mat M = cv::getRotationMatrix2D(_center, angle, 1);
    cv::warpAffine(_plastron_img, _plastron_img, M,
						_plastron_img.size());
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
    // assigns to interest points
    for (int ii = 0; ii < 7; ++ii) {
	_l_juncs[ii] = _r_juncs[ii] = interest_pts[ii];
    }
    _left_side = interest_pts[7];
    _right_side = interest_pts[8];
    _center = interest_pts[9];

    _angle = general_hough.get_best_angle();
}


void Tortoise::locate_junctions()
{
    std::cout << "locate_junctions entered" << std::endl;
    this->rotate_img_and_pts(_angle);
    this->display_points_on_plastron();
    this->locate_Abd_juncs();
}


void Tortoise::locate_Abd_juncs()
{
//~ //using edges on a wider central-seam stripe
//~ const int nAbdSeamLength = sJunctions.s5AbdToFem.yCoor - sJunctions.s4PecToAbd.yCoor;
//~ CImg<int> imgPlastron(imgRotatedOriginal.get_crop(
  //~ min(m_imgRotatedOriginal.width()-1,max(sPlastronCentreOnRotatedImg.xCoor-nAbdSeamLength,0)),sJunctions.s1Head.yCoor,
  //~ min(m_imgRotatedOriginal.width()-1,max(sPlastronCentreOnRotatedImg.xCoor+nAbdSeamLength,0)),sJunctions.s7Tail.yCoor));
//~ imgPlastron = GetUniformlyResizedImg(imgPlastron,500);
//~ int nCentreYCoor = (sPlastronCentreOnRotatedImg.yCoor - sJunctions.s1Head.yCoor)*imgPlastron.height()/(sJunctions.s7Tail.yCoor-sJunctions.s1Head.yCoor);
//~ CImg<int> imgEdges(GetEdgeImg(imgPlastron,0,0,0));
//~ //imgEdges.display();
      //~ const int nImgEdgesWidth = imgEdges.width();
      //~ const int nImgEdgesHeight = imgEdges.height();

//~ //every edge pix is part of a some path - all edge pixs are assigned an index of the path they belong to
      //~ CImg<int> sidePathsIndexes(imgEdges.get_fill(0));
//~ //pixs on the cetral path have 2-column index array, beacause there meet paths from left and right side  - they meet in the middle: x = nImgEdgesWidth/2
      //~ CImg<int> indexesOnTheCentralPath(2,nImgEdgesHeight,1,1, 0);
//~ //prefilling - side paths to the value 10000, for future easier processing
//~ const int nMaxPathValue = 10000;
      //~ CImg<int> sidePathsFromLeft(imgEdges.get_fill(nMaxPathValue)), sidePathsFromRight(sidePathsFromLeft);

//~ for(int y = 0; y < nImgEdgesHeight; ++y)
      //~ {
  //~ if(imgEdges(nImgEdgesWidth/2,y)==0) imgEdges(nImgEdgesWidth/2,y) = 1;
  //~ if(imgEdges(4,y)==0) imgEdges(4,y) = 1;
  //~ if(imgEdges(nImgEdgesWidth-5,y)==0) imgEdges(nImgEdgesWidth-5,y) = 1;
      //~ }
//~ //prefilling - the sides (5 the most left columns) of sidePathsFromLeft and sidePathsIndexes
      //~ int nPathIndex = 2;
//~ for(int y = 0; y < nImgEdgesHeight; ++y)
      //~ {
	      //~ int x = 4;
	      //~ {
//~ //prefilling - for all edge pixs in the area
		      //~ if(imgEdges(x,y))
		      //~ {
//~ //prefilling - initial value (distance) of sidePathsFromLeft is set to 2
			      //~ sidePathsFromLeft(x,y) = 2;
//~ //prefilling - edge pixs in the area are all assigned an index in sidePathsIndexes, which they will pass on to pixs that will be in the same path
			      //~ sidePathsIndexes(x,y) = nPathIndex;
		      //~ }
	      //~ }
	      //~ nPathIndex++;
      //~ }
//~ //prefilling - the sides (5 the most right columns) of sidePathsFromRight and sidePathsIndexes
      //~ nPathIndex += 1000;
      //~ for(int y = 0; y < nImgEdgesHeight; ++y)
      //~ {
  //~ int x = nImgEdgesWidth-5;
	      //~ {
		      //~ if(imgEdges(x,y))
		      //~ {
//~ //prefilling - initial value (distance) of sidePathsFromRight is set to 2
			      //~ sidePathsFromRight(x,y) = 2;
//~ //prefilling - edge pixs in the area are all assigned an index in sidePathsIndexes, which they will pass on to pixs that will be in the same path
			      //~ sidePathsIndexes(x,y) = nPathIndex;
		      //~ }
	      //~ }
  //~ nPathIndex++;
      //~ }
//~ //(sidePathsIndexes,sidePathsFromLeft,sidePathsFromRight).display();
//~ //prefilling - distanceMask is used to determined the distance to a pix from pixs near it (3 pix up and down and 5 pixs backwards, i.e. 7x5)
      //~ CImg<int> distanceMask(7,5,1,1, 0);
//~ //prefilling - the rows of distanceMask are as follows [4,1,1,1,1,1,4 ; 7,4,4,4,4,4,7 ; 12,9,9,9,9,9,12 ; 19,16,16,16,16,16,19 ; 28,25,25,25,25,25,28]
      //~ for(int x = 0; x < 7; x++) for(int y = 0; y < 5; y++) distanceMask(x,y) = static_cast<int>(pow((y+1),2.0));
      //~ for(int y = 0; y < 5; y++)
//~ {
  //~ distanceMask(0,y) += 3;
  //~ distanceMask(6,y) += 3;
  //~ distanceMask(1,y) += 2;
  //~ distanceMask(5,y) += 2;
  //~ distanceMask(2,y) += 1;
  //~ distanceMask(4,y) += 1;
//~ }
      //~ //distanceMask.display();
//~ //filling left side paths - rotating distanceMask
      //~ distanceMask.rotate(90);
//~ //distanceMask.display();
//~ //filling left side paths - going though all img pixs (rows only those where the central path is)
//~ for(int x = 5; x < nImgEdgesWidth; x++)
      //~ {
  //~ for(int y = 5; y < nImgEdgesHeight-5; ++y)
	      //~ {
//~ //filling left side paths - evaluating only edge pixs on the left of the central path
    //~ if(imgEdges(x,y) && x<=nImgEdgesWidth/2)
		      //~ {
//~ //filling left side paths - for an edge pixs an area next to it is cropped and distance is added so cropAreaDistances contains total distances from the side to the current pix though nearby pixs
			      //~ CImg<int> cropAreaDistances(sidePathsFromLeft.get_crop(x-1,y-3,x-5,y+3) + distanceMask);
//~ //filling left side paths - minimum distance is saved
			      //~ int min = cropAreaDistances.min();
//~ //filling left side paths - if in the area there is any edge pixs then min < 10000, and we will change the distance to the current pix from 10000 to min
      //~ if(min<nMaxPathValue)
			      //~ {
//~ //filling left side paths - minimum distance is saved in the current pix
				      //~ sidePathsFromLeft(x,y) = min;
//~ //filling left side paths - index of the pix with minimum distace is passed on to the current pix
				//~ sidePathsIndexes(x,y) = sidePathsIndexes(x-5+static_cast<int>(cropAreaDistances.get_stats()(4)),y-3+static_cast<int>(cropAreaDistances.get_stats()(5)));
//~ //filling left side paths - if the pix is on the central path, the index is saved in indexesOnTheBestCentralPath
	//~ if(x==nImgEdgesWidth/2) indexesOnTheCentralPath(0,y) = sidePathsIndexes(x,y);
			      //~ }
		      //~ }
	      //~ }
      //~ }
//~ //  (sidePathsFromLeft,sidePathsIndexes).display();
//~ //filling right side paths - rotating distanceMask
     //~ distanceMask.rotate(180);
//~ //filling right side paths - the same as is in the left side
//~ for(int x = nImgEdgesWidth-6; x >= 0; --x)
      //~ {
  //~ for(int y = 5; y < nImgEdgesHeight-5; ++y)
	      //~ {
    //~ if(imgEdges(x,y) && x >= nImgEdgesWidth/2)
		      //~ {
			      //~ CImg<int> cropAreaDistances(sidePathsFromRight.get_crop(x+1,y-3,x+5,y+3) + distanceMask );
			      //~ int min = cropAreaDistances.min();
      //~ if(min<nMaxPathValue)
			      //~ {
				      //~ sidePathsFromRight(x,y) = min;
				      //~ sidePathsIndexes(x,y) = sidePathsIndexes(x+1+static_cast<int>(cropAreaDistances.get_stats()(4)),y-3+static_cast<int>(cropAreaDistances.get_stats()(5)));
	//~ if(x==nImgEdgesWidth/2) indexesOnTheCentralPath(1,y) = sidePathsIndexes(x,y);
			      //~ }
		      //~ }
	      //~ }
      //~ }
      //~ //(sidePathsFromLeft,sidePathsFromRight).display();
//~ //postprocessing - set values 10000 to 0 in sidePathsFromLeft and sidePathsFromRight
//~ for(int y = 0; y < nImgEdgesHeight; y++)
      //~ {
  //~ for(int x = 0; x < nImgEdgesWidth; x++)
	      //~ {
    //~ if(sidePathsFromLeft(x,y)  == nMaxPathValue) sidePathsFromLeft(x,y) = 0;
    //~ if(sidePathsFromRight(x,y) == nMaxPathValue) sidePathsFromRight(x,y) = 0;
	      //~ }
      //~ }
     //~ //(sidePathsFromLeft+sidePathsFromRight).display("sidePaths");
//~ //postprocessing - CentralSideCrosses will contain sum of 2x values of distances of pixs on the central seam of the left and right paths
//~ CImg<int> CentralSideCrosses(1,nImgEdgesHeight,1,1, nMaxPathValue);
//~ //postprocessing - for all central path pixs near the tortoise center
//~ for(int y = 6; y < nImgEdgesHeight-6; ++y)
      //~ {
  //~ int x = nImgEdgesWidth/2;
//~ //postprocessing - if there are paths from left and right that end in the central-path pix
	      //~ if(sidePathsFromLeft(x,y)>0 && sidePathsFromRight(x,y)>0)
	      //~ {
//~ //postprocessing - the values are added to CentralSideCrosses
		      //~ CentralSideCrosses(y) = sidePathsFromLeft(x,y) + sidePathsFromRight(x,y);
	      //~ }
//~ //postprocessing - else, left and right indexes of the central-path pix are set to 0
	      //~ else
	      //~ {
		      //~ indexesOnTheCentralPath(0,y) = 0;
		      //~ indexesOnTheCentralPath(1,y) = 0;
	      //~ }
      //~ }
//~ // (CentralSideCrosses,indexesOnTheCentralPath,(sidePathsFromLeft+sidePathsFromRight)).display();
//~ //finding best side paths - setting number of pixs that are to be discarded, being next (based on index) to a pix where we will find the local minimum (of distances in CentralSideCrosses)
//~ //                        - with 3000pixs height img time constant 0.025 (based on real images) so that we get ~ 70pixs
//~ //                        - (it has to be so much that we eliminate all neighbouring paths with similar distances, but so that we do not eliminate other correct side paths)
//~ int halfIntervalOfIndexes = static_cast<int>(nImgEdgesHeight*0.1);
//~ //finding best side paths - countOfCrossing counts side paths and is used to index the paths
      //~ int countOfCrossing = 0;
//~ //finding best side paths - we want no more than 9 paths (5-side paths, 2-end and beggining of the plastron, 2-ends of tail and head/ends of the shell)
      //~ const int permittedNumberOfPaths = 9;
      //~ CImg<int> yCoorOfTheJunctions(1,permittedNumberOfPaths,1,1, 0);
//~ bool bEnd = 0;
//~ int nNumOfUpJunctions = 0;
//~ int nNumOfDownJunctions = 0;
//~ //finding best side paths - going though CentralSideCrosses until there is no pix with smaller distance than upperSupremum or the is permittedNumberOfPaths reached
//~ // (CentralSideCrosses,indexesOnTheCentralPath,yCoorOfTheJunctions,valuesOfTheJunctions).display();
      //~ while(CentralSideCrosses.min()<nMaxPathValue && countOfCrossing < permittedNumberOfPaths && !bEnd)
      //~ {
//~ //finding best side paths - finding the current minimum in CentralSideCrosses (distances of left and right paths)
	      //~ int yMin = static_cast<int>(CentralSideCrosses.get_stats()(5));
  //~ int currentMin =  CentralSideCrosses.min();
  //~ if(yMin<nCentreYCoor) nNumOfUpJunctions++;
  //~ else nNumOfDownJunctions++;
  //~ if(nNumOfUpJunctions==2 && nNumOfDownJunctions==2) bEnd=1;
  //~ if((nNumOfUpJunctions<=2&&yMin<nCentreYCoor) || (nNumOfDownJunctions<=2&&yMin>nCentreYCoor)) // spatne neco jinak || &&
  //~ {
		//~ yCoorOfTheJunctions(countOfCrossing) = yMin;
  //~ }
//~ //finding best side paths - getting left index of the pix with the minimum and its y coordinate
		//~ int leftIndexOfMinPix = indexesOnTheCentralPath(0,yMin);
//~ //finding best side paths - getting right index of the pix with the minimum and its y coordinate
		//~ int rightIndexOfMinPix = indexesOnTheCentralPath(1,yMin);
//~ //finding best side paths - deleting indexes close to index of the pix with the current min value
    //~ for(int y = 0; y < nImgEdgesHeight; y++)
		//~ {
      //~ if( (indexesOnTheCentralPath(0,y) <= leftIndexOfMinPix + halfIntervalOfIndexes) &&
				//~ (indexesOnTheCentralPath(0,y) >= max(1,leftIndexOfMinPix - halfIntervalOfIndexes)))
			//~ {
//~ //finding best side paths - index to be deleted from indexesOnTheBestCentralPath
				//~ indexesOnTheCentralPath(0,y) = 0;
	//~ indexesOnTheCentralPath(1,y) = 0;
//~ //finding best side paths - coresponding distances in CentralSideCrosses set to upperSupremum (cannot become minimum in the future)
	//~ CentralSideCrosses(y) = nMaxPathValue;
			//~ }
//~ //finding best side paths - the same for the right side indexes
      //~ else if( indexesOnTheCentralPath(1,y) <= rightIndexOfMinPix + halfIntervalOfIndexes &&
				     //~ indexesOnTheCentralPath(1,y) >= max(1,rightIndexOfMinPix - halfIntervalOfIndexes) )
			//~ {
	//~ indexesOnTheCentralPath(0,y) = 0;
				//~ indexesOnTheCentralPath(1,y) = 0;
				//~ CentralSideCrosses(y) = nMaxPathValue;
			//~ }
		//~ }
//~ //finding best side paths - numbering side paths
		//~ CentralSideCrosses(yMin) = nMaxPathValue + currentMin + 1;

//~ //finding best side paths - rising the count of side paths
    //~ //(CentralSideCrosses,indexesOnTheCentralPath).display();
   //~ //(CentralSideCrosses,indexesOnTheCentralPath,(sidePathsFromLeft+sidePathsFromRight)).display();
  //~ countOfCrossing++;
      //~ }
     //~ //(sidePathsFromLeft+sidePathsFromRight).display("sidePaths");

//~ //  (CentralSideCrosses,indexesOnTheCentralPath,(sidePathsFromLeft+sidePathsFromRight),yCoorOfTheJunctions).display();
//~ for(int ii = 0; ii < permittedNumberOfPaths; ++ii)
//~ {
  //~ if(yCoorOfTheJunctions(ii)>0)
  //~ {
    //~ yCoorOfTheJunctions(ii) = sJunctions.s1Head.yCoor + yCoorOfTheJunctions(ii)
      //~ *(sJunctions.s7Tail.yCoor - sJunctions.s1Head.yCoor)/nImgEdgesHeight;
  //~ }
//~ }
//~ const int nPositionOfPecToAbdInTheJunctions = static_cast<int>( (yCoorOfTheJunctions.get_mul((-yCoorOfTheJunctions).get_threshold(-sPlastronCentreOnRotatedImg.yCoor))-sPlastronCentreOnRotatedImg.yCoor).abs().get_stats()(5) );
//~ const int nPositionOfAbdToFemInTheJunctions = static_cast<int>( (yCoorOfTheJunctions.get_mul((yCoorOfTheJunctions).get_threshold(sPlastronCentreOnRotatedImg.yCoor))-sPlastronCentreOnRotatedImg.yCoor).abs().get_stats()(5) );
//~ m_sJunctions.s4PecToAbd.yCoor = yCoorOfTheJunctions(nPositionOfPecToAbdInTheJunctions);
//~ m_sJunctions.s5AbdToFem.yCoor = yCoorOfTheJunctions(nPositionOfAbdToFemInTheJunctions);
//~ m_sPlastronCentreOnRotatedImg.yCoor = (m_sJunctions.s5AbdToFem.yCoor + m_sJunctions.s4PecToAbd.yCoor)/2;
}

#endif /* _TORTOISE_HPP_ */
