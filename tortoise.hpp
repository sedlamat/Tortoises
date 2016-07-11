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
//#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

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

    struct JunctionIndices {
	const int j0Head = 0;
	const int j1GulHum = 1;
	const int j2HumPec = 2;
	const int j3PecAbd = 3;
	const int j4AbdFem = 4;
	const int j5FemAna = 5;
	const int j6Tail = 6;
    } _idx;

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
    void locate_central_seam();
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
	    _color(),
	    _idx()
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
    this->locate_central_seam();
    this->display_points_on_plastron();

}


void Tortoise::locate_central_seam()
{
    std::cout << "Abdominal junctions localization" << std::endl;

    // calculates the to-be half width of the central seam stripe
    const int abd_seam_length = _l_juncs[_idx.j4AbdFem].y -
				_l_juncs[_idx.j3PecAbd].y;

    // prepares the area rectangle of the central seam stripe
    cv::Rect stripe_rect(
		cv::Point(_l_juncs[_idx.j0Head].x - abd_seam_length/2,
		          _l_juncs[_idx.j0Head].y ),
		cv::Point(_l_juncs[_idx.j6Tail].x + abd_seam_length/2,
		          _l_juncs[_idx.j6Tail].y ));

    // gets the area of the central seam stripe
    cv::Mat stripe_img = _plastron_img(stripe_rect).clone();


    // resizes the area
    int stripe_w = stripe_img.cols;
    int stripe_h = stripe_img.rows;
    const double resize_coeff = 800.0 / std::max(stripe_w,stripe_h);
    cv::resize(stripe_img, stripe_img, cv::Size(0,0), resize_coeff,
						resize_coeff);
    stripe_w = stripe_img.cols;
    stripe_h = stripe_img.rows;

    // gets the edge image of the area
    cv::Mat stripe_edges;
    cv::Canny(stripe_img, stripe_edges, 0, 0);
    stripe_edges.convertTo(stripe_edges, CV_32F);
    stripe_edges /= 255;


    cv::Mat stripe_values(stripe_edges.size(), CV_32F,
						    cv::Scalar(10000));

    // prepares distance mask 5rows x 7columns
    cv::Mat dist_mask(5, 7, CV_32F, cv::Scalar(0));

    /* fills the rows of dist_mask as follows
     * [ 4, 3, 2, 1, 2, 3, 4;
     *   7, 6, 5, 4, 5, 6, 7;
     *  12,11,10, 9,10,11,12;
     *  19,18,17,16,17,18,19;
     *  28,27,26,25,26,27,28 ]
    **/

    for (int x = 0; x < 7; ++x)  {
	for (int y = 0; y < 5; ++y) {
	    dist_mask.at<float>(y,x) = pow((y+1),2.0);
	}
    }
    for (int y = 0; y < 5; ++y) {
      dist_mask.at<float>(y,0) += 3;
      dist_mask.at<float>(y,6) += 3;
      dist_mask.at<float>(y,1) += 2;
      dist_mask.at<float>(y,5) += 2;
      dist_mask.at<float>(y,2) += 1;
      dist_mask.at<float>(y,4) += 1;
    }

    cv::Rect mask_rect(cv::Point(0,0), cv::Size(7,5));
    std::cout << mask_rect << std::endl;
    //~ cv::flip(dist_mask, dist_mask, 0);

    const int low_line_y = stripe_h*1/6;
    const int mid_line_y = stripe_h*3/6;
    const int high_line_y = stripe_h*5/6;

    const cv::Rect low_line_rect(cv::Point(0,low_line_y),
				 cv::Point(stripe_w-1,low_line_y+1));
    const cv::Rect mid_line_rect(cv::Point(0,mid_line_y),
				 cv::Point(stripe_w-1,mid_line_y+1));
    const cv::Rect high_line_rect(cv::Point(0,high_line_y),
				 cv::Point(stripe_w-1,high_line_y+1));

    stripe_edges(low_line_rect) = 1.0;
    stripe_edges(mid_line_rect) = 1.0;
    stripe_edges(high_line_rect) = 1.0;

    // locates the central seam at the low_line position

    stripe_values(mid_line_rect) = 1.0;

    for (int y = mid_line_y-1; y >= 6 ; --y) {
	for (int x = 3; x < stripe_w - 3; ++x) {
	    if (stripe_edges.at<int>(y,x)) {
		double min;
		mask_rect.x = x-3;
		mask_rect.y = y+1;
		//~ std::cout << x << std::endl;
		//~ std::cout << y << std::endl;
		//~ std::cout << mask_rect << std::endl;
		//~ sedlamat::display(stripe_values(mask_rect));
		//~ sedlamat::display(stripe_edges(mask_rect));
		//~ sedlamat::display(dist_mask);
		//~ sedlamat::display((stripe_edges(mask_rect).mul(dist_mask) +
			      //~ stripe_values(mask_rect)));
		cv::minMaxLoc(stripe_edges(mask_rect).mul(dist_mask) +
			      stripe_values(mask_rect), &min);
		stripe_values.at<float>(y,x) = min;
	    }
	}
    }
    cv::Mat high_values;
    cv::Mat(stripe_values < 9000).convertTo(high_values, CV_32F);
    sedlamat::display(stripe_values.mul(high_values));
    //~ cv::Mat non_zero_coordinates;
    //~ findNonZero(img, nonZeroCoordinates);
    //~ for (int i = 0; i < nonZeroCoordinates.total(); i++ ) {
        //~ cout << "Zero#" << i << ": " << nonZeroCoordinates.at<Point>(i).x << ", " << nonZeroCoordinates.at<Point>(i).y << endl;
    //~ }
    //~ return 0;
    //~ stripe_values
    //~ abd_values(cv::Range(search_start_y_idx,abd_edges.size().height), cv::Range(0,3)) -= 10000;
    //~ abd_values(cv::Range(search_start_y_idx,abd_edges.size().height), cv::Range(abd_edges.size().width - 3,abd_edges.size().width)) -= 10000;
    //~ abd_values(cv::Range(0,search_start_y_idx-1), cv::Range::all()) -= 10000;
    //~ abd_values(cv::Range(search_start_y_idx-1,abd_edges.size().height), cv::Range::all()) = abd_values(cv::Range(search_start_y_idx-1,abd_edges.size().height), cv::Range::all()).mul(abd_edges(cv::Range(search_start_y_idx-1,abd_edges.size().height), cv::Range::all()));
    //abd_values = abd_values.mul(abd_edges);
    //std::cout << dist_mask << std::endl;
    //~ sedlamat::display(abd_values);
    //sedlamat::display(abd_edges);
      //~ //distance

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


cv::Mat Tortoise::get_skeletonized_min_dist(const cv::Mat &edge_img)
{
    cv::Mat edges = edge_img.clone();
    this->eliminate_one_pix_edges(edges);
    this->eliminate_two_pix_edges(edges);
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

void Tortoise::eliminate_one_pix_edges(cv::Mat &edges)
{
    int edges_w = edges.size().width;
    int edges_h = edges.size().height;
    for (int y = 1; y < edges_h-1; ++y) {
	for (int x = 1; x < edges_w-1; ++x) {
	    if(edges.at<int>(y,x)) {
		if(!edges.at<int>(y,x-1) &&
		   !edges.at<int>(y-1,x-1) &&
		   !edges.at<int>(y+1,x-1) &&
		   !edges.at<int>(y-1,x) &&
		   !edges.at<int>(y,x+1) &&
		   !edges.at<int>(y-1,x+1) &&
		   !edges.at<int>(y+1,x+1) &&
		   !edges.at<int>(y+1,x)) {

		    edges.at<int>(y,x) = 0;
		}
	    }
	}
    }
}

void Tortoise::eliminate_two_pix_edges(cv::Mat &imgToDoItOn)
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


}

#endif /* _TORTOISE_HPP_ */
