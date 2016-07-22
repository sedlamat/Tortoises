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
    std::vector<cv::Point> get_best_path(const cv::Mat &edges);

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
	std::cout << "File Tg#####_info.txt does not exist. It will";
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
    //~ this->display_points_on_plastron();
}


void Tortoise::locate_central_seam()
{
    // central seam localization

    // measures the to-be half width of the central seam stripe
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

    const int new_size = 800;

    const double resize_coeff = new_size * 1.0 / std::max(stripe_w,
							  stripe_h);
    cv::resize(stripe_img, stripe_img, cv::Size(0,0), resize_coeff,
						      resize_coeff);

    stripe_w = stripe_img.cols;
    stripe_h = stripe_img.rows;

    // gets the edge image of the area
    cv::Mat stripe_edges, stripe_gray;
    cv::cvtColor(stripe_img, stripe_gray, CV_BGR2GRAY);
    cv::Canny(stripe_gray, stripe_edges, 0, 0);
    stripe_edges.convertTo(stripe_edges, CV_32F);
    stripe_edges /= 255;

    cv::Mat stripe_edges_orig = stripe_edges.clone();

    sedlamat::display(stripe_edges);

    std::vector<cv::Point> central_seam;

    central_seam = get_best_path(stripe_edges);

    for (auto &pt : central_seam) {
	stripe_img.at<cv::Vec3b>(pt) = _color.red;
    }

    // locates side seams

    stripe_edges = stripe_edges_orig.t();

    std::vector<cv::Point> side_seam;

    side_seam = get_best_path(stripe_edges);

    stripe_img = stripe_img.t();

    for (auto &pt : side_seam) {
	stripe_img.at<cv::Vec3b>(pt) = _color.red;
    }
    stripe_img = stripe_img.t();

    sedlamat::display(stripe_img);
}


std::vector<cv::Point> Tortoise::get_best_path(const cv::Mat &edges)
{
    // edges in 0/1 float matrix

    int edges_w = edges.cols;
    int edges_h = edges.rows;

    const float supreme_value = 10000;

    cv::Mat values(edges.size(), CV_32F, cv::Scalar(supreme_value));

    // prepares distance mask 5rows x 7columns
    cv::Mat dist_mask(5, 7, CV_32F, cv::Scalar(0));

    // fills the dist_mask
    for (int x = 0; x < 7; ++x)  {
	dist_mask.at<float>(0,x) = 1.0;
	dist_mask.at<float>(1,x) = 3.0;
	dist_mask.at<float>(2,x) = 7.0;
	dist_mask.at<float>(3,x) = 12.0;
	dist_mask.at<float>(4,x) = 18.0;
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

    const int low_line_y = edges_h*1/8;
    const int high_line_y = edges_h*7/8;

    const cv::Rect low_line_rect(cv::Point(0,low_line_y),
				 cv::Point(edges_w-1,low_line_y+1));
    const cv::Rect high_line_rect(cv::Point(0,high_line_y),
				  cv::Point(edges_w-1,high_line_y+1));

    edges(low_line_rect) = 1.0;
    edges(high_line_rect) = 1.0;

    // locates the central seam at the low_line position

    values(high_line_rect) = 1.0;

    for (int y = high_line_y-1; y >= 6; --y) {
	for (int x = 3; x < edges_w - 3; ++x) {
	    if (edges.at<float>(y,x)) {
		double min;
		mask_rect.x = x-3;
		mask_rect.y = y+1;
		cv::minMaxLoc(edges(mask_rect).mul(dist_mask) +
			      values(mask_rect), &min);
		values.at<float>(y,x) = min;
	    }
	}
    }

    double low_min, low_max;
    cv::Point low_min_pt;
    cv::minMaxLoc(values(low_line_rect), &low_min, &low_max,
							&low_min_pt);

    std::vector<cv::Point> best_path;

    // locates central seam
    int y = low_line_y;
    int x = low_min_pt.x;

    while (y != high_line_y) {
	double loc_min, loc_max;
	cv::Point loc_min_pt;
	mask_rect.x = x-3;
	mask_rect.y = y+1;
	cv::minMaxLoc(dist_mask + values(mask_rect), &loc_min,
						&loc_max, &loc_min_pt);
	best_path.push_back(cv::Point(x,y));
	x += loc_min_pt.x - 3;
	y += loc_min_pt.y + 1;
    }

    return best_path;
}

#endif /* _TORTOISE_HPP_ */
