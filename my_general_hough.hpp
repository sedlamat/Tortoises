/**
    my_general_hough.hpp

    Implementation of the generalized Hough transform, using OpenCV
    library.

    Do not forget that images are displayed with y-axis
    pointing down, it might fool you if you try to visualize
    individual steps of the algorithms below.

    The code below is dependent on the "my_img_proc.hpp", but it can
    be easily undone (only few procedures from "my_img_proc.hpp" are
    called).

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
#include <ctime>
#include <limits>
#include <deque>

/* THIRD PARTY LIBRARIES */
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"

/*
    Beware:
	cv::Point is alias for cv::Point_<int>
	cv::Size is alias for cv::Size_<int>
	cv::Rect is alias for cv::Rect_<int>
*/

namespace sedlamat
{
    /* deque with negative indexing and rotation of elements */
    template <typename T> class DequeNegIdx: public std::deque<T> {
    public:
	DequeNegIdx(): std::deque<T>() {}
	DequeNegIdx(size_t size): std::deque<T>(size) {}
	~DequeNegIdx() {};

	T &operator[](int idx)
	{
	    size_t size = this->size();
	    idx = idx % size;
	    if (idx < 0) idx += size;
	    return this->at(idx);
	}
	// only prefix increment and decrement operators are defined
	void operator++()
	{
	    this->push_front(this->back());
	    this->pop_back();
	}
	void operator--()
	{
	    this->push_back(this->front());
	    this->pop_front();
	}
	void shift(int num_shifts)
	{
	    if (num_shifts < 0) {
		for (; num_shifts < 0; ++num_shifts) {
		    operator--();
		}
	    } else {
		for (; num_shifts > 0; --num_shifts) {
		    operator++();
		}
	    }
	}
    };

    typedef DequeNegIdx<std::vector<cv::Point>> HoughTable;


    class GeneralHough {
	/*** constructor input parameters: ***/
	// source(with objects to be detected) and template images
	const cv::Mat src_img_orig, tmpl_img;
	cv::Mat src_img, src_edges, tmpl_edges;
	// template reference point
	const cv::Point ref_pt;

	std::vector<int> angles;
	std::vector<float> scales;
	// number of possible scales of the object
	const int num_scales;
	// size to which the src_img is resized
	const int max_img_size;
	// thresholds for creating edge images(lower values->more edges)
	const int canny_low_thresh = 50;
	const int canny_high_thresh = 150;



	/*** other parameters set up in contructor: ***/
	// kernel for smoothing the accumulator
	cv::Mat gauss;
	// Size and Rectangle are of the src/tmpl_img
	cv::Size src_img_size, tmpl_img_size;
	cv::Rect src_img_rect, tmpl_img_rect;
	// scaling step in the accumulator
	float scale_step;
	// number of neighbouring quants for pooling in the accumulator
	int num_quant_neighbours;
	// template size in pixels (area with non-zero pixels)
	cv::Rect tmpl_bound_rect;
	cv::Size tmpl_size;
	// resize coefficient from original source image
	double resize_coeff;
	// min/max_scaling of the template
	double max_scale;
	double min_scale;
	// flag - if 1 then all accumulators will be displayed
	const bool display_accum;
	// interest point to be located after GeneralHough run
	std::map<std::string, cv::Point> tmpl_interest_pts;

	// cv::Points grouped by edge directions
	HoughTable r_table, src_hough_table;

	/*** constants with int and double minima: ***/
	const double MIN_DOUBLE = -std::numeric_limits<double>::max();
	const int MIN_INT = -std::numeric_limits<int>::max();

	/*** Variables for the best fit of the hough transform: ***/
	double best_accum_val;
	double best_scale;
	double best_angle;
	cv::Point best_ref_pt;

    public:
	GeneralHough(const cv::Mat &source_image,
	     const cv::Mat &template_image,
	     const cv::Point reference_point,
	     const std::vector<int> &angles_vec,
	     const int number_of_scales = 20,
	     const int maximum_image_size = 150,
	     const double maximum_template_scale = 1.0,
	     const double minimum_template_scale = 0.3,
	     const bool display_accumulator = 0,
	     std::map<std::string, cv::Point> template_interest_points
			= std::map<std::string, cv::Point>());

	~GeneralHough() {}
	void detect();
	cv::Mat get_result_img() const;
	double get_best_accum_val() const { return best_accum_val; }
	double get_best_angle() const { return best_angle; }
	double get_best_scale() const { return best_scale; }
	cv::Point get_best_ref_pt() const { return best_ref_pt; }
	cv::Mat get_src_img() const { return src_img; }
	cv::Mat get_tmpl_img() const { return tmpl_img; }
	cv::Mat get_tmpl_edges() const { return tmpl_edges;}
	cv::Mat get_src_edges() const { return src_edges; }

    private:
	void fill_r_table();
	void fill_src_hough_table();
	HoughTable get_hough_table(const cv::Mat &img,
				   const cv::Mat &img_edges);
	HoughTable get_rotated_r_table(const int angle);
	void accumulate(HoughTable &rotated_r_table, const int angle);
	HoughTable get_quanted_table(HoughTable &table,
				     const int angle,
				     const int num_quants = 4);
    };

    /**
	GeneralHough constructor. Sets various parameters, see
	class declaration.

	@param see class declaration.
	@return void.
    */
    GeneralHough::GeneralHough(const cv::Mat &source_image,
	    const cv::Mat &template_image,
	    const cv::Point reference_point,
	    const std::vector<int> &angles_vec,
	    const int number_of_scales,
	    const int maximum_image_size,
	    const double maximum_template_scale,
	    const double minimum_template_scale,
	    const bool display_accumulator,
	    std::map<std::string, cv::Point> template_interest_points):
	    src_img_orig(source_image),
	    tmpl_img(template_image),
	    ref_pt(reference_point),
	    angles(angles_vec),
	    num_scales(number_of_scales),
	    max_img_size(maximum_image_size),
	    display_accum(display_accumulator),
	    tmpl_interest_pts(template_interest_points)
    {
	// checks max_img_size
	if ( ! (1 <= max_img_size && max_img_size <= 1000) ) {
	    throw "Wrong: max_img_size must be in [1,1000]";
	}

	/* checks src/tempalate_img AND resizes src_img_orig
	   AND sets src/templ_size/rect */
	if (!src_img_orig.empty() && !tmpl_img.empty()) {
	    tmpl_img_size = tmpl_img.size();
	    tmpl_img_rect = cv::Rect(cv::Point(0,0), tmpl_img_size);

	    // Warns if ref_pt outside template_img
	    if (!ref_pt.inside(tmpl_img_rect)) {
		throw "Wrong: Reference point outside template image.";
	    }

	    // Resizes the original source image
	    int src_w = src_img_orig.cols;
	    int src_h = src_img_orig.rows;
	    resize_coeff = max_img_size * 1.0 / std::max(src_w,src_h);
	    cv::resize(src_img_orig, src_img, cv::Size(0,0),
					    resize_coeff, resize_coeff);

	    src_img_size = src_img.size();
	    src_img_rect = cv::Rect(cv::Point(0,0), src_img_size);
	    cv::GaussianBlur(src_img, src_img, cv::Size(0,0), 0.5);
	} else {
	    throw "Source and/or template image is empty!";
	}

	// sets the src_edges
	CV_Assert(src_img.channels() == 3 || src_img.channels() == 1);
	if (src_img.channels() == 3) {
	    cv::cvtColor(src_img, src_edges, cv::COLOR_BGR2GRAY);
	} else {
	    src_img.copyTo(src_edges);
	}
	cv::Canny(src_edges, src_edges,	canny_low_thresh,
						    canny_high_thresh);

	// sets the tmpl_edges
	CV_Assert(tmpl_img.channels() == 3 || tmpl_img.channels() == 1);
	if (tmpl_img.channels() == 3) {
	    cv::cvtColor(tmpl_img, tmpl_edges, cv::COLOR_BGR2GRAY);
	} else {
	    tmpl_img.copyTo(tmpl_edges);
	}
	cv::Canny(tmpl_edges, tmpl_edges, canny_low_thresh,
						    canny_high_thresh);


	// checks min/max_tmpl_scale AND sets min/max_scale
	if (minimum_template_scale < maximum_template_scale &&
				0 < minimum_template_scale &&
					maximum_template_scale <= 1) {

	    cv::Mat tmpl_pts;
	    cv::findNonZero(tmpl_img, tmpl_pts);
	    tmpl_bound_rect = cv::boundingRect(tmpl_pts);
	    tmpl_size = tmpl_bound_rect.size();
	    double src_templ_size_koef = max_img_size*1.0 /
			    std::max(tmpl_size.width,tmpl_size.height);

	    //~ for (int scale_idx = 0; scale_idx < num_scales; ++scale_idx) {
	    //~ cv::Mat accum = cv::Mat(src_img_size, CV_32F,
					    //~ cv::Scalar_<float>(0.0));
	    //~ double s = min_scale + scale_idx * scale_step;
	    max_scale = maximum_template_scale * src_templ_size_koef;
	    min_scale = minimum_template_scale * src_templ_size_koef;
	} else {
	    throw "Wrong: min_tmpl_scale param must be < then \
			  max_tmpl_scale and both in (0,1].";
	}

	// checks num_scales AND sets scale_step
	double scale_step;
	if (1 < num_scales && num_scales <= 1000) {
	    scale_step = (max_scale - min_scale) / (num_scales - 1);
	} else if (num_scales == 1) {
	    scale_step = 0;
	    min_scale = max_scale;
	} else {
	    throw "Wrong: num_scales must be in [1,1000].";
	}

	for (int ii = 0; ii < num_scales; ++ii) {
	    scales.push_back(min_scale + ii * scale_step);
	}

	// sets the values of the best fit
	best_accum_val = MIN_DOUBLE;
	best_scale = MIN_DOUBLE;
	best_angle = MIN_DOUBLE;
	best_ref_pt = cv::Point(MIN_INT, MIN_INT);

	// prepares gaussian kernel for smothing of the accumulator
	int gauss_size = 9.0*max_img_size/100;
	double gauss_sigma = 1.0*max_img_size/100;
	gauss_size = (gauss_size % 2) ? gauss_size : gauss_size + 1;
	gauss = cv::getGaussianKernel(gauss_size, gauss_sigma, CV_32F);
	gauss = gauss * gauss.t();



    }




    /**
	Fills in the hough_points container with Hough points (edge
	points with their orientation.

	@param hough_points - HoughTable container to be filled.
	@param img - Image which will provide edge points and
		     orientations.
	@param img_edges - Edge image of img parameter.
	@return void.
    */
    //~ void GeneralHough::set_hough_points(HoughTable &hough_points,
					    //~ const cv::Mat &img,
					    //~ const cv::Mat &img_edges)
    //~ {
	//~ hough_points = HoughTable(num_quant_directions);
//~
	//~ // changes gradient orientations to gradient directions
	//~ cv::Mat orient, orient_adjust, directions;
	//~ orient = sedlamat::get_gradient_orientation(img);
	//~ orient_adjust = (orient >= 180)/255;
	//~ orient_adjust.convertTo(orient_adjust, orient.depth());
	//~ directions = orient + orient_adjust * -180;
	//~ directions.convertTo(directions, CV_8UC1);
	//~ sedlamat::print(directions);
	//~ sedlamat::display(directions);
	//~ float quant_width = 180.0 / num_quant_directions;
	//~ //prt(quant_width/2);
	//~ // goes through all edge pixels
	//~ for (int yy = 0; yy < img_edges.rows; yy++) {
	    //~ const uchar *ptr_src_edges_irow = img_edges.ptr<uchar>(yy);
	    //~ for (int xx = 0; xx < img_edges.cols; xx++) {
		//~ if (ptr_src_edges_irow[xx]) {
		    //~ cv::Point pt(xx, yy);
		    //~ // quantizes the edge direction
		    //~ float phi = directions.at<float>(pt);
		    //~ int quant_idx = -1;
		    //~ int idx = 0;
		    //~ if (phi >= 180 - quant_width/2.0 ||
			//~ phi < quant_width/2.0) {
			//~ quant_idx = idx;
		    //~ } else {
			//~ ++idx;
			//~ for (float quant_low_bound = quant_width/2.0;
			     //~ quant_low_bound < 180.0 - quant_width;
			     //~ quant_low_bound += quant_width) {
			    //~ //prt(quant_low_bound);
			    //~ if (phi >= quant_low_bound &&
				//~ phi < quant_low_bound + quant_width) {
				//~ quant_idx = idx;
				//~ break;
			    //~ }
			    //~ ++idx;
			//~ }
		    //~ }
		    //~ // fills the HoughTable
		    //~ hough_points[quant_idx].push_back(pt);
		//~ }
	    //~ }
	//~ }
    //~ }


    /**
	Fills in the hough_points container with Hough points (edge
	points with their orientation.

	@param hough_points - HoughTable container to be filled.
	@param img - Image which will provide edge points and
		     orientations.
	@param img_edges - Edge image of img parameter.
	@return void.
    */
    HoughTable GeneralHough::get_hough_table(const cv::Mat &img,
					      const cv::Mat &img_edges)
    {
	HoughTable hough_table(180); // 0...179
//~ //~
	// changes gradient orientations to gradient directions
	cv::Mat orient, orient_adjust, directions;
	orient = sedlamat::get_gradient_orientation(img);
	orient_adjust = (orient >= 180)/255;
	orient_adjust.convertTo(orient_adjust, orient.depth());
	directions = orient + orient_adjust * -180;
	directions.convertTo(directions, CV_8UC1);
	// goes through all edge pixels
	for (int yy = 0; yy < img_edges.rows; ++yy) {
	    for (int xx = 0; xx < img_edges.cols; ++xx) {
		if (img_edges.at<uchar>(yy, xx)) {
		    cv::Point pt(xx, yy);
		    int phi = directions.at<uchar>(pt);
		    hough_table[phi].push_back(pt);
		}
	    }
	}
	//~ for (int ii = 0; ii < 180; ++ii) {
	//~ sedlamat::visualize_points(hough_table[ii], cv::Size(500,500), cv::Point(150,150));
	//~ }
	//sedlamat::display(directions);
	return hough_table;
    }

    /**
	Fills the R-table. First fills the r_table with hough points
	from tmpl_img and tmpl_edges and then shifts them by
	the reference point.

	@param void.
	@return void.
    */
    void GeneralHough::fill_r_table()
    {
	r_table = get_hough_table(tmpl_img, tmpl_edges);

	// shifts all points in HoughTable from the reference point
	for(auto &pts : r_table) {
	    for(auto &pt : pts) {
		pt -= ref_pt;
	    }
	}
    }


    /**
	Fills the src_hough_pts using src_img and src_edges.

	@param void.
	@return void.
    */
    void GeneralHough::fill_src_hough_table()
    {
	src_hough_table = get_hough_table(src_img, src_edges);
    }


    /**
	Sets rotated_r_table by shifting the r_table by num_table_shift
	and rotating the points by rot_step_rad.

	@param rot_step_rad - Angle in radians by which points are
			   rotated.
	@param num_table_shift - Number of times the R-Table quants
				 are shifted.
	@return void.
    */
    HoughTable GeneralHough::get_rotated_r_table(const int angle)
    {
	HoughTable rotated_r_table = r_table;

	rotated_r_table.shift(angle);

	// rotates the R-Table points by angle in rad(counter-clockwise)
	double cs = std::cos(angle * M_PI / 180.0);
	double sn = std::sin(angle * M_PI / 180.0);

	for(auto &pts : rotated_r_table) {
	    for(auto &pt : pts) {
		int x = pt.x;
		pt.x = cs*pt.x - sn*pt.y;
		pt.y = sn*x + cs*pt.y;
	    }
	}
	return rotated_r_table;
    }


    /**
	Gets the hough accumulator for current rot_step_rad and
	rotated_r_table. Best results are stored in best_... variables.

	@param void.
	@return void.
    */
    void GeneralHough::accumulate(HoughTable &rotated_r_table,
				  const int angle)
    {

	HoughTable quanted_rotated_r_table, quanted_src_hough_table;
	quanted_rotated_r_table = get_quanted_table(rotated_r_table,
						    angle);
	quanted_src_hough_table = get_quanted_table(src_hough_table,
						    angle);
	int num_quant = quanted_src_hough_table.size();

	for (auto &s : scales) {
	    cv::Mat accum = cv::Mat(src_img_size, CV_32F,
					    cv::Scalar_<float>(0.0));
	    for (int quant = 0; quant < num_quant; ++quant) {
		//sedlamat::visualize_points(quanted_src_hough_table[quant], cv::Size(200,200));
		//sedlamat::visualize_points(quanted_rotated_r_table[quant], cv::Size(200,200));
		//sedlamat::visualize_points(quanted_rotated_r_table[quant], cv::Size(200,200));

		cv::Mat quant_accum = cv::Mat(src_img_size, CV_32F,
					      cv::Scalar_<float>(0.0));
		double num_quant_pts
				= quanted_rotated_r_table[quant].size();
		for(const auto &pt_diff : quanted_rotated_r_table[quant]) {
		    for(const auto &src_pt : quanted_src_hough_table[quant]) {
			cv::Point refer_pt(src_pt - (pt_diff*s));
			if (src_img_rect.contains(refer_pt)) {
			    float &accum_pt =
					quant_accum.at<float>(refer_pt);
			    // should not be more that points in r_table
			    if (accum_pt < num_quant_pts) {
				accum_pt += 1.0;
			    }
			}
		    }
		}
		// divide by the total number of points in the quant
		// big quants are not more important
		if (num_quant_pts > 0) {
		    accum += quant_accum*1.0/num_quant_pts;
		}
	    }

	    // smooths the accumulator
	    cv::filter2D(accum, accum, CV_32F, gauss);

	    if (display_accum) sedlamat::display(accum);

	    // finds accumulator maximum
	    double local_max, local_min;
	    cv::Point local_max_pt, local_min_pt;
	    cv::minMaxLoc(accum, &local_min, &local_max,
					&local_min_pt, &local_max_pt);

	    if (display_accum) sedlamat::print(local_max);

	    // if accum_max greatest so far, then sets the best_ params
	    if (local_max > best_accum_val) {
		best_accum_val = local_max;
		best_scale = s;
		best_angle = angle;
		best_ref_pt.x = local_max_pt.x;
		best_ref_pt.y = local_max_pt.y;
	    }
	}
    }

    HoughTable GeneralHough::get_quanted_table(HoughTable &table,
					       const int angle,
					       const int num_quants)
    {
	//~ for (int ii = 0; ii < 180; ++ii) {
	    //~ sedlamat::visualize_points(table[ii], cv::Size(500,500), cv::Point(150,150));
	//~ }
	HoughTable quanted_table(num_quants);
	int quant_idx = 0;
	const int quant_size = 180 / num_quants;
	int idx = angle < 0 ? angle % 180 + 180 : angle % 180;
	idx -= quant_size / 2;
	for (int ii = 1; ii <= 180; ++ii) {
	    //sedlamat::visualize_points(table[ii], cv::Size(500,500), cv::Point(150,150));
	    //std::cout << idx << std::endl;
	    if (!(ii % quant_size)) {
		++quant_idx;
		//std::cout << quanted_table[quant_idx-1] << std::endl;
 		//sedlamat::visualize_points(quanted_table[quant_idx-1], cv::Size(500,500), cv::Point(150,150));

	    }
	    //std::cout << table[idx] << std::endl;
	    quanted_table[quant_idx].insert(
				    quanted_table[quant_idx].end(),
				    table[idx].begin(),
				    table[idx].end());
	    //std::cout << quanted_table[quant_idx] << std::endl;
	    ++idx;
	}
	//~ for (int ii = 0; ii < 180; ++ii) {
	    //~ sedlamat::visualize_points(quanted_table[ii], cv::Size(500,500), cv::Point(150,150));
	//~ }
	//std::cout << quant_size << std::endl;
	return quanted_table;
    }

    /**
	Finds and sets the best_ params of the Hough transform
	by rotating the r_table and filling the Hough accumulator.
	After execution, best_accum_val/angle/ref_pt/scale will
	hold the information of the estimated location of the
	target object(s).

	@param void.
	@return void.
    */
    void GeneralHough::detect()
    {
	this->fill_r_table();
	this->fill_src_hough_table();
	// accumulate for all rotated r-tables
	for (auto angle : angles) { //all angles in degrees integers!!
	    HoughTable rotated_r_table(get_rotated_r_table(angle));
	    accumulate(rotated_r_table, angle);
	}

	//~ const int num_of_rot = num_quant_directions * 2;
	//~ const double rot_step_rad = 2.0 * M_PI / num_of_rot;
	//~ for (int rot_idx = 0; rot_idx < num_of_rot; ++rot_idx) {
	    //~ rot_rad = rot_idx * rot_step_rad;
	//~ }

	double inv_resize_coeff = 1.0 / resize_coeff;
	best_scale *= inv_resize_coeff;
	best_ref_pt *= inv_resize_coeff;

	if (!tmpl_interest_pts.empty()) {
	    double cs = std::cos(best_angle*M_PI/180.0);
	    double sn = std::sin(best_angle*M_PI/180.0);

	    std::map<std::string, cv::Point>::iterator it, it_end;
	    it = tmpl_interest_pts.begin();
	    it_end = tmpl_interest_pts.end();
	    while (it != it_end) {
		it->second -= ref_pt;
		int x = it->second.x;
		it->second.x = cs*it->second.x - sn*it->second.y;
		it->second.y = sn*x + cs*it->second.y;
		it->second *= best_scale;
		it->second += best_ref_pt;
		++it;
	    }
	}
    }

    /**
	Gets the source image with the best template match marked on it.

	@param void.
	@return Image with best template fit.
    */
    cv::Mat GeneralHough::get_result_img() const
    {
	cv::Mat dst;
	//make a deep copy
	src_img_orig.copyTo(dst);
	cv::Rect src_img_orig_rect(cv::Point(0,0), src_img_orig.size());

	// rotates and shift the R-Table points
	double cs = std::cos(best_angle*M_PI/180.0);
	double sn = std::sin(best_angle*M_PI/180.0);

	double inv_resize_coeff = 1.0 / resize_coeff;
	cv::Point shift_pt(inv_resize_coeff,inv_resize_coeff);

	for(auto table_quant : r_table) {
	    for(auto pt : table_quant) {
		int x = pt.x;
		pt.x = cs*pt.x - sn*pt.y;
		pt.y = sn*x + cs*pt.y;
		pt = pt * best_scale;
		pt = pt + best_ref_pt;
		if (src_img_orig_rect.contains(pt)) {
		    cv::rectangle(dst,
				  pt-shift_pt,
				  pt+shift_pt,
				  cv::Scalar(0,0,255), CV_FILLED);
		}
	    }
	}
	if (!tmpl_interest_pts.empty()) {
	    std::map<std::string, cv::Point>::const_iterator it, it_end;
	    it = tmpl_interest_pts.cbegin();
	    it_end = tmpl_interest_pts.cend();
	    while (it != it_end) {
		if (src_img_orig_rect.contains(it->second)) {
		    cv::rectangle(dst,
				  it->second-shift_pt,
				  it->second+shift_pt,
				  cv::Scalar(0,255,0), CV_FILLED);
		}
		++it;
	    }
	}
	return dst;
    }

} /* namespace sedlamat */

#endif /* _MY_GENERAL_HOUGH_HPP_ */
