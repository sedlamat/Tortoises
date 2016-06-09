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
    class GeneralHough {
	/*** constructor input parameters: ***/
	// source(with objects to be detected) and template images
	const cv::Mat src_img_orig, tmpl_img;
	cv::Mat src_img, src_edges, tmpl_edges;
	// template reference point
	const cv::Point ref_pt;
	// number of possible orientations of the object in src_img
	const int num_quant_directions;
	// number of possible scales of the object
	const int num_scales;
	// size to which the src_img is resized
	const int max_img_size;
	// min/max size of the object in percents of the src_img (0..1]
	const double max_tmpl_scale;
	const double min_tmpl_scale;
	// thresholds for creating edge images(lower values->more edges)
	const int canny_low_thresh = 25;
	const int canny_high_thresh = 50;

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
	// rotation step of the rotated_r_table in radians
	double rot_rad;
	// flag - if 1 then all accumulators will be displayed
	const bool display_accum;
	// interest point to be located after GeneralHough run
	std::map<std::string, cv::Point> tmpl_interest_pts;

	/*** Declarations for Hough R-table and src Hough points: ***/
	// alias used to store points grouped by gradient direction
	typedef std::vector<std::vector<cv::Point> > HoughTable;
	// cv::Points grouped by orientation
	HoughTable r_table, rotated_r_table, src_hough_pts;

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
	     const int number_of_directions = 4,
	     const int number_of_scales = 20,
	     const int maximum_image_size = 150,
	     const double maximum_template_scale = 1.0,
	     const double minimum_template_scale = 0.3,
	     const bool display_accumulator = 0,
	     std::map<std::string, cv::Point> template_interest_points
			= std::map<std::string, cv::Point>());

	~GeneralHough() {}
	void run();
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
	void fill_src_hough_pts();
	void set_hough_points(HoughTable &hough_points,
				const cv::Mat &img,
				const cv::Mat &img_edges);
	void set_rotated_r_table(double rot_step_rad,
					int num_table_shift);
	bool tmpl_inside_src(cv::Point refer_pt) {return 1;}
	void accumulate();
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
	    const int number_of_directions,
	    const int number_of_scales,
	    const int maximum_image_size,
	    const double maximum_template_scale,
	    const double minimum_template_scale,
	    const bool display_accumulator,
	    std::map<std::string, cv::Point> template_interest_points):
	    src_img_orig(source_image),
	    tmpl_img(template_image),
	    ref_pt(reference_point),
	    num_quant_directions(number_of_directions),
	    num_scales(number_of_scales),
	    max_img_size(maximum_image_size),
	    max_tmpl_scale(maximum_template_scale),
	    min_tmpl_scale(minimum_template_scale),
	    display_accum(display_accumulator),
	    tmpl_interest_pts(template_interest_points)
    {
	// checks num_quant_directions AND sets num_quant_neighbours
	switch (num_quant_directions) {
	    case 4:
		num_quant_neighbours = 0;
		break;
	    case 12:
		num_quant_neighbours = 1;
		break;
	    case 36:
		num_quant_neighbours = 4;
		break;
	    default:
		throw "Wrong num_quant_directions. \
						Can be 4, 12 or 36.";
	}

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
	if (min_tmpl_scale < max_tmpl_scale && 0 < min_tmpl_scale
					     && max_tmpl_scale <= 1) {

	    cv::Mat tmpl_pts;
	    cv::findNonZero(tmpl_img, tmpl_pts);
	    tmpl_bound_rect = cv::boundingRect(tmpl_pts);
	    tmpl_size = tmpl_bound_rect.size();
	    double src_templ_size_koef = max_img_size*1.0 /
			    std::max(tmpl_size.width,tmpl_size.height);
	    max_scale = max_tmpl_scale * src_templ_size_koef;
	    min_scale = min_tmpl_scale * src_templ_size_koef;
	} else {
	    throw "Wrong: min_tmpl_scale param must be < then \
			  max_tmpl_scale and both in (0,1].";
	}
	// checks num_scales AND sets scale_step
	if (1 < num_scales && num_scales <= 1000) {
	    scale_step = (max_scale - min_scale) / (num_scales - 1);
	} else if (num_scales == 1) {
	    scale_step = 0;
	    min_scale = max_scale;
	} else {
	    throw "Wrong: num_scales must be in [1,1000].";
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

	// initialize the rotation step
	rot_rad = 0;
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
    void GeneralHough::set_hough_points(HoughTable &hough_points,
					    const cv::Mat &img,
					    const cv::Mat &img_edges)
    {
	hough_points = HoughTable(num_quant_directions);

	// changes gradient orientations to gradient directions
	cv::Mat orient, orient_adjust, directions;
	orient = sedlamat::get_gradient_orientation(img);
	orient_adjust = (orient >= 180)/255;
	orient_adjust.convertTo(orient_adjust, orient.depth());
	directions = orient + orient_adjust * -180;
	//my::display(directions);
	float quant_width = 180.0 / num_quant_directions;
	//prt(quant_width/2);
	// goes through all edge pixels
	for (int yy = 0; yy < img_edges.rows; yy++) {
	    const uchar *ptr_src_edges_irow = img_edges.ptr<uchar>(yy);
	    for (int xx = 0; xx < img_edges.cols; xx++) {
		if (ptr_src_edges_irow[xx]) {
		    cv::Point pt(xx, yy);
		    // quantizes the edge direction
		    float phi = directions.at<float>(pt);
		    int quant_idx = -1;
		    int idx = 0;
		    if (phi >= 180 - quant_width/2.0 ||
			phi < quant_width/2.0) {
			quant_idx = idx;
		    } else {
			++idx;
			for (float quant_low_bound = quant_width/2.0;
			     quant_low_bound < 180.0 - quant_width;
			     quant_low_bound += quant_width) {
			    //prt(quant_low_bound);
			    if (phi >= quant_low_bound &&
				phi < quant_low_bound + quant_width) {
				quant_idx = idx;
				break;
			    }
			    ++idx;
			}
		    }
		    // fills the HoughTable
		    hough_points[quant_idx].push_back(pt);
		}
	    }
	}
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
	set_hough_points(r_table, tmpl_img, tmpl_edges);

	// shifts all points in HoughTable from the reference point
	for(auto & pts : r_table) {
	    for(auto & pt : pts) {
		pt -= ref_pt;
	    }
	}
    }


    /**
	Fills the src_hough_pts using src_img and src_edges.

	@param void.
	@return void.
    */
    void GeneralHough::fill_src_hough_pts()
    {
	set_hough_points(src_hough_pts, src_img, src_edges);
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
    void GeneralHough::set_rotated_r_table(double rot_step_rad,
					    int num_table_shift)
    {
	rotated_r_table = r_table;

	// shifts the R-Table quants, last-becomes-first shifting
	for (int shift = 0; shift < num_table_shift; ++shift) {
	    std::vector<cv::Point> table_quant;
	    table_quant = rotated_r_table.back();
	    rotated_r_table.pop_back();
	    HoughTable::iterator it = rotated_r_table.begin();
	    rotated_r_table.insert(it, table_quant);
	}

	// rotates the R-Table points by angle_rad (counter-clockwise)
	double cs = std::cos(rot_step_rad);
	double sn = std::sin(rot_step_rad);

	for(auto & table_quant : rotated_r_table) {
	    for(auto & pt : table_quant) {
		int x = pt.x;
		pt.x = cs*pt.x - sn*pt.y;
		pt.y = sn*x + cs*pt.y;
	    }
	}
    }


    /**
	Gets the hough accumulator for current rot_step_rad and
	rotated_r_table. Best results are stored in best_... variables.

	@param void.
	@return void.
    */
    void GeneralHough::accumulate()
    {
	// for all scales
	for (int scale_idx = 0; scale_idx < num_scales; ++scale_idx) {
	    cv::Mat accum = cv::Mat(src_img_size, CV_32F,
					    cv::Scalar_<float>(0.0));
	    double s = min_scale + scale_idx * scale_step;
	    // for all quantified directions
	    for(int quant_idx = 0; quant_idx < num_quant_directions;
							++quant_idx) {
		cv::Mat quant_accum = cv::Mat(src_img_size, CV_32F,
					      cv::Scalar_<float>(0.0));
		std::vector<cv::Point> r_table_pts, src_pts;

		// pooling with neighbouring quants
		for (int quant_shift = -num_quant_neighbours;
				quant_shift <= num_quant_neighbours;
						      ++quant_shift) {
		    int idx = quant_idx + quant_shift;
		    if (idx < 0) {
			idx = idx + num_quant_directions;
		    }
		    if (idx >= num_quant_directions) {
			idx = idx - num_quant_directions;
		    }
		    r_table_pts.insert(r_table_pts.end(),
				    rotated_r_table[idx].begin(),
				    rotated_r_table[idx].end());
		    src_pts.insert(src_pts.end(),
				   src_hough_pts[idx].begin(),
				   src_hough_pts[idx].end());
		}

		float num_quant_pts = r_table_pts.size();
		// for all points in the r_table quant
		for(auto & pt_diff : r_table_pts) {
		    // for all points in the src_hough_points quant
		    for(auto & src_pt : src_pts) {
			// increments accumulator
			cv::Point refer_pt(src_pt - (pt_diff*s));
			// only if inside source image area
			if (src_img_rect.contains(refer_pt) &&
					tmpl_inside_src(refer_pt) ) {
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
		    accum += quant_accum/num_quant_pts;
		}
	    }

	    // smooths the accumulator
	    cv::filter2D(accum, accum, CV_32F, gauss);

	    if (display_accum) sedlamat::display(accum);

	    // finds accumulator maximum
	    double local_max = 0, local_min = 0;
	    cv::Point local_max_pt(0,0), local_min_pt(0,0);
	    cv::minMaxLoc(accum, &local_min, &local_max,
					&local_min_pt, &local_max_pt);

	    if (display_accum) sedlamat::print(local_max);

	    // if accum_max greatest so far, then sets the best_ params
	    if (local_max > best_accum_val) {
		best_accum_val = local_max;
		best_scale = s;
		best_angle = rot_rad;
		best_ref_pt.x = local_max_pt.x;
		best_ref_pt.y = local_max_pt.y;
	    }
	}
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
    void GeneralHough::run()
    {
	this->fill_r_table();
	this->fill_src_hough_pts();
	// accumulate for all rotated r-tables
	const int num_of_rot = num_quant_directions * 2;
	const double rot_step_rad = 2.0 * M_PI / num_of_rot;
	for (int rot_idx = 0; rot_idx < num_of_rot; ++rot_idx) {
	    rot_rad = rot_idx * rot_step_rad;
	    set_rotated_r_table(rot_rad, rot_idx);
	    accumulate();
	}

	double inv_resize_coeff = 1.0 / resize_coeff;
	best_scale *= inv_resize_coeff;
	best_ref_pt *= inv_resize_coeff;

	double cs = std::cos(best_angle);
	double sn = std::sin(best_angle);

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
	double cs = std::cos(best_angle);
	double sn = std::sin(best_angle);

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
	}
	return dst;
    }

} /* namespace sedlamat */

#endif /* _MY_GENERAL_HOUGH_HPP_ */
