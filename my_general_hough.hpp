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
	cv::Mat src_img, tmpl_img, src_edges, tmpl_edges;
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
	const int canny_low_thresh;
	const int canny_high_thresh;

	/*** other parameters set up in contructor: ***/
	// kernel for smoothing the accumulator
	cv::Mat gauss;
	// Size and Rectangle are of the src/tmpl_img
	cv::Size src_size, tmpl_img_size;
	cv::Rect src_img_rect, tmpl_img_rect;
	// scaling step in the accumulator
	float scale_step;
	// number of neighbouring quants for pooling in the accumulator
	int num_quant_neighbours;
	// template size in pixels (area with non-zero pixels)
	cv::Rect tmpl_bound_rect;
	cv::Size tmpl_size;
	// min/max_scaling of the template
	double max_scale;
	double min_scale;
	// rotation step of the rotated_r_table in radians
	double rot_rad;
	// flag - if 1 then all accumulators will be displayed
	const bool display_accum;

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
		     const cv::Mat &source_edges = cv::Mat(),
		     const cv::Mat &template_edges = cv::Mat(),
		     const int number_of_directions = 4,
		     const int number_of_scales = 20,
		     const int maximum_image_size = 150,
		     const double maximum_template_scale = 1.0,
		     const double minimum_template_scale = 0.3,
		     const int canny_lower_threshold = 20,
		     const int canny_higher_threshold = 50,
		     const bool display_accumulator = 0);

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
			    const cv::Mat &source_edges,
			    const cv::Mat &template_edges,
			    const int number_of_directions,
			    const int number_of_scales,
			    const int maximum_image_size,
			    const double maximum_template_scale,
			    const double minimum_template_scale,
			    const int canny_lower_threshold,
			    const int canny_higher_threshold,
			    const bool display_accumulator):
			    src_img(source_image),
			    tmpl_img(template_image),
			    src_edges(source_edges),
			    tmpl_edges(template_edges),
			    ref_pt(reference_point),
			    num_quant_directions(number_of_directions),
			    num_scales(number_of_scales),
			    max_img_size(maximum_image_size),
			    max_tmpl_scale(maximum_template_scale),
			    min_tmpl_scale(minimum_template_scale),
			    canny_low_thresh(canny_lower_threshold),
			    canny_high_thresh(canny_higher_threshold),
			    display_accum(display_accumulator)
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
	if ( ! (1 <= >max_img_size && max_img_size <= 1000) ) {
	    throw "Wrong: max_img_size must be in [1,1000]";
	}

	/* checks src/tempalate_img AND resizes src_img
	   AND sets src/templ_size/area */
	if (!src_img.empty() && !tmpl_img.empty()) {
	    tmpl_img_size = tmpl_img.size();
	    tmpl_img_rect = cv::Rect(cv::Point(0,0), tmpl_img_size);

	    // Warns if ref_pt outside template_img
	    if (!ref_pt.inside(tmpl_img_rect)) {
		throw "Wrong: Reference point outside template image. \
		       If it has to be outside, make changes in the  \
		       GeneralHough constructor and in the \
		       accumulator() method.";
	    }

	    // Resizes the source image
	    int src_w = src_img.cols;
	    int src_h = src_img.rows;
	    double resize_koef = max_img_size * 1.0 /
					         std::max(src_w,src_h);
	    cv::resize(src_img, src_img, cv::Size(0,0), resize_koef,
							  resize_koef);

	    src_size = src_img.size();
	    src_img_rect = cv::Rect(cv::Point(0,0), src_size);
	} else {
	    throw "Source and/or template image is empty!";
	}

	// checks canny_low/high_thresh
	if ( ! (this->canny_low_thresh <= this->canny_high_thresh
		&& this->canny_low_thresh >= 0
		&& this->canny_high_thresh <= 1000) ) {
	    throw "Wrong: canny_low_thresh param must be <= then \
			  canny_high_thresh and both in [0,1000].";
	}

	// checks and sets the src_edges
	if (this->src_edges.empty()) {
	    cv::cvtColor(this->src_img, this->src_edges,
						cv::COLOR_BGR2GRAY);
	    cv::Canny(this->src_edges, this->src_edges,
				canny_low_thresh, canny_high_thresh);
	} else if (this->src_edges.size() != this->src_size) {
	    throw "Size of src image and src edges is not equal";
	}

	// checks and sets the tmpl_edges
	if (this->tmpl_edges.empty()) {
	    if (this->tmpl_img.channels() == 3) {
		cv::cvtColor(this->tmpl_img, this->tmpl_edges,
						cv::COLOR_BGR2GRAY);
	    } else {
		this->tmpl_img.copyTo(this->tmpl_edges);
	    }
	    cv::Canny(this->tmpl_edges, this->tmpl_edges,
				canny_low_thresh, canny_high_thresh);
	} else if (this->tmpl_edges.size()
				    != this->tmpl_img_size) {
	    throw "Size of src image and src edges is not equal";
	}

	// checks min/max_tmpl_scale AND sets min/max_scale
	if (this->min_tmpl_scale < this->max_tmpl_scale
	    && 0 < this->min_tmpl_scale
	    && this->max_tmpl_scale <= 1)
	{
	    cv::Mat tmpl_pts;
	    cv::findNonZero(this->tmpl_img, tmpl_pts);
	    this->tmpl_bound_rect = cv::boundingRect(tmpl_pts);
	    this->tmpl_size = this->tmpl_bound_rect.size();
	    double src_templ_size_koef = this->max_img_size*1.0 /
			std::max(this->tmpl_size.width,this->tmpl_size.height);
	    this->max_scale = max_tmpl_scale * src_templ_size_koef;
	    this->min_scale = min_tmpl_scale * src_templ_size_koef;
	} else {
	    throw "Wrong: min_tmpl_scale param must be < then \
			  max_tmpl_scale and both in (0,1].";
	}
	// checks num_scales AND sets scale_step
	if (1 < this->num_scales && this->num_scales <= 1000) {
	    this->scale_step = (this->max_scale - this->min_scale) /
						(this->num_scales - 1);
	} else if (this->num_scales == 1) {
	    this->scale_step = 0;
	    this->min_scale = this->max_scale;
	} else {
	    throw "Wrong: num_scales must be in [1,1000].";
	}

	// sets the values of the best fit
	this->best_accum_val = MIN_DOUBLE;
	this->best_scale = MIN_DOUBLE;
	this->best_angle = MIN_DOUBLE;
	this->best_ref_pt = cv::Point(MIN_INT, MIN_INT);

	// prepares gaussian kernel for smothing of the accumulator
	int gauss_size = 9.0*max_img_size/100;
	double gauss_sigma = 1.0*max_img_size/100;
	gauss_size = (gauss_size % 2) ? gauss_size : gauss_size + 1;
	gauss = cv::getGaussianKernel(gauss_size, gauss_sigma, CV_32F);
	gauss = gauss * gauss.t();

	// initialize the rotation step
	this->rot_rad = 0;
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
	hough_points = HoughTable(this->num_quant_directions);

	// changes gradient orientations to gradient directions
	cv::Mat orient, orient_adjust, directions;
	orient = sedlamat::get_gradient_orientation(img);
	orient_adjust = (orient >= 180)/255;
	orient_adjust.convertTo(orient_adjust, orient.depth());
	directions = orient + orient_adjust * -180;
	//my::display(directions);
	float quant_width = 180.0 / this->num_quant_directions;
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
	this->set_hough_points(this->r_table, this->tmpl_img,
						this->tmpl_edges);

	// shifts all points in HoughTable from the reference point
	for(auto & pts : this->r_table) {
	    for(auto & pt : pts) {
		pt -= this->ref_pt;
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
	this->set_hough_points(this->src_hough_pts, this->src_img,
						this->src_edges);
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
	this->rotated_r_table = this->r_table;

	// shifts the R-Table quants, last-becomes-first shifting
	for (int shift = 0; shift < num_table_shift; ++shift) {
	    std::vector<cv::Point> table_quant;
	    table_quant = this->rotated_r_table.back();
	    this->rotated_r_table.pop_back();
	    HoughTable::iterator it = this->rotated_r_table.begin();
	    this->rotated_r_table.insert(it, table_quant);
	}

	// rotates the R-Table points by angle_rad (counter-clockwise)
	double cs = std::cos(rot_step_rad);
	double sn = std::sin(rot_step_rad);

	for(auto & table_quant : this->rotated_r_table) {
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
	for (int scale_idx = 0; scale_idx < this->num_scales;
						    ++scale_idx) {
	    cv::Mat accum = cv::Mat(this->src_size, CV_32F,
					    cv::Scalar_<float>(0.0));
	    double s = this->min_scale + scale_idx * this->scale_step;
	    // for all quantified directions
	    for(int quant_idx = 0;
		quant_idx < this->num_quant_directions;
		++quant_idx) {
		cv::Mat quant_accum = cv::Mat(this->src_size, CV_32F,
					    cv::Scalar_<float>(0.0));
		std::vector<cv::Point> r_table_pts, src_pts;

		// pooling with neighbouring quants
		for (int quant_shift = -this->num_quant_neighbours;
		     quant_shift <= this->num_quant_neighbours;
		     ++quant_shift) {
		    int idx = quant_idx + quant_shift;
		    if (idx < 0) {
			idx = idx + this->num_quant_directions;
		    }
		    if (idx >= this->num_quant_directions) {
			idx = idx - this->num_quant_directions;
		    }
		    r_table_pts.insert(r_table_pts.end(),
				    this->rotated_r_table[idx].begin(),
				    this->rotated_r_table[idx].end());
		    src_pts.insert(src_pts.end(),
				   this->src_hough_pts[idx].begin(),
				   this->src_hough_pts[idx].end());
		}

		float num_quant_pts = r_table_pts.size();
		// for all points in the r_table quant
		for(auto & pt_diff : r_table_pts) {
		    // for all points in the src_hough_points quant
		    for(auto & src_pt : src_pts) {
			// increments accumulator
			cv::Point refer_pt(src_pt - (pt_diff*s));
			// only if inside source image area
			if (this->src_img_rect.contains(refer_pt)) {
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
	    //if (display_accum) sedlamat::display(accum);

	    // smooths the accumulator
	    cv::filter2D(accum, accum, CV_32F, this->gauss);

	    if (display_accum) sedlamat::display(accum);

	    // finds accumulator maximum
	    double local_max = 0, local_min = 0;
	    cv::Point local_max_pt(0,0), local_min_pt(0,0);
	    cv::minMaxLoc(accum, &local_min, &local_max,
					&local_min_pt, &local_max_pt);

	    if (display_accum) sedlamat::print(local_max);

	    // if accum_max greatest so far, then sets the best_ params
	    if (local_max > this->best_accum_val) {
		this->best_accum_val = local_max;
		this->best_scale = s;
		this->best_angle = this->rot_rad;
		this->best_ref_pt.x = local_max_pt.x;
		this->best_ref_pt.y = local_max_pt.y;
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
	const int num_of_rot = this->num_quant_directions * 2;
	const double rot_step_rad = 2.0 * M_PI / num_of_rot;
	for (int rot_idx = 0; rot_idx < num_of_rot; ++rot_idx) {
	    this->rot_rad = rot_idx * rot_step_rad;
	    this->set_rotated_r_table(rot_rad, rot_idx);
	    this->accumulate();
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
	this->src_img.copyTo(dst);

	// rotates and shift the R-Table points
	double cs = std::cos(this->best_angle);
	double sn = std::sin(this->best_angle);

	for(auto table_quant : this->r_table) {
	    for(auto pt : table_quant) {
		int x = pt.x;
		pt.x = cs*pt.x - sn*pt.y;
		pt.y = sn*x + cs*pt.y;
		pt = pt * this->best_scale;
		pt = pt + this->best_ref_pt;
		if (this->src_img_rect.contains(pt)) {
		    dst.at<cv::Vec3b>(pt) = cv::Vec3b(0,0,255);
		}
	    }
	}
	cv::resize(dst, dst, cv::Size(0,0), 4.0, 4.0);
	return dst;
    }

} /* namespace sedlamat */

#endif /* _MY_GENERAL_HOUGH_HPP_ */
