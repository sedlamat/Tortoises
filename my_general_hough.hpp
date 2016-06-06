/**
    my_general_hough.hpp

    Implementation of the generalized Hough transform, using OpenCV
    library.

    Do not forget that images are displayed with y-axis
    pointing down, it might fool you if you try to visualize
    individual steps of the algorithms below.

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

/* THIRD PARTY LIBRARIES */
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"

namespace sedlamat
{
    /*
	Beware:
	    cv::Point is alias for cv::Point_<int>
	    cv::Size is alias for cv::Size_<int>
	    cv::Rect is alias for cv::Rect_<int>
    */
    class GeneralHough {

	// alias used to store points grouped by gradient direction
	typedef std::vector<std::vector<cv::Point> > HoughTable;

	cv::Mat src_img, src_edges, template_img, template_edges;
	HoughTable r_table, src_hough_points;

	cv::Size src_size, template_img_size;
	cv::Rect src_area, template_img_area;

	cv::Point ref_pt; // reference point of the template

	const int num_quant_directions; // 4 * 3 [* 3] ...
	const int num_scales;
	const int max_img_size;
	const double max_template_scale;
	const double min_template_scale;

	int template_size;
	double max_scale; //= 1.1*MAX_IMG_SIZE/MAX_TEMPLATE_SIZE;
	double min_scale; //= 0.3 * max_template_scale;

	// parameters of the best fit
	double best_scale;
	double best_angle;
	cv::Point best_ref_pt;

    public:
	GeneralHough(const cv::Mat &src_img,
		     const cv::Mat &template_img,
		     const cv::Point ref_pt,
		     const cv::Mat &src_edges = cv::Mat(),
		     const cv::Mat &template_edges = cv::Mat(),
		     const int num_quant_directions = 4,
		     const int num_scales = 50,
		     const int max_img_size = 150,
		     const double max_template_scale = 1.0,
		     const double min_template_scale = 0.0); // 0.3

	~GeneralHough() {}

	void run();
	void display() const;
	cv::Point get_ref_pt() const { return ref_pt; };

    private:
	void fill_r_table();
	HoughTable get_hough_points();
    };

    GeneralHough::GeneralHough(const cv::Mat &src_img,
			    const cv::Mat &template_img,
			    const cv::Point ref_pt,
			    const cv::Mat &src_edges,
			    const cv::Mat &template_edges,
			    const int num_quant_directions,
			    const int num_scales,
			    const int max_img_size,
			    const double max_template_scale,
			    const double min_template_scale):
			    src_img(src_img),
			    src_edges(src_edges),
			    template_img(template_img),
			    template_edges(template_edges),
			    ref_pt(ref_pt),
			    num_quant_directions(num_quant_directions),
			    num_scales(num_scales),
			    max_img_size(max_img_size),
			    max_template_scale(max_template_scale),
			    min_template_scale(min_template_scale)
    {
	if (this->src_img.empty() || this->template_img.empty()) {
	    throw "Source and/or template image is empty!";
	} else {
	    this->template_img_size = this->template_img.size();
	    this->template_img_area = cv::Rect(cv::Point(0,0),
					    this->template_img_size);
	    if ( ! this->ref_pt.inside(this->template_img_area)) {
		throw "Reference point outside template image!";
	    } else {
		int src_w = this->src_img.cols;
		int src_h = this->src_img.rows;
		double resize_koef = this->max_img_size * 1.0 /
						std::max(src_w,src_h);
		cv::resize(this->src_img, this->src_img, cv::Size(0,0),
					resize_koef, resize_koef);

		this->src_size = this->src_img.size();
		this->src_area = cv::Rect(cv::Point(0,0),
						    this->src_size);
	    }
	}





    }



    // alias used to store points grouped by gradient direction
    typedef std::vector<std::vector<cv::Point_<int> > > HoughTable;

    // global constants for general hough transform
    const int NUM_OF_QUANT_DIRECTIONS = 4; // 4 * 3 [* 3] ...
    const int NUM_OF_SCALES = 50;
    const int MAX_IMG_SIZE = 150;
    const int MAX_TEMPLATE_SIZE = 150;
    const float MAX_SCALE = 1.1*MAX_IMG_SIZE/MAX_TEMPLATE_SIZE;
    const float SCALES_LOW_BOUND = 0.3 * MAX_SCALE;

    /**
	Gets HoughTable containing coordinates of edges grouped by
	quantized gradient directions of the edges.

	Gradient orientations in [180, 360) are shifted to [0, 180) as
	only the directions of the edges are used, e.g., 45deg is
	the same direction as 225deg. This fixes the problem of a
	gradient with 45deg orientation becoming 225deg when only the
	intensity of the object/background is inversed.

	@param src - Color (BGR) image.
	@param src_edges -  Edge image of the src image.
	@return HoughTable generated from src and src_edges.
    */
    HoughTable get_hough_points(const cv::Mat &src,
				const cv::Mat &src_edges)
    {
	HoughTable hough_points(NUM_OF_QUANT_DIRECTIONS);

	// changes gradient orientations to gradient directions
	cv::Mat orient, orient_adjust, directions;
	orient = my::get_gradient_orientation(src);
	orient_adjust = (orient >= 180)/255;
	orient_adjust.convertTo(orient_adjust, orient.depth());
	directions = orient + orient_adjust * -180;
	//my::display(directions);
	float quant_width = 180.0 / NUM_OF_QUANT_DIRECTIONS;
	//prt(quant_width/2);
	// goes through all edge pixels
	for (int yy = 0; yy < src_edges.rows; yy++) {
	    const uchar *ptr_src_edges_irow = src_edges.ptr<uchar>(yy);
	    for (int xx = 0; xx < src_edges.cols; xx++) {
		if (ptr_src_edges_irow[xx]) {
		    cv::Point_<int> pt(xx, yy);
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
		    // checks if quant for a given direction was found
		    CV_Assert(quant_idx != -1);
		    // fills the HoughTable
		    hough_points[quant_idx].push_back(pt);
		}
	    }
	}
	return hough_points;
    }

    /**
	Gets Hough R-Table, i.e., HoughTable containing coordinates of
	edges of a predefined template shifted by a given reference
	point on the template and grouped by quantized gradient
	directions of the edges.

	The R-Table holds diffences of edges points from the reference
	point, grouped by directions of the edge points.

	@param src_templ - Template image for the hough tranform.
	@param src_edges - Edge image of the template image.
	@param ref_point - Reference point of the template.
	@return HoughTable of a template shifted by a reference point.
    */
    HoughTable get_r_table(const cv::Mat &src_templ,
			   const cv::Mat &src_templ_edges,
			   const cv::Point_<int> &ref_point)
    {
	HoughTable r_table;
	r_table = get_hough_points(src_templ, src_templ_edges);

	// shifts all points in HoughTable from the reference point
	for(auto & points : r_table) {
	    for(auto & point : points) {
		point = point - ref_point;
	    }
	}
	return r_table;
    }

    /**
	Gets rotated R-Table. BEWARE: rotating counter-clockwise,
	parameter angle_rad is negative, the shifting the R-Table
	quants is then easier to implement.

	@param r_table - Hough R-Table.
	@param angle_rad - Angle in radians by which points are
			   rotated.
	@param num_table_shift - Number of times the R-Table quants
				 are shifted.
	@return Rotated R-Table.
    */

    HoughTable get_rotated_r_table(const HoughTable &r_table,
				   double angle_rad,
				   int num_table_shift)
    {
	HoughTable rotated_r_table = r_table;

	// shifts the R-Table quants, last-becomes-first shifting
	for (int shift = 0; shift < num_table_shift; ++shift) {
	    //my::visualize_points(r_table[shift], cv::Size_<int>(500,500), cv::Point_<int>(250,250));
	    std::vector<cv::Point> table_quant;
	    table_quant = rotated_r_table.back();
	    rotated_r_table.pop_back();
	    HoughTable::iterator it = rotated_r_table.begin();
	    rotated_r_table.insert(it, table_quant);
	    //my::visualize_points(rotated_r_table[2], cv::Size_<int>(500,500), cv::Point_<int>(250,250));
	}

	// rotates the R-Table points by angle_rad (counter-clockwise)
	double cs = std::cos(angle_rad);
	double sn = std::sin(angle_rad);

	for(auto & table_quant : rotated_r_table) {
	    for(auto & pt : table_quant) {
		int x = pt.x;
		pt.x = cs*pt.x - sn*pt.y;
		pt.y = sn*x + cs*pt.y;
	    }
	}
	//my::visualize_points(rotated_r_table[1], cv::Size_<int>(500,500), cv::Point_<int>(250,250));
	return rotated_r_table;
    }

    /**
	Gets the hough accumulator source feature points and reference
	points from the hough r-table.

	@param r_table - Reference point shifts grouped by their gradient
			 orientations.
	@param src_hough_points - Feature points grouped by their gradient
			      orientations.
	@param size - Size of the accumulator/source_image.
	@return void.
    */
    void get_accumulator(HoughTable &r_table,
			 HoughTable &src_hough_points,
			 cv::Size_<int> &size,
			 double &accum_max,
			 double &scale_max,
			 cv::Point_<int> &ref_point_found)
    {
	//std::clock_t t_start = std::clock();

	accum_max = 0;
	float multiplier = 1.0;
	int quant_neighbour = NUM_OF_QUANT_DIRECTIONS/4;
	switch (quant_neighbour) {
	    case 1:
		quant_neighbour = 0;
		break;
	    case 3:
		quant_neighbour = 1;
		break;
	    case 9:
		quant_neighbour = 4;
		break;
	    default:
		throw "Unexpected quant_neighbour";
	}
	//quant_neighbour = quant_neighbour == 1 ? 0 //NUM_OF_QUANT_DIRECTIONS/12;

	cv::Size_<int> sizeM(size.width*multiplier, size.height*multiplier);
	cv::Rect_<int> src_rect(cv::Point_<int>(0,0),size);
	// prepare 2 kernels - plain and gauss
	int plain_size = 2*multiplier;
	int gauss_size = 11*multiplier;
	cv::Mat plain(plain_size,plain_size, CV_32F, cv::Scalar_<float>(1.0));
	cv::Mat gauss = cv::getGaussianKernel(gauss_size,1, CV_32F);
	gauss = gauss * gauss.t();
	//~ cv::Mat minus(gauss_size, gauss_size, CV_32F,
			//~ cv::Scalar_<float>(-my::maxMat(gauss)/1000.0));
	//~ gauss += minus;


	float scale = (MAX_SCALE-SCALES_LOW_BOUND)/(NUM_OF_SCALES-1);

	for (int scale_idx = 0; scale_idx < NUM_OF_SCALES;
							++scale_idx) {
	    cv::Mat accum;
	    accum = cv::Mat(sizeM, CV_32F, cv::Scalar_<float>(0.0));
	    float s = SCALES_LOW_BOUND + scale_idx * scale;
	    for(int quant_idx = 0; quant_idx < NUM_OF_QUANT_DIRECTIONS;
							++quant_idx) {
		cv::Mat quant_accum =
			cv::Mat(sizeM, CV_32F, cv::Scalar_<float>(0.0));
		std::vector<cv::Point_<int> > r_table_pts, src_pts;
		//prt(r_table[quant_idx].size());
		//my::visualize_points(r_table[quant_idx], cv::Size_<int>(500,500), cv::Point_<int>(250,250));
		for (int quant_shift = -quant_neighbour;
		     quant_shift <= quant_neighbour; ++quant_shift) {
		    int idx = quant_idx + quant_shift;
		    if (idx < 0) {
			idx = idx + NUM_OF_QUANT_DIRECTIONS;
		    }
		    if (idx >= NUM_OF_QUANT_DIRECTIONS) {
			idx = idx - NUM_OF_QUANT_DIRECTIONS;
		    }
		    r_table_pts.insert(r_table_pts.end(),
				       r_table[idx].begin(),
				       r_table[idx].end());
		    src_pts.insert(src_pts.end(),
				   src_hough_points[idx].begin(),
				   src_hough_points[idx].end());
		}
		float num_quant_pts = r_table_pts.size();
		//prt(num_quant_pts);
		for(auto & pt_diff : r_table_pts) {
		    for(auto & src_pt : src_pts) {
			cv::Point_<int> ref_pt(src_pt - (pt_diff*s));
			if (src_rect.contains(ref_pt)) {
			    float &accum_pt =
			      quant_accum.at<float>(ref_pt*multiplier);
			    if (accum_pt < num_quant_pts) {
				accum_pt += 1.0;
			    }
			}
		    }
		}
		if (num_quant_pts > 0) {
		    accum += quant_accum/num_quant_pts;
		}

		//prt(accum);
	    }
	    cv::filter2D(accum, accum, CV_32F, plain);
	    cv::filter2D(accum, accum, CV_32F, gauss);

	    //my::display(accum);

	    double local_max = 0, local_min = 0;
	    cv::Point_<int> local_max_pt(0,0), local_min_pt(0,0);
	    cv::minMaxLoc(accum, &local_min, &local_max,
					&local_min_pt, &local_max_pt);
	    //prt(local_max);
	    if (local_max > accum_max) {
		accum_max = local_max;
		scale_max = s;
		ref_point_found.x = local_max_pt.x/multiplier;
		ref_point_found.y = local_max_pt.y/multiplier;
	    }
	}
	//std::clock_t t_middle = std::clock();
	//prt((t_middle-t_start)*1.0/CLOCKS_PER_SEC);
    }

    /**
	Computes the generalized Hough transform. The info about the
	best match is saved into non-const parameters of the function,
	see the declaration.

	@param templ - Template image for the hough tranform.
	@param templ_edges - Edges for the template image.
	@param ref_point - Reference point of the template.
	@param src - Source image.
	@param src_edges - Edges of the source image.
	@return void.
    */
    void general_hough(const cv::Mat &templ,
		       const cv::Mat &templ_edges,
		       const cv::Point_<int> &ref_point,
		       const cv::Mat &src,
		       const cv::Mat &src_edges,
		       double &accu_max,
		       double &rot_max,
		       double &scale_max,
		       cv::Point_<int> &ref_point_found)
    {
	HoughTable r_table, src_hough_points;
	// sets accu_max -1 to check at the end if it changed
	accu_max = -1;

	r_table = get_r_table(templ, templ_edges, ref_point);
	src_hough_points = get_hough_points(src, src_edges);

	cv::Size_<int> size = src.size();

	// rotates the R-Table
	const int num_of_rot = NUM_OF_QUANT_DIRECTIONS * 2;
	const double rot_step_rad = 2.0 * M_PI / num_of_rot;
	for (int rot_idx = 0; rot_idx < num_of_rot; ++rot_idx) {
	    // BEWARE: rotating counter-clockwise (see the minus sign)
	    double angle_rad = rot_idx * rot_step_rad;
	    HoughTable rotated_r_table;
	    rotated_r_table = get_rotated_r_table(r_table, angle_rad,
						  rot_idx);
	    /* gets accumulator maximum at one of the scales for a
	       given rotated R-Table and saves it into rot_accu_max
	       and rot_scale_max */
	    double rot_accu_max = -1;
	    double rot_scale_max = -1;
	    cv::Point_<int> rot_max_ref_point;
	    get_accumulator(rotated_r_table,
				src_hough_points,
				size,
				rot_accu_max,
				rot_scale_max,
				rot_max_ref_point);
	    // keeps track of the global accumulator
	    if (rot_accu_max > accu_max) {
		accu_max = rot_accu_max;
		rot_max = angle_rad;
		scale_max = rot_scale_max;
		ref_point_found = rot_max_ref_point;
	    }
	}
	CV_Assert(accu_max != -1);
    }

    /**
	Computes the generalized Hough transform. The best match is
	displayed on the source image with the template edge image.

	@param templ - Template image for the Hough tranform.
	@param templ_edges - Edges for the template image.
	@param ref_point - Reference point of the template.
	@param src - Source image.
	@param src_edges - Edges of the source image.
	@return Best fit.
    */
    cv::Mat general_hough_fit_on_img(const cv::Mat &templ,
			          const cv::Mat &templ_edges,
				  const cv::Point_<int> &ref_point,
				  const cv::Mat &src,
				  const cv::Mat &src_edges)
    {
	cv::Mat dst;
	src.copyTo(dst);
	//prt("dst type"); my::prt(dst.type());
	double accu_max = 0, rot_max = 0, scale_max = 0;
	cv::Point_<int> ref_point_found(0,0);

	general_hough(templ, templ_edges, ref_point, src, src_edges,
		      accu_max, rot_max, scale_max, ref_point_found);

	//prt(accu_max); prt(rot_max); prt(scale_max); prt(ref_point_found);

	HoughTable r_table;
	r_table = get_r_table(templ, templ_edges, ref_point);

	cv::Rect_<int> src_rect(cv::Point_<int>(0,0), src.size());

	// rotates and shift the R-Table points
	double cs = std::cos(rot_max);
	double sn = std::sin(rot_max);

	for(auto & table_quant : r_table) {
	    for(auto & pt : table_quant) {
		int x = pt.x;
		pt.x = cs*pt.x - sn*pt.y;
		pt.y = sn*x + cs*pt.y;
		pt = pt * scale_max;
		pt = pt + ref_point_found;
		if (src_rect.contains(pt)) {
		    dst.at<cv::Vec3b>(pt) = cv::Vec3b(0,0,255);
		}
	    }
	}
	cv::resize(dst, dst, cv::Size_<int>(0,0), 4.0, 4.0);
	return dst;
    }

} /* namespace my */

#endif /* _MY_GENERAL_HOUGH_HPP_ */
