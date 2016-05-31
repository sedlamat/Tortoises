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

namespace my
{
    // alias used to store points grouped by gradient direction
    typedef std::vector<std::vector<cv::Point_<int> > > HoughTable;

    // global constants for general hough transform
    const int NUM_OF_QUANT_DIRECTIONS = 16;
    const int NUM_OF_SCALES = 50;
    const int MAX_IMG_SIZE = 256;
    const int MAX_TEMPLATE_SIZE = 150;
    float MAX_SCALE = MAX_IMG_SIZE*1.0/MAX_TEMPLATE_SIZE;
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

	float quant_width = 180.0 / NUM_OF_QUANT_DIRECTIONS;
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
		    if (phi >= 180 - quant_width/2 ||
			phi < quant_width/2) {
			quant_idx = idx;
		    } else {
			++idx;
			for (float quant_low_bound = quant_width/2;
			     quant_low_bound < 180 - quant_width;
			     quant_low_bound += quant_width) {
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
	@return The r-table (vector of vector of [points - ref_point]).
    */
    cv::Mat get_accumulator(HoughTable &r_table,
			    HoughTable &src_hough_points,
			    cv::Size_<int> &size,
			    double &accum_max,
			    double &scale_max,
			    cv::Point_<int> &ref_point_found)
    {
	std::clock_t t_start = std::clock();

	const int sizes[3] = {size.height, size.width, NUM_OF_SCALES};
	cv::Mat accum[NUM_OF_QUANT_DIRECTIONS];
	for (int ii = 0; ii < NUM_OF_QUANT_DIRECTIONS; ++ii) {
	    accum[ii] = cv::Mat(3, sizes, CV_32F, cv::Scalar_<float>(0.0));
	}

	cv::Rect_<int> src_rect(cv::Point_<int>(0,0),size);

	float scale = MAX_SCALE/NUM_OF_SCALES;
	int ii = 0;
	for (float s = SCALES_LOW_BOUND; s < MAX_SCALE + scale; s += scale, ii++) {

	    //std::cout << s << std::endl;
	    for(int quant_idx = 0; quant_idx < NUM_OF_QUANT_DIRECTIONS; ++quant_idx) {

		//std::cout << quant_idx << std::endl;
		//my::visualize_points(src_hough_points[quant_idx], cv::Size_<int>(500,500), cv::Point_<int>(250,250));
		//my::visualize_points(r_table[quant_idx], cv::Size_<int>(500,500), cv::Point_<int>(250,250));
		//cv::Mat accum_quant(size, CV_32FC1, cv::Scalar_<float>(0.0));
		float quant_num = r_table[quant_idx].size();
		//std::cout << quant_num << std::endl;
		for(auto & pt_diff : r_table[quant_idx]) {
		    //std::vector<cv::Point_<int> > ref_pts;

		    for(auto & src_pt : src_hough_points[quant_idx]) {

			// OPTIONAL - ref_point of the template inside img
			cv::Point_<int> ref_pt(src_pt - (pt_diff*s) );
			//std::cout << pt_diff*s << std::endl;
			if (src_rect.contains(ref_pt)) {
			    //std::cout << ii << std::endl;
			    if (accum[quant_idx].at<float>(ref_pt.y, ref_pt.x, ii) < quant_num) {
				accum[quant_idx].at<float>(ref_pt.y, ref_pt.x, ii) += 1.0;
			    }
			    //std::cout << accum.at<float>(ref_pt.y, ref_pt.x) << std::endl;//+= 1.0;
			}
			//std::cout << pt_diff*s << std::endl;
			    //ref_pts.push_back(ref_pt);
		    }
		    //cv::Mat accum_single = my::get_accum_layer(ref_pts, src_pt,
							    //  size);
		    //std::cout << my::maxMat(accum_layer) << std::endl;
		    //my::display(accum_single);
		    //exit(1);
		    //my::display(accum_quant);
		    //cv::bitwise_or(accum_quant,accum_single,accum_quant);
		    //accum_quant += accum_single;
		    //my::display(accum_quant);
		}
		//cv::threshold(accum_quant, accum_quant, 1, 1, cv::THRESH_BINARY);
		//return accum;
		//my::display(accum_quant);
		//accum += accum_quant;

	    }
	}
	cv::Mat dst;

	accum_max = 0;

	for (int ii = 0; ii < NUM_OF_SCALES; ++ii) {
	    cv::Mat channel(size, CV_32FC1, cv::Scalar_<float>(1.0));
	    for (int xx = 0; xx < size.width; ++xx) {
		for (int yy = 0; yy < size.height; ++yy) {
		    for (int quant_ii = 0; quant_ii < NUM_OF_QUANT_DIRECTIONS; ++quant_ii) {
			if (r_table[quant_ii].size() > 0) {
			    channel.at<float>(yy, xx) += accum[quant_ii].at<float>(yy, xx, ii)*1.0/r_table[quant_ii].size();
			}
		    }
		}
	    }
	    cv::Mat kernel(2,2, CV_32F, cv::Scalar_<float>(1.0));
	    cv::Mat gauss = cv::getGaussianKernel(11,1, CV_32F);
	    gauss = gauss * gauss.t();
	    cv::Mat minus(11, 11, CV_32F, cv::Scalar_<float>(-my::maxMat(gauss)/100.0));

	    gauss += minus;
	    //prt(gauss);
	    //my::display(gauss);
	    //cv::rectangle(kernel,cv::Point_<int>(24,24),cv::Point_<int>(28,28), cv::Scalar_<float>(1.0), CV_FILLED);
	    cv::filter2D(channel, channel, CV_32F, kernel);
	    cv::filter2D(channel, channel, CV_32F, gauss);
	    double channel_max = 0, channel_min = 0;
	    cv::Point_<int> channel_pt_max(0,0), channel_pt_min(0,0);
	    cv::minMaxLoc(channel, &channel_min, &channel_max, &channel_pt_min, &channel_pt_max);
	    //double channel_max = my::maxMat(channel);
	    if (channel_max > accum_max) {
		accum_max = channel_max;
		scale_max = SCALES_LOW_BOUND + ii*scale;
		ref_point_found = channel_pt_max;
	    }
	    //std::cout << my::maxMat(channel) << std::endl;
	    //my::display(channel);
	}

	std::clock_t t_end = std::clock();
	prt((t_end-t_start)*1.0/CLOCKS_PER_SEC);
	return accum[0];
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

	r_table = my::get_r_table(templ, templ_edges, ref_point);
	src_hough_points = my::get_hough_points(src, src_edges);

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
	    cv::Mat accumulator;
	    accumulator = my::get_accumulator(rotated_r_table,
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
	@return void.
    */
    void general_hough_fit_on_img(const cv::Mat &templ,
			          const cv::Mat &templ_edges,
				  const cv::Point_<int> &ref_point,
				  const cv::Mat &src,
				  const cv::Mat &src_edges)
    {
	cv::Mat dst;
	src.copyTo(dst);
	prt("dst type"); my::prt(dst.type());
	double accu_max = 0, rot_max = 0, scale_max = 0;
	cv::Point_<int> ref_point_found(0,0);

	general_hough(templ, templ_edges, ref_point, src, src_edges,
		      accu_max, rot_max, scale_max, ref_point_found);

	prt(accu_max); prt(rot_max); prt(scale_max); prt(ref_point_found);

	HoughTable r_table;
	r_table = my::get_r_table(templ, templ_edges, ref_point);

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
	my::display(dst);
    }

} /* namespace my */

#endif /* _MY_GENERAL_HOUGH_HPP_ */
