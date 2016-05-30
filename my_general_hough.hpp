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
    const int NUM_OF_QUANT_DIRECTIONS = 4;
    const int NUM_OF_SCALES = 100;
    const float SCALES_LOW_BOUND = 0.3;

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
	orient_adjust.convertTo(orient_adjust,orient.depth());
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
	@param angle_rad - Angle in radians by which points are rotated.
	@param num_table_shift - Number of times the R-Table quants are
				 shifted.
	@return Rotated R-Table.
    */

    HoughTable get_rotated_r_table(const HoughTable &r_table,
				   double angle_rad,
				   int num_table_shift)
    {
	HoughTable rotated_r_table = r_table;

	// shifts the R-Table quants, last-becomes-first shifting
	for (int shift = 0; shift < num_table_shift; ++shift) {
	    std::vector<cv::Point> table_quant;
	    table_quant = rotated_r_table.back();
	    rotated_r_table.pop_back();
	    HoughTable::iterator it = rotated_r_table.begin();
	    rotated_r_table.insert(it, table_quant);
	}

	// rotating the R-Table points by angle_rad (counter-clockwise)
	double cs = std::cos(angle_rad);
	double sn = std::sin(angle_rad);

	for(auto & table_quant : rotated_r_table) {
	    for(auto & pt : table_quant) {
		int x = pt.x;
		pt.x = cs*pt.x - sn*pt.y;
		pt.y = sn*x + cs*pt.y;
	    }
	}
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

    cv::Mat get_accumulator(
	std::vector<std::vector<cv::Point_<int> > > r_table,
	std::vector<std::vector<cv::Point_<int> > > src_hough_points,
	cv::Size_<int> size,
	cv::Mat ref_pt_mask,
	double &accu_max,
	double &scale_max)
    {
	CV_Assert(ref_pt_mask.size() == size);
	const int sizes[3] = {size.height, size.width, NUM_OF_SCALES};
	cv::Mat accum[NUM_OF_QUANT_DIRECTIONS];
	for (int ii = 0; ii < NUM_OF_QUANT_DIRECTIONS; ++ii) {
	    accum[ii] = cv::Mat(3, sizes, CV_32F, cv::Scalar_<float>(0.0));
	}

	cv::Rect_<int> src_rect(cv::Point_<int>(0,0),size);
	float scale = 1.0/NUM_OF_SCALES;
	int ii = 0;
	for (float s = SCALES_LOW_BOUND; s < 1 + scale; s += scale, ii++) {

	    //std::cout << s << std::endl;
	    for(int quant_idx = 0; quant_idx < NUM_OF_QUANT_DIRECTIONS; ++quant_idx) {

		//std::cout << quant_idx << std::endl;
		//my::visualize_points(src_hough_points[quant_idx], size);
		//cv::Mat accum_quant(size, CV_32FC1, cv::Scalar_<float>(0.0));
		float quant_num = r_table[quant_idx].size();
		//std::cout << quant_num << std::endl;
		for(auto & pt_diff : r_table[quant_idx]) {
		    //std::vector<cv::Point_<int> > ref_pts;

		    for(auto & src_pt : src_hough_points[quant_idx]) {

			// OPTIONAL - ref_point of the template inside img
			cv::Point_<int> ref_pt(src_pt - pt_diff*s);
			//std::cout << pt_diff*s << std::endl;
			if (src_rect.contains(ref_pt) &&
			    ref_pt_mask.at<int>(ref_pt)) {
			    //std::cout << ii << std::endl;
			    if (accum[quant_idx].at<float>(ref_pt.y, ref_pt.x, ii) < quant_num)
				accum[quant_idx].at<float>(ref_pt.y, ref_pt.x, ii) += 1.0;
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

	accu_max = 0;

	for (int ii = 0; ii < NUM_OF_SCALES; ++ii) {
	    cv::Mat channel(size, CV_32FC1, cv::Scalar_<float>(1.0));
	    for (int xx = 0; xx < size.width; ++xx) {
		for (int yy = 0; yy < size.height; ++yy) {
		    for (int quant_ii = 0; quant_ii < NUM_OF_QUANT_DIRECTIONS; ++quant_ii) {
			channel.at<float>(yy, xx) += accum[quant_ii].at<float>(yy, xx, ii)*1.0/r_table[quant_ii].size();
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
	    double channel_max = my::maxMat(channel);
	    if (channel_max > accu_max) {
		accu_max = channel_max;
		scale_max = SCALES_LOW_BOUND + ii*scale;
	    }
	    //std::cout << my::maxMat(channel) << std::endl;
	    //my::display(channel);
	}
	return accum[0];
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

    const int num_of_rot = NUM_OF_QUANT_DIRECTIONS * 2;
    const double rot_step_rad = 2 * M_PI / num_of_rot;
    double accu_max = -1;
    double rot_max = -1;
    double scale_max = -1;
    for (int rot_idx = 0; rot_idx < num_of_rot; ++rot_idx) {
	// BEWARE: rotating counter-clockwise (see the minus sign below)
	double angle_rad = -rot_idx * rot_step_rad;
	//std::cout << angle_rad << std::endl;
	HoughTable rotated_r_table;
	rotated_r_table = get_rotated_r_table(r_table, angle_rad, rot_idx, size, ref_point);

	cv::Mat ref_pt_mask(size, CV_8UC1, cv::Scalar_<uchar>(255));

    //my::display(ref_pt_mask);
	double rot_accu_max = -1;
	double rot_scale_max = -1;
	cv::Mat accumulator = my::get_accumulator(rotated_r_table,
						src_hough_points,
						size, ref_pt_mask,
						rot_accu_max,
						rot_scale_max);
	if (rot_accu_max > accu_max) {
	    accu_max = rot_accu_max;
	    rot_max = angle_rad * 180 / M_PI;
	    scale_max = rot_scale_max;
	}
	prt(accu_max); prt(" where current "); prt(rot_accu_max);
	prt(rot_idx);
	prt(scale_max);
    }
    prt(accu_max);
    prt(rot_max);
    prt(scale_max);
    //my::display(accumulator);

    //cv::threshold(accumulator, ref_pt_mask, my::maxMat(accumulator)*0.8,
	//	  255, cv::THRESH_BINARY);
    //my::display(ref_pt_mask);

    //~ for(int ii = 0; ii < 10; ii++) {
	//~ accumulator = my::get_accumulator(r_table, src_hough_points,
				      //~ size, ref_pt_mask);
	//~ cv::threshold(accumulator, ref_pt_mask, my::maxMat(accumulator)*0.8,
		  //~ 255, cv::THRESH_BINARY);
	//~ my::display(accumulator);
    //~ }
 /*   cv::threshold(accumulator, ref_pt_mask, my::maxMat(accumulator)*0.8,
		  255, cv::THRESH_BINARY);
    my::display(ref_pt_mask);

    accumulator = my::get_accumulator(r_table, src_hough_points,
				      size, ref_pt_mask);
*/
    //~ my::display(accumulator);
//~
    //~ my::display(src);
    //~ my::display(src_edges);
    //~ my::display(src_templ);
    //~ my::display(src_templ_edges);
  }



} /* namespace my */


#endif /* _MY_GENERAL_HOUGH_HPP_ */
