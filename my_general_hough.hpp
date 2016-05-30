/**
  my_general_hough.hpp

  Implementation of the generalized Hough transform, using OpenCV
  library.

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
    // cv::Point is alias (typedef) for cv::Point_<int>
    typedef std::vector<cv::Point> HoughTableQuant;
    typedef std::vector<HoughTableQuant> HoughTable;

    const int NUM_OF_QUANT_DIRECTIONS = 4;
    const int NUM_OF_SCALES = 100;
    const float SCALES_LOW_BOUND = 0.3;

    /**
	Gets the hough points, i.e., a vector (indexed by quantized
	gradient orientations - 0 to 3) of vector of points (the
	coordinates of edges).


	@param src - A color image.
	@param src_edges -  An edge image of the color image.
	@return Vector of vector of points generated from src and
		src_edges.
    */
    HoughTable get_hough_points(
					       const cv::Mat& src,
					       const cv::Mat& src_edges)
    {
	cv::Mat orient = my::get_gradient_orientation(src);
	std::vector<std::vector<cv::Point_<int> > > hough_points(
					    NUM_OF_QUANT_DIRECTIONS);

	// orientations in [0, 360) => make [180, 360) to [0, 180), only
	// direction needed (45deg is the same direction as 225deg)
	// this prevents 45deg being 225deg with inverse intensity.
	// inverse intesity changes orientation, do not want that.
	// direction changes only if the edge changes.

	cv::Mat orient_adjust = (orient >= 180)/255;
	orient_adjust.convertTo(orient_adjust,orient.depth());
	orient += orient_adjust * -180;
	//my::display(orient);
	int max = 0;
	float quant_width = 180.0 / NUM_OF_QUANT_DIRECTIONS;
	// points are quantized from [0,180) into 0,1,2,3 indices
	for (int yy = 0; yy < src_edges.rows; yy++) {
	    const uchar *ptr_src_edges_irow = src_edges.ptr<uchar>(yy);
	    for (int xx = 0; xx < src_edges.cols; xx++) {
		if (ptr_src_edges_irow[xx]) {
		    cv::Point_<int> pt(xx, yy);
		    float phi = orient.at<float>(pt);
		    int quant_idx = -1;
		    if (phi >= 180 - quant_width/2 ||
			phi < quant_width/2)
			quant_idx = 0;
		    else {
			int idx = 1;
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
		    CV_Assert(quant_idx != -1&&"Undefined orientation");
		    if (quant_idx > max) max = quant_idx;
		    hough_points[quant_idx].push_back(pt);
		}
	    }
	}
	std::cout << max << std::endl;
	return hough_points;
    }

  /**
    Gets the hough r-table, i.e., the hough points minus reference
    point of the template.

    @param src_templ - A template image for the hough tranform.
    @param src_edges - An edge image of the template image.
    @param ref_point - Reference point of the template.
    @return The r-table (vector of vector of [points - ref_point]).
  */

  std::vector<std::vector<cv::Point_<int> > >  get_r_table(
				      const cv::Mat& src_templ,
				      const cv::Mat& src_templ_edges,
				      const cv::Point_<int>& ref_point)
  {
    std::vector<std::vector<cv::Point_<int> > >
	  r_table = get_hough_points(src_templ, src_templ_edges);

    for(auto & points : r_table) {
      //my::visualize_points(points, src_templ_edges.size());
      for(auto & point : points) {
	point = ref_point - point;
      }
    }
    return r_table;
  }

  /**
    Rotates r-table.

    @param p - Something.
    @param src_edges - An sdfaassssssssthe template image.
    @param ref_point - Referessssssssssss template.
    @return The sssssssssssssssssssss).
  */

    HoughTable get_rotated_r_table(const HoughTable &r_table, double angle_rad,
			       int num_table_shift, cv::Size_<int> size, cv::Point ref_point)
  {
    HoughTable rotated_r_table = r_table;
    // shift the table
    for (int shift = 0; shift < num_table_shift; ++shift) {
	HoughTableQuant table_quant = rotated_r_table.back();
	rotated_r_table.pop_back();
	HoughTable::iterator it = rotated_r_table.begin();
	rotated_r_table.insert(it, table_quant);
    }
    //my::visualize_points(r_table[0], cv::Size(1000,1000), -cv::Point(400,400));
    //my::visualize_points(rotated_r_table[1], cv::Size(1000,1000), -cv::Point(400,400));
    // rotate the points
    //angle_rad = 1.0*M_PI/2;

    double cs = std::cos(angle_rad);
    double sn = std::sin(angle_rad);
    //std::cout << cs << " " << sn << std::endl;
    //rotated_r_table[1][0] = cv::Point(10,10);
    for(auto & table_quant : rotated_r_table) {
	for(auto & pt : table_quant) {
	    int x = pt.x;
	    pt.x = (int) (cs*pt.x - sn*pt.y);
	    pt.y = (int) (sn*x + cs*pt.y);
	}
    }
    //std::cout << rotated_r_table[1][0] << std::endl;
    //my::visualize_points(r_table[0], cv::Size(1000,1000), -cv::Point(400,400));
    //my::visualize_points(rotated_r_table[1], cv::Size(1000,1000), -cv::Point(400,400));
    return rotated_r_table;
  }


    /**
	Gets a filled accumulator for one point.

	@param ref_points - Possible reference points for a given src_pt.
	@param src_pt - One feature point in the src_image.
	@param size - Size of the src_image/accumulator.
	@return One-layer accumulator for ref_points and src_pt.
    */

    cv::Mat get_accum_layer(std::vector<cv::Point_<int> > ref_points,
			  cv::Point_<int> src_pt,
			  cv::Size_<int> size)
    {
	cv::Point_<int> pt_shift(3,3);
	cv::Mat accum_layer(size, CV_32FC1, cv::Scalar_<float>(0.0));
	for(auto const& ref_pt : ref_points) {
	    //std::cout << ref_pt << " " << src_pt << std::endl;
	    cv::Point_<double> pt1, pt2;
	    double a, b;
	    if (src_pt.x == ref_pt.x) {
		pt1.x = pt2.x = src_pt.x;
		pt1.y = 0;
		pt2.y = size.height;
	    } else {
		a = (src_pt.y - ref_pt.y)*1.0/(src_pt.x - ref_pt.x);
		b = (ref_pt.y*src_pt.x - src_pt.y*ref_pt.x)*1.0/
							(src_pt.x - ref_pt.x);
		if (a == 0) {
		    pt1.y = pt2.y = src_pt.y;
		    pt1.x = 0;
		    pt2.x = size.width;
		} else if (std::abs(a) < 1) { // then use y = a*x + b
		    //continue;
		    pt1.x = 0;
		    pt1.y = a*pt1.x + b;
		    pt2.x = size.width;
		    pt2.y = a*pt2.x + b;
		} else { // then use x = (y-b)/a
		    //continue;
		    pt1.y = 0;
		    pt1.x = (pt1.y - b) / a;
		    pt2.y = size.height;
		    pt2.x = (pt2.y - b) / a;
		}
	    }
	    cv::line(accum_layer, pt1, pt2, cv::Scalar_<float>(1.0), 1);
	}
	//~ cv::rectangle(accum_layer, src_pt-pt_shift, src_pt+pt_shift,
		      //~ cv::Scalar_<float>(0.0), CV_FILLED);
	return accum_layer;
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
			cv::Point_<int> ref_pt(pt_diff*s + src_pt);
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
