/**
    general_hough.hpp

    Implementation of the generalized Hough transform, using OpenCV
    library.

    Usage:
    *

    @author Matej Sedlacek
    @version 0.0
*/

#ifndef _GENERAL_HOUGH_HPP_
#define _GENERAL_HOUGH_HPP_

/* STANDARD C++ LIBRARIES */
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <cmath>
#include <ctime>
#include <limits>

/* THIRD PARTY LIBRARIES */
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

/**
    Class GeneralHough
*/
class GeneralHough {

    /**
	Class DequeNegIdx - For easier implementation of the Hough
	 		    R-table.
			  - Double ended queue with negative indexing.
			  - Inheriting deque from standard library.
    */
    /* deque with negative indexing and rotation of elements */
    template <typename T> struct DequeNegIdx: public std::deque<T> {
	DequeNegIdx(): std::deque<T>() {}
	DequeNegIdx(size_t size): std::deque<T>(size) {}
	virtual ~DequeNegIdx() {};

	T & operator[](int idx);
	void operator++();
	void operator--();
	void shift(int num_shifts);

    };

    /**
	Class HoughTable - DequeNegIdx (deque) of vectors of points
    */
    struct HoughTable: public DequeNegIdx<std::vector<cv::Point>> {
	HoughTable(): DequeNegIdx<std::vector<cv::Point>>() {}
	HoughTable(size_t size):
			    DequeNegIdx<std::vector<cv::Point>>(size){}
	virtual ~HoughTable() {};

	void rotate_points(const int angle);
    };

    /** GeneralHough data members **/

    // source image
    const cv::Mat _src_img;

    // Hough table for the _src_img
    HoughTable _src_hough_table;

    // Hough R-table
    HoughTable _r_table;

    // template reference point
    const cv::Point _ref_pt;

    // angles and scales of the template to be tried
    std::vector<int> _angles;
    std::vector<float> _scales;

    // kernel for smoothing the accumulator
    cv::Mat _gauss;

    // resize coefficient from original source image
    double _resize_coeff;

    // boundary rectangle of the resized source image
    cv::Rect _resized_src_img_rect;

    // if 1 then all accumulators will be displayed and only one thread
    const bool _display_accum;

    // interest point to be located after GeneralHough detection
    std::map<std::string, cv::Point> _tmpl_interest_pts;

    // object for safe mutual execution of threads over parts of code
    std::mutex _mutualexec;

    // variables for the best fit of the hough transform
    double _best_accum_val;
    double _best_scale;
    int _best_angle;
    cv::Point _best_ref_pt;

public: /** GeneralHough public member functions */

    GeneralHough(
	    const cv::Mat &src_img,
	    const cv::Mat &tmpl_img,
	    const cv::Point ref_pt,
	    const std::vector<int> &angles,
	    const int num_scales = 20,
	    const int max_img_size = 150,
	    const double max_scale = 1.0,
	    const double min_scale = 0.3,
	    const int canny_low_thresh = 50,
	    const int canny_high_thresh = 100,
	    const double blur_sigma = 0.5,
	    const bool display_accum = 0,
	    const std::map<std::string, cv::Point> tmpl_interest_pts
				= std::map<std::string, cv::Point>());
    virtual ~GeneralHough() {}
    void detect();
    cv::Mat get_result_img() const;
    double get_best_accum_val() const { return _best_accum_val; }
    double get_best_angle() const { return _best_angle; }
    double get_best_scale() const { return _best_scale; }
    cv::Point get_best_ref_pt() const { return _best_ref_pt; }
    cv::Mat get_src_img() const { return _src_img; }

private: /** GeneralHough private member functions */

    HoughTable get_hough_table(const cv::Mat &img,
			       const cv::Mat &img_edges) const;
    HoughTable get_rotated_r_table(const int angle) const;
    void accumulate(const int angle);
    HoughTable get_quanted_table(HoughTable &table,
				 const int angle,
				 const int num_quants = 4);
    cv::Mat get_gradient_orientation(const cv::Mat& src) const;
};


/********* class DequeNegIdx members' declaration *******************/

/**
    Allowing negative indexing. E.g., index -1 equals index
    size-1 of the deque.
*/
template <typename T>
T &GeneralHough::DequeNegIdx<T>::operator[](int idx)
{
    size_t size = this->size();
    idx = idx % size;
    if (idx < 0) idx += size;
    return this->at(idx);
}

/**
    Rotating the elements of the deque to the right.
*/
template <typename T>
void GeneralHough::DequeNegIdx<T>::operator++()
{
    this->push_front(this->back());
    this->pop_back();
}
/**
    Rotating the elements of the deque to the left.
*/
template <typename T>
void GeneralHough::DequeNegIdx<T>::operator--()
{
    this->push_back(this->front());
    this->pop_front();
}

/**
    Rotating the elements of the deque by num_shifts, to the
    left or to the right depending on the sign of num_shifts.
*/
template <typename T>
void GeneralHough::DequeNegIdx<T>::shift(int num_shifts)
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

/********* class HoughTable members' declaration ********************/

/**
    Rotating the table points by angle in degrees.
*/
void GeneralHough::HoughTable::rotate_points(const int angle)
{
    double cs = std::cos(angle * M_PI / 180.0);
    double sn = std::sin(angle * M_PI / 180.0);

    const auto it_start = this->begin(), it_end = this->end();
    for (auto it = it_start; it != it_end; ++it) {
	for(auto &pt : *it) {
	    int x = pt.x;
	    pt.x = cs*pt.x - sn*pt.y;
	    pt.y = sn*x + cs*pt.y;
	}
    }
}

/********* class GeneralHough members' declaration ******************/

/**
    GeneralHough constructor. Sets various parameters, see
    class declaration.

    @param see class declaration.
    @return void.
*/
GeneralHough::GeneralHough(
	const cv::Mat &src_img,
	const cv::Mat &tmpl_img,
	const cv::Point ref_pt,
	const std::vector<int> &angles,
	const int num_scales,
	const int max_img_size,
	double max_scale,
	double min_scale,
	const int canny_low_thresh,
	const int canny_high_thresh,
	const double blur_sigma,
	const bool display_accum,
	const std::map<std::string, cv::Point> tmpl_interest_pts)
	:
	_src_img(src_img),
	_ref_pt(ref_pt),
	_angles(angles),
	_display_accum(display_accum),
	_tmpl_interest_pts(tmpl_interest_pts),
	_best_accum_val(-1.0),
	_best_scale(-1.0),
	_best_angle(-1),
	_best_ref_pt(cv::Point(-1,-1))
{
    // checks max_img_size
    if ( ! (1 <= max_img_size && max_img_size <= 1000) ) {
	throw "Wrong: max_img_size must be in [1,1000]";
    }

    // checks tmpl_img
    if (tmpl_img.empty()) {
	throw "Template image is empty!";
    }

    //checks _ref_pt
    if ( ! _ref_pt.inside(cv::Rect(cv::Point(0,0), tmpl_img.size()))) {
	    throw "Wrong: Reference point outside template image.";
    }

     // checks _src_img
    if (_src_img.empty()) {
	throw "Source image is empty!";
    }

    // resizes _src_img into src_img_resized
    cv::Mat src_img_resized;

    int src_w = _src_img.cols;
    int src_h = _src_img.rows;

    _resize_coeff = max_img_size * 1.0 / std::max(src_w,src_h);
    cv::resize(_src_img, src_img_resized, cv::Size(0,0),
				    _resize_coeff, _resize_coeff);

    // sets _resized_src_img_rect
    _resized_src_img_rect = cv::Rect(cv::Point(0,0),
					src_img_resized.size());

    // blurs the src_img_resized
    cv::GaussianBlur(src_img_resized, src_img_resized, cv::Size(0,0),
							blur_sigma);

    // sets source image edges
    cv::Mat src_edges;

    if (src_img_resized.channels() == 3) {
	cv::cvtColor(src_img_resized, src_edges, cv::COLOR_BGR2GRAY);
    } else if (src_img_resized.channels() == 1) {
	src_img_resized.copyTo(src_edges);
    } else {
	throw "Wrong: Source image has to have 1 or 3 channels.";
    }

    cv::Canny(src_edges, src_edges, canny_low_thresh,
						    canny_high_thresh);

    // fills source Hough table
    _src_hough_table = this->get_hough_table(src_img_resized,
							    src_edges);

    // sets template edges
    cv::Mat tmpl_edges;

    if (tmpl_img.channels() == 3) {
	cv::cvtColor(tmpl_img, tmpl_edges, cv::COLOR_BGR2GRAY);
    } else if (tmpl_img.channels() == 1) {
	tmpl_img.copyTo(tmpl_edges);
    } else {
	throw "Wrong: Template image has to have 1 or 3 channels.";
    }
    cv::Canny(tmpl_edges, tmpl_edges, canny_low_thresh,
						    canny_high_thresh);

    // fills Hough R-table
    _r_table = this->get_hough_table(tmpl_img, tmpl_edges);

    // shifts all points in the R-table by the reference point
    for(auto &pts : _r_table) {
	for(auto &pt : pts) {
	    pt -= _ref_pt;
	}
    }

    // checks min/max_scale
    if (min_scale < max_scale && 0 < min_scale && max_scale <= 1) {
    	throw "Wrong: min_scale must be < then max_scale and both \
							    in (0,1].";
    }

    // sets proper min/max_scale based on size of tmpl_img/src_img_res.
    cv::Mat tmpl_pts;
    cv::findNonZero(tmpl_img, tmpl_pts);
    cv::Size tmpl_size = cv::boundingRect(tmpl_pts).size();
    const double src_templ_size_koef = max_img_size*1.0 /
			std::max(tmpl_size.width,tmpl_size.height);

    max_scale *= src_templ_size_koef;
    min_scale *= src_templ_size_koef;

    // checks num_scales
    if ( ! (1 <= num_scales && num_scales <= 1000) ) {
	throw "Wrong: num_scales must be in [1,1000].";
    }

    // sets scale_step
    float scale_step;

    if (num_scales == 1) {
	scale_step = 0;
	min_scale = max_scale;
    } else {
	scale_step = (max_scale - min_scale) / (num_scales - 1);
    }

    // fills _scales
    for (int ii = 0; ii < num_scales; ++ii) {
	_scales.push_back(min_scale + ii * scale_step);
    }

    // fills gaussian kernel
    int gauss_size = 2.0*max_img_size/100;
    double gauss_sigma = 0.8*max_img_size/100;
    gauss_size = (gauss_size % 2) ? gauss_size : gauss_size + 1;
    _gauss = cv::getGaussianKernel(gauss_size, gauss_sigma, CV_32F);
    _gauss = _gauss * _gauss.t();

}

/**
    Fills in the HoughTable container with Hough points (edge
    points with their direction).

    @param hough_points - HoughTable container to be filled.
    @param img - image which will provide edge points and
		 directions.
    @param img_edges - edge image of img.
    @return filled HoughTable container.
*/
GeneralHough::HoughTable GeneralHough::get_hough_table(
					const cv::Mat &img,
					const cv::Mat &img_edges) const
{
    HoughTable hough_table(180); // 0...179

    // changes gradient orientations to gradient directions
    cv::Mat orient, orient_adjust, directions;
    orient = this->get_gradient_orientation(img);
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
    return hough_table;
}

cv::Mat GeneralHough::get_gradient_orientation(const cv::Mat& src) const
{
    cv::Mat gray, grad_orient;
    if (src.channels() == 3) {
	cv::cvtColor(src, gray, CV_BGR2GRAY);
    } else if (src.channels() == 1) {
	gray =  src;
    } else {
	throw "Incompatible number of channels";
    }
    cv::Mat dx, dy;
    cv::Mat kernel = (cv::Mat_<int>(3,3) << 3, 0, -3,
					   10, 0,-10,
					    3, 0, -3);
    cv::filter2D(gray, dx, CV_32F, kernel);
    cv::transpose(kernel, kernel);
    cv::filter2D(gray, dy, CV_32F, kernel);
    cv::phase(dx, dy, grad_orient, 1);
    return grad_orient;
}


/**
    Sets rotated_r_table by shifting the r_table by num_table_shift
    and rotating the points by rot_step_rad.

    @param angle - Angle in degrees by which points are rotated.
    @return void.
*/
GeneralHough::HoughTable GeneralHough::get_rotated_r_table(
						const int angle) const
{
    HoughTable rotated_r_table = _r_table;

    rotated_r_table.shift(angle);

    rotated_r_table.rotate_points(angle);

    return rotated_r_table;
}


/**
    Computes the hough accumulator for a given angle (rotation of the
    template/R-table. Best results are stored in best_... variables.

    @param angle - rotation angle of the R-table in degrees.
    @return void.
*/
void GeneralHough::accumulate(const int angle)
{
    int down_size = 2;
    cv::Size src_img_size(_resized_src_img_rect.size());
    cv::Size accum_size(src_img_size.width/down_size + 1,
			src_img_size.height/down_size + 1 );
    HoughTable rotated_r_table(get_rotated_r_table(angle));
    HoughTable quanted_rotated_r_table, quanted_src_hough_table;
    quanted_rotated_r_table = get_quanted_table(rotated_r_table,
						angle);
    quanted_src_hough_table = get_quanted_table(_src_hough_table,
						angle);
    int num_quant = quanted_src_hough_table.size();
    double total_num_r_table_pts = 0;
    for (auto &s : _scales) {
	cv::Mat accum = cv::Mat(accum_size, CV_32F,
					cv::Scalar_<float>(0.0));
	for (int quant = 0; quant < num_quant; ++quant) {
	    cv::Mat quant_accum = cv::Mat(accum_size, CV_32F,
					  cv::Scalar_<float>(0.0));
	    double num_quant_pts
			    = quanted_rotated_r_table[quant].size();
	    total_num_r_table_pts += num_quant_pts;
	    for(const auto &pt_diff : quanted_rotated_r_table[quant]) {
		for(const auto &src_pt : quanted_src_hough_table[quant]) {
		    cv::Point refer_pt(src_pt - (pt_diff*s));
		    if (_resized_src_img_rect.contains(refer_pt)) {
			float &accum_pt =
				quant_accum.at<float>(cv::Point(refer_pt.x/down_size, refer_pt.y/down_size));
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
		accum += quant_accum*1.0/(num_quant_pts);
	    }
	}

	// smooths the accumulator
	cv::filter2D(accum, accum, CV_32F, _gauss);

	if (_display_accum) sedlamat::display(accum);

	double local_max, local_min;
	cv::Point local_max_pt, local_min_pt;
	cv::minMaxLoc(accum, &local_min, &local_max,
				    &local_min_pt, &local_max_pt);

	cv::Mat dst;
	double min_val, max_val;
	cv::minMaxLoc(accum, &min_val, &max_val);
	double alpha = 1.0 / (max_val - min_val);
	double beta = - min_val * alpha;
	accum.convertTo(dst, CV_32F, alpha, beta);

	cv::Point lpt, rpt;
	lpt = cv::Point(_gauss.size().width/2, _gauss.size().height/2);
	rpt = cv::Point(accum.size().width, accum.size().height);
	rpt -= lpt;
	cv::Rect accum_insides(lpt, rpt);
	cv::Rect gauss_rect(cv::Point(0,0), cv::Size(lpt.x,lpt.y));

	//sedlamat::print(lpt);
	//sedlamat::print(rpt);
	//cv::Mat accum_mask(accum.size(), CV_32F, cv::Scalar_<float>(1));
	//cv::rectangle(accum_mask, gauss_rect,
					//~ cv::Scalar_<float>(0),
					//~ CV_FILLED);
    //    sedlamat::display(accum_mask);
	// finds accumulator maximum
	cv::matchTemplate(dst, _gauss, dst, CV_TM_SQDIFF);
	//sedlamat::display(dst);
	//double minimum = dst.at<float>(local_max_pt)/max_val;
	//std::cout << minimum << " " << dst.at<float>(local_max_pt) << " " << max_val << std::endl;
	cv::Mat accum_mask(dst.size(), CV_32F, cv::Scalar_<float>(1));
	gauss_rect += local_max_pt;
	gauss_rect -= lpt;
	gauss_rect -= cv::Point(lpt.x/2, lpt.y/2);
	cv::rectangle(accum_mask, gauss_rect,
					cv::Scalar_<float>(0),
					CV_FILLED);
	//sedlamat::display(accum);
	cv::Mat accum_mask_inverse;
	accum_mask.copyTo(accum_mask_inverse);
	accum_mask_inverse -= 1;
	//sedlamat::display(accum_mask);
	accum_mask_inverse = cv::Mat(cv::abs(accum_mask_inverse));
	//sedlamat::display(accum_mask);
	//sedlamat::display(accum_mask);
	//sedlamat::display(accum_mask_inverse);
	double sum_outer = cv::sum(accum_mask.mul(dst))[0];
	//double sum_inner = cv::sum(accum_mask_inverse.mul(dst))[0];
	double minimum;
	if (local_max)
	minimum = sum_outer/1000.0/local_max;
	if(_display_accum) sedlamat::print(minimum);
	//~ sedlamat::print(sum_inner);
	    std::lock_guard<std::mutex> guarding(_mutualexec);

	    //if (_display_accum) sedlamat::print(local_max);
	{
	    // if accum_max greatest so far, then sets the best_ params
	    if ( minimum < _best_accum_val) {

		_best_accum_val = minimum;
		_best_scale = s;
		_best_angle = angle;
		_best_ref_pt.x = local_max_pt.x*down_size;
		_best_ref_pt.y = local_max_pt.y*down_size;
	    }
	}
    }
}

GeneralHough::HoughTable GeneralHough::get_quanted_table(
					    HoughTable &table,
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
    After execution, _best_accum_val/angle/_ref_pt/scale will
    hold the information of the estimated location of the
    target object(s).

    @param void.
    @return void.
*/
void GeneralHough::detect()
{
    // accumulate for all rotated r-tables
    if (!_display_accum) {
	std::vector<std::thread> threads;
	for (auto angle : _angles) { //all angles in degrees integers!!
	    threads.push_back(std::thread(&GeneralHough::accumulate,this,angle));
	    //accumulate(angle);
	}
	for (auto &t : threads) {
	    t.join();
	}
    } else {
	for (auto angle : _angles) {
	    accumulate(angle);
	}
    }
    //~ const int num_of_rot = num_quant_directions * 2;
    //~ const double rot_step_rad = 2.0 * M_PI / num_of_rot;
    //~ for (int rot_idx = 0; rot_idx < num_of_rot; ++rot_idx) {
	//~ rot_rad = rot_idx * rot_step_rad;
    //~ }

    double inv_resize_coeff = 1.0 / _resize_coeff;
    _best_scale *= inv_resize_coeff;
    _best_ref_pt *= inv_resize_coeff;

    if (!_tmpl_interest_pts.empty()) {
	double cs = std::cos(_best_angle*M_PI/180.0);
	double sn = std::sin(_best_angle*M_PI/180.0);

	std::map<std::string, cv::Point>::iterator it, it_end;
	it = _tmpl_interest_pts.begin();
	it_end = _tmpl_interest_pts.end();
	while (it != it_end) {
	    it->second -= _ref_pt;
	    int x = it->second.x;
	    it->second.x = cs*it->second.x - sn*it->second.y;
	    it->second.y = sn*x + cs*it->second.y;
	    it->second *= _best_scale;
	    it->second += _best_ref_pt;
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
    _src_img.copyTo(dst);
    cv::Rect src_img_orig_rect(cv::Point(0,0), _src_img.size());

    // rotates and shift the R-Table points
    double cs = std::cos(_best_angle*M_PI/180.0);
    double sn = std::sin(_best_angle*M_PI/180.0);

    double inv_resize_coeff = 1.0 / _resize_coeff;
    cv::Point shift_pt(inv_resize_coeff,inv_resize_coeff);

    for(auto table_quant : _r_table) {
	for(auto pt : table_quant) {
	    int x = pt.x;
	    pt.x = cs*pt.x - sn*pt.y;
	    pt.y = sn*x + cs*pt.y;
	    pt = pt * _best_scale;
	    pt = pt + _best_ref_pt;
	    if (src_img_orig_rect.contains(pt)) {
		cv::rectangle(dst,
			      pt-shift_pt,
			      pt+shift_pt,
			      cv::Scalar(0,0,255), CV_FILLED);
	    }
	}
    }
    if (!_tmpl_interest_pts.empty()) {
	std::map<std::string, cv::Point>::const_iterator it, it_end;
	it = _tmpl_interest_pts.cbegin();
	it_end = _tmpl_interest_pts.cend();
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

#endif /* _GENERAL_HOUGH_HPP_ */
