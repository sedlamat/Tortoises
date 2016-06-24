/**
    general_hough.hpp

    Implementation of the generalized Hough transform, using OpenCV
    library. Finds the best (one) instance of the searched object.

    Usage:

    // GeneralHough object declaration and initialization
    GeneralHough general_hough(source_img,	//non-optional param
			       template_img,	//non-optional param
			       reference_point,	//non-optional param
			       angles,		//non-optional param
			       20,	//number of scales
			       200,	//max source img size
			       1.0, 	//max template/source scale
			       0.3, 	//min template/source scale
			       50, 	//canny lower threshold
			       100, 	//canny higher threshold
			       0.5, 	//sigma of gaussinan blur
			       0,	//if accumulator shall be shown
			       ptr_interest_pts); //points to be found
    // -> the detection takes place in the constructor of the object

    //recoved detected data using public getters
    cv::Mat result_img = general_hough.get_result_img();
    double best_accum_value = general_hough.get_best_accum_val();
    double best_angle = general_hough.get_best_angle();
    double best_scale = get_best_scale();
    cv::Point best_reference_point = general_hough.get_best_ref_pt();

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
	DequeNegIdx(std::size_t size): std::deque<T>(size) {}
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
	HoughTable(std::size_t size):
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
    std::vector<cv::Point> *_ptr_interest_pts;

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
	    const int max_img_size = 200,
	    const double max_scale = 1.0,
	    const double min_scale = 0.3,
	    const int canny_low_thresh = 50,
	    const int canny_high_thresh = 100,
	    const double blur_sigma = 0.5,
	    const bool display_accum = 0,
	    std::vector<cv::Point> *ptr_interest_pts
				= nullptr);
    virtual ~GeneralHough() {}

    cv::Mat get_result_img() const;
    double get_best_accum_val() const { return _best_accum_val; }
    double get_best_angle() const { return _best_angle; }
    double get_best_scale() const { return _best_scale; }
    cv::Point get_best_ref_pt() const { return _best_ref_pt; }

private: /** GeneralHough private member functions */
    void detect();
    HoughTable get_hough_table(const cv::Mat &img,
			       const cv::Mat &img_edges) const;
    HoughTable get_rotated_r_table(const int angle) const;
    void accumulate(const int angle);
    HoughTable get_quanted_table(HoughTable &table,
				 const int angle,
				 const int num_quants = 4);
    cv::Mat get_gradient_direction(const cv::Mat& src) const;
    void display(const cv::Mat &src,
		    const std::string &window_name="image") const;
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
    GeneralHough constructor. Sets various parameters, see class
    declaration.

    @param See class declaration.
    @return Void.
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
	std::vector<cv::Point> *ptr_interest_pts)
	:
	_src_img(src_img),
	_src_hough_table(),
	_r_table(),
	_ref_pt(ref_pt),
	_angles(angles),
	_scales(),
	_gauss(),
	_resize_coeff(0.0),
	_resized_src_img_rect(),
	_display_accum(display_accum),
	_ptr_interest_pts(ptr_interest_pts),
	_mutualexec(),
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

    const int src_w = _src_img.cols;
    const int src_h = _src_img.rows;

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
    if ( !(min_scale < max_scale && 0 < min_scale && max_scale <= 1)) {
    	throw "Wrong: min_scale must be < max_scale " \
						  "and both in (0,1].";
    }

    // sets proper min/max_scale based on size of tmpl_img/src_img_res.
    cv::Mat tmpl_pts;
    cv::findNonZero(tmpl_img, tmpl_pts);
    const cv::Size tmpl_size = cv::boundingRect(tmpl_pts).size();
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
    const double gauss_sigma = 0.8*max_img_size/100;
    gauss_size = (gauss_size % 2) ? gauss_size : gauss_size + 1;
    _gauss = cv::getGaussianKernel(gauss_size, gauss_sigma, CV_32F);
    _gauss = _gauss * _gauss.t();

    this->detect();
}

/**
    Gets a filled HoughTable container with Hough points (edge
    points with their direction).

    @param img - Image which will provide edge points and directions.
    @param img_edges - Edge image of img.
    @return Filled HoughTable container.
*/
GeneralHough::HoughTable GeneralHough::get_hough_table(
					const cv::Mat &img,
					const cv::Mat &img_edges) const
{
    HoughTable hough_table(180); // 0...179

    // transforms gradient orientations to gradient directions
    cv::Mat grad_directions = this->get_gradient_direction(img);

    // goes through all edge pixels
    for (int yy = 0; yy < img_edges.rows; ++yy) {
	for (int xx = 0; xx < img_edges.cols; ++xx) {
	    if (img_edges.at<uchar>(yy, xx)) {
		cv::Point pt(xx, yy);
		int phi = grad_directions.at<uchar>(pt);
		hough_table[phi].push_back(pt);
	    }
	}
    }
    return hough_table;
}

/**
    Gets gradient direction image of the source image.

    @param src - Source image.
    @return Gradient direction image.
*/
cv::Mat GeneralHough::get_gradient_direction(const cv::Mat& src) const
{
    // gets gradient orientation
    cv::Mat gray, grad_orient;
    if (src.channels() == 3) {
	cv::cvtColor(src, gray, CV_BGR2GRAY);
    } else if (src.channels() == 1) {
	gray =  src;
    } else {
	throw "Incompatible number of channels";
    }
    cv::Mat dx, dy;
    const cv::Mat kernel = (cv::Mat_<int>(3,3) << 3, 0, -3,
						  10, 0,-10,
						  3, 0, -3);
    cv::filter2D(gray, dx, CV_32F, kernel);
    cv::transpose(kernel, kernel);
    cv::filter2D(gray, dy, CV_32F, kernel);
    cv::phase(dx, dy, grad_orient, 1);

    // transforms the gradient orientation to direction, orientation
    // angle above 180 is shifted by -180 to get direction
    cv::Mat orient_adjust, grad_dir;
    orient_adjust = (grad_orient >= 180)/255;
    orient_adjust.convertTo(orient_adjust, grad_orient.depth());
    grad_dir = grad_orient + orient_adjust * -180;
    grad_dir.convertTo(grad_dir, CV_8UC1);

    return grad_dir;
}


/**
    Gets rotated R-table by shifting and rotating the points of
    R-table by a given angle.

    @param angle - Angle in degrees by which points are rotated and
		   shifted.
    @return Rotated R-table (in HoughTable container).
*/
GeneralHough::HoughTable GeneralHough::get_rotated_r_table(
						const int angle) const
{
    HoughTable rotated_r_table(_r_table);

    rotated_r_table.shift(angle);
    rotated_r_table.rotate_points(angle);

    return rotated_r_table;
}


/**
    Computes the hough accumulator for a given angle (rotation of the
    template/R-table. Best results are stored in _best_... variables.
    Safe for multithreaded access.

    @param angle - Rotation angle of the R-table in degrees.
    @return Void.
*/
void GeneralHough::accumulate(const int angle)
{
    // scales down the accumulator size
    const float down_size = 0.5;
    const cv::Size src_img_size(_resized_src_img_rect.size());
    const cv::Size accum_size(src_img_size.width * down_size + 1,
			      src_img_size.height * down_size + 1);

    // prepares the rotated R-table
    HoughTable rotated_r_table(this->get_rotated_r_table(angle));

    // quantizes the HoughTable containers with respect to the angle
    HoughTable quanted_rotated_r_table, quanted_src_hough_table;
    quanted_rotated_r_table = this->get_quanted_table(rotated_r_table,
								angle);
    quanted_src_hough_table = this->get_quanted_table(_src_hough_table,
								angle);

    // fills the accumulator for all quants and scales
    const int num_quant = quanted_src_hough_table.size();

    for (auto &s : _scales) {
	cv::Mat accum = cv::Mat(accum_size, CV_32F,
					    cv::Scalar_<float>(0.0));
	for (int quant = 0; quant < num_quant; ++quant) {
	    cv::Mat quant_accum = cv::Mat(accum_size, CV_32F,
					  cv::Scalar_<float>(0.0));
	    const float num_quant_pts
			    = quanted_rotated_r_table[quant].size();
	    for(auto &pt_diff : quanted_rotated_r_table[quant]) {
		for(auto &src_pt : quanted_src_hough_table[quant]) {
		    cv::Point refer_pt(src_pt - (pt_diff*s));
		    if (_resized_src_img_rect.contains(refer_pt)) {
			float &accum_pt =
			    quant_accum.at<float>(refer_pt * down_size);

			// accum value should not become more than
			// the number of points in the r_table quant
			if (accum_pt < num_quant_pts) {
			    accum_pt += 1.0;
			}
		    }
		}
	    }
	    // divide by the total number of points in the quant
	    // big quants are not more important
	    if (num_quant_pts > 0) {
		accum += quant_accum / num_quant_pts;
	    }
	}

	// smooths the accumulator
	cv::filter2D(accum, accum, CV_32F, _gauss);

	if (_display_accum) this->display(accum);

	// finds local minimum and maximum
	double local_max, local_min;
	cv::Point local_max_pt, local_min_pt;
	cv::minMaxLoc(accum, &local_min, &local_max,
				    &local_min_pt, &local_max_pt);

	// gets scaled accumulator for analysis of its distribution
	cv::Mat dst;
	const double alpha = 1.0 / (local_max - local_min);
	const double beta = - local_min * alpha;
	accum.convertTo(dst, CV_32F, alpha, beta);
	cv::matchTemplate(dst, _gauss, dst, CV_TM_SQDIFF);

	// prepares mask that will cover the area of local maximum
	cv::Mat dst_mask(dst.size(), CV_32F, cv::Scalar_<float>(1));
	cv::Rect gauss_rect(cv::Point(0,0), _gauss.size()*2
						    + cv::Size(1,1));
	auto size_diff = accum.size() - dst.size();
	gauss_rect += local_max_pt;
	gauss_rect -= cv::Point(_gauss.size().width,
						_gauss.size().height);
	gauss_rect -= cv::Point(size_diff.width/2, size_diff.height/2);
	cv::rectangle(dst_mask, gauss_rect, cv::Scalar_<float>(0),
							    CV_FILLED);

	// sums pixel values around the area of local maximum
	const double sum_outer = cv::sum(dst_mask.mul(dst))[0];
	const double local_accum_max = local_max*100000.0/sum_outer;

	if(_display_accum) std::cout << local_accum_max << std::endl;

	// safe multithreading acces shared variables _best_...
	std::lock_guard<std::mutex> guarding(_mutualexec);

	// saves the best-so-far values
	if (local_accum_max > _best_accum_val) {
	    _best_accum_val = local_accum_max;
	    _best_scale = s;
	    _best_angle = angle;
	    _best_ref_pt.x = local_max_pt.x / down_size;
	    _best_ref_pt.y = local_max_pt.y / down_size;
	}
    }
}


/**
    Displays an unsigned 8bit (CV_8U) image, values [0...255], the
    number of channels allowed is 1 and 3. The image is scaled and
    casted to CV_8U if it is of other type.

    @param src - Image to be displayed.
    @param window_name - Displayed window name (optional).
    @return Void.
*/
void GeneralHough::display(const cv::Mat &src,
				const std::string &window_name) const
{
    CV_Assert(src.channels() == 1 || src.channels() == 3);
    cv::Mat dst;
    // checks if src is CV_8U, if not make it
    if (src.type() != CV_8UC1 && src.type() != CV_8UC3) {
	// scales src image into [0,...,255] and converts to CV_8U
	double min_val, max_val;
	cv::minMaxLoc(src, &min_val, &max_val);
	const double alpha = 255 / (max_val - min_val);
	const double beta = - min_val * alpha;
	src.convertTo(dst, CV_8U, alpha, beta);
    } else {
	dst = src;
    }
    // displays the image, exit with ESC key
    cv::namedWindow(window_name, cv::WINDOW_NORMAL);
    int k = 1;
    while (k != 27) {
	cv::imshow(window_name, dst);
	k = cv::waitKey(0);
    }
}


/**
    Gets HoughTable container with num_quants quanted from [0,...,179]
    original quants of the r/src_hough_table and the new quants are
    oriented based on a given angle.

    @param table - A HoughTable to be requanted.
    @param angle - Angle in degrees defining the offset of the new
								quants.
    @param num_quants - Number of new quants.
    @return Requanted HougTable container.
*/
GeneralHough::HoughTable GeneralHough::get_quanted_table(
					    HoughTable &table,
					    const int angle,
					    const int num_quants)
{
    HoughTable quanted_table(num_quants);
    int quant_idx = 0;
    const int quant_size = 180 / num_quants;
    int idx = angle < 0 ? angle % 180 + 180 : angle % 180;
    idx -= quant_size / 2;

    for (int ii = 1; ii <= 180; ++ii, ++idx) {
	if (!(ii % quant_size)) {
	    ++quant_idx;
	}
	quanted_table[quant_idx].insert(
				quanted_table[quant_idx].end(),
				table[idx].begin(),
				table[idx].end());
    }

    return quanted_table;
}


/**
    Finds and sets the _best_... params of the GeneralHough class
    by callling the accumulate member function with all the different
    angles. If template interest points were specified, then their
    position is computed here.

    @param Void.
    @return Void.
*/
void GeneralHough::detect()
{
    // if no displaying, the processing is done in parallel
    if (!_display_accum) {
	std::vector<std::thread> threads;
	for (auto angle : _angles) {
	    threads.push_back(std::thread(&GeneralHough::accumulate,
							this, angle));
	}
	for (auto &t : threads) {
	    t.join();
	}
    } else {
	for (auto angle : _angles) {
	    accumulate(angle);
	}
    }

    // rescales found _best_... params to the original size
    const double inv_resize_coeff = 1.0 / _resize_coeff;
    _best_scale *= inv_resize_coeff;
    _best_ref_pt *= inv_resize_coeff;

    // computes template interest points
    if (_ptr_interest_pts != nullptr) {
	const double cs = std::cos(_best_angle * M_PI / 180.0);
	const double sn = std::sin(_best_angle * M_PI / 180.0);

	auto it = _ptr_interest_pts->begin();
	const auto it_end = _ptr_interest_pts->end();
	while (it != it_end) {
	    *it -= _ref_pt;
	    int x = it->x;
	    it->x = cs * it->x - sn * it->y;
	    it->y = sn * x + cs * it->y;
	    *it *= _best_scale;
	    *it += _best_ref_pt;
	    ++it;
	}
    }
}


/**
    Gets the source image with the best template match marked on it.
    It also draws the template interest points if any were specified.

    @param Void.
    @return Image with the best template fit.
*/
cv::Mat GeneralHough::get_result_img() const
{
    //make a deep copy of the source image - it would be affected
    cv::Mat dst =  _src_img.clone();

    cv::Rect src_img_orig_rect(cv::Point(0,0), _src_img.size());

    // rotates and shift the R-Table points
    const double cs = std::cos(_best_angle * M_PI / 180.0);
    const double sn = std::sin(_best_angle * M_PI / 180.0);

    // points from template will be rectangles in the source image in
    // order to get resolution-invariant drawing result
    const double inv_resize_coeff = 1.0 / _resize_coeff;
    cv::Point shift_pt(inv_resize_coeff,inv_resize_coeff);

    // draws the template
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

    // draws the position of the template interest points
    if (_ptr_interest_pts != nullptr) {
	auto it = _ptr_interest_pts->cbegin();
	const auto it_end = _ptr_interest_pts->cend();
	while (it != it_end) {
	    if (src_img_orig_rect.contains(*it)) {
		cv::rectangle(dst,
			      *it - shift_pt,
			      *it + shift_pt,
			      cv::Scalar(0,255,0), CV_FILLED);
	    }
	    ++it;
	}
    }

    return dst;
}

#endif /* _GENERAL_HOUGH_HPP_ */
