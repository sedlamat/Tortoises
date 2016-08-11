/**
    GeneralHough.hpp

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

    // Detection takes place in the constructor of the object

    // User interface - recoved detected data using public getters
    //                |
    //                V
    cv::Mat result_img = general_hough.get_result_img();
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

    typedef std::deque<std::vector<cv::Point>> TableType;
    /**
	Class HoughTable - double-ended queue of vectors of points
    */
    class HoughTable {
	TableType _table;
    public:
	HoughTable() : _table() {}
	HoughTable(std::size_t size) : _table(size) {}
	~HoughTable() {};

	HoughTable &operator=(const HoughTable &) = default;
	HoughTable(const HoughTable &) = default;

	std::vector<cv::Point> &operator[](int idx);
	void operator++();
	void operator--();
	void shift(int num_shifts);
	void rotate_points(const int angle);
	TableType::iterator begin() { return _table.begin(); }
	TableType::iterator end() { return _table.end(); }
	TableType::const_iterator begin() const {
					    return _table.cbegin(); }
	TableType::const_iterator end() const { return _table.cend(); }
	std::size_t get_size() const { return _table.size(); }
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
    ~GeneralHough() {} // virtual only if the class is to be used
		       // as a base class and polymorfic (with other
		       // virtual functions

    /** user interface */
    cv::Mat get_result_img() const;
    double get_best_angle() const { return _best_angle; }
    double get_best_scale() const { return _best_scale; }
    cv::Point get_best_ref_pt() const { return _best_ref_pt; }

private: /** GeneralHough private member functions */
    GeneralHough(const GeneralHough &) = delete;
    GeneralHough &operator=(const GeneralHough &) = delete;

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


#endif /* _GENERAL_HOUGH_HPP_ */
