/**
    tortoise.cpp

    Implementation of the Tortoise class member functions.

    @author Matej Sedlacek
    @version 0.0
*/

#include "tortoise.hpp"

/**
    Constructor, where the feature detection is executed.

    @param plastron_image - Reference to a cv::Mat plastron image.
    @param tortoise_name - String with the toroise name/tag.
    @return Void.
*/
Tortoise::Tortoise(const cv::Mat &plastron_image,
		   const std::string &tortoise_name)
	    :
	    _plastron_img(plastron_image),
	    _tortoise_name(tortoise_name),
	    _info_file_name(tortoise_name.substr(0,7) + "_info.txt"),
	    _angle(0),
	    _l_juncs(std::vector<cv::Point>(7)),
	    _r_juncs(std::vector<cv::Point>(7)),
	    _center(),
	    _left_side(),
	    _right_side(),
	    _plastron_found(0),
	    _junctions_found(0),
	    _color(),
	    _idx()
{
    this->read_info();

    if (!_plastron_found) {
	this->locate_plastron();
	_plastron_found = 1;
	this->write_info();
    }

    if (_plastron_found && !_junctions_found) {
	this->locate_junctions();
	//_junctions_found = 1; enable only if localization is complete
	//this->write_info();
    }
}


/**
    Overloading operator >> for reading cv::Point from a stream.

    @param stream - Reference to a stream to read from.
    @param pt - Reference to a point to write.
    @return Returns the stream.
*/
std::istream &operator>>(std::istream &stream, cv::Point &pt)
{
    // cv::Point is printed in format [pt.x, pt.y]
    int a = 0;
    // moves through the stream until an integer is ahead
    while (!std::isdigit(a = stream.peek()) && a != '-') stream.get();
    // loads pt.x
    stream >> pt.x;
    // again moves through the stream until an integer is ahead
    while (!std::isdigit(a = stream.peek()) && a != '-') stream.get();
    // loads pt.y
    stream >> pt.y;
    // gets beyond the cv::Point expression
    while (stream.get() != ' ') ;
    return stream;
}


/**
    Read already computed information from a file if it exists.

    @param Void.
    @return Void.
*/
void Tortoise::read_info()
{
    std::ifstream info(_info_file_name);
    if (info.good()) {
	info >> _plastron_found;
	info >> _junctions_found;
	info >> _angle;
	info >> _center;
	info >> _left_side;
	info >> _right_side;
	for (auto &junc : _l_juncs) {
	    info >> junc;
	}
	for (auto &junc : _r_juncs) {
	    info >> junc;
	}
	std::cout << _plastron_found << '\n';
	std::cout << _junctions_found << '\n';
	std::cout << _angle << '\n';
	std::cout << _center << '\n';
	std::cout << _left_side << '\n';
	std::cout << _right_side << '\n';

	this->print_juncs();

	if (!info.good()) { // error when reading -> recompute data
	    _plastron_found = 0;
	    _junctions_found = 0;
	}
    } else {
	std::cout << "File Tg#####_info.txt does not exist. It will";
	std::cout << " be created." << std::endl;
    }

    info.close();
}


/**
    Write already computed information to a file.

    @param Void.
    @return Void.
*/
void Tortoise::write_info() const
{
    std::ofstream info(_info_file_name, std::fstream::out);
    if (info.good()) {
	info << _plastron_found << " ";
	info << _junctions_found << " ";
	info << _angle << " ";
	info << _center << " ";
	info << _left_side << " ";
	info << _right_side << " ";
	for (const auto &junc : _l_juncs) {
	    info << junc  << " ";
	}
	for (const auto &junc : _r_juncs) {
	    info << junc << " ";
	}

	info.close();
    } else {
	info.close();
	std::cout << "Error when opening/reading Tg#####_info.txt";
	std::cout << " file."<< std::endl;
	std::exit(1);
    }
}


/**
    Print values of junctions to standard output.

    @param Void.
    @return Void.
*/
void Tortoise::print_juncs() const
{
    std::cout << "Left junctions:" << std::endl;
    for (const auto &junc : _l_juncs) {
	std::cout << junc << std::endl;
    }

    std::cout << "Right junctions:" << std::endl;
    for (const auto &junc : _r_juncs) {
	std::cout << junc << std::endl;
    }
}

void Tortoise::put_point_on_img(cv::Mat &img,
				const cv::Rect &img_rect,
				const cv::Point &pt,
				const cv::Vec3b &color,
				const cv::Point &shift)
{
    if (img_rect.contains(pt)) {
	cv::rectangle(img,
		      pt-shift,
		      pt+shift,
		      cv::Scalar(color), CV_FILLED);
    }
}


/**
    Displays junction or feature points on the corresponding plastron
    image.

    @param Void.
    @return Void.
*/
void Tortoise::display_points_on_plastron()
{
    cv::Mat img;
    _plastron_img.copyTo(img);
    cv::Rect img_rect(cv::Point(0,0), img.size());

    const double inv_resize_coeff = std::max(img.cols, img.rows) / 200;
    cv::Point shift(inv_resize_coeff,inv_resize_coeff);

    for (const auto &junc : _l_juncs) {
	this->put_point_on_img(img, img_rect, junc, _color.red, shift);
    }
    for (const auto &junc : _r_juncs) {
	this->put_point_on_img(img, img_rect, junc, _color.green,
								shift);
    }

    this->put_point_on_img(img, img_rect, _center, _color.blue,
								shift);
    this->put_point_on_img(img, img_rect, _left_side, _color.purple,
								shift);
    this->put_point_on_img(img, img_rect, _right_side, _color.cyan,
								shift);
    this->put_point_on_img(img, img_rect, _r_juncs[0], _color.red,
								shift);

    sedlamat::display(img);
}

/**
    Rotates point (cv::Point) by an angle around a given center point.

    @param pt - The point to be rotated.
    @param center - The center of rotation.
    @param angle - The angle of rotation.
    @return Void.
*/
void Tortoise::rotate_pt(cv::Point &pt, const cv::Point &center,
			 const int angle)
{
    const double cs = std::cos(angle * M_PI / 180.0);
    const double sn = std::sin(angle * M_PI / 180.0);

    pt -= center;

    const int x = pt.x;
    const int y = pt.y;
    pt.x = cs*x - sn*y;
    pt.y = sn*x + cs*y;

    pt += center;
}


/**
    Rotates the plaston image and all interest points by an angle
    around the center of the plastron.

    @param angle - The angle of rotation.
    @return Void.
*/
void Tortoise::rotate_img_and_pts(const int angle)
{

    // makes _plastron_img square matrix with wide borders
    const int w = _plastron_img.size().width;
    const int h = _plastron_img.size().height;
    const int max_dist = std::max(std::max(_center.x, _center.y),
			          std::max(w-_center.x, h-_center.y));
    const int top = max_dist - _center.y;
    const int bottom = max_dist - (h-_center.y);
    const int left = max_dist - _center.x;
    const int right = max_dist - (w-_center.x);
    cv::copyMakeBorder(_plastron_img, _plastron_img, top, bottom, left,
		       right, cv::BORDER_CONSTANT);

    // rotates points and _plastron_img around _center
    const cv::Point shift = cv::Point(max_dist, max_dist) - _center;

    for (auto &junc : _l_juncs) {
	this->rotate_pt(junc, _center, -angle);
	junc += shift;
    }
    for (auto &junc : _r_juncs) {
	this->rotate_pt(junc, _center, -angle);
	junc += shift;
    }
    this->rotate_pt(_left_side, _center, -angle);
    _left_side += shift;
    this->rotate_pt(_right_side, _center, -angle);
    _right_side += shift;

    _center += shift;

    cv::Mat M = cv::getRotationMatrix2D(_center, angle, 1);
    cv::warpAffine(_plastron_img, _plastron_img, M,
						_plastron_img.size());
}


/**
    Locates the plastron position in the plastron image using
    generalized Hough transform.

    @param Void.
    @return Void.
*/
void Tortoise::locate_plastron()
{
    const cv::Mat plastron_template
			    = cv::imread("plastron_template.bmp", 0);

    if (plastron_template.empty()) {
	std::cout << "Error: plastron_template.bmp image is not" \
	" inside the directory of the executable file."
	<< std::endl;
	exit(1);
    }

    const std::vector<int> angles = {-5,0,5,-85,-90,-95,180,85,90,95};

    const cv::Point reference_point(65,125);

    std::vector<cv::Point> interest_pts;
    interest_pts.push_back(cv::Point(64,58)); 	// junction 0 - Head
    interest_pts.push_back(cv::Point(64,74));
    interest_pts.push_back(cv::Point(64,97));
    interest_pts.push_back(cv::Point(64,105));
    interest_pts.push_back(cv::Point(64,155));
    interest_pts.push_back(cv::Point(64,168));
    interest_pts.push_back(cv::Point(64,185)); 	// junction 6 - Tail
    interest_pts.push_back(cv::Point(7,131));  	// left side border
    interest_pts.push_back(cv::Point(123,131)); // right side border
    interest_pts.push_back(cv::Point(65,125));  // plastron center

    GeneralHough general_hough(_plastron_img, plastron_template,
			       reference_point,	angles, 20, 200, 1.0,
			       0.3, 50, 100, 0.5, 0,
			       &interest_pts);
    // assigns to interest points
    for (int ii = 0; ii < 7; ++ii) {
	_l_juncs[ii] = _r_juncs[ii] = interest_pts[ii];
    }
    _left_side = interest_pts[7];
    _right_side = interest_pts[8];
    _center = interest_pts[9];

    _angle = general_hough.get_best_angle();
}


/**


    @param Void.
    @return Void.
*/
void Tortoise::locate_junctions()
{
    std::cout << "locate_junctions entered" << std::endl;
    this->rotate_img_and_pts(_angle);
    //~ this->locate_central_seam();
    //~ this->display_points_on_plastron();
    this->locate_seams();
}


/**


    @param Void.
    @return Void.
*/
void Tortoise::locate_seams()
{
    // measures the to-be half width of the central seam stripe
    const int abd_seam_length = _l_juncs[_idx.j4AbdFem].y -
				_l_juncs[_idx.j3PecAbd].y;

    // prepares the area rectangle of the central seam stripe
    cv::Rect stripe_rect(
		cv::Point(_l_juncs[_idx.j0Head].x - abd_seam_length/2,
		          _l_juncs[_idx.j0Head].y ),
		cv::Point(_l_juncs[_idx.j6Tail].x + abd_seam_length/2,
		          _l_juncs[_idx.j6Tail].y ));

    // gets the area of the central seam stripe
    cv::Mat stripe_img = _plastron_img(stripe_rect).clone();

    // resizes the area
    int stripe_w = stripe_img.cols;
    int stripe_h = stripe_img.rows;

    const int new_size = 1200;

    const double resize_coeff = new_size * 1.0 / std::max(stripe_w,
							  stripe_h);
    cv::resize(stripe_img, stripe_img, cv::Size(0,0), resize_coeff,
						      resize_coeff);

    cv::cvtColor(stripe_img, stripe_img, CV_BGR2GRAY);

    stripe_w = stripe_img.cols;
    stripe_h = stripe_img.rows;

    //~ sedlamat::display(stripe_img);

    cv::Mat gauss_x1, gauss_x2, gauss_y, gauss, out;
    int size_x = 50;
    int size_y = 50;
    const double sigma_x1 = 1.0;
    const double sigma_x2 = 1.5;
    const double sigma_y = 10.0;
    size_x = (size_x % 2) ? size_x : size_x + 1;
    size_y = (size_y % 2) ? size_y : size_y + 1;
    gauss_x1 = cv::getGaussianKernel(size_x, sigma_x1, CV_32F);
    gauss_x2 = cv::getGaussianKernel(size_x, sigma_x2, CV_32F);
    gauss_y = cv::getGaussianKernel(size_y, sigma_y, CV_32F);
    gauss = gauss_y * gauss_x1.t();
    gauss -= gauss_y * gauss_x2.t();

    gauss(cv::Range(size_y/2 + 1, size_y), cv::Range::all()) = cv::Scalar(0);

    cv::Mat M = cv::getRotationMatrix2D(cv::Point(size_x/2, size_y/2), -90, 1);
    cv::warpAffine(gauss, gauss, M, gauss.size());
    //~ M = cv::getRotationMatrix2D(cv::Point(size_x/2, size_y/2), -283, 1);
    //~ cv::warpAffine(gauss, gauss, M, gauss.size());
    std::cout << gauss << std::endl;

    // We need 4 times cv::Mat - Two maxima, Two minima images

    sedlamat::display(gauss);

    cv::filter2D(stripe_img, out, CV_32F, gauss);
    std::cout << out.size() << std::endl;
    std::cout << stripe_img.size() << std::endl;


    sedlamat::display(out);
    sedlamat::display(stripe_img);
}



void Tortoise::locate_central_seam()
{
    // central seam localization

    // measures the to-be half width of the central seam stripe
    const int abd_seam_length = _l_juncs[_idx.j4AbdFem].y -
				_l_juncs[_idx.j3PecAbd].y;

    // prepares the area rectangle of the central seam stripe
    cv::Rect stripe_rect(
		cv::Point(_l_juncs[_idx.j0Head].x - abd_seam_length/2,
		          _l_juncs[_idx.j0Head].y ),
		cv::Point(_l_juncs[_idx.j6Tail].x + abd_seam_length/2,
		          _l_juncs[_idx.j6Tail].y ));

    // gets the area of the central seam stripe
    cv::Mat stripe_img = _plastron_img(stripe_rect).clone();

    // resizes the area
    int stripe_w = stripe_img.cols;
    int stripe_h = stripe_img.rows;

    const int new_size = 400;

    const double resize_coeff = new_size * 1.0 / std::max(stripe_w,
							  stripe_h);
    cv::resize(stripe_img, stripe_img, cv::Size(0,0), resize_coeff,
						      resize_coeff);

    stripe_w = stripe_img.cols;
    stripe_h = stripe_img.rows;

    cv::flip(stripe_img, stripe_img, 1);

    // gets the edge image of the area
    cv::Mat stripe_edges, stripe_gray;
    cv::cvtColor(stripe_img, stripe_gray, CV_BGR2GRAY);
    cv::Canny(stripe_gray, stripe_edges, 0, 0);
    stripe_edges.convertTo(stripe_edges, CV_32F);
    stripe_edges /= 255;

    cv::Mat stripe_edges_orig = stripe_edges.clone();

    sedlamat::display(stripe_edges);

    std::vector<cv::Point> central_seam;

    central_seam = get_best_path(stripe_edges);

    for (auto &pt : central_seam) {
	stripe_img.at<cv::Vec3b>(pt) = _color.red;
    }

    // locates side seams

    stripe_edges = stripe_edges_orig.t();

    stripe_edges.colRange(0, stripe_h/6).setTo(0);
    stripe_edges.colRange(stripe_h*5/6, stripe_h).setTo(0);

    //~ sedlamat::display(stripe_edges);

    for (int ii = 0; ii < 4; ++ii) {
	std::vector<cv::Point> side_seam;

	side_seam = get_best_path(stripe_edges);

	int best_path_x_coor = 0;
	for (auto &pt : side_seam) {
	    best_path_x_coor += pt.x;
	}
	if (side_seam.size()) {
	    best_path_x_coor /= side_seam.size();
	} else {
	    std::cout << "seam not found" << std::endl;
	    break;
	}

	std::cout << best_path_x_coor << std::endl;

	stripe_edges.colRange(best_path_x_coor-stripe_h*1/24,
			      best_path_x_coor+stripe_h*1/24).setTo(0);

	sedlamat::display(stripe_edges);

	stripe_img = stripe_img.t();

	for (auto &pt : side_seam) {
	    stripe_img.at<cv::Vec3b>(pt) = _color.red;
	}
	stripe_img = stripe_img.t();

	sedlamat::display(stripe_img);
    }

    // TODO: first approximate localization of the seven junctions
}


std::vector<cv::Point> Tortoise::get_best_path(const cv::Mat &edges)
{
    // edges in 0/1 float matrix

    int edges_w = edges.cols;
    int edges_h = edges.rows;

    const float supreme_value = 10000;

    cv::Mat values(edges.size(), CV_32F, cv::Scalar(supreme_value));

    // prepares distance mask 5rows x 7columns
    cv::Mat dist_mask(5, 7, CV_32F, cv::Scalar(0));

    // fills the dist_mask
    for (int x = 0; x < 7; ++x)  {
	dist_mask.at<float>(0,x) = 1.0;
	dist_mask.at<float>(1,x) = 3.0;
	dist_mask.at<float>(2,x) = 7.0;
	dist_mask.at<float>(3,x) = 12.0;
	dist_mask.at<float>(4,x) = 18.0;
    }

    for (int y = 0; y < 5; ++y) {
	dist_mask.at<float>(y,0) += 3;
	dist_mask.at<float>(y,6) += 3;
	dist_mask.at<float>(y,1) += 2;
	dist_mask.at<float>(y,5) += 2;
	dist_mask.at<float>(y,2) += 1;
	dist_mask.at<float>(y,4) += 1;
    }

    cv::Rect mask_rect(cv::Point(0,0), cv::Size(7,5));

    const int low_line_y = edges_h*1/8;
    const int high_line_y = edges_h*7/8;

    const cv::Rect low_line_rect(cv::Point(0,low_line_y),
				 cv::Point(edges_w-1,low_line_y+1));
    const cv::Rect high_line_rect(cv::Point(0,high_line_y),
				  cv::Point(edges_w-1,high_line_y+1));

    edges(low_line_rect) = 1.0;
    edges(high_line_rect) = 1.0;

    // locates the central seam at the low_line position

    values(high_line_rect) = 1.0;

    for (int y = high_line_y-1; y >= 6; --y) {
	for (int x = 3; x < edges_w - 3; ++x) {
	    if (edges.at<float>(y,x)) {
		double min;
		mask_rect.x = x-3;
		mask_rect.y = y+1;
		cv::minMaxLoc(edges(mask_rect).mul(dist_mask) +
			      values(mask_rect), &min);
		values.at<float>(y,x) = min;
	    }
	}
    }

    double low_min, low_max;
    cv::Point low_min_pt;
    cv::minMaxLoc(values(low_line_rect), &low_min, &low_max,
							&low_min_pt);

    std::vector<cv::Point> best_path;

    // locates central seam
    int y = low_line_y;
    int x = low_min_pt.x;

    while (y != high_line_y) {
	double loc_min, loc_max;
	cv::Point loc_min_pt;
	mask_rect.x = x-3;
	mask_rect.y = y+1;
	cv::minMaxLoc(dist_mask + values(mask_rect), &loc_min,
						&loc_max, &loc_min_pt);
	best_path.push_back(cv::Point(x,y));
	x += loc_min_pt.x - 3;
	y += loc_min_pt.y + 1;
    }

    return best_path;
}
