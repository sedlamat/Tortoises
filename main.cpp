﻿/* STANDARD C++ LIBRARIES */
#include <iostream>
#include <string>
#include <limits>
#include <cmath>
#include <iomanip>
#include <ctime>

/* THIRD PARTY LIBRARIES */
#include <pwd.h>
#include <unistd.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

/* FIRST PARTY LIBRARIES */
#include "tortoise.hpp"

int main(int argc, char *argv[])
{
    try {
	time_t t0, t1;
	std::time(&t0);
	passwd *pw = getpwuid(getuid());
	std::string path(pw->pw_dir);
//~ //~
	std::string file_name = R"(Tg00301.jpg)";
	std::string imgs_path = path + R"(/Images/Tortoises/)";
	std::string img_path = imgs_path + file_name;
	//~ std::string img_path = argv[1];
	cv::Mat plastron_img = cv::imread(img_path,1);


	std::string templ_path = "plastron_template.bmp";
	cv::Mat plastron_template = cv::imread(templ_path, 0);

	cv::Point reference_point(65,125);

	Tortoise tg(plastron_img, plastron_template, reference_point);

	//std::vector<int> angles = {-5,0,5,-85,-90,-95,180,85,90,95};
	std::vector<int> angles = {-5,0,5};
	//std::vector<int> angles = {85,90,95};
	//std::vector<int> angles = {-85,-90,-95};
	//std::vector<int> angles = {0};


	//~ std::map<std::string, cv::Point> map_template_junctions;
	//~ map_template_junctions["1"] = cv::Point(64,58);
	//~ map_template_junctions["2"] = cv::Point(64,74);
	//~ map_template_junctions["3"] = cv::Point(64,97);
	//~ map_template_junctions["4"] = cv::Point(64,105);
	//~ map_template_junctions["5"] = cv::Point(64,155);
	//~ map_template_junctions["6"] = cv::Point(64,168);
	//~ map_template_junctions["7"] = cv::Point(64,185);
	//sedlamat::display(img);
    // GeneralHough object declaration and initialization
    GeneralHough general_hough(plastron_img,	//non-optional param
			       plastron_template,	//non-optional param
			       reference_point,	//non-optional param
			       angles,		//non-optional param
			       20,	//number of scales
			       200,	//max source img size
			       1.0, 	//max template/source scale
			       0.3, 	//min template/source scale
			       50, 	//canny lower threshold
			       100, 	//canny higher threshold
			       0.5, 	//sigma of gaussinan blur
			       0	//if accumulator shall be shown
			       ); //points to be found

    // runs the detection
    //general_hough.detect();

    //recoved detected data using public getters
    cv::Mat result_img = general_hough.get_result_img();
    double best_accum_value = general_hough.get_best_accum_val();
    double best_angle = general_hough.get_best_angle();
    double best_scale = general_hough.get_best_scale();
    cv::Point best_reference_point = general_hough.get_best_ref_pt();
	//~ GeneralHough general_hough(plastron_img,
					     //~ plastron_template,
					     //~ reference_point,
					     //~ angles,
					     //~ 20, 200, 1.0, 0.3, 50, 100, 0.5, 0);
					     //map_template_junctions);
	//sedlamat::display(img);
	//sedlamat::print("running");
	//sedlamat::display(general_hough.get_src_edges());
	//general_hough.detect();
	//~ time_t t1 = clock();
	//std::cout << general_hough.get_best_accum_val() << std::endl;

	//~ sedlamat::display(general_hough.get_template_edges());
	//sedlamat::display(general_hough.get_src_edges());
	//~ sedlamat::display(general_hough.get_src_img());
	//~ sedlamat::display(general_hough.get_template_img());
	//~ sedlamat::print("hotovo");
	//cv::imwrite(argv[2],general_hough.get_result_img());
	std::time(&t1);
	std::cout << difftime(t1,t0) <<std::endl;
	//sedlamat::display(general_hough.get_result_img());
    } catch (const char *e) {
	std::cout << e << std::endl;
	exit(1);
    }
    catch (...) {
	std::cout << "Unexpected error." << std::endl;
    }
    return 0;
}
