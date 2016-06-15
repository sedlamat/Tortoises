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
	time_t t0 = clock();
	passwd *pw = getpwuid(getuid());
	std::string path(pw->pw_dir);
//~ //~
	std::string file_name = "Tg58400.jpg";
	std::string imgs_path = path + "/Images/Tortoises/";
	std::string img_path = imgs_path + file_name;
	//~ std::string img_path = argv[1];
	cv::Mat plastron_img = cv::imread(img_path,1);


	std::string templ_path = "plastron_template.bmp";
	cv::Mat plastron_template = cv::imread(templ_path, 0);

	cv::Point reference_point(65,125);

	std::vector<int> angles = {0,84,180,-96};


	//~ std::map<std::string, cv::Point> map_template_junctions;
	//~ map_template_junctions["1"] = cv::Point(64,58);
	//~ map_template_junctions["2"] = cv::Point(64,74);
	//~ map_template_junctions["3"] = cv::Point(64,97);
	//~ map_template_junctions["4"] = cv::Point(64,105);
	//~ map_template_junctions["5"] = cv::Point(64,155);
	//~ map_template_junctions["6"] = cv::Point(64,168);
	//~ map_template_junctions["7"] = cv::Point(64,185);
	//sedlamat::display(img);
	sedlamat::GeneralHough general_hough(plastron_img,
					     plastron_template,
					     reference_point,
					     angles,
					     20, 200, 1.0, 0.3, 0);
					     //map_template_junctions);
	//sedlamat::display(img);
	//sedlamat::print("running");
	//sedlamat::display(general_hough.get_src_edges());
	general_hough.detect();
	//~ time_t t1 = clock();
	std::cout << general_hough.get_best_accum_val() << std::endl;

	//~ sedlamat::display(general_hough.get_template_edges());
	sedlamat::display(general_hough.get_src_edges());
	//~ sedlamat::display(general_hough.get_src_img());
	//~ sedlamat::display(general_hough.get_template_img());
	//~ sedlamat::print("hotovo");
	//~ cv::imwrite(argv[2],general_hough.get_result_img());
	time_t t1 = clock();
	std::cout << (t1-t0)*1.0/CLOCKS_PER_SEC << std::endl;
	sedlamat::display(general_hough.get_result_img());
    } catch (std::string e) {
	std::cout << "Error when processing " << argv[1] << std::endl;
	std::cout << " " << e << std::endl;
	exit(1);
    }
    return 0;
}
