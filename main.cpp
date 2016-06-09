/* STANDARD C++ LIBRARIES */
#include <iostream>
#include <string>
#include <limits>
#include <cmath>
#include <iomanip>

/* THIRD PARTY LIBRARIES */
#include <pwd.h>
#include <unistd.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"
#include "my_general_hough.hpp"
//#include "tortoise.hpp"

int main(int argc, char *argv[])
{
    try {
	passwd *pw = getpwuid(getuid());
	std::string path(pw->pw_dir);

	std::string ght_path = path + "/Images/Generalized_Hough_Transform/";
	std::string templ_path = ght_path + "tortoise_template.bmp";
	std::string templ_edges_path = ght_path + "ght_template_edges1.bmp";
	cv::Mat templ = cv::imread(templ_path, 0);
	cv::Mat templ_edges = cv::imread(templ_edges_path, 0);
	cv::Mat ref_point_loc;
	cv::findNonZero(templ_edges==127,ref_point_loc);
	cv::Point ref_point(ref_point_loc.at<int>(0,0),
			ref_point_loc.at<int>(0,1));

	std::string file_name = "Tg59800.jpg";
	std::string imgs_path = path + "/Images/Tortoises/";
	std::string img_path = imgs_path + file_name;
	//std::string img_path = argv[1];
	cv::Mat img = cv::imread(img_path,1);

	//sedlamat::display(img);
	sedlamat::GeneralHough general_hough(img, templ, ref_point,
					 4, 15, 100, 1.0, 0.3, 1);
	//sedlamat::display(img);
	//sedlamat::print("running");
	general_hough.run();
	sedlamat::display(general_hough.get_result_img());
	//~ sedlamat::display(general_hough.get_template_edges());
	//~ sedlamat::display(general_hough.get_src_edges());
	//~ sedlamat::display(general_hough.get_src_img());
	//~ sedlamat::display(general_hough.get_template_img());
	//~ sedlamat::print("hotovo");
	//cv::imwrite(argv[2],general_hough.get_result_img());
    } catch (std::string e) {
	std::cout << "Error when processing " << argv[1] << std::endl;
	std::cout << " " << e << std::endl;
	exit(1);
    }
    return 0;
}
