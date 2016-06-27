/* STANDARD C++ LIBRARIES */
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

	const passwd * const pw = getpwuid(getuid());
	const std::string path(pw->pw_dir);
	const std::string file_name = R"(Tg30000.jpg)";
	const std::string imgs_path = path + R"(/Images/Tortoises/)";
	const std::string img_path = imgs_path + file_name;

	//const std::string img_path = argv[1];

	const cv::Mat plastron_img = cv::imread(img_path,1);

	Tortoise tg(plastron_img, file_name);

	//cv::imwrite(argv[2],general_hough.get_result_img());

	std::time(&t1);
	std::cout << difftime(t1,t0) <<std::endl;

    } catch (const char *e) {
	std::cout << e << std::endl;
	exit(1);
    }
    catch (...) {
	std::cout << "Unexpected error." << std::endl;
    }
    return 0;
}
