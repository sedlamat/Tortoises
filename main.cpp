#include <iostream>

class Something
{
public:
    static int s_value;
public:
    static int getValue();
};


int Something::getValue() { return s_value; } // static member function


int Something::s_value = 1; // initializer

int main()
{
    std::cout << Something::getValue() << '\n';
}


//~ /* STANDARD C++ LIBRARIES */
//~ #include <iostream>
//~ #include <string>
//~ #include <limits>
//~ #include <cmath>
//~ #include <iomanip>
//~ #include <ctime>

//~ /* THIRD PARTY LIBRARIES */
//~ #include <pwd.h>
//~ #include <unistd.h>
//~ #include <opencv2/core/core.hpp>
//~ #include <opencv2/highgui/highgui.hpp>
//~ #include <opencv2/imgproc/imgproc.hpp>

//~ /* FIRST PARTY LIBRARIES */
//~ #include "tortoise.hpp"

//~ int main(int argc, char *argv[])
//~ {
    //~ try {
	//~ time_t t0, t1;
	//~ std::time(&t0);
	//~ passwd *pw = getpwuid(getuid());
	//~ std::string path(pw->pw_dir);
	//~ std::string file_name = R"(Tg30000.jpg)";
	//~ std::string imgs_path = path + R"(/Images/Tortoises/)";
	//~ std::string img_path = imgs_path + file_name;

	//~ //std::string img_path = argv[1];

	//~ cv::Mat plastron_img = cv::imread(img_path,1);

	//~ Tortoise tg(plastron_img);

	//~ //cv::imwrite(argv[2],general_hough.get_result_img());

	//~ std::time(&t1);
	//~ std::cout << difftime(t1,t0) <<std::endl;
    //~ } catch (const char *e) {
	//~ std::cout << e << std::endl;
	//~ exit(1);
    //~ }
    //~ catch (...) {
	//~ std::cout << "Unexpected error." << std::endl;
    //~ }
    //~ return 0;
//~ }
