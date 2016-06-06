/* STANDARD C/C++ LIBRARIES */
#include <iostream>
#include <string>

/* THIRD PARTY LIBRARIES */
#include <pwd.h>
#include <unistd.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"
#include "my_general_hough.hpp"
//#include "tortoise.h"

using namespace std;
using namespace cv;
using namespace my;


int main(int argc, char *argv[])
{
    try {
	passwd *pw = getpwuid(getuid());
	string path(pw->pw_dir);

	string ght_path = path + "/Images/Generalized_Hough_Transform/";
	string templ_path = ght_path + "ght_template.bmp";
	string templ_edges_path = ght_path + "ght_template_edges1.bmp";
	Mat templ = imread(templ_path, 0);
	Mat templ_edges = imread(templ_edges_path, 0);
	Mat ref_point_loc;
	findNonZero(templ_edges==127,ref_point_loc);
	Point ref_point(ref_point_loc.at<int>(0,0),
			ref_point_loc.at<int>(0,1));

	string file_name = "Tg59800.jpg";
	string imgs_path = path + "/Images/Tortoises/";
	string img_path = imgs_path + file_name;
	//string img_path = argv[1];
	Mat img = imread(img_path,1);

	my::display(img);
	sedlamat::GeneralHough general_hough(img, templ, ref_point);
	my::display(general_hough.get_template_edges());
	my::display(general_hough.get_src_edges());
	my::display(general_hough.get_src_img());
	my::display(general_hough.get_template_img());
	prt(general_hough.get_template_max_size());
	prt("hotovo");
	exit(0);

	double resize_koef = 150.0/max(img.rows,img.cols);
	resize(img, img, Size(0,0), resize_koef, resize_koef);
	//my::display(img);
	Mat edges = get_edges_color_based(img);

	//my::display(edges);
	//exit(0);
	Mat dst = sedlamat::general_hough_fit_on_img(templ, templ_edges,
						ref_point, img, edges);
	//my::display(dst);
	//exit(0);

	//imwrite(argv[2], dst);
	//cout << "image written" << endl;
    } catch (string e) {
	cout << "Error when processing " << argv[1];
	cout << " " << e << endl;
	exit(0);
    }
    return 0;
}
