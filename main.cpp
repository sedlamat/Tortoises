/* STANDARD C/C++ LIBRARIES */
#include <iostream>
#include <string>
#include <algorithm>
#include <pwd.h>
#include <unistd.h>
#include <bitset>

/* THIRD PARTY LIBRARIES */
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

/* FIRST PARTY LIBRARIES */
#include "my_img_proc.hpp"
#include "my_general_hough.hpp"
//#include "tortoise.h"

using namespace std;
using namespace cv;

int main(int argc, char **argv)

{ 
  try
  {
    passwd* pw = getpwuid(getuid());
    string path(pw->pw_dir);
    
    string ght_path = path + "/Images/Generalized_Hough_Transform/";
    string templ_path = ght_path + "ght_template.bmp";
    string templ_edges_path = ght_path + "ght_template_edges.bmp";
    Mat templ = imread(templ_path, 0);
    Mat templ_edges = imread(templ_edges_path, 0);
    Mat ref_point_loc;
    findNonZero(templ_edges==127,ref_point_loc);
    Point ref_point(ref_point_loc.at<int>(0,0),
		    ref_point_loc.at<int>(0,1));
		    
    cout << ref_point<< endl;
    my::display(templ);
    my::display(templ_edges);
    
    string file_name = "Tg00301.pnm";
    string imgs_path = path + "/Images/Tortoises/";
    string img_path = imgs_path + file_name;
    Mat img = imread(img_path,1);
    double resize_koef = 256.0/max(img.rows,img.cols);
    resize(img, img, Size(0,0), resize_koef, resize_koef);
    
    Mat edges = my::get_edges_color_based(img);
    cout << "here" << endl;
    my::get_hough_points(img, edges);
    my::display(img);
  }
  catch( cv::Exception& e )
  {
    const char * err_msg = e.what();
    std::cout << "exception caught: " << err_msg << std::endl;
  }
    return 1;
}


/*
void oneTortoiseRecognition(string imgDirectory, string tortoiseName, string processID)
{
	try
	{
		Tortoise testudoGraeca(imgDirectory, tortoiseName, processID);

		cout << tortoiseName << " ... STARTS." << endl;

    testudoGraeca.Recognition();	
    //testudoGraeca.Classification();

		cout << tortoiseName << " ... DONE." << endl;
	}
	catch(...)
	{
		cout << "Processing of " << tortoiseName << " ended with error!" << endl;
	}
}

int main(int argc, char *argv[])
{

  for(int groupIndex = 0; groupIndex <= 0; ++groupIndex)
  {
    int isUsingDatabaseWithGoodTGsWithoutTheYoungest = groupIndex;
    string imgDirectory;
    if (isUsingDatabaseWithGoodTGsWithoutTheYoungest)
    {
      imgDirectory.assign("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGoodNoYoungTo6cm\\");
    }
    else
    {
      imgDirectory.assign("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGood\\");
    }
	  
    Classification();
    return 1;
    //kNNclassifier();
    //AccuracyEvaluation(imgDirectory);

    if(argc > 1)
	  {
		  int tortoiseProcess = 0;
		       if (*argv[1] == '1') tortoiseProcess = 1;
		  else if (*argv[1] == '2') tortoiseProcess = 2;
		  else cout << "Wrong arguments!! Enter zelvy.exe [para] where para = 1 or 2"<< endl;

      string numbersOfTGsFileName (imgDirectory + "numbersOfTGs.cimg");
      CImg<int> numberOfTGs(numbersOfTGsFileName.c_str());
      string processID;
      if(tortoiseProcess==1)
      {
        processID.assign("1");
        for( int tg = numberOfTGs.width()/2-1; tg >= 0; --tg )
        {
          stringstream tortoiseNumber;
          tortoiseNumber << setfill('0') << setw(5) << to_string( (long double)numberOfTGs(tg) );
          string tortoiseName ("Tg"+ tortoiseNumber.str());
          oneTortoiseRecognition(imgDirectory, tortoiseName,processID);
        }
      }
      else if(tortoiseProcess==2)
      {
        processID.assign("2");
        for( int tg = numberOfTGs.width()/2; tg < numberOfTGs.width(); ++tg )
        {
          stringstream tortoiseNumber;
          tortoiseNumber << setfill('0') << setw(5) << to_string( (long double)numberOfTGs(tg) );
          string tortoiseName ("Tg"+ tortoiseNumber.str());
          oneTortoiseRecognition(imgDirectory, tortoiseName,processID);
        }
      }
	  }
	  else
	  {
		  bool DoOnlyOne = 1;
		  bool DoMore = 1;
		  bool DoAll = 1;

      string processID("0");
		  if(DoOnlyOne)
		  {
       //char tortoiseName[10] = {'T' , 'g' , 0+48,  7+48, 0+48, 1+48 , 6+48 };
       //char tortoiseName[10] = {'T' , 'g' , 1+48,  9+48, 0+48, 0+48 , 6+48 };
       //char tortoiseName[10] = {'T' , 'g' , 5+48,  9+48, 6+48, 0+48 , 0+48 };
			 //char tortoiseName[10] = {'T' , 'g' , 3+48, 6+48,7+48, 0+48 , 0+48 };
			 char tortoiseName[10] = {'T' , 'g' , 0+48, 2+48, 7+48, 0+48 ,1+48 }; 
        oneTortoiseRecognition(imgDirectory, tortoiseName,processID);
		  }
		  else if (DoMore)
		  {
			  for(int i = 1; i <= 16; ++i)
			  {
				       if (i== 1) { char tortoiseName[10] = {'T' , 'g' , 5+48, 3+48, 1+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i== 2) { char tortoiseName[10] = {'T' , 'g' , 5+48, 9+48, 1+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i== 3) { char tortoiseName[10] = {'T' , 'g' , 5+48, 2+48, 9+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i== 4) { char tortoiseName[10] = {'T' , 'g' , 5+48, 9+48, 8+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  //else if (i== 5) { char tortoiseName[10] = {'T' , 'g' , 5+48, 4+48, 0+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i== 6) { char tortoiseName[10] = {'T' , 'g' , 0+48, 1+48, 5+48, 0+48 , 5+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i== 7) { char tortoiseName[10] = {'T' , 'g' , 0+48, 5+48, 1+48, 0+48 , 4+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i== 8) { char tortoiseName[10] = {'T' , 'g' , 3+48, 9+48, 6+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i== 9) { char tortoiseName[10] = {'T' , 'g' , 4+48, 5+48, 2+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i==10) { char tortoiseName[10] = {'T' , 'g' , 3+48, 0+48, 0+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i==11) { char tortoiseName[10] = {'T' , 'g' , 3+48, 0+48, 7+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i==12) { char tortoiseName[10] = {'T' , 'g' , 5+48, 3+48, 5+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i==13) { char tortoiseName[10] = {'T' , 'g' , 5+48, 9+48, 6+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i==14) { char tortoiseName[10] = {'T' , 'g' , 5+48, 2+48, 3+48, 0+48 , 0+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i==15) { char tortoiseName[10] = {'T' , 'g' , 1+48, 9+48, 0+48, 0+48 , 6+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
				  else if (i==16) { char tortoiseName[10] = {'T' , 'g' , 1+48, 1+48, 5+48, 0+48 , 4+48 }; oneTortoiseRecognition(imgDirectory, tortoiseName, processID); }
			  }
		  }
 		  else if(DoAll)
		  {
        string numbersOfTGsFileName (imgDirectory + "numbersOfTGs.cimg");
        CImg<int> numberOfTGs(numbersOfTGsFileName.c_str()); // numbers of all TGs in the given database directory
	      //for( int tg = numberOfTGs.width()-1; tg >= 0; --tg )
	      for( int tg = 0; tg <  numberOfTGs.width(); ++tg )
        {
          stringstream tortoiseNumber;
          tortoiseNumber << setfill('0') << setw(5) << to_string( (long double)numberOfTGs(tg) ); // puts together 5 digits TG number - e.g.: 00310
          string tortoiseName ("Tg"+ tortoiseNumber.str());
          oneTortoiseRecognition(imgDirectory, tortoiseName, processID);
        }
		  }
	  }
  }

 return 1;
}
//*/
