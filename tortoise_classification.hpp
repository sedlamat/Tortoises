





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////                                                                                                                   ///////////////////////////////////////////////////
////////////////////////////  Methods for measuring the seam segments by hand, their evaluation and a tortoise classification based on them.   ///////////////////////////////////////////////////
////////////////////////////                                                                                                                   ///////////////////////////////////////////////////
////////////////////////////                   The methods and functions are not part of the Tortoise class.                                   ///////////////////////////////////////////////////
////////////////////////////                                                                                                                   ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  //! event Call Back Function used for recording mouse clicks -> for hand measuring 
  /**
      \note 
      - Procedure that retrieves x and y coordinate of a mouse click and saves them into userdata array
      \outline
      - When right mouse button is pressed, x and y coordinate of the mouse poiter is saved
        into userdata array. Last element in userdata stores place of next free elements where
        next point will be saved.
        If CTRL key is pressed when pressing right mouse click, then -1 is assign to userdata[0] 
        (a flag, which will cause skipping the image without saving the measurements).
        NumberOfPointsToBeMeasured has to be passed here manually.
  **/
  void CallBackFunc(int event, int x, int y, int flags, void* userdata)
  {
    const int numberOfPointsToBeMeasured = 12;
	  if ( event == cv::EVENT_RBUTTONDOWN)
	  {
		  ((int*)userdata)[((int*)userdata)[2*numberOfPointsToBeMeasured]] = x;
		  ((int*)userdata)[((int*)userdata)[2*numberOfPointsToBeMeasured]+1] = y;
		  ((int*)userdata)[2*numberOfPointsToBeMeasured] += 2;
		  cout << "Right mouse button is released - position (" << x << ", " << y << ")" << endl;
	  }
	  if( flags == (cv::EVENT_FLAG_CTRLKEY + cv::EVENT_FLAG_RBUTTON))
	  {
		  ((int*)userdata)[0] = -1;
	  }
  }

  //! records position of the mouse pointer when clicked on a specified image
  /**
      \note 
      - Function opens an image and by calling CallBackFunction records coordinates of 
        the right mouse button click.
      \outline
      - Loads image.
        XYcoordinatesOfLeftAndRightJuncitons[24], which contains the current position in the 
        XYcoordinatesOfLeftAndRightJuncitons array, where recorded coordinates are saved, 
        is set to 0, so the saving of coordinates starts at the beginning of the array.
        setMouseCallback -  ensures that while mouse keys are pressed, the CallBackFunc
        save mouse pointer coordinates into XYcoordinatesOfLeftAndRightJuncitons array.
        When any other key (i.e. keyboard key) is pressed, the function ends and returns
        XYcoordinatesOfLeftAndRightJuncitons array.
  **/
  void gettingXYcoordinatesOfMouseClick(string loadFileName, int* XYcoordinatesOfLeftAndRightJunctions, const int numberOfPointsToBeMeasured)
  {
	  // Read image from file 
	  cv::Mat img = cv::imread(loadFileName);
	  //if fail to read the image
	  if( img.empty() ) 
	  { 
		  cout << "Error loading the image" << endl;
	  }
	  //Create a window
	  cv::namedWindow("My window", CV_WINDOW_NORMAL);
	  XYcoordinatesOfLeftAndRightJunctions[2*numberOfPointsToBeMeasured] = 0;
	  //set the callback function for any mouse event
	  cv::setMouseCallback("My window", CallBackFunc, XYcoordinatesOfLeftAndRightJunctions);
	  //show the image
	  cv::imshow("My window", img);
	  // Wait until user press some key
	  cv::waitKey(0);
  }

   //! records position of the central seam junctions for all tortoises in loadFileDirectory
  /**
      \note 
      - Calls gettingXYcoordinatesOfMouseClick to record all 12 points.
      \outline
      - All tortoises in loadFileDirectory must have their tortoiseNumber, e.g. Tg00301, writen
        in TortoisesNumbers.txt which has to be in loadFileDirectory
        If all 12 points are not recorded, it calls gettingXYcoordinatesOfMouseClickto again.
        Sequence in which points are to be recorded: 
          1. centralSeamEndHead
          2.-6. leftJunctions
          7. centralSeamEndTail
          8.-12. rightJunctinos
	      For each point x and y coordinates are recorded, so the array contains 24 numbers in total.
        If any of the images are not to be measured(recorded), then press CTRL key and right button key. It 
        assigns XYcoordinatesOfLeftAndRightJuncitons[0] = -1, a flag that ensures data will not be
        recorded. To go to anothe image, press any keyboard key.
        Recorded date is saved in .cimg file measuredPointsByHandOfSeamSegments, in saveFileDirectory.
  **/
  void getAndSaveCoordinatesOfCentralSeamJunctionsForAllBenderTortoises(const int numberOfPointsToBeMeasured)
  {
    ifstream tortoisesNumbers ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnm\\tortoisesNumbers.txt");
    string loadFileDirectory ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnm\\");
    string saveFileDirectory ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizp\\");
    string tortoiseNumber;
    for(int i = 0; i < 1304; ++i) // 1304 is the number of tortoises in tortoisesNumbers.txt
    {
      getline( tortoisesNumbers, tortoiseNumber);
      if(i >= 1005) // condition if you want to measure only some of the tortoises
      {
        string loadFileName (loadFileDirectory + tortoiseNumber + ".pnm");
        int* XYcoordinatesOfLeftAndRightJunctions = new int[2*numberOfPointsToBeMeasured+1]; // 2(for x,y)*numberOfPointsToBeMeasured+1(for the count of measured points)
        gettingXYcoordinatesOfMouseClick( loadFileName, XYcoordinatesOfLeftAndRightJunctions, numberOfPointsToBeMeasured);
	      while(XYcoordinatesOfLeftAndRightJunctions[2*numberOfPointsToBeMeasured] != 2*numberOfPointsToBeMeasured && XYcoordinatesOfLeftAndRightJunctions[0] != -1)
	      {
		      cout << "NEW TRY" << endl;
		      gettingXYcoordinatesOfMouseClick( loadFileName, XYcoordinatesOfLeftAndRightJunctions, numberOfPointsToBeMeasured);
	      }
	      if(XYcoordinatesOfLeftAndRightJunctions[0] != -1) // if == -1 it skips saving for current image
	      {
		      CImg<double> measuredPointsByHandOfSeamSegments(1,2*numberOfPointsToBeMeasured,1,1, 0);
		      for(int i = 0; i < 2*numberOfPointsToBeMeasured; i++)
		      {
			      measuredPointsByHandOfSeamSegments(0,i) = XYcoordinatesOfLeftAndRightJunctions[i]; // copying data to CImg that will be saved
		      }
          string saveFileName(saveFileDirectory + tortoiseNumber + "measuredPointsByHandOfSeamSegments.cimg");
          measuredPointsByHandOfSeamSegments.save(saveFileName.c_str());
          cout << "Measured data for " << tortoiseNumber << " was saved." << endl;
	      }
      }
    }
  }
  
  //! calculate seam segments lengths and inner junction differences from measured points
  /**
      \note 
      - Load the meaured points and return seam segments lengths and inner junction differences
      \outline
      - First the measured points are loaded from loadFileDirectoryOfMeasuredPoints, the files
        have form of for example "Tg03301measuredPointsByHandOfSeamSegments.cimg".
        The directory has to be changed depending on where the measured points are stored.
        First the loaded points are ordered in a new array MP of size 2(left,right)x7(junctions)x2(x,y-coordinates).
        The first and last junctions on left and right will have the same cooordinates, for the head and tail 
        junctions had only one measurement.
        Seam segments are saved in two columns of 6 elements, where in left column are left segments and in right column 
        are right segments.
        Segment is not a mere L2 distance between cosecutive measured point on the left or the right. Because the central
        seam is not straight always, espectially in inner-junction space.
        So the segments lenghts are made of the distance between the nearest left and/or right connections of two consecutive 
        juntions and then the left-right/right-left inner-junction distances are added.
        In the same time, the leftRight inner junction differences are computed and saved in leftRightJunctionDifferences.
        Both arrays are then normalized by the central seam lenght and data is saved in the saveFileDirectory directory.
  **/
  void calculateSeamLengthsAndJunctionDifferencesFromMeasuredPoints()
  {
    ifstream tortoisesNumbers ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnm\\tortoisesNumbers.txt");
    string loadFileDirectoryOfMeasuredPoints ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizp\\");
    string saveFileDirectory ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\data\\");
    string tortoiseNumber;
    for(int i = 0; i < 1304; ++i)
    {
      getline( tortoisesNumbers, tortoiseNumber);
      if(i>=1005)
      {
        string loadNameOfMeasuredPoints(loadFileDirectoryOfMeasuredPoints + tortoiseNumber + "measuredPointsByHandOfSeamSegments.cimg");
        CImg<double> measuredPoints(loadNameOfMeasuredPoints.c_str()), MP(2,7,1,2,0);
        for(int i = 0; i <= 1; i++)
        {
          MP(i,0,0) = measuredPoints(0, 0); // left/right head junction x-coord
          MP(i,0,1) = measuredPoints(0, 1); // left/right head junction y-coord
          MP(i,6,0) = measuredPoints(0,12); // left/right tail junction x-coord
          MP(i,6,1) = measuredPoints(0,13); // left/right tail junction y-coord   
        }
        for(int i = 1; i <= 5; i++)
        {
          MP(0,i,0) = measuredPoints(0,i*2   ); // left  i-th junction x-coord
          MP(0,i,1) = measuredPoints(0,i*2+ 1); // left  i-th junction y-coord
          MP(1,i,0) = measuredPoints(0,i*2+12); // right i-th junction x-coord
          MP(1,i,1) = measuredPoints(0,i*2+13); // right i-th junction y-coord       
        }
        //calculate segments lenghts and junction differences at the same time
	      CImg<double> seamSegments(2,6,1,1, 0), leftRightJunctionDifferences(1,5,1,1, 0);
        double lastNearestX = MP(0,0,0);
        double lastNearestY = MP(0,0,1);
        for(int i = 0; i < 6; i++)
        {
          double distanceFromLastNearestPointOnLeft = sqrt( (double) (lastNearestX-MP(0,i+1,0))*(lastNearestX-MP(0,i+1,0)) 
                                                                   + (lastNearestY-MP(0,i+1,1))*(lastNearestY-MP(0,i+1,1)) );
          double distanceFromLastNearestPointOnRight = sqrt( (double) (lastNearestX-MP(1,i+1,0))*(lastNearestX-MP(1,i+1,0)) 
                                                                    + (lastNearestY-MP(1,i+1,1))*(lastNearestY-MP(1,i+1,1)) );
          double differenceBetweenLeftRightJunctionPoints = sqrt( (double) (MP(0,i+1,0)-MP(1,i+1,0))*(MP(0,i+1,0)-MP(1,i+1,0)) 
                                                                         + (MP(0,i+1,1)-MP(1,i+1,1))*(MP(0,i+1,1)-MP(1,i+1,1)) );
          if (distanceFromLastNearestPointOnLeft < distanceFromLastNearestPointOnRight)
          {
            seamSegments(0,i) += distanceFromLastNearestPointOnLeft;
            seamSegments(1,i) += distanceFromLastNearestPointOnLeft + differenceBetweenLeftRightJunctionPoints;
            lastNearestX = MP(1,i+1,0);
            lastNearestY = MP(1,i+1,1);
            if (i!=5)
            {
             seamSegments(0,i+1) += differenceBetweenLeftRightJunctionPoints;
             leftRightJunctionDifferences(0,i) = - differenceBetweenLeftRightJunctionPoints;
            }      
          }
          else
          {
            seamSegments(1,i) += distanceFromLastNearestPointOnLeft;
            seamSegments(0,i) += distanceFromLastNearestPointOnLeft + differenceBetweenLeftRightJunctionPoints;
            lastNearestX = MP(0,i+1,0);
            lastNearestY = MP(0,i+1,1);
            if (i!=5)
            {
             seamSegments(1,i+1) += differenceBetweenLeftRightJunctionPoints;
             leftRightJunctionDifferences(0,i) = differenceBetweenLeftRightJunctionPoints;
            }  
          }
        }
        //normalizing to the size of the plastron
	      double lengthOfCentralSeam = seamSegments.get_column(0).sum();
	      for(int i = 0; i < 6; i++) 
        {
          seamSegments(0,i) /= lengthOfCentralSeam;
          seamSegments(1,i) /= lengthOfCentralSeam;
          if (i!=5) leftRightJunctionDifferences(0,i) /= lengthOfCentralSeam;
        }
        //genarate save directory
        string saveNameForSeamSegments(saveFileDirectory + tortoiseNumber + "seamSegments.cimg");
        string saveNameForLeftRightJunctionDifferences(saveFileDirectory + tortoiseNumber + "leftRightJunctionDifferences.cimg");
        //save in C-string directory
        seamSegments.save(saveNameForSeamSegments.c_str());
        leftRightJunctionDifferences.save(saveNameForLeftRightJunctionDifferences.c_str());
      }
    }
  }

  //! a procedure that puts all measured seam segments lengths and inner junction differences into one file
  /**
      \note 
      - We want all seam segments lenghts and all inner junction distances in one file 
        for all 1304 Bender and Cizp tortoises. And we want left and right seam segments lengths and inner
        junction differences in one column.
      \outline
      - With use of numbering of Bender and Cizp Tortoises in tortoisesNumbers we go though all of them (segments lengths 
        and inner junction distances) and append them to the file MeasuredFeaturesForAllTortoises which are saved at the end.
        Also seamSegments that are in two columns, are put in one column (column 0 is appended below column 1) and
        below that the inner junction differences are appended.
  **/
  void putTheCalculatedSeamDistancesForAllTortoisesIntoOneFile()
  {
    ifstream tortoisesNumbers ("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnm\\tortoisesNumbers.txt");
    CImg<double> MeasuredFeaturesForAllTortoises;
    string loadFileDirectoryOfMeasuredFeatures ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\data\\");
    string tortoiseNumber;
    for(int i = 0; i < 1304; ++i)
    {
      getline( tortoisesNumbers, tortoiseNumber);

      string loadSeamSegments(loadFileDirectoryOfMeasuredFeatures + tortoiseNumber + "seamSegments.cimg");
      string loadJunctionDifferences(loadFileDirectoryOfMeasuredFeatures + tortoiseNumber + "leftRightJunctionDifferences.cimg");
      CImg<double> seamSegments(loadSeamSegments.c_str()), junctionDifferences(loadJunctionDifferences.c_str());
      seamSegments.column(0).append(seamSegments.get_column(1),'y').append(junctionDifferences,'y');
      MeasuredFeaturesForAllTortoises.append(seamSegments,'x');
    }
    MeasuredFeaturesForAllTortoises.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\MeasuredFeaturesForAllTortoises.cimg");
  }

  //! save TG numbers, read from file names in a given directory, into a cimg file 
  /**
      \outline 
      - We want it to read file names, in our case pnm files fo chosen selection of TGs, from a given directory.
        Then, for each TG file name, extract the TG number and convert it to integer (Tg00301 -> 00301 -> 301) and then save it
        into CImg variable, which is saved at the end.
  **/
  void getAndSaveTGnumbersOfTGfileNamesFromGivenDirectory()
  {
    CImg<double> SomeTgNumbers;
    _chdir("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGoodNoYoungTo6cm\\"); //changes the directory
    struct _finddata_t c_file;
    long hFile;
    if ( (hFile = _findfirst("*.pnm", &c_file)) == -1L ) //checks if there is any pnm file in the directory
      printf("No *.pnm files in current directory");
    else
    {
      do
      {
        int tgNumberINT;
        string tgNumber( c_file.name); //get the file name
        istringstream (tgNumber.substr(2,5) ) >> tgNumberINT; //get the number from the name and convert it to integer
        CImg<int> TGnumber(1,1,1,1, tgNumberINT);
        SomeTgNumbers.append(TGnumber,'x');
      } while ( _findnext(hFile, &c_file) == 0 ); //find next file name
      _findclose(hFile);
    } 
    SomeTgNumbers.display();
    SomeTgNumbers.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\TgNumbersInINTonlyGoodTGsAndNoYoungTGs.cimg");
  }

  //! given a subset of TG numbers ( < 1304), it gets corresponding subset of  measured features
  /**
      \outline 
      - Based on the set of TG numbers of all TGs and its subset of some TG numbers, it will select
        the corresponding measured features for the subset of TGs and will save it.
  **/
  void getAndSaveTgFeatureSubSetForAnyChosenSubSetOfTGnumbers()
  {
    CImg<double> MeasuredFeaturesForAllTortoises("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\MeasuredFeaturesForAllTortoises.cimg");
    CImg<int> AllTgNumbers("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\AllTgNumbersInINT.cimg");
    CImg<int> SomeTgNumbers("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\TgNumbersInINTonlyGoodTGs.cimg");
    int numOfFeatures = MeasuredFeaturesForAllTortoises.height();
    for(int i = 0; i < AllTgNumbers.width(); ++i)
    {
      if( (SomeTgNumbers - AllTgNumbers(i)).abs().min() > 0 ) //checks if a TG number is not in the subset of TGs
      {
        for(int j = 0; j < MeasuredFeaturesForAllTortoises.height(); ++j) //for such a TG all its measurements are set to 1000
        {
          MeasuredFeaturesForAllTortoises(i,j) = 1000;
        }
      }
    }
    for(int i = 0; i < MeasuredFeaturesForAllTortoises.width(); ++i) //eliminates measured features of TG not belonging to the subset
    {
      if( MeasuredFeaturesForAllTortoises.get_column(i).sum() == numOfFeatures*1000 )
      {
        if (i==0)
        {
          MeasuredFeaturesForAllTortoises = MeasuredFeaturesForAllTortoises.get_columns(1,MeasuredFeaturesForAllTortoises.width()-1);
        }
        else if (i==MeasuredFeaturesForAllTortoises.width()-1)
        {
          MeasuredFeaturesForAllTortoises = MeasuredFeaturesForAllTortoises.get_columns(0,MeasuredFeaturesForAllTortoises.width()-2);
        }
        else
        {
          MeasuredFeaturesForAllTortoises = MeasuredFeaturesForAllTortoises.get_columns(0,i-1).append(MeasuredFeaturesForAllTortoises.get_columns(i+1,MeasuredFeaturesForAllTortoises.width()-1),'x');
        }
        --i;
      }
    }
    MeasuredFeaturesForAllTortoises.display();
    MeasuredFeaturesForAllTortoises.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\MeasuredFeaturesForOnlyGoodTGs.cimg");
  }

  //! get and save groups of features and groups of numbers of same tortoises and of different tortoises into a CImgList file
  /**
      \outline 
      - For a given selection of TG measured features and coresponding TG numbers, it compose all possible groups of TGs,
        where in a group there are multiple measurements/names (of different images) of the same individual. In the last 
        group there are measurements/names for all individual tortoises, for each only one measurement.
  **/
  void getAndSaveGroupsOfFeaturesAndGroupsOfNumbersOfSameAndDifferentTGs()
  {
    CImg<double> MeasuredFeaturesForSelectedTortoises ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\MeasuredFeaturesForOnlyGoodTGsAndNoYoungTGs.cimg");
    CImg<int> NumbersOfSelectedTortoises ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\numbersOfTGsInINTonlyGoodTGsAndNoYoungTGs.cimg");
    CImgList<double> groupsOfFeatures;
    CImgList<int> groupsOfTGnumbers;
    int indexOfCurrentTGgroup = 0;  // numbers/measurements index to the TG that already belong to the next group
    while (indexOfCurrentTGgroup < NumbersOfSelectedTortoises.width()-1)
    {
      int currentTGnumber = floor(NumbersOfSelectedTortoises(indexOfCurrentTGgroup)*1.0/100); //gets the number of current TG individual (301 -> 3, 12603 -> 126)
      int indexOfNextTGgroup = indexOfCurrentTGgroup + 1;
      int nextTGnumber = floor(NumbersOfSelectedTortoises(indexOfNextTGgroup)*1.0/100); //gets the number of next TG individual
      if(currentTGnumber==nextTGnumber) //if next measurement belongs to the same individual, than it will forms a new group
      {
        while(currentTGnumber==nextTGnumber && indexOfNextTGgroup < NumbersOfSelectedTortoises.width())//finds the last number/measurement that belong to the same individual
        {
          ++indexOfNextTGgroup;
          nextTGnumber = floor(NumbersOfSelectedTortoises(indexOfNextTGgroup)*1.0/100);
        }
        groupsOfFeatures.push_back(MeasuredFeaturesForSelectedTortoises.get_columns(indexOfCurrentTGgroup,indexOfNextTGgroup-1)); //stors the new group of measurements
        groupsOfTGnumbers.push_back(NumbersOfSelectedTortoises.get_columns(indexOfCurrentTGgroup,indexOfNextTGgroup-1)); //stors the new group of TG numbers
        indexOfCurrentTGgroup = indexOfNextTGgroup; //set index to the next TG individual
      }
      else
      {
        indexOfCurrentTGgroup = indexOfNextTGgroup; //set index to the next TG individual
      }
    }
    indexOfCurrentTGgroup = 0;
    CImg<double> groupOfDifferentTortoises; //each individual will be in the group only once
    CImg<int> numbersOfDifferentTortoises;
    while (indexOfCurrentTGgroup < NumbersOfSelectedTortoises.width())
    {
      if(indexOfCurrentTGgroup==NumbersOfSelectedTortoises.width()-1) //if index is on the last TG it has therefore only single measurement TG so it is stored
      {
        groupOfDifferentTortoises.append(MeasuredFeaturesForSelectedTortoises.get_column(indexOfCurrentTGgroup),'x');
        numbersOfDifferentTortoises.append(NumbersOfSelectedTortoises.get_column(indexOfCurrentTGgroup),'x');
        ++indexOfCurrentTGgroup;
      }
      else //for any other index it checks if the individual has one or more measurements
      {
        int currentTGnumber = floor(NumbersOfSelectedTortoises(indexOfCurrentTGgroup)*1.0/100);
        int indexOfNextTGgroup = indexOfCurrentTGgroup + 1;
        int nextTGnumber = floor(NumbersOfSelectedTortoises(indexOfNextTGgroup)*1.0/100);
        if(currentTGnumber==nextTGnumber) //if the individual has more measurements
        {
          while(currentTGnumber==nextTGnumber && indexOfNextTGgroup < NumbersOfSelectedTortoises.width())
          {
            ++indexOfNextTGgroup;
            nextTGnumber = floor(NumbersOfSelectedTortoises(indexOfNextTGgroup)*1.0/100);
          }
          groupOfDifferentTortoises.append(MeasuredFeaturesForSelectedTortoises.get_column(floor((indexOfCurrentTGgroup+indexOfNextTGgroup-1)/2.0)+1),'x'); //it selects one of the measurements (in the middle)
          numbersOfDifferentTortoises.append(NumbersOfSelectedTortoises.get_column(floor((indexOfCurrentTGgroup+indexOfNextTGgroup-1)/2.0)+1),'x'); //it selects one of the measurements (in the middle)
          indexOfCurrentTGgroup = indexOfNextTGgroup;
        }
        else//if the individual doesnt have more measurements it is stored 
        {
          groupOfDifferentTortoises.append(MeasuredFeaturesForSelectedTortoises.get_column(indexOfCurrentTGgroup),'x');
          numbersOfDifferentTortoises.append(NumbersOfSelectedTortoises.get_column(indexOfCurrentTGgroup),'x');
          indexOfCurrentTGgroup = indexOfNextTGgroup;
        }
      }
    }
    groupsOfFeatures.push_back(groupOfDifferentTortoises);
    groupsOfTGnumbers.push_back(numbersOfDifferentTortoises);
    groupsOfTGnumbers.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfTGnumbersForOnlyGoodTGsAndNoYoungTGs.cimg");
    groupsOfFeatures.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfMeasuredFeaturesForOnlyGoodTGsAndNoYoungTGs.cimg");
  }

  //! PCA - get eigen values and eigen vectors of correlation matrix of hand measured lengths for all groups (of same and of different TGs) and return transformed feature vectors
  /**
      \outline 
      - For all the measurements (the segments lengths and the junction difference), that are first standardized it will compute the correlation 
        matrix and that compute and save its eigenvectors and eigenvalues. All is saved in a txt file, plus the correlation
        matrix is computed and saved there as well.
        The method return tranformed measurements.

  **/
  CImg<double> PCAmethod(CImg<double> measurementsStandardized, double PCAthreshold, int whichGroups)
  {
    //measurementsStandardized.display();
	  unsigned int numOfFeatures = measurementsStandardized.width();
	  unsigned int numOfMeasurements = measurementsStandardized.height();
	  CImg<double> measurementsCovarianceMatrix ( (measurementsStandardized.get_transpose()*measurementsStandardized)/(numOfMeasurements-1) ); //computing the covariance matrix
    //measurementsCovarianceMatrix.display();
	  CImg<double> eigenValues, eigenVectors;
	  measurementsCovarianceMatrix.symmetric_eigen(eigenValues,eigenVectors); //computing the eigen vals and vecs
    eigenValues /= eigenValues.sum(); //normalizing the eigen vals
    //(eigenValues,eigenVectors,measurementsCovarianceMatrix).display();
    ofstream eigenvaluesLog("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\eigenvalues.txt",ios_base::app);
    eigenvaluesLog << "Group: " << whichGroups << endl;
    for(int i = 0; i < numOfFeatures; i++) eigenvaluesLog << right << setw(20) << eigenValues(i) << endl;
    eigenvaluesLog << endl << endl;
    //eigenValues.display_graph();

    int numberOfUsefulEigenVectors = 0;
    for(int i = 0; i < numOfFeatures; ++i)
    {
      if(eigenValues(i) > PCAthreshold) // threshold for the eigenvalues
      {
        numberOfUsefulEigenVectors = i;
      }
    }
    if(numberOfUsefulEigenVectors>0) eigenVectors.columns(0,numberOfUsefulEigenVectors);
    else eigenVectors.column(numberOfUsefulEigenVectors);
    CImg<double> transformedMeasurementsStandardized;
    transformedMeasurementsStandardized = measurementsStandardized*eigenVectors;
    //(measurementsStandardized,transformedMeasurementsStandardized).display();
    return transformedMeasurementsStandardized;
  }

  //! get the training/test set, the training/test labels and the training/test TG numbers from a given group a selected features
  /**
      \outline 
      - For a given group of feature we want to get a training set and a test set that will contain differences of features of same
        individuals and of different individuals, in both sets there will be equally differences of pairs of same TGs and different TGs.
        Then we want to get training labels and test labels that will indicate whether a difference of a pair belong to same TGs or
        different TGs, will have values 1 or -1. And finally we want to know TG numbers of all pairs that are is training and test sets.
  **/
  CImgList<double> getTrainingAndTestSetsLabelsAndTGnumbersForAGivenGroupOfFeatures(int whichGroups, bool doStandardization, bool doPCA, double PCAthreshold)
  {
    ofstream log ("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\kNNResults_log.txt",ios_base::app);
    CImgList<double> groupsOfFeatures;
    CImgList<int> groupsOfTGnumbers;
    if(whichGroups==0) // choose which group of features will be used
    {
      groupsOfFeatures.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfMeasuredFeaturesForAllTortoises.cimg");
      groupsOfTGnumbers.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfTGnumbersForAllTGs.cimg");
    }
    else if (whichGroups==1)
    {
      groupsOfFeatures.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfMeasuredFeaturesForOnlyGoodTGs.cimg");
      groupsOfTGnumbers.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfTGnumbersForOnlyGoodTGs.cimg");
    }
    else if (whichGroups==2)
    {
      groupsOfFeatures.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfMeasuredFeaturesForOnlyGoodTGsAndNoYoungTGs.cimg");
      groupsOfTGnumbers.assign("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\groupsOfTGnumbersForOnlyGoodTGsAndNoYoungTGs.cimg");
    }
    log << left <<setw(10) << "Groups:" << setw(20) << whichGroups << setw(15) << "withStandrd:" << setw(10) << doStandardization << setw(10) << "withPCA:" << setw(10) << doPCA << setw(10) << "PCAthres:" << setw(20) << PCAthreshold << endl;

    unsigned int numOfGroups = groupsOfFeatures.size();
    log << setw(30) << "NumberOfGroups:" << setw(20) << numOfGroups << endl;

    unsigned int numOfPairsInTheSetOfSameTGs = 0;
    
    (groupsOfFeatures,groupsOfTGnumbers).display();
    for(unsigned int groupIndex = 0; groupIndex < numOfGroups; ++groupIndex) 
    {
      groupsOfFeatures(groupIndex).transpose(); //transpose of the groups
      if(groupIndex < numOfGroups-1) 
      {
        int numOfTGsInGroup = groupsOfFeatures(groupIndex).height();
        numOfPairsInTheSetOfSameTGs += (numOfTGsInGroup*(numOfTGsInGroup-1)/2);
      }
    }
    //groupsOfFeatures.display();
    unsigned int numOfFeatures = groupsOfFeatures(0).width();
    log << setw(30) << "NumOfFeatures:" << setw(20) << numOfFeatures << endl;
    log << setw(30) << "NumOfSameAndDiffPairs:" << setw(20) << numOfPairsInTheSetOfSameTGs << endl;
    cout << numOfFeatures << " " << numOfPairsInTheSetOfSameTGs << endl;

    cout << 1 << endl;
    CImg<double> setOfDifferencesForSameTGs, setOfTGnumbersForSameTGs; // (i,0) - segment differences, (i,1) - junction differences, (i,2) - label
    for(int groupIndex = 0; groupIndex < numOfGroups -1 ; ++groupIndex) // prepare the set of same TGs and their differences
    {
      CImg<double> groupOfTGs(groupsOfFeatures(groupIndex)), groupOfTGnumbers(groupsOfTGnumbers(groupIndex));
      for(int i = 0; i < groupOfTGs.height(); ++i)
      {
        for(int j = i+1; j < groupOfTGs.height(); ++j)
        {
          CImg<float> class1 (1,1,1,1, 1.0);
          setOfDifferencesForSameTGs.append( (groupOfTGs.get_row(i)-groupOfTGs.get_row(j)).abs().append(class1,'x'), 'y');
          setOfTGnumbersForSameTGs.append(groupOfTGnumbers.get_crop(i,0,i,0).append(groupOfTGnumbers.get_crop(j,0,j,0),'x'),'y');
        }
      }
    }

    cout << 2 << endl;
    CImg<double> groupOfDifferentTGs(groupsOfFeatures(numOfGroups-1)), groupOfNumbersOfDifferentTGs(groupsOfTGnumbers(numOfGroups-1)), setOfDifferencesForDifferentTGs, setOfTGnumbersForDifferentTGs;
    unsigned int numOfDifferentTGs = groupOfDifferentTGs.height();
    unsigned int numOfPairsOfDifferentTGsAllCombinations = numOfDifferentTGs*(numOfDifferentTGs-1)/2;
    unsigned int step = numOfPairsOfDifferentTGsAllCombinations/numOfPairsInTheSetOfSameTGs;
    unsigned int counterOfPairs = 0;
    for(int i = 0; i < numOfDifferentTGs; ++i) // prepare the set of different TGs and their differences
    {
      for(int j = i+step; j < numOfDifferentTGs; j+=step)
      {
        if(counterOfPairs < numOfPairsInTheSetOfSameTGs)
        {
          CImg<float> class2 (1,1,1,1, 0);
          setOfDifferencesForDifferentTGs.append( (groupOfDifferentTGs.get_row(i)-groupOfDifferentTGs.get_row(j)).abs().append(class2,'x'), 'y');
          setOfTGnumbersForDifferentTGs.append(groupOfNumbersOfDifferentTGs.get_crop(i,0,i,0).append(groupOfNumbersOfDifferentTGs.get_crop(j,0,j,0),'x'),'y');
          ++counterOfPairs;
        }
      }
    }

    cout << 3 << endl;

    if(doStandardization)
    {
      CImg<double> measurementsStandardized;
      CImg<double> setOfDifferencesForSameAndDifferentTGs(setOfDifferencesForSameTGs.get_append(setOfDifferencesForDifferentTGs,'y'));
      //(setOfDifferencesForSameTGs,setOfDifferencesForDifferentTGs).display();
	    unsigned int numOfMeasurements = setOfDifferencesForSameAndDifferentTGs.height();
	    CImg<double> meansOfFeatures(numOfFeatures,1);
      CImg<double> stdsOfFeatures(numOfFeatures,1);
	    for(int i = 0; i < numOfFeatures; i++) meansOfFeatures(i) = setOfDifferencesForSameAndDifferentTGs.get_column(i).mean();
      for(int i = 0; i < numOfFeatures; i++) stdsOfFeatures(i) = sqrt(setOfDifferencesForSameAndDifferentTGs.get_column(i).variance());
	    
      for(int i = 0; i < numOfFeatures; i++) measurementsStandardized.append((setOfDifferencesForSameAndDifferentTGs.get_column(i)-meansOfFeatures(i))/stdsOfFeatures(i),'x'); //standardizing the features
      //(setOfDifferencesForSameAndDifferentTGs,measurementsStandardized).display("standardized");
      if(doPCA)
      {
        CImg<double> PCAresult( PCAmethod( measurementsStandardized.get_columns(0,numOfFeatures-1), PCAthreshold, whichGroups ) );
        setOfDifferencesForSameAndDifferentTGs = PCAresult.append(setOfDifferencesForSameAndDifferentTGs.get_column(numOfFeatures),'x');
      }
      else
      {
        setOfDifferencesForSameAndDifferentTGs = measurementsStandardized.append(setOfDifferencesForSameAndDifferentTGs.get_column(numOfFeatures),'x');
      }
      setOfDifferencesForSameTGs = setOfDifferencesForSameAndDifferentTGs.get_rows(0,numOfPairsInTheSetOfSameTGs-1);
      setOfDifferencesForDifferentTGs = setOfDifferencesForSameAndDifferentTGs.get_rows(numOfPairsInTheSetOfSameTGs,2*numOfPairsInTheSetOfSameTGs-1);

    }
    numOfFeatures = setOfDifferencesForDifferentTGs.width()-1;
    log << setw(30) << "NewNumOfFeatures:" << setw(20) << numOfFeatures << endl;
    cout << numOfFeatures << " " << setOfDifferencesForSameTGs.height() << " " << setOfDifferencesForDifferentTGs.height() << endl;
    //(setOfDifferencesForSameTGs,setOfTGnumbersForSameTGs,setOfDifferencesForDifferentTGs,setOfTGnumbersForDifferentTGs).display("sets ater trans");
    //setOfDifferencesForSameTGs


    CImg<double> trainingSet, testSet, trainingLabels, testLabels, trainingSetTGnumbers, testSetTGnumbers;
    for(int i = 0; i < numOfPairsInTheSetOfSameTGs; ++i) // divide sets into training set and test set - each second is taken to test set from sets with differences of same and different TGs 
    {
      if(i%2)
      {
        trainingSet.append(setOfDifferencesForSameTGs.get_row(i),'y');
        trainingSet.append(setOfDifferencesForDifferentTGs.get_row(i),'y');
        trainingSetTGnumbers.append(setOfTGnumbersForSameTGs.get_row(i),'y');
        trainingSetTGnumbers.append(setOfTGnumbersForDifferentTGs.get_row(i),'y');
      }
      else
      {
        testSet.append(setOfDifferencesForSameTGs.get_row(i),'y');
        testSet.append(setOfDifferencesForDifferentTGs.get_row(i),'y');
        testSetTGnumbers.append(setOfTGnumbersForSameTGs.get_row(i),'y');
        testSetTGnumbers.append(setOfTGnumbersForDifferentTGs.get_row(i),'y');
      }
    }
    cout << testSet.height() << " " << trainingSetTGnumbers.height() << endl;
    log << setw(30) << "numOfPairsInTrainingSet:" << setw(20) << trainingSet.height() << endl;
    log << setw(30) << "numOfPairsInTestSet:" << setw(20) << testSet.height() << endl;
    trainingLabels = trainingSet.get_column(numOfFeatures);
    testLabels = testSet.get_column(numOfFeatures);
    trainingSet.columns(0,numOfFeatures-1);
    testSet.columns(0,numOfFeatures-1);
    log << endl << endl;

    //(trainingSet, testSet).display();
    //for(int i = 0; i < trainingSet.height(); ++i) // here you can choose to sum the features
    //{
    //  trainingSet(0,i) = trainingSet.get_row(i).sum(); 
    //}
    //trainingSet.column(0);
    //for(int i = 0; i < testSet.height(); ++i)
    //{
    //  testSet(0,i) = testSet.get_row(i).sum(); 
    //}
    //testSet.column(0);
    //(trainingSet,trainingLabels,testSet,testLabels,trainingSetTGnumbers,testSetTGnumbers);
    return trainingSet,trainingLabels,testSet,testLabels,trainingSetTGnumbers,testSetTGnumbers;
  }

    //void getTestSetForAutoMeasuredFeatures()
  
  //! k - nearest neigbour classifier - OpenCV implementation
  /**
      \outline 
      - Classifier designed to classify pairs of TGs based on their differences of features. 
        As input it takes training set, training labels, test set, test labels, training numbers and test numbers.
        Successfulness of a classification is save into txt file, where faulse positive and faulse negative errors are.
        Also total errors for each k are also saved in extra txt file with suffix _Errors.
        Also TG numbers of outliers are save in extra txt file with suffix _Outliers.
  **/  
  void kNNclassifier()
  {

    for(int grp = 1; grp <=2; ++grp)
    {
      string grps;
      if(grp==0) grps.assign("_group0");
      if(grp==1) grps.assign("_group1");
      if(grp==2) grps.assign("_group2");

      string classificationResultSuffix(grps+"ONLY11FEATURES");
      string classificationResultsDirectory("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\");
      string kNNclassificationResultsFileName (classificationResultsDirectory + "kNNResults" + classificationResultSuffix);
      CImgList<double> featureSet(getTrainingAndTestSetsLabelsAndTGnumbersForAGivenGroupOfFeatures(grp,0,0,0));
      CImg<double> trainingSet(featureSet(0)), trainingLabels(featureSet(1)),
                  testSet(featureSet(2)), testLabels(featureSet(3)),trainingSetTGnumbers(featureSet(4)),testSetTGnumbers(featureSet(5));
      //trainingSet.columns(6,16);
      //testSet.columns(6,16);
      featureSet.display();

      const int K = 100; // max k allowed
      int numberOfTrainingData = trainingSet.height();
      int numberOfDataDimensions = trainingSet.width();
      int numberOfTestData = testSet.height();
      cv::Mat trainData ( numberOfTrainingData, numberOfDataDimensions, CV_32FC1 );
      cv::Mat trainClasses ( numberOfTrainingData, 1, CV_32FC1 );
      for(int i = 0; i < numberOfTrainingData; ++i)
      {
        for(int j = 0; j < numberOfDataDimensions; ++j)
        {
          trainData.at<float>(i,j) = trainingSet(j,i);
        }
      }
      for(int i = 0; i < numberOfTrainingData; ++i)
      {
        trainClasses.at<float>(i,0) = trainingLabels(0,i);
      }
      CvKNearest knn( trainData, trainClasses, cv::Mat(), false, K ); // learning of the classifier, using training set and traning labels
    
      ofstream kNNClassificationResults (kNNclassificationResultsFileName + ".txt");
      ofstream kAndErrors (kNNclassificationResultsFileName + "_Errors" + ".txt");
      ofstream outliers (kNNclassificationResultsFileName + "_Outliers" + ".txt");
      kAndErrors << setw(20) << "Number k - neighbor" << setw(30) << "totalOptimisticError in %" << setw(30) << "totalPesimisticError in %" << setw(30) << "averageTotalError in %" << setw(30) << "ratting of k" << '\n';
      outliers << "Class 0 is the class of different tortoises and Class 1 is the class of same tortoises." << endl;
      for(int k = 1; k < K; k+=1)
      {
        if(k>30 && k < K-10) k+=10;
        kAndErrors << setw(20) << k;
        cout << k << endl;
        //int imgSize = 5000;
        //CImg<int> resultsOnImg(imgSize,3,1,2, 10);  // if feature space is 2d then results can be save onto an image
        //double sumOf1DfeaturesInClass0 = 0;
        //double sumOfSquaresOf1DfeaturesInClass0 = 0;
        //double countOfSamplesInClass0 = 0;
        //double sumOf1DfeaturesInClass1 = 0;
        //double sumOfSquaresOf1DfeaturesInClass1 = 0;
        //double countOfSamplesInClass1 = 0;
        //for (int i =  0; i < numberOfTrainingData; ++i) // OPTIMISTIC CLASSIFICATION - using training set
        //{
        //  cv::Mat sampleMat (1,1,CV_32FC1);
        //  sampleMat.at<float>(0,0) = trainingSet(0,i);
        //  float response = knn.find_nearest(sampleMat,k,0,0,0,0);
        //  if (response == 1) //TGs are same
        //  {
        //    sumOf1DfeaturesInClass1 +=  trainingSet(0,i);
        //    sumOfSquaresOf1DfeaturesInClass1 +=  trainingSet(0,i)*trainingSet(0,i);
        //    ++countOfSamplesInClass1;
        //  }
        //  else if (response == 0) //TGs are not same
        //  {
        //    sumOf1DfeaturesInClass0 +=  trainingSet(0,i);
        //    sumOfSquaresOf1DfeaturesInClass0 +=  trainingSet(0,i)*trainingSet(0,i);
        //    ++countOfSamplesInClass0;
        //  }
        //  else
        //  {
        //    cout << "Error in kNN classification" << endl;
        //  }
        ////}
        //double meanOfClass0 = sumOf1DfeaturesInClass0/countOfSamplesInClass0;
        //double stdOfClass0 = sqrt( ( sumOfSquaresOf1DfeaturesInClass0 - (sumOf1DfeaturesInClass0)*(sumOf1DfeaturesInClass0)/countOfSamplesInClass0 ) / (countOfSamplesInClass0-1) );
        //double meanOfClass1 = sumOf1DfeaturesInClass1/countOfSamplesInClass1;
        //double stdOfClass1 = sqrt( ( sumOfSquaresOf1DfeaturesInClass1 - (sumOf1DfeaturesInClass1)*(sumOf1DfeaturesInClass1)/countOfSamplesInClass1 ) / (countOfSamplesInClass1-1) );
        //cout << meanOfClass0 << '\t' << stdOfClass0 << '\t' << meanOfClass1 << '\t' << stdOfClass1 << '\t' << endl;
        //{
        //  for(int x = meanOfClass0*imgSize/trainingSet.max()+1; x > (meanOfClass0 - stdOfClass0)*imgSize/trainingSet.max()+1; --x)
        //  {
        //    resultsOnImg(x,0,0) =  5;
        //    resultsOnImg(x,0,1) =  5;
        //  }
        //  for(int x = meanOfClass1*imgSize/trainingSet.max()+1; x < (meanOfClass1 + stdOfClass1)*imgSize/trainingSet.max()+1; ++x)
        //  {
        //    resultsOnImg(x,0,0) = 15;
        //    resultsOnImg(x,0,1) = 15;
        //  }
        //}
        int numberOfCorrectlyClassifiedAsSameTortoises = 0;
        int numberOfCorrectlyClassifiedAsDifferentTortoises = 0;
        int numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises = 0;
        int numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises = 0;
        outliers << '\n' << k << " - nearest neigbours" << '\n';
        //outliers << setw(20) << "meanClass1" << setw(20) << "stdClass1" << setw(20) << "meanClass0" << setw(20) << "stdClass0"<< '\n';
        //outliers << setw(20) << meanOfClass1 << setw(20) << stdOfClass1 << setw(20) << meanOfClass0 << setw(20) << stdOfClass0 << '\n';
        outliers << "TrainingSet Classification" << '\n';
        outliers << setw(20) <<"TgNumber1" << setw(20) << "TgNumber2" << setw(20) << "sumOfFeatures" << setw(20) << "GroundTruthClass" << setw(20) << "ClassifiedClass" << '\n';
      
        for (int i =  0; i < numberOfTrainingData; ++i) // OPTIMISTIC CLASSIFICATION - using training set
        {
          cv::Mat sampleMat (1,numberOfDataDimensions,CV_32FC1);
          for(int j = 0; j < numberOfDataDimensions; ++j)
          {
            sampleMat.at<float>(0,j) = trainingSet(j,i);
          }
          //int x = trainingSet(0,i)*imgSize/trainingSet.max()+1;
          //int y = 1;
          //int y = trainingSet(1,i)*1000/trainingSet.max()+1;
          float response = knn.find_nearest(sampleMat,k,0,0,0,0);
          outliers << setw(20) << trainingSetTGnumbers(0,i) << setw(20) << trainingSetTGnumbers(1,i) << setw(20) << trainingSet.get_row(i).sum() << setw(20) << trainingLabels(0,i) << setw(20) << response << '\n';

          if (response == 1 && trainingLabels(0,i)==1)
          {
            ++numberOfCorrectlyClassifiedAsSameTortoises;
            //if(x < imgSize && y < imgSize)
            //{
            //  resultsOnImg(x,y,0) =  15;
            //  resultsOnImg(x,y,1) =  15;
            //}
          }
          else if (response == 0 && trainingLabels(0,i)==0)
          {
            ++numberOfCorrectlyClassifiedAsDifferentTortoises;
            // if(x < imgSize && y < imgSize)
            // {
            //   resultsOnImg(x,y,0) =  5;
            //   resultsOnImg(x,y,1) =  5;
            //}
          }
          else if (response == 1 && trainingLabels(0,i)==0) // false positive
          {
            ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises;
            //if(x < imgSize && y < imgSize)
            //{
            //  resultsOnImg(x,y,0) =  5;
            //  resultsOnImg(x,y,1) =  15;
            //}
          }
          else if (response ==0 && trainingLabels(0,i)==1) //false negative
          {
            ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises;
            //outliers << setw(20) << trainingSetTGnumbers(0,i) << setw(20) << trainingSetTGnumbers(1,i) << setw(20) << " " << setw(20) << (trainingSet(0,i)-meanOfClass0)/stdOfClass0 << '\n';
            //if(x < imgSize && y < imgSize)
            //{
            //  resultsOnImg(x,y,0) =  15;
            //  resultsOnImg(x,y,1) =  5;
            //}
          }
          else
          {
            cout << response << " " << trainingLabels(0,i) << " neco se pokazilo" << endl;
          }
        }
        //if(k==51) resultsOnImg.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\classificationImgwithoutPow2.pnm");
        kNNClassificationResults << "  " <<  k << " - NN" << '\n';
        kNNClassificationResults << "OPTIMISTIC CLASSIFICATION ERROR" << '\n';
        kNNClassificationResults << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalOptimisticSuccess in %" << '\n';
        kNNClassificationResults << left << setw(30) << numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTrainingData << setw(40) << numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTrainingData;
        kNNClassificationResults << left << setw(40) << (numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0 << '\n';
        kNNClassificationResults << left << setw(30) << "incorrectlyClassifiedAsSame" << setw(40) << "incorrectlyClassifiedAsDifferent" << setw(40) << "totalOptimisticError in %" << '\n';
        kNNClassificationResults << left << setw(30) << numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData << setw(40) << numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData;
        kNNClassificationResults << left << setw(40) << (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0 << '\n';
        double optimisticTotalError = (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTrainingData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTrainingData)*100.0;
        kAndErrors << setw(30) << optimisticTotalError;
        kNNClassificationResults << '\n';
      
        numberOfCorrectlyClassifiedAsSameTortoises = 0;
        numberOfCorrectlyClassifiedAsDifferentTortoises = 0;
        numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises = 0;
        numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises = 0;
        outliers << '\n' << "TestSet Classification" << '\n';
        outliers << setw(20) << "TgNumber1" << setw(20) << "TgNumber2" << setw(20) << "sumOfFeatures" << setw(20) << "GroundTruthClass" << setw(20) << "ClassifiedClass" << '\n';
        for (int i =  0; i < numberOfTestData; ++i) //PESSIMISTIC CLASSIFICATION - using test set
        {
          cv::Mat sampleMat (1,numberOfDataDimensions,CV_32FC1);
          for(int j = 0; j < numberOfDataDimensions; ++j)
          {
            sampleMat.at<float>(0,j) = testSet(j,i);
          }
          //int x = testSet(0,i)*imgSize/testSet.max()+1;
          //int y = 2;
          float response = knn.find_nearest(sampleMat,k,0,0,0,0);
          outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << testSet.get_row(i).sum() << setw(20) << testLabels(0,i) << setw(20) << response << '\n';
          if (response == 1 && testLabels(0,i)==1)
          {
            ++numberOfCorrectlyClassifiedAsSameTortoises;
            //if(x < imgSize && y < imgSize)
            //{
            //  resultsOnImg(x,y,0) =  15;
            //  resultsOnImg(x,y,1) =  15;
            //}
          }
          else if (response ==0 && testLabels(0,i)==0)
          {
            ++numberOfCorrectlyClassifiedAsDifferentTortoises;
            //if(x < imgSize && y < imgSize)
            //{
            //  resultsOnImg(x,y,0) =  5;
            //  resultsOnImg(x,y,1) =  5;
            //}
          }
          else if (response == 1 && testLabels(0,i)==0)
          {
            ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises;
            //outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << 1 << setw(20) << " " << '\n';
            //if(x < imgSize && y < imgSize)
            //{
            //  resultsOnImg(x,y,0) =  5;
            //  resultsOnImg(x,y,1) =  15;
            //}
          }
          else if (response == 0 && testLabels(0,i)==1)
          {
            ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises;
            //outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << " " << setw(20) << 1 << '\n';
            //if(x < imgSize && y < imgSize)
            //{
            //  resultsOnImg(x,y,0) =  15;
            //  resultsOnImg(x,y,1) =  5;
            //}
          }
          else
          {
            cout << response << " " << testLabels(0,i) << " neco se pokazilo" << endl;
          }
        }
        //if(k%10==1) resultsOnImg.display();

        kNNClassificationResults << "PESSIMISTIC CLASSIFICATION ERROR" << '\n';
        kNNClassificationResults << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalPesimisticSuccess in %" << '\n';
        kNNClassificationResults << left << setw(30) << numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTestData << setw(40) << numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTestData;
        kNNClassificationResults << left << setw(40) << (numberOfCorrectlyClassifiedAsSameTortoises*1.0/numberOfTestData + numberOfCorrectlyClassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0 << '\n';
        kNNClassificationResults << left << setw(30) << "incorrectlyClassifiedAsSame" << setw(40) << "incorrectlyClassifiedAsDifferent" << setw(40) << "totalPesimisticError in %" << '\n';
        kNNClassificationResults << left << setw(30) << numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData << setw(40) << numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData;
        kNNClassificationResults << left << setw(40) << (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0 << '\n';
        double pessimisticTotalError = (numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*1.0/numberOfTestData + numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*1.0/numberOfTestData)*100.0;
        kAndErrors << setw(30) << pessimisticTotalError <<  setw(30) << (optimisticTotalError + pessimisticTotalError)/2 <<  setw(30) << (optimisticTotalError + pessimisticTotalError)*pow(pessimisticTotalError,2.0)/2 << endl;
        kNNClassificationResults << '\n' << endl;
        outliers << endl;
      }
    }
  }

  void removeIthColumnFromIntMatrix(CImg<int> &Matrix, int indexI)
  {
    if(indexI==Matrix.width()-1)
    {
      Matrix = Matrix.get_columns(0,indexI-1);
    }
    else if(indexI==0)
    {
      Matrix = Matrix.get_columns(indexI+1,Matrix.width()-1);
    }
    else
    {
      Matrix = Matrix.get_columns(0,indexI-1).append(Matrix.get_columns(indexI+1,Matrix.width()-1),'x');
    }
  }

  void removeIthColumnFromDoubleMatrix(CImg<double> &Matrix, int indexI)
  {
    if(indexI==Matrix.width()-1)
    {
      Matrix = Matrix.get_columns(0,indexI-1);
    }
    else if(indexI==0)
    {
      Matrix = Matrix.get_columns(indexI+1,Matrix.width()-1);
    }
    else
    {
      Matrix = Matrix.get_columns(0,indexI-1).append(Matrix.get_columns(indexI+1,Matrix.width()-1),'x');
    }
  }

  void standardizationOfFeatures(CImg<double> &featureMatrixTestSet,CImg<double> &featureMatrixTrainingSet) // columns ... feature vectors, rows ... measurements of individual features
  {
      const int sizeTestSet = featureMatrixTestSet.width();
      CImg<double> featureMatrix(featureMatrixTestSet.get_append(featureMatrixTrainingSet,'x'));
      CImg<double> standardizedMatrix;
	    const int nNumOfFeatures = featureMatrix.height();
	    CImg<double> meansOfFeatures(1,nNumOfFeatures);
      CImg<double> stdsOfFeatures(1,nNumOfFeatures);
	    for(int i = 0; i < nNumOfFeatures; i++) meansOfFeatures(i) = featureMatrix.get_row(i).mean();
      for(int i = 0; i < nNumOfFeatures; i++) stdsOfFeatures(i) = sqrt(featureMatrix.get_row(i).variance());
      for(int i = 0; i < nNumOfFeatures; i++) standardizedMatrix.append((featureMatrix.get_row(i)-meansOfFeatures(i))/stdsOfFeatures(i),'y'); //standardizing the features
      featureMatrix = standardizedMatrix;
      featureMatrixTestSet = featureMatrix.get_columns(0,sizeTestSet-1);
      featureMatrixTrainingSet = featureMatrix.get_columns(sizeTestSet,featureMatrix.width()-1);
  }

  void AccuracyEvaluation(string strLoadDirectory)
  {
    string kNNclassificationResultsFileName (strLoadDirectory + "AccuracyAllResults772vs1030.txt");
    ofstream xtxFileResults(kNNclassificationResultsFileName, ios::app);
    string strFeaturesMeasuredByComputerFileName(strLoadDirectory + "FeaturesMeasuredByComputer.cimg");
    string strFeaturesMeasuredByHandFileName(strLoadDirectory + "FeaturesMeasuredByHand.cimg");
    string strNumbersOfTGsFileName(strLoadDirectory + "numbersOfTGs.cimg");
    CImg<double> featuresMeasuredByComputer(strFeaturesMeasuredByComputerFileName.c_str());
    CImg<double> featuresMeasuredByHand(strFeaturesMeasuredByHandFileName.c_str());
    CImg<int>    tgNumbers(strNumbersOfTGsFileName.c_str());
    for(int i = 852; i < tgNumbers.width(); ++i) 
    {
      featuresMeasuredByComputer(i,0) = 999;
      //featuresMeasuredByHand.columns(0,851);
      //tgNumbers.columns(0,851);
    }
    //featuresMeasuredByComputer.columns(852,tgNumbers.width()-1);
    //featuresMeasuredByHand.columns(852,tgNumbers.width()-1); 107 772
    //tgNumbers.columns(852,tgNumbers.width()-1); 0 981 982 width-1
    //featuresMeasuredByComputer.columns(0,851);
    //featuresMeasuredByHand.columns(0,851);
    //tgNumbers.columns(0,851);
    //featuresMeasuredByComputer.rows(6,16);
    //featuresMeasuredByHand.rows(6,16);

    (featuresMeasuredByComputer,featuresMeasuredByHand,tgNumbers).display();
    const int nNumOfFirstGood = tgNumbers.width();
    (featuresMeasuredByComputer,featuresMeasuredByHand,tgNumbers).display();
    CImg<int> arrayOfHits(nNumOfFirstGood,1,1,1, 0);
	  for( int tg = 0; tg <  tgNumbers.width(); ++tg )
    {
      int actualTgNumber = tgNumbers(tg);
      CImg<double> featureDifference(tgNumbers.width(),1,1,1, 0);
      if(featuresMeasuredByComputer(tg,0) < 999)
      {
        for(int ii = 0; ii < tgNumbers.width(); ++ii)
        {
          featureDifference(ii,0) =  (featuresMeasuredByHand.get_column(ii) - featuresMeasuredByComputer.get_column(tg)).abs().sum();
        }
        int kk = 0;
        while( kk < nNumOfFirstGood )
        {
          const int nCoorOfMinFeatDist = featureDifference.get_stats()(4);
          int tgNumberWithMinDifference = tgNumbers(nCoorOfMinFeatDist,0);
          if(actualTgNumber/100 == tgNumberWithMinDifference/100) 
          {
            arrayOfHits(kk)++;
            kk=nNumOfFirstGood+1;
          }
          //else if(actualTgNumber/100 == tgNumberWithMinDifference/100)
          //{
          //  arrayOfAlmostHits[kk]++;
          //  kk=nNumOfFirstGood*2;
          //}
          else
          {
            featureDifference(nCoorOfMinFeatDist,0) = 999;
          }
          ++kk;
        }
        if(kk==nNumOfFirstGood) cout << "Something wrong - no tg number match" << endl;;
      }
    }
    arrayOfHits.save((strLoadDirectory + "AccuracyAllResults772vs1030.cimg").c_str());
    //for(int i = 0; i < nNumOfFirstGood; ++i) xtxFileResults << i+1 << "." << '\t';
    //xtxFileResults << endl;
    const int totalNum = arrayOfHits.sum();
    int sum = 0;
    xtxFileResults << totalNum << " auto-measurements tested against " << nNumOfFirstGood << " hand-measurements (CEI + Bender images) "<< endl << endl;
    for(int i = 0; i < nNumOfFirstGood; ++i) 
    {
      sum += arrayOfHits[i];
      xtxFileResults << i+1 << "." << '\t' << arrayOfHits[i] << '\t'<< " " << sum << "     " << '\t' << 100.0*sum/totalNum << '\t' << '\t'<< i+1 << "." <<'\n';
    }
    xtxFileResults << endl;
  }

  CImgList<double> GetClassificationSet(const CImg<double> &features, const CImg<int> &tgNames)
  {
    CImg<double> featureDifferences, labels;
    CImg<int> pairsOfTgNames;
    for(int ii = 0; ii < features.width(); ++ii)
    {
      if(features(ii,0) < 999) // value 999 indicates, that for the given tg, the features were not successfully measured 
      {
        int step = 1;
        for(int jj = ii+1; jj < features.width(); jj+=step)
        {
          if(features(jj,0) < 999)
          {
            //cout << ii << " " << jj << endl;
            CImg<double> classLabel(1,1,1,1);
            if(tgNames(ii,0)/100==tgNames(jj,0)/100) classLabel(0,0) = 1.0;
            else { classLabel(0,0) = 0.0; step = 50;}
            featureDifferences.append( (features.get_column(ii)-features.get_column(jj)).abs() ,'x');
            labels.append(classLabel,'x');
            pairsOfTgNames.append(tgNames.get_column(ii).append(tgNames.get_column(jj),'y'),'x');
          }
        }
      }
    }
    //equalizing the proportions of the classes
    const int nNumOfClass1Features = labels.sum();
    const int nNumOfClass0Features = labels.width() - nNumOfClass1Features;
    const float stepForRemovingClass0 = nNumOfClass0Features*1.0f/(nNumOfClass0Features-nNumOfClass1Features);
    cout << nNumOfClass1Features << " " <<  nNumOfClass0Features << " " << stepForRemovingClass0 << endl;
    int nCountOfClass0 = 0;
    int boltNumber = labels.width();
    for(int kk = 0; kk < boltNumber; ++kk)
    {
      if(labels(kk,0) == 0)
      {
        nCountOfClass0++;
        //cout << nCountOfClass0/stepForRemovingClass0 << " " << nCountOfClass0 - (int)(nCountOfClass0/stepForRemovingClass0)*stepForRemovingClass0 << endl;
        if(nCountOfClass0 - (int)(nCountOfClass0/stepForRemovingClass0)*stepForRemovingClass0 < 1)
        {
          removeIthColumnFromIntMatrix(pairsOfTgNames,kk);
          removeIthColumnFromDoubleMatrix(labels,kk);
          removeIthColumnFromDoubleMatrix(featureDifferences,kk);
          boltNumber = labels.width();
          --kk;
        }
      }
    }
    return featureDifferences, labels, pairsOfTgNames;
  }

  void performPCA(CImg<double> &featureMatrixTestSet, CImg<double> &featureMatrixTrainingSet, const int nNumOfPCs, const string strPCAEigenValues)
  {
    const int sizeTestSet = featureMatrixTestSet.width();
    CImg<double> featuresMatrix(featureMatrixTestSet.get_append(featureMatrixTrainingSet,'x'));
	  const int numOfFeatures = featuresMatrix.height(); //rows = different features
	  const int numOfMeasurements = featuresMatrix.width(); //columns = different measurements of features
	  CImg<double> measurementsCovarianceMatrix ( (featuresMatrix*featuresMatrix.get_transpose())/(numOfMeasurements-1) ); //computing the covariance matrix
    //measurementsCovarianceMatrix.display();
	  CImg<double> eigenValues, eigenVectors;
	  measurementsCovarianceMatrix.symmetric_eigen(eigenValues,eigenVectors); //computing the eigen vals and vecs
    eigenValues /= eigenValues.sum(); //normalizing the eigen vals
    //eigenValues.display();
    ofstream EigenValuesStream (strPCAEigenValues,ios_base::app);
    EigenValuesStream << strPCAEigenValues << endl;
    for(int ii = 0; ii < eigenValues.height(); ++ii)
    {
      EigenValuesStream << ii << ". " << '\t' << eigenValues(ii) << '\n';
    }
    EigenValuesStream << endl << endl;
    //(eigenValues,eigenVectors,measurementsCovarianceMatrix).display();
    //ofstream eigenvaluesLog("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\eigenvalues.txt",ios_base::app);
    //eigenvaluesLog << "Group: " << whichGroups << endl;
    //for(int i = 0; i < numOfFeatures; i++) eigenvaluesLog << right << setw(20) << eigenValues(i) << endl;
    //eigenvaluesLog << endl << endl;
    //eigenValues.display_graph();
    //cout << "enter threshold" << endl;
    //float PCAthreshold = 0;
    //cin >> PCAthreshold;
    //int numberOfUsefulEigenVectors = 0;
    //for(int i = 0; i < numOfFeatures; ++i)
    //{
    //  if(eigenValues(i) > PCAthreshold) // threshold for the eigenvalues
    //  {
    //    numberOfUsefulEigenVectors = i;
    //  }
    //}
    if(nNumOfPCs>0) eigenVectors.columns(0,nNumOfPCs-1);
    //else eigenVectors.column(numberOfUsefulEigenVectors);
    featuresMatrix = (featuresMatrix.get_transpose()*eigenVectors).get_transpose();
    featureMatrixTestSet = featuresMatrix.get_columns(0,sizeTestSet-1);
    featureMatrixTrainingSet = featuresMatrix.get_columns(sizeTestSet,featuresMatrix.width()-1);
  }

  void Classification()
  {
    for(int p1 = 0; p1 <= 0; ++p1)
    {
    for(int p2 = 0; p2 <= 0; ++p2)
    {
    for(int nPCA = 0; nPCA <= 0; ++nPCA)
    {
    cout << p1 << p2 << nPCA << endl;
    const int  database           = p1; // 0 ... database A, 1 ... database B;
    const bool usingAutoMeasured  = p2; // 0 ... hand-measured features divided into Test and Training set, 1 ... auto(computer)-measured form Test set and all hand-measured form the Training set;
    const bool doStandardization  = 0;  
    const int doPCA               = nPCA; //doPCA == 1 ... 11PCs, doPCA == 2 ... 6PCs;

    string strDirectory;
    string strAutoMeasuredFeatures;
    string strHandMeasuredFeatures;
    string strTgNumbers;
    string strStandardization;
    string strPCA;

    if (database == 0)
    {
      strDirectory.assign("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGood\\");
    }
    else
    {
      strDirectory.assign("C:\\Users\\Matej\\Documents\\Zelvy\\ZelvyFotky\\TestudoGraecaBenderCizp1304pnmOnlyGoodNoYoungTo6cm\\");
    }

    strAutoMeasuredFeatures.assign(strDirectory + "FeaturesMeasuredByComputer.cimg");
    strHandMeasuredFeatures.assign(strDirectory + "FeaturesMeasuredByHand.cimg");
    strTgNumbers.assign(strDirectory + "numbersOfTGs.cimg");
    strStandardization.assign("_Standrd");
    strPCA.assign("_PCA");
    if(doPCA==1) strPCA.append("11PCs");
    else if(doPCA==2) strPCA.append("6PCs");
    string strPCAEigenValues;
    string strKNNClassificationResults;
    if(usingAutoMeasured)
    {
      strPCAEigenValues.assign(strDirectory + "AutoClassResults" + strStandardization + strPCA + "_EigenValues.txt");
      strKNNClassificationResults.assign(strDirectory + "AutoClassResults");
    }
    else
    {
      strPCAEigenValues.assign(strDirectory + "HandClassResults" + strStandardization + strPCA + "_EigenValues.txt");
      strKNNClassificationResults.assign(strDirectory + "HandClassResults");
    }
    if(doStandardization) strKNNClassificationResults.append(strStandardization);
    if(doPCA>0) strKNNClassificationResults.append(strPCA);
  
    CImgList<double> testSet, trainingSet;

    //if(usingAutoMeasured)
    //{
    //  CImg<double> featuresForTest((strAutoMeasuredFeatures).c_str()), featuresForTraining((strHandMeasuredFeatures).c_str());
    //  CImg<int>    tgNamesForTest((strTgNumbers).c_str()), tgNamesForTraining((strTgNumbers).c_str());

    //  cout << "Forming of test set started ... " << endl;
    //  testSet.assign(GetClassificationSet(featuresForTest, tgNamesForTest));
    //  cout << "Forming of training set started ... " << endl;
    //  trainingSet.assign(GetClassificationSet(featuresForTraining, tgNamesForTraining));
    //}
    //else // if there are only hand measurements
    //{
    //  CImg<double> features((strHandMeasuredFeatures).c_str()), tgNames((strTgNumbers).c_str());
    //  cout << "Forming of classification set started ... " << endl;
    //  CImgList<double> classificationSet(GetClassificationSet(features, tgNames));
    //  cout << "Dividing of the classification set into trainig and test sets ... " << endl;
    //  testSet.assign(classificationSet);
    //  trainingSet.assign(classificationSet);
    //  CImg<int> labelsTotal(classificationSet(1).threshold(0.5)); 
    //  const int nTotalNumOfMeasurements = labelsTotal.width();
    //  int nCountClass1 = 0;
    //  int nCountClass0 = 0;
    //  for(int ii = 1; ii < nTotalNumOfMeasurements; ++ii)
    //  {
    //    if(labelsTotal(ii,0)==1) 
    //    {
    //      nCountClass1++;
    //      if(nCountClass1%2)
    //      {
    //        trainingSet(1)(ii,0) = 999;
    //      }
    //      else
    //      {
    //        testSet(1)(ii,0) = 999;
    //      }
    //    }
    //    else
    //    {
    //      nCountClass0++;
    //      if(nCountClass0%2)
    //      {
    //        trainingSet(1)(ii,0) = 999;
    //      }
    //      else
    //      {
    //        testSet(1)(ii,0) = 999;
    //      }
    //    }
    //  }
    //  //(trainingSet(1)+testSet(1)).display();
    //  for(int jj = 0; jj < trainingSet(1).width(); ++jj)
    //  {
    //    if(trainingSet(1)(jj,0)==999)
    //    {
    //      removeIthColumnFromDoubleMatrix(trainingSet(0),jj);
    //      removeIthColumnFromDoubleMatrix(trainingSet(1),jj);
    //      removeIthColumnFromDoubleMatrix(trainingSet(2),jj);
    //      --jj;
    //    }
    //  }
    //  for(int jj = 0; jj < testSet(1).width(); ++jj)
    //  {
    //    if(testSet(1)(jj,0)==999)
    //    {
    //      removeIthColumnFromDoubleMatrix(testSet(0),jj);
    //      removeIthColumnFromDoubleMatrix(testSet(1),jj);
    //      removeIthColumnFromDoubleMatrix(testSet(2),jj);
    //      --jj;
    //    }
    //  }
    //  (trainingSet,testSet).display();
    //}
    //if(usingAutoMeasured)
    //{
    //  trainingSet.save( (strDirectory + "ClassifAutoTrainingSet.cimg").c_str() );
    //  testSet.save( (strDirectory + "ClassifAutoTestSet.cimg").c_str() );
    //}
    //else
    //{
    //  trainingSet.save( (strDirectory + "ClassifHandTrainingSet.cimg").c_str() );
    //  testSet.save( (strDirectory + "ClassifHandTestSet.cimg").c_str() );
    //}

    if(usingAutoMeasured)
    {
      trainingSet.assign( (strDirectory + "ClassifAutoTrainingSet.cimg").c_str() );
      testSet.assign( (strDirectory + "ClassifAutoTestSet.cimg").c_str() );
    }
    else
    {
      trainingSet.assign( (strDirectory + "ClassifHandTestSet.cimg").c_str() );
      testSet.assign( (strDirectory + "ClassifHandTrainingSet.cimg").c_str() );
    }
    //features have to be ordered in columns one next to each other, columns of features in same order with the columns of corresponding tg names

    CImg<double> featureDifferencesForTest(testSet(0)), featureDifferencesForTraining(trainingSet(0));
    CImg<float>  labelsForTest(testSet(1)), labelsForTraining(trainingSet(1));
    CImg<int>    pairsOfTgNamesForTest(testSet(2).round()), pairsOfTgNamesForTraining(trainingSet(2).round());

    //(testSet,trainingSet).display();
      cout << "Forming of test set started ... " << endl;
    if(doStandardization)
    {
      standardizationOfFeatures(featureDifferencesForTest,featureDifferencesForTraining);
      if(doPCA>0)
      {
        int nNumOfPCs = 17;
        if (doPCA==1) nNumOfPCs = 11;
        else if (doPCA==2) nNumOfPCs = 6;
        performPCA(featureDifferencesForTest,featureDifferencesForTraining,nNumOfPCs,strPCAEigenValues);
      }
    }
    //(featureDifferencesForTest, labelsForTest, pairsOfTgNamesForTest, featureDifferencesForTraining, labelsForTraining, pairsOfTgNamesForTraining).display();
    
    ///*


    featureDifferencesForTest.transpose();
    labelsForTest.transpose(); 
    pairsOfTgNamesForTest.transpose();
    featureDifferencesForTraining.transpose();
    labelsForTraining.transpose();
    pairsOfTgNamesForTraining.transpose();

    const int K = 100; // max k allowed
    const int numberOfTrainingData = featureDifferencesForTraining.height();
    const int numberOfDataDimensions = featureDifferencesForTraining.width();
    const int numberOfTestData = featureDifferencesForTest.height();
    const int numberOfTestClass1Data = labelsForTest.get_round().sum();
    const int numberOfTestClass0Data = numberOfTestData - numberOfTestClass1Data;
    const int numberOfTrainingClass1Data = labelsForTraining.get_round().sum();
    const int numberOfTrainingClass0Data = numberOfTrainingData - numberOfTrainingClass1Data;
    //cout << labelsForTest.get_round().sum() << " " << numberOfTestClass1Data << endl;
    //labelsForTest.get_round().display();
    cv::Mat trainData ( numberOfTrainingData, numberOfDataDimensions, CV_32FC1 );
    cv::Mat trainClasses ( numberOfTrainingData, 1, CV_32FC1 );
    for(int i = 0; i < numberOfTrainingData; ++i)
    {
      for(int j = 0; j < numberOfDataDimensions; ++j)
      {
        trainData.at<float>(i,j) = featureDifferencesForTraining(j,i);
      }
    }
    for(int i = 0; i < numberOfTrainingData; ++i)
    {
      trainClasses.at<float>(i,0) = labelsForTraining(0,i);
    }
    CvKNearest knn( trainData, trainClasses, cv::Mat(), false, K ); // learning of the classifier, using training set and traning labels
    
    ofstream kNNClassificationResults (strKNNClassificationResults + ".txt");
    ofstream kAndErrors (strKNNClassificationResults + "_Errors" + ".txt");
    ofstream CeiBender (strKNNClassificationResults + "_HandDatasets" + ".txt");
    //ofstream outliers (strKNNClassificationResults + "_Outliers" + ".txt");
    kAndErrors << setw(20) << "Number k - neighbor" << setw(30) << "FalPosOptimisticError" << setw(30) << "FalNegOptimisticError" << setw(30) << "FalPosPesimisticError" << setw(30) << "FalNegPesimisticError" << setw(30) << "averageTotalError in %" << '\n';
    //outliers << "Class 0 is the class of different tortoises and Class 1 is the class of same tortoises." << endl;
    for(int k = 30; k < 31; k+=1)
    {
      //if(k>30 && k < K-10) k+=10;
      kAndErrors << setw(20) << k;
      cout << k << endl;
      //int imgSize = 5000;
      //CImg<int> resultsOnImg(imgSize,3,1,2, 10);  // if feature space is 2d then results can be save onto an image
      //double sumOf1DfeaturesInClass0 = 0;
      //double sumOfSquaresOf1DfeaturesInClass0 = 0;
      //double countOfSamplesInClass0 = 0;
      //double sumOf1DfeaturesInClass1 = 0;
      //double sumOfSquaresOf1DfeaturesInClass1 = 0;
      //double countOfSamplesInClass1 = 0;
      //for (int i =  0; i < numberOfTrainingData; ++i) // OPTIMISTIC CLASSIFICATION - using training set
      //{
      //  cv::Mat sampleMat (1,1,CV_32FC1);
      //  sampleMat.at<float>(0,0) = trainingSet(0,i);
      //  float response = knn.find_nearest(sampleMat,k,0,0,0,0);
      //  if (response == 1) //TGs are same
      //  {
      //    sumOf1DfeaturesInClass1 +=  trainingSet(0,i);
      //    sumOfSquaresOf1DfeaturesInClass1 +=  trainingSet(0,i)*trainingSet(0,i);
      //    ++countOfSamplesInClass1;
      //  }
      //  else if (response == 0) //TGs are not same
      //  {
      //    sumOf1DfeaturesInClass0 +=  trainingSet(0,i);
      //    sumOfSquaresOf1DfeaturesInClass0 +=  trainingSet(0,i)*trainingSet(0,i);
      //    ++countOfSamplesInClass0;
      //  }
      //  else
      //  {
      //    cout << "Error in kNN classification" << endl;
      //  }
      ////}
      //double meanOfClass0 = sumOf1DfeaturesInClass0/countOfSamplesInClass0;
      //double stdOfClass0 = sqrt( ( sumOfSquaresOf1DfeaturesInClass0 - (sumOf1DfeaturesInClass0)*(sumOf1DfeaturesInClass0)/countOfSamplesInClass0 ) / (countOfSamplesInClass0-1) );
      //double meanOfClass1 = sumOf1DfeaturesInClass1/countOfSamplesInClass1;
      //double stdOfClass1 = sqrt( ( sumOfSquaresOf1DfeaturesInClass1 - (sumOf1DfeaturesInClass1)*(sumOf1DfeaturesInClass1)/countOfSamplesInClass1 ) / (countOfSamplesInClass1-1) );
      //cout << meanOfClass0 << '\t' << stdOfClass0 << '\t' << meanOfClass1 << '\t' << stdOfClass1 << '\t' << endl;
      //{
      //  for(int x = meanOfClass0*imgSize/trainingSet.max()+1; x > (meanOfClass0 - stdOfClass0)*imgSize/trainingSet.max()+1; --x)
      //  {
      //    resultsOnImg(x,0,0) =  5;
      //    resultsOnImg(x,0,1) =  5;
      //  }
      //  for(int x = meanOfClass1*imgSize/trainingSet.max()+1; x < (meanOfClass1 + stdOfClass1)*imgSize/trainingSet.max()+1; ++x)
      //  {
      //    resultsOnImg(x,0,0) = 15;
      //    resultsOnImg(x,0,1) = 15;
      //  }
      //}
      int numberOfCorrectlyClassifiedAsSameTortoises = 0;
      int numberOfCorrectlyClassifiedAsDifferentTortoises = 0;
      int numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises = 0;
      int numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises = 0;
      //outliers << '\n' << k << " - nearest neigbours" << '\n';
      //outliers << setw(20) << "meanClass1" << setw(20) << "stdClass1" << setw(20) << "meanClass0" << setw(20) << "stdClass0"<< '\n';
      //outliers << setw(20) << meanOfClass1 << setw(20) << stdOfClass1 << setw(20) << meanOfClass0 << setw(20) << stdOfClass0 << '\n';
      //outliers << "TrainingSet Classification" << '\n';
      //outliers << setw(20) <<"TgNumber1" << setw(20) << "TgNumber2" << setw(20) << "sumOfFeatures" << setw(20) << "GroundTruthClass" << setw(20) << "ClassifiedClass" << '\n';
      
      for (int i =  0; i < numberOfTrainingData; ++i) // OPTIMISTIC CLASSIFICATION - using training set
      {
        cv::Mat sampleMat (1,numberOfDataDimensions,CV_32FC1);
        for(int j = 0; j < numberOfDataDimensions; ++j)
        {
          sampleMat.at<float>(0,j) = featureDifferencesForTraining(j,i);
        }
        //int x = trainingSet(0,i)*imgSize/trainingSet.max()+1;
        //int y = 1;
        //int y = trainingSet(1,i)*1000/trainingSet.max()+1;
        float response = knn.find_nearest(sampleMat,k,0,0,0,0);
        //outliers << setw(20) << pairsOfTgNamesForTraining(0,i) << setw(20) << pairsOfTgNamesForTraining(1,i) << setw(20) << featureDifferencesForTraining.get_row(i).sum() << setw(20) << labelsForTraining(0,i) << setw(20) << response << '\n';

        if (response == 1 && labelsForTraining(0,i)==1)
        {
          ++numberOfCorrectlyClassifiedAsSameTortoises;

          //if(x < imgSize && y < imgSize)
          //{
          //  resultsOnImg(x,y,0) =  15;
          //  resultsOnImg(x,y,1) =  15;
          //}
        }
        else if (response == 0 && labelsForTraining(0,i)==0)
        {
          ++numberOfCorrectlyClassifiedAsDifferentTortoises;
          // if(x < imgSize && y < imgSize)
          // {
          //   resultsOnImg(x,y,0) =  5;
          //   resultsOnImg(x,y,1) =  5;
          //}
        }
        else if (response == 1 && labelsForTraining(0,i)==0) // false positive
        {
          ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises;
          //if(x < imgSize && y < imgSize)
          //{
          //  resultsOnImg(x,y,0) =  5;
          //  resultsOnImg(x,y,1) =  15;
          //}
        }
        else if (response ==0 && labelsForTraining(0,i)==1) //false negative
        {
          ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises;
          //outliers << setw(20) << trainingSetTGnumbers(0,i) << setw(20) << trainingSetTGnumbers(1,i) << setw(20) << " " << setw(20) << (trainingSet(0,i)-meanOfClass0)/stdOfClass0 << '\n';
          //if(x < imgSize && y < imgSize)
          //{
          //  resultsOnImg(x,y,0) =  15;
          //  resultsOnImg(x,y,1) =  5;
          //}
        }
        else
        {
          cout << response << " " << labelsForTraining(0,i) << " neco se pokazilo" << endl;
        }
      }
      //if(k==51) resultsOnImg.save("C:\\Users\\Matej\\Documents\\Zelvy\\Vysledky\\measuredPointsByHandOfSeamSegmentsBenderAndCizpSegmentsLenghtsAndDifferences\\classificationImgwithoutPow2.pnm");
      float fPercentageCorrectlyClassifiedAsSame      = numberOfCorrectlyClassifiedAsSameTortoises*100.0f/numberOfTrainingClass1Data;
      float fPercentageCorrectlyClassifiedAsDifferent = numberOfCorrectlyClassifiedAsDifferentTortoises*100.0f/numberOfTrainingClass0Data;
      float fPercentageTotalOptimisticSuccess         = (fPercentageCorrectlyClassifiedAsSame + fPercentageCorrectlyClassifiedAsDifferent)/2.0f;
      float fPercentageIncorrectClassifAsDifferent    = numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*100.0f/numberOfTrainingClass1Data;
      float fPercentageIncorrectClassifAsSame         = numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*100.0f/numberOfTrainingClass0Data;
      float fPercentageTotalOptimisticError           = (fPercentageIncorrectClassifAsDifferent + fPercentageIncorrectClassifAsSame)/2.0f;

      kNNClassificationResults << "  " <<  k << " - NN" << '\n';
      kNNClassificationResults << "OPTIMISTIC CLASSIFICATION ERROR" << '\n';
      kNNClassificationResults << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalOptimisticSuccess in %" << '\n';
      kNNClassificationResults << left << setw(30) << fPercentageCorrectlyClassifiedAsSame << setw(40) << fPercentageCorrectlyClassifiedAsDifferent << setw(40) << fPercentageTotalOptimisticSuccess << '\n';
      kNNClassificationResults << left << setw(30) << "incorrectClassifAsDifferent" << setw(40) << "incorrectClassifAsSame" << setw(40) << "totalOptimisticError in %" << '\n';
      kNNClassificationResults << left << setw(30) << fPercentageIncorrectClassifAsDifferent << setw(40) << fPercentageIncorrectClassifAsSame << setw(40) << fPercentageTotalOptimisticError << '\n';
      kAndErrors << setw(30) << fPercentageIncorrectClassifAsSame << setw(30) << fPercentageIncorrectClassifAsDifferent;
      
      numberOfCorrectlyClassifiedAsSameTortoises = 0;
      numberOfCorrectlyClassifiedAsDifferentTortoises = 0;
      numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises = 0;
      numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises = 0;
      int numberOfCorrectlyClassifiedAsSameTortoisesOneDataset = 0;
      int numberOfCorrectlyClassifiedAsDifferentTortoisesOneDataset = 0;
      int numberOfFalsePositiveErrorsMisclassifiedAsSameTortoisesOneDataset = 0;
      int numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoisesOneDataset = 0;
      int totalClass0OneDataset = 0;
      int totalClass1OneDataset = 0;
      //outliers << '\n' << "TestSet Classification" << '\n';
      //outliers << setw(20) << "TgNumber1" << setw(20) << "TgNumber2" << setw(20) << "sumOfFeatures" << setw(20) << "GroundTruthClass" << setw(20) << "ClassifiedClass" << '\n';
      for (int i =  0; i < numberOfTestData; ++i) //PESSIMISTIC CLASSIFICATION - using test set
      {
        cv::Mat sampleMat (1,numberOfDataDimensions,CV_32FC1);
        for(int j = 0; j < numberOfDataDimensions; ++j)
        {
          sampleMat.at<float>(0,j) = featureDifferencesForTest(j,i);
        }
        //int x = testSet(0,i)*imgSize/testSet.max()+1;
        //int y = 2;
        float response = knn.find_nearest(sampleMat,k,0,0,0,0);          
        if (pairsOfTgNamesForTest(0,i)/100 < 300 || pairsOfTgNamesForTest(1,i)/100 < 300 ) 
        {
           if(labelsForTest(0,i)==1)
           {
             ++totalClass1OneDataset;
           }
           else
           {
             ++totalClass0OneDataset;
           }
        }
        //outliers << setw(20) << pairsOfTgNamesForTest(0,i) << setw(20) << pairsOfTgNamesForTest(1,i) << setw(20) << featureDifferencesForTest.get_row(i).sum() << setw(20) << labelsForTest(0,i) << setw(20) << response << '\n';

        if (response == 1 && labelsForTest(0,i)==1)
        {
          ++numberOfCorrectlyClassifiedAsSameTortoises;
          if (pairsOfTgNamesForTest(0,i)/100 < 300 || pairsOfTgNamesForTest(1,i)/100 < 300 ) 
          {
            ++numberOfCorrectlyClassifiedAsSameTortoisesOneDataset;
          }
          //if(x < imgSize && y < imgSize)
          //{
          //  resultsOnImg(x,y,0) =  15;
          //  resultsOnImg(x,y,1) =  15;
          //}
        }
        else if (response ==0 && labelsForTest(0,i)==0)
        {
          ++numberOfCorrectlyClassifiedAsDifferentTortoises;
          if (pairsOfTgNamesForTest(0,i)/100 < 300 || pairsOfTgNamesForTest(1,i)/100 < 300 ) 
          {
            ++numberOfCorrectlyClassifiedAsDifferentTortoisesOneDataset;
          }
          //if(x < imgSize && y < imgSize)
          //{
          //  resultsOnImg(x,y,0) =  5;
          //  resultsOnImg(x,y,1) =  5;
          //}
        }
        else if (response == 1 && labelsForTest(0,i)==0)
        {
          ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises;
          if (pairsOfTgNamesForTest(0,i)/100 < 300 || pairsOfTgNamesForTest(1,i)/100 < 300 ) 
          {
            ++numberOfFalsePositiveErrorsMisclassifiedAsSameTortoisesOneDataset;
          }
          //outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << 1 << setw(20) << " " << '\n';
          //if(x < imgSize && y < imgSize)
          //{
          //  resultsOnImg(x,y,0) =  5;
          //  resultsOnImg(x,y,1) =  15;
          //}
        }
        else if (response == 0 && labelsForTest(0,i)==1)
        {
          ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises;
          if (pairsOfTgNamesForTest(0,i)/100 < 300 || pairsOfTgNamesForTest(1,i)/100 < 300 ) 
          {
            ++numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoisesOneDataset;
          }
          //outliers << setw(20) << testSetTGnumbers(0,i) << setw(20) << testSetTGnumbers(1,i) << setw(20) << " " << setw(20) << 1 << '\n';
          //if(x < imgSize && y < imgSize)
          //{
          //  resultsOnImg(x,y,0) =  15;
          //  resultsOnImg(x,y,1) =  5;
          //}
        }
        else
        {
          cout << response << " " << labelsForTest(0,i) << " neco se pokazilo" << endl;
        }
      }
      //if(k%10==1) resultsOnImg.display();
      float fPercentageTotalOptimisticErrorOptimistic = fPercentageTotalOptimisticError;
      fPercentageCorrectlyClassifiedAsSame      = numberOfCorrectlyClassifiedAsSameTortoises*100.0f/numberOfTestClass1Data;
      fPercentageCorrectlyClassifiedAsDifferent = numberOfCorrectlyClassifiedAsDifferentTortoises*100.0f/numberOfTestClass0Data;
      fPercentageTotalOptimisticSuccess         = (fPercentageCorrectlyClassifiedAsSame + fPercentageCorrectlyClassifiedAsDifferent)/2.0f;
      fPercentageIncorrectClassifAsDifferent    = numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoises*100.0f/numberOfTestClass1Data;
      fPercentageIncorrectClassifAsSame         = numberOfFalsePositiveErrorsMisclassifiedAsSameTortoises*100.0f/numberOfTestClass0Data;
      fPercentageTotalOptimisticError           = (fPercentageIncorrectClassifAsDifferent + fPercentageIncorrectClassifAsSame)/2.0f;

      float fPercentageCorrectlyClassifiedAsSameOneDataset      = numberOfCorrectlyClassifiedAsSameTortoisesOneDataset*100.0f/totalClass1OneDataset;
      float fPercentageCorrectlyClassifiedAsDifferentOneDataset = numberOfCorrectlyClassifiedAsDifferentTortoisesOneDataset*100.0f/totalClass0OneDataset;
      float fPercentageIncorrectClassifAsDifferentOneDataset    = numberOfFalseNegativeErrorsMisclassifiedAsDifferentTortoisesOneDataset*100.0f/totalClass1OneDataset;
      float fPercentageIncorrectClassifAsSameOneDataset         = numberOfFalsePositiveErrorsMisclassifiedAsSameTortoisesOneDataset*100.0f/totalClass0OneDataset;

      kNNClassificationResults << "  " <<  k << " - NN" << '\n';
      kNNClassificationResults << "PESSIMISTIC CLASSIFICATION ERROR" << '\n';
      kNNClassificationResults << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalOptimisticSuccess in %" << '\n';
      kNNClassificationResults << left << setw(30) << fPercentageCorrectlyClassifiedAsSame << setw(40) << fPercentageCorrectlyClassifiedAsDifferent << setw(40) << fPercentageTotalOptimisticSuccess << '\n';
      kNNClassificationResults << left << setw(30) << "incorrectClassifAsDifferent" << setw(40) << "incorrectClassifAsSame" << setw(40) << "totalOptimisticError in %" << '\n';
      kNNClassificationResults << left << setw(30) << fPercentageIncorrectClassifAsDifferent << setw(40) << fPercentageIncorrectClassifAsSame << setw(40) << fPercentageTotalOptimisticError << '\n';
      kAndErrors << setw(30) << fPercentageIncorrectClassifAsSame <<  setw(30) << fPercentageIncorrectClassifAsDifferent << setw(30) << (fPercentageTotalOptimisticErrorOptimistic + fPercentageTotalOptimisticError)/2 << endl;
      kNNClassificationResults << '\n' << endl;

      CeiBender << "  " <<  k << " - NN" << '\n';
      CeiBender << "PESSIMISTIC CLASSIFICATION ERROR" << '\n';
      CeiBender << left << setw(30) << "correctlyClassifiedAsSame" << setw(40) << "correctlyClassifiedAsDifferent" << setw(40) << "totalOptimisticSuccess in %" << '\n';
      CeiBender << left << setw(30) << fPercentageCorrectlyClassifiedAsSameOneDataset << setw(40) << fPercentageCorrectlyClassifiedAsDifferentOneDataset << '\n';
      CeiBender << left << setw(30) << "incorrectClassifAsDifferent" << setw(40) << "incorrectClassifAsSame" << setw(40) << "totalOptimisticError in %" << '\n';  
      CeiBender << left << setw(30) << fPercentageIncorrectClassifAsDifferentOneDataset << setw(40) << fPercentageIncorrectClassifAsSameOneDataset << '\n';
      CeiBender << '\n' << endl;
      CeiBender << "Number of Class1 test features one set" << '\t' << "Number of Class0 test features one set" << endl;
      CeiBender << numberOfTestClass1Data <<  '\t' << numberOfTestClass0Data << endl;
      CeiBender << "Number of Class1 test features" << '\t' << "Number of Class0 test features" << endl;
      CeiBender << totalClass1OneDataset <<  '\t' << totalClass0OneDataset<< endl;
      CeiBender << "Number of Class1 training features" << '\t' << "Number of Class0 traning features" << endl;
      CeiBender << numberOfTrainingClass1Data <<  '\t' << numberOfTrainingClass0Data<< endl;
      //outliers << endl;
    }
    //*/
    }
    }
    }
  }


