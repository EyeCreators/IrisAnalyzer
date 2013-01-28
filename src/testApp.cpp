#include "testApp.h"

//--------------------------------------------------------------
void testApp::exit(){
	gui->saveSettings("GUI/guiSettings.xml");
    delete gui; 
}

//--------------------------------------------------------------
void testApp::guiEvent(ofxUIEventArgs &e)
{
    if(e.widget->getName() == "SLIDER1"){
        ofxUISlider *slider = (ofxUISlider *) e.widget;    
        sliderFloat01 = slider->getScaledValue();
    }
    else if(e.widget->getName() == "FULLSCREEN")
    {
        ofxUIToggle *toggle = (ofxUIToggle *) e.widget;
        ofSetFullscreen(toggle->getValue());
    }    
}


//--------------------------------------------------------------
void testApp::setup(){
    propertiesXML.loadFile("properties.xml");
    bExportedInformaticsFile = false;
    ofSetVerticalSync(false);
    ofSetFrameRate(999);
    
    bShowThings = true;
    int showThingsVal = propertiesXML.getValue("SHOW_THINGS", 1);
    bShowThings = (bool)(showThingsVal > 0);
    
    bReportDeadLabels = true;
    int deadLabelsVal = propertiesXML.getValue("REPORT_DEAD_LABELS", 1);
    bReportDeadLabels = (bool)(deadLabelsVal > 0);
    
	bVerbose = true;
    bExitWhenFinished = true;
    bAllocatedSomeImages = false;
    bIsTheEyeThere = false;
    bDrawHalfScale = true;
    
	pupilContourAnalyzer = new PupilContourAnalyzer(); 
	irisContourAnalyzer  = new PupilContourAnalyzer();
	
	ofEnableSmoothing();
	gui = new ofxUICanvas(600,100,320,320);
	gui->addWidgetDown(new ofxUILabel("IRIS", OFX_UI_FONT_LARGE)); 
	
	//(string _name, float _min, float _max, float _value, float w, float h, float x = 0, float y = 0)
	//ofxUISlider(string _name, float _min, float _max, float _value, float w, float h, float x = 0, float y = 0)
	gui->addWidgetDown(new ofxUISlider("SLIDER1", 0.0, 1.0, 0.5, 300,16,0,0));
	gui->addWidgetDown(new ofxUIToggle(32, 32, false, "FULLSCREEN"));
	ofAddListener(gui->newGUIEvent, this, &testApp::guiEvent); 
	gui->loadSettings("GUI/guiSettings.xml");
    gui->toggleVisible();
    
    unwarper = new threesixtyUnwarp();
	
	
	ipp8u_incoming			= ippiMalloc_8u_C3 (PROCESS_W, PROCESS_H_MAX, &processS3);	// storage for incoming frame
	ipp8u_Y 				= ippiMalloc_8u_C1 (PROCESS_W, PROCESS_H_MAX, &processS1);	// Y of incoming
	ipp8u_U 				= ippiMalloc_8u_C1 (PROCESS_W, PROCESS_H_MAX, &processS1);	// U of incoming
	ipp8u_V 				= ippiMalloc_8u_C1 (PROCESS_W, PROCESS_H_MAX, &processS1);	// V of incoming
	ipp8u_Y2				= ippiMalloc_8u_C1 (PROCESS_W, PROCESS_H_MAX, &processS1);	// temp;
	ipp8u_HandyC1			= ippiMalloc_8u_C1 (PROCESS_W, PROCESS_H_MAX, &processS1);	// temp;
	
	iplImage_ipp8uY			= cvCreateImage(cvSize(PROCESS_W, PROCESS_H_MAX), 8,1);
	iplImage_ipp8uY2		= cvCreateImage(cvSize(PROCESS_W, PROCESS_H_MAX), 8,1);
	
	histEqer				= new tmVisThresholderC1 (PROCESS_W, PROCESS_H_MAX);
	pupilThresholder		= new tmVisThresholderC1 (PROCESS_W, PROCESS_H_MAX);
    


    int unwarpedW = unwarper->unwarpedW;
    int unwarpedH = unwarper->unwarpedH;
    int unwarpedS = unwarpedW;
    featureSearchW = unwarpedW;
    featureSearchH = (int)(unwarpedH * 0.625);
    // we only examine the top two-thirds of the iris image
    
    unwarpedImageOfxCvC3.allocate(unwarpedW, unwarpedH);
    unwarpedImageOfxCvC1.allocate(unwarpedW, unwarpedH);
    unwarpedImageHistEqC1.allocate(unwarpedW, unwarpedH);
    ipp8uUnwarpedHistEqC1 = ippiMalloc_8u_C1 (unwarpedW, unwarpedH, &unwarpedS);
    unwarpEqer = new tmVisThresholderC1 (unwarpedW, unwarpedH);
    HistEqC1Pix.allocate(featureSearchW, featureSearchH, 1); 
    surfFeatureAccumImage.allocate(featureSearchW, featureSearchH);
    

    
    pupilPrevPx = pupilPx = processW/2;
    pupilPrevPy = pupilPy = processH/2;
    maxNSurfTrackerLabels = 0;
    maxNGfttTrackerLabels = 0; 
    
    /*
    IppiSize processRoi = {processW, processH};
    ippiSet_8u_C1R(0, ipp8u_Y, processS1, processRoi);
    pupilSearchImageC1.setFromPixels(ipp8u_Y, processW, processH);
    pupilSearchImagePrevC1.setFromPixels(ipp8u_Y, processW, processH);
    pupilSearchImageDiffC1.setFromPixels(ipp8u_Y, processW, processH);
     */
    
    
    
    int     GFTT_maxCorners = 8;
    double  GFTT_qualityLevel = 0.01;
    double  GFTT_minDistance = 8.0;
    int     GFTT_blockSize = 9;
    bool    GFTT_useHarrisDetector = true; // true
    double  GFTT_k = 0.04;
    goodFeatureTracker = new GoodFeaturesToTrackDetector (GFTT_maxCorners,
                                                          GFTT_qualityLevel,
                                                          GFTT_minDistance,
                                                          GFTT_blockSize,
                                                          GFTT_useHarrisDetector,
                                                          GFTT_k);
    SFD = new cv::SurfFeatureDetector();
    nSurfPointsToKeep   = 40;
    nSurfPointsToReport = 16;
    surfFeatureAccumImage.set(0);
    bComputeFeatures = true;
    
    // These will change (get larger, slightly ) when we include the Tracker. 
    maxXmlSurfFeatures = nSurfPointsToReport;
    maxXmlGfttFeatures = GFTT_maxCorners;
    
    
    
    bUseMovie = true;
    if (bUseMovie){
        loadEyeMovie();
        getNextFrameImage();
    } else {
        loadImage();
    }
    
}


//--------------------------------------------------------------
void testApp::loadEyeMovie(){
	
	lastTagNumber   = 0;
	analysisScaleRatio   = 1.0;
	
    

    string videoFilenameToLoad = propertiesXML.getValue("VIDEO_FILENAME_TO_LOAD", "default.mov");
	eyeMovie.loadMovie(videoFilenameToLoad);
	eyeMovie.play();
	eyeMovie.setPaused(true);
    
	analysisScaleRatio = (float)eyeMovie.getWidth() / (float)PROCESS_W;
    
    
    // Resize this frame to ensure width is PROCESS_W
	int oldw = eyeMovie.getWidth();
	int oldh = eyeMovie.getHeight();
	float oldAspect = (float)oldh/(float)oldw;
	processW = PROCESS_W; //64*(oldw/64);
	processH = (int)(processW * oldAspect);
    
    unwarper->setup(processW, processH);
    pupilSearchImageC1.allocate    (processW, processH);
    pupilSearchImagePrevC1.allocate(processW, processH);
    pupilSearchImageDiffC1.allocate(processW, processH);
    
    bRecordingXML   = true;
    lastTagNumber	= XML.addTag("RECORDING");
}



//--------------------------------------------------------------
void testApp::getNextFrameImage(){
	
	ofImage	inputImage;
	eyeMovie.update();
	inputImage.setFromPixels(eyeMovie.getPixelsRef()); 

	// Resize this frame to ensure width is PROCESS_W
    /*
	int oldw = inputImage.getWidth();
	int oldh = inputImage.getHeight();
	float oldAspect = (float)oldh/(float)oldw;
	processW = PROCESS_W; //64*(oldw/64);
	processH = (int)(processW * oldAspect);
    */
	inputImage.resize(processW, processH);
    
    
	if (!bAllocatedSomeImages){
        pupilSearchMaskPix.allocate(processW, processH, 1);
        pupilSearchMaskC1.allocate(processW, processH);
        pupilSearchImageThreshedC1.allocate(processW, processH);
        
        irisSearchMaskPix.allocate(processW, processH, 1);
        irisSearchMaskC1.allocate(processW, processH);
    }
    bAllocatedSomeImages = true;
	
	rawImageC3.setFromPixels(inputImage.getPixelsRef());
	
	int nFrames = eyeMovie.getTotalNumFrames() -1;
	if (eyeMovie.getCurrentFrame() < nFrames){
		eyeMovie.nextFrame();
        
	} else {
        int nAnalyzedFrames = eyeMovie.getCurrentFrame();
		eyeMovie.setFrame(0);
        
        if (bExportedInformaticsFile == false){
            
            if (XML.pushTag("RECORDING", lastTagNumber) ){
                int fcTagnum = XML.addTag("FEATURE_COUNTS");
                XML.setValue("FEATURE_COUNTS:MAX_N_PUPILS", 1, fcTagnum);
                XML.setValue("FEATURE_COUNTS:MAX_N_IRISES", 1, fcTagnum);
                XML.setValue("FEATURE_COUNTS:MAX_N_HIGHLIGHTS", maxNHighlights, fcTagnum);
                
                XML.setValue("FEATURE_COUNTS:MAX_N_S_FEATURES", maxNSurfTrackerLabels /*maxXmlSurfFeatures*/, fcTagnum);
                XML.setValue("FEATURE_COUNTS:MAX_N_G_FEATURES", maxNGfttTrackerLabels /*maxXmlGfttFeatures*/, fcTagnum);

                XML.popTag();
            }
            
            string videoFilenameToLoad = propertiesXML.getValue("VIDEO_FILENAME_TO_LOAD", "default.mov");
            sprintf(handyString, "%s_Informatics.xml", videoFilenameToLoad.c_str());
            XML.saveFile(handyString);
            
            bRecordingXML = false;
            bExportedInformaticsFile = true;
            printf("Completed Recording with nFrames = %d\n", nAnalyzedFrames);
            printf("Data Recording saved to xml: %s\n", handyString);
            
            
            int exitWhenFinishedVal = propertiesXML.getValue("EXIT_WHEN_FINISHED", 0);
            bExitWhenFinished = (bool)(exitWhenFinishedVal > 0);
            if (bExitWhenFinished){
                ofExit(0); 
            }
        }
	}
	
    
    
   
}



//--------------------------------------------------------------
void testApp::loadImage(){
	
	lastTagNumber        = 0;
	bRecordingXML        = false;
	analysisScaleRatio   = 1.0; 
	
	ofImage	inputImage;
	
	int r = (int) ofRandom(0, 8.9); 
	sprintf(handyString, "images/eye-%d.jpg", r); 
	inputImage.loadImage(handyString); 

	// Resize image to ensure width is a multiple of 64 pixels. 
	// This is helpful later with certain image processing libraries, like IPP.
	int oldw = inputImage.getWidth();
	int oldh = inputImage.getHeight();
	float oldAspect = (float)oldh/(float)oldw;
	processW = PROCESS_W; //64*(oldw/64);
	processH = (int)(processW * oldAspect);
	inputImage.resize(processW, processH); 
	analysisScaleRatio = (float)oldw / (float)processW;
	
	pupilSearchMaskPix.allocate(processW, processH, 1);
	pupilSearchMaskC1.allocate(processW, processH);
	pupilSearchImageThreshedC1.allocate(processW, processH);
    irisSearchMaskPix.allocate(processW, processH, 1);
    irisSearchMaskC1.allocate(processW, processH);

    
    bAllocatedSomeImages = true;
	
	rawImageC3.setFromPixels(inputImage.getPixelsRef());
}


//--------------------------------------------------------------
void testApp::findHighlights(){
    
    // Create another grayscale copy of the raw image. 
    highlightSearchImageC1.setFromColorImage(rawImageC3);
    
    // assume the mask spot has already been created.
    // Subtract it, creating an image which is dark on the outside. 
    highlightSearchImageC1 -= pupilSearchMaskC1;
    
    // Apply a nonlinear LUT to the highlightSearchImageC1
    unsigned char *highlightSearchPix = highlightSearchImageC1.getPixels();
    for (int y=0; y<processH; y++){
        for (int x=0; x<processW; x++){
            int ind = y*processW + x;
            float val = ofMap(highlightSearchPix[ind], 0,255, 0,1, true);
            if (val < 0.75){ val = powf(val, 2.0); }
            highlightSearchPix[ind] = (unsigned char) ((int)ofMap(val,0,1, 0,255));
        }
    }
    
    // Threshold the highlightSearchImageC1 to binarize the highlights.
    // We'll use an empirically-determined threshold (0.70), instead of automatic.
    float highlightThreshold01 = 0.60;
    int highlightThreshold = ofMap(highlightThreshold01, 0,1, 0,255); //sliderFloat01
    highlightSearchImageC1.threshold(highlightThreshold);
    
    //Apply the contour finders to locate the highlights.
    int maxHighlightArea = 0.00400 * processW * processH;
	int minHighlightArea = 0.00005 * processW * processH;
    highlightContourFinder.findContours (highlightSearchImageC1, minHighlightArea,maxHighlightArea, 2, false, false);
    maxNHighlights = max(maxNHighlights, highlightContourFinder.nBlobs); 
    
}

//--------------------------------------------------------------
void testApp::process (){
	
    bIsTheEyeThere = true;  // initially, we believe it's there.
	
	// printf("sliderFloat01 = %f\n", sliderFloat01); 
	
	//--------------------------
	// Take the RGB rawImageC3, and histogram-equalize it. 
	histogramNormalizeC3(); 
	
	
	//--------------------------
	// CREATE SEARCH MASK SPOT
	// Assume the pupil is a dark spot, close to the center of the image. 
	// Generate a probability mask, with the same dimensions as our input image, 
	// which priviliges the center of the image as a search place for the pupil. 
	int pupilSearchMaskPixN = pupilSearchMaskPix.size(); 
	int whalf = processW/2; 
	int hhalf = processH/2;
	float farDist = (whalf+hhalf)/2;   // distance to what will be the far edge.  
	float nearDist = 0.30 * min(whalf,hhalf);
	for (int y=0; y<processH; y++){
		for (int x=0; x<processW; x++){
			
			float dx = (x - whalf)*1.0; 
			float dy = (y - hhalf)*1.5; 
			float dh = sqrtf(dx*dx + dy*dy); 
			dh = ofMap(dh, nearDist,farDist, 0,1, true); 
			dh = powf(dh, 2.0); 
			dh = ofMap(dh,0,1, 0,255, true); 
			char c = (char)(dh); 
			
			int i = (y*processW) + x;
			pupilSearchMaskPix[i] = c;
		}
	}
	pupilSearchMaskC1.setFromPixels(pupilSearchMaskPix.getPixels(), processW,processH);
    
   
	//--------------------------
	// CREATE IMAGE IN WHICH TO SEARCH FOR PUPIL
	// Create a grayscale pupilSearchImageC1, within which we will search for the pupil. 
	// This is the histogram equalized luminance of the color image, plus the mask.
    pupilSearchImagePrevC1 = pupilSearchImageC1;
	pupilSearchImageC1.setFromPixels(ipp8u_Y, processW, processH);
	pupilSearchImageC1 += pupilSearchMaskC1;
    
    //--------------------------
    // DO FRAME DIFFERENCING to detect blinking.
    // If the frame difference is too high -- it's a blink or something.
    pupilSearchImageDiffC1.absDiff(pupilSearchImagePrevC1, pupilSearchImageC1);
    IplImage*  diffCvImg = pupilSearchImageDiffC1.getCvImage();
    double diffSum = cvSumPixels(diffCvImg);
    float differenceFactor = (float) diffSum / (float)(processW * processH);
    //printf("differenceFactor = %f\n", differenceFactor);
    if ((differenceFactor < 0.001) || (differenceFactor > 3.0)){ // empirical
        bIsTheEyeThere = false;
    }
	
	//--------------------------
	// COMPUTE OPTIMAL THRESHOLD FOR PUPIL
	// Compute a threshold, just above "black" (pupil-color), which we will use to threshold. 
	pupilThreshold = (int)(255 * 0.25);
	float pupilThresholdLo01 = 0.06;
	float pupilThresholdHi01 = 0.34;
	int nThresholdSweepSteps = (int)(255.0 * (pupilThresholdHi01 - pupilThresholdLo01));
	
	float leastCompactness = 999; 
	int idealPupilThreshold = pupilThreshold;
	
	bool bInvertThreshedPupilSearchImage = true;
	int maxPupilArea = 0.160 * processW * processH;
	int minPupilArea = 0.001 * processW * processH;
	
	for (int s=0; s<nThresholdSweepSteps; s++){
		pupilThreshold = (int) round(255.0 * ofMap(s, 0,nThresholdSweepSteps-1, pupilThresholdLo01,pupilThresholdHi01)); 
		
		// Threshold the pupilSearchImageC1 so as to find the pupil blob. 
		// Obtain the blob with the cv contour finder, and obtain its centroid. 
		pupilSearchImageThreshedC1.setFromPixels(pupilSearchImageC1.getPixels(), processW,processH); 
		pupilSearchImageThreshedC1.threshold(pupilThreshold, bInvertThreshedPupilSearchImage);
		pupilContourFinder.findContours(pupilSearchImageThreshedC1, minPupilArea, maxPupilArea, 1, false, false);
		
		// Compute the compactness of the obtained blob. 
		// Store the pupilThreshold which produces the most compact pupil blob. 
		if (pupilContourFinder.nBlobs == 1){
			ofxCvBlob pupilBlob = pupilContourFinder.blobs.at(0);
			int nPupilBlobPts = pupilBlob.nPts;
			if (nPupilBlobPts > 1){ // pathology check
				
				// compute compactness of blob.
				float blobArea      = pupilBlob.area;
				float blobPerimeter = pupilBlob.length;
				float compactness = ((blobPerimeter*blobPerimeter)/blobArea) / (4.0 * PI);
				
				if (compactness < leastCompactness){
					leastCompactness = compactness; 
					idealPupilThreshold = pupilThreshold;
				}
			}
		}
	}
	
	// now that we know the idealPupilThreshold, use it once more: 
	// Threshold the pupilSearchImageC1 so as to find the pupil blob. 
	// Obtain the blob with the cv contour finder, and obtain its centroid. 
	pupilThreshold = idealPupilThreshold;
	pupilSearchImageThreshedC1.setFromPixels(pupilSearchImageC1.getPixels(), processW,processH); 
	pupilSearchImageThreshedC1.threshold(pupilThreshold, bInvertThreshedPupilSearchImage);
	pupilContourFinder.findContours(pupilSearchImageThreshedC1, minPupilArea, maxPupilArea, 1, false, false);
	
	
	//--------------------------
	// FIT PUPIL WITH RANSAC ELLIPSE
	// Use ransac ellipse fitting to obtain the pupil location
    pupilPrevPx = pupilPx;
    pupilPrevPy = pupilPy;
	irisPx = pupilPx = processW/2;
	irisPy = pupilPy = processH/2;
	pupilPr = processW/30;
    irisPr = pupilPr * 3;// a guess
    
    // while we're at it, reckon the pupil motion as a way of determining eye presence
    float pupilDx = pupilPrevPx - pupilPx;
    float pupilDy = pupilPrevPy - pupilPy;
    float pupilDh = sqrtf(pupilDx*pupilDx + pupilDy*pupilDy);
    float pupilMovementFactor = 1000.0 * pupilDh / processW;
    // printf("pupilMovementFactor = %f\n", pupilMovementFactor);
    if ((pupilMovementFactor < 0.001) || (pupilMovementFactor > 100.0)){
        bIsTheEyeThere = false;
    }
    
	
	if (pupilContourFinder.nBlobs == 1){
		ofxCvBlob pupilBlob = pupilContourFinder.blobs.at(0);
		int nPupilBlobPts = pupilBlob.nPts;
		if (nPupilBlobPts > 1){ // pathology check
			pupilContourAnalyzer->fitEllipseToProvidedBlob(pupilBlob);
            ofVec3f pupilReport = pupilContourAnalyzer->getCircleReport();
            pupilPx = pupilReport.x;
            pupilPy = pupilReport.y;
            pupilPr = pupilReport.z;
        }
    }
    
		
	

	
}

//--------------------------------------------------------------
void testApp::addDataToXMLInformaticReport(){
    if (bRecordingXML){
        if (XML.pushTag("RECORDING", lastTagNumber) ){
            
            currentFrameTag = XML.addTag("FRAME");
            XML.pushTag("FRAME", currentFrameTag);
            
            float currentMovieTime = ofGetElapsedTimef();
            if (bUseMovie){
                currentMovieTime = eyeMovie.getPosition()* eyeMovie.getDuration();
            }
            XML.setValue("TIME", currentMovieTime, currentFrameTag);
            XML.setValue("EYE_THERE", bIsTheEyeThere, currentFrameTag);
            
            int tagNum = XML.addTag("PUPIL");
            XML.setValue("PUPIL:ID", "Pupil-0", tagNum);
            XML.setValue("PUPIL:X", pupilPx * analysisScaleRatio, tagNum);
            XML.setValue("PUPIL:Y", pupilPy * analysisScaleRatio, tagNum);
            XML.setValue("PUPIL:R", pupilPr * analysisScaleRatio, tagNum);
            // printf("%f %f %f\n", pupilPx, pupilPy, pupilPr);
            
            int irisTagNum = XML.addTag("IRIS");
            XML.setValue("IRIS:ID", "Iris-0", tagNum);
            XML.setValue("IRIS:X", irisPx    * analysisScaleRatio, irisTagNum); // using irisPx now instead of pupilPx, OKWITCHU?
            XML.setValue("IRIS:Y", irisPy    * analysisScaleRatio, irisTagNum);
            XML.setValue("IRIS:R", irisRadius * analysisScaleRatio, irisTagNum);
            
            
            // add highlights to XML
            int nHighlights = highlightContourFinder.nBlobs;
            if (nHighlights > 0){
                for (int i=0; i<nHighlights; i++){
                    ofxCvBlob highlightBlob = highlightContourFinder.blobs.at(i);
                    int highlightTagNum = XML.addTag("HIGHLIGHT");
                    sprintf(handyString, "Highlight-%d", i);
                    float hiPx = highlightBlob.centroid.x;
                    float hiPy = highlightBlob.centroid.y;
                    float hiPw = highlightBlob.boundingRect.getWidth();
                    float hiPh = highlightBlob.boundingRect.getHeight();
                    float hiPr = max(hiPw, hiPh);
                    
                    // compute color at the 4 corners
                    ofVec3f c0 = getColorAt(hiPx - hiPw, hiPy - hiPh, rawImageC3);
                    ofVec3f c1 = getColorAt(hiPx + hiPw, hiPy - hiPh, rawImageC3);
                    ofVec3f c2 = getColorAt(hiPx + hiPw, hiPy + hiPh, rawImageC3);
                    ofVec3f c3 = getColorAt(hiPx - hiPw, hiPy + hiPh, rawImageC3);
                    float hiR = (c0.x + c1.x + c2.x + c3.x)/4.0;
                    float hiG = (c0.y + c1.y + c2.y + c3.y)/4.0;
                    float hiB = (c0.z + c1.z + c2.z + c3.z)/4.0;
                    
                    XML.setValue("HIGHLIGHT:ID", handyString, highlightTagNum);
                    XML.setValue("HIGHLIGHT:X",   hiPx * analysisScaleRatio, highlightTagNum);
                    XML.setValue("HIGHLIGHT:Y",   hiPy * analysisScaleRatio, highlightTagNum);
                    XML.setValue("HIGHLIGHT:R",   hiPr * analysisScaleRatio, highlightTagNum);
                    
                    XML.setValue("HIGHLIGHT:RED", hiR, highlightTagNum);
                    XML.setValue("HIGHLIGHT:GRN", hiG, highlightTagNum);
                    XML.setValue("HIGHLIGHT:BLU", hiB, highlightTagNum);
                }
            }
            
            bool bUseOldNonTrackingMethodForExportingFeatures = false;
            if (bUseOldNonTrackingMethodForExportingFeatures){
        
                // Export the SURF Features.  // OLD OLD OLD, NOT USED
                int nCartesianSurfPts = cartesianSurfFeatureCvPoints.size();
                for (int i=0; i<nCartesianSurfPts; i++){
                    Point2f pt = cartesianSurfFeatureCvPoints[i];
                    int surfFeatureTagNum = XML.addTag("S_FEATURE");
                    sprintf(handyString, "S-Feature-%d", i);
                    XML.setValue("S_FEATURE:ID",  handyString, surfFeatureTagNum);
                    XML.setValue("S_FEATURE:X",   pt.x * analysisScaleRatio, surfFeatureTagNum);
                    XML.setValue("S_FEATURE:Y",   pt.y * analysisScaleRatio, surfFeatureTagNum);
                }
                maxNSurfTrackerLabels = max(maxNSurfTrackerLabels, nCartesianSurfPts);
                
            } else {
             
                //NEW 
                vector<unsigned int> trackedSPointLabels = surfTracker.getCurrentLabels();
                int nTrackedSPointLabels = trackedSPointLabels.size();
                int sCount = 0; 
                for (int i=0; i<nTrackedSPointLabels; i++){
                    unsigned int aLabel = trackedSPointLabels[i];
                    cv::Point2f pt = surfTracker.getCurrent(aLabel); // a tracked point
                    int surfFeatureTagNum = XML.addTag("S_FEATURE");
                    sprintf(handyString, "S-Feature-%d", aLabel);
                    XML.setValue("S_FEATURE:ID",  handyString, surfFeatureTagNum);
                    XML.setValue("S_FEATURE:X",   pt.x * analysisScaleRatio, surfFeatureTagNum);
                    XML.setValue("S_FEATURE:Y",   pt.y * analysisScaleRatio, surfFeatureTagNum);
                    sCount++;
                }
                if (bReportDeadLabels){
                    const vector<unsigned int>& deadLabels = surfTracker.getDeadLabels();
                    for(int i = 0; i < deadLabels.size(); i++) {
                        int aDeadLabel = deadLabels[i];
                        cv::Point2f pt = surfTracker.getPrevious(aDeadLabel); // a tracked point
                        int surfFeatureTagNum = XML.addTag("S_FEATURE");
                        sprintf(handyString, "S-Feature-%d", aDeadLabel);
                        XML.setValue("S_FEATURE:ID",  handyString, surfFeatureTagNum);
                        XML.setValue("S_FEATURE:X",   pt.x * analysisScaleRatio, surfFeatureTagNum);
                        XML.setValue("S_FEATURE:Y",   pt.y * analysisScaleRatio, surfFeatureTagNum);
                        sCount++;
                    }
                }
                if (sCount > maxNSurfTrackerLabels){
                    maxNSurfTrackerLabels = sCount;
                }
            }
            
            
            if (bUseOldNonTrackingMethodForExportingFeatures){
           
                // Export the GFTT Features. OLD OLD OLD NOT USED!
                int nCartesianGfttPts = cartesianGfttFeatureCvPoints.size();
                for (int i=0; i<nCartesianGfttPts; i++){
                    Point2f pt = cartesianGfttFeatureCvPoints[i];
                    int gfttFeatureTagNum = XML.addTag("G_FEATURE");
                    sprintf(handyString, "G-Feature-%d", i);
                    XML.setValue("G_FEATURE:ID",  handyString, gfttFeatureTagNum);
                    XML.setValue("G_FEATURE:X",   pt.x * analysisScaleRatio, gfttFeatureTagNum);
                    XML.setValue("G_FEATURE:Y",   pt.y * analysisScaleRatio, gfttFeatureTagNum);
                }
                maxNGfttTrackerLabels = max (maxNGfttTrackerLabels, nCartesianGfttPts);
            
            } else {
            
                vector<unsigned int> trackedGPointLabels = gfttTracker.getCurrentLabels();
                int nTrackedGPointLabels = trackedGPointLabels.size();
                int gCount = 0; 
                for (int i=0; i<nTrackedGPointLabels; i++){
                    unsigned int aLabel = trackedGPointLabels[i];
                    cv::Point2f pt = gfttTracker.getCurrent(aLabel); // a tracked point
                    int gfttFeatureTagNum = XML.addTag("G_FEATURE");
                    sprintf(handyString, "G-Feature-%d", aLabel);
                    XML.setValue("G_FEATURE:ID",  handyString, gfttFeatureTagNum);
                    XML.setValue("G_FEATURE:X",   pt.x * analysisScaleRatio, gfttFeatureTagNum);
                    XML.setValue("G_FEATURE:Y",   pt.y * analysisScaleRatio, gfttFeatureTagNum);
                    gCount++; 
                }
                if (bReportDeadLabels){
                    const vector<unsigned int>& deadLabels = gfttTracker.getDeadLabels();
                    for(int i = 0; i < deadLabels.size(); i++) {
                        int aDeadLabel = deadLabels[i];
                        cv::Point2f pt = gfttTracker.getPrevious(aDeadLabel); // a tracked point
                        int gfttFeatureTagNum = XML.addTag("G_FEATURE");
                        sprintf(handyString, "G-Feature-%d", aDeadLabel);
                        XML.setValue("G_FEATURE:ID",  handyString, gfttFeatureTagNum);
                        XML.setValue("G_FEATURE:X",   pt.x * analysisScaleRatio, gfttFeatureTagNum);
                        XML.setValue("G_FEATURE:Y",   pt.y * analysisScaleRatio, gfttFeatureTagNum);
                        gCount++; 
                    }
                }
                if (gCount > maxNGfttTrackerLabels){
                    maxNGfttTrackerLabels = gCount;
                }
           
            }
            
            
            XML.popTag(); // frame
            XML.popTag(); // recording
        }
    }

}

//--------------------------------------------------------------
ofVec3f testApp::getColorAt (float x, float y, ofxCvColorImage colImg){
    
    int px = (int)round(ofClamp(x, 0,processW-1));
    int py = (int)round(ofClamp(y, 0,processH-1));
    int index = (py * processW * 3) + (px * 3);
    unsigned char *colImgPix = colImg.getPixels();
    
    int r = (int) colImgPix[index  ];
    int g = (int) colImgPix[index+1];
    int b = (int) colImgPix[index+2];
    
    ofVec3f out;
    out.x = r;
    out.y = g;
    out.z = b;
    
    // printf(" - color %d %d // %f %f %f\n",  px, py, out.x, out.y, out.z);
    return out;
}


//--------------------------------------------------------------
void testApp::clobberHighlights(){
    
    // add highlights to XML
    int nHighlights = highlightContourFinder.nBlobs;
    for (int i=0; i<nHighlights; i++){
        ofxCvBlob highlightBlob = highlightContourFinder.blobs.at(i);
        sprintf(handyString, "Highlight-%d", i);
        float hiPx = highlightBlob.centroid.x;
        float hiPy = highlightBlob.centroid.y;
        float hiPw = highlightBlob.boundingRect.getWidth();
        float hiPh = highlightBlob.boundingRect.getHeight();
        float hiPr = max(hiPw, hiPh);
        
        // compute color at the 4 corners
        ofVec3f c0 = getColorAt(hiPx - hiPw, hiPy - hiPh, rawImageC3);
        ofVec3f c1 = getColorAt(hiPx + hiPw, hiPy - hiPh, rawImageC3);
        ofVec3f c2 = getColorAt(hiPx + hiPw, hiPy + hiPh, rawImageC3);
        ofVec3f c3 = getColorAt(hiPx - hiPw, hiPy + hiPh, rawImageC3);
        float hiR = (c0.x + c1.x + c2.x + c3.x)/4.0;
        float hiG = (c0.y + c1.y + c2.y + c3.y)/4.0;
        float hiB = (c0.z + c1.z + c2.z + c3.z)/4.0;
        
        
        int T = (int) (hiPy - hiPh*1.25);
        int B = (int) (hiPy + hiPh*1.25);
        int L = (int) (hiPx - hiPw*1.25);
        int R = (int) (hiPx + hiPw*1.25);
        T = (int)ofClamp(T, 0, processH-1);
        B = (int)ofClamp(B, 0, processH-1);
        L = (int)ofClamp(L, 0, processW-1);
        R = (int)ofClamp(R, 0, processW-1);
        
        int cx = (int) hiPx;
        int cy = (int) hiPy; 
        
        float hiPr07 = hiPr*0.7;
        float hiPr12 = hiPr*1.2;
        
        unsigned char *colImgPix = rawImageC3.getPixels();
        for (int y=T; y<B; y++){
            for (int x=L; x<R; x++){
                
                int dx = x - cx;
                int dy = y - cy;
                float dh = sqrtf(dx*dx + dy*dy);
                int index = (y * processW * 3) + (x * 3);
                
                if (dh < hiPr07){
                    colImgPix[index  ] = hiR;
                    colImgPix[index+1] = hiG;
                    colImgPix[index+2] = hiB;
                } else if (dh < hiPr12){
                    float A = ofMap(dh, hiPr07, hiPr12, 1.0, 0.0);
                    float B = 1.0-A;
                    
                    colImgPix[index  ] = (int)(B*colImgPix[index  ] + A*hiR);
                    colImgPix[index+1] = (int)(B*colImgPix[index+1] + A*hiG);
                    colImgPix[index+2] = (int)(B*colImgPix[index+2] + A*hiB);
                }
            
            }
        }
        
    }

}


//--------------------------------------------------------------
void testApp::analyzeIris(){
  
    //--------------------------
	// CREATE IRIS SEARCH MASK SPOT
	// Assume the iris is centered on the pupil.
	// Generate a probability mask, with the same dimensions as our input image,
	// which priviliges the pupil as the center from which to search for the iris.
	int irisSearchMaskPixN = irisSearchMaskPix.size();
    
    float farDist  = processH * 0.8; 
	float nearDist = pupilPr * 5.0;
    float farDist2  = farDist * farDist;
    float nearDist2 = nearDist * nearDist;
    float irisPow = 2.0;
    
    int procRow; 
	for (int y=0; y<processH; y++){
        procRow = (y*processW);
		for (int x=0; x<processW; x++){
			
			float dx = (x - pupilPx);
			float dy = (y - pupilPy);
			float dh2 = (dx*dx + dy*dy*2.25);
            int i = procRow + x;
			
            if (dh2 > farDist2){
                irisSearchMaskPix[i] = (char) 255;
            } else if (dh2 < nearDist2) {
                irisSearchMaskPix[i] = (char) 0;
            } else {
                float dh = ofMap(dh2, nearDist2,farDist2, 0,1, true);
                dh = dh*dh;    //powf(dh, irisPow); same as for irisPow = 2.0
                dh = dh*255.0; //ofMap(dh,0,1, 0,255, true);
                irisSearchMaskPix[i] = (char)dh;
            }
		}
	}
	irisSearchMaskC1.setFromPixels(irisSearchMaskPix.getPixels(), processW,processH);
	
	//--------------------------
	// CREATE IMAGE IN WHICH TO SEARCH FOR IRIS
	// Create a grayscale irisSearchImageC1, within which we will search for the iris.
	// This is the histogram equalized luminance of the color image, plus the iris mask.
	irisSearchImageC1.setFromPixels(ipp8u_Y, processW, processH);
	irisSearchImageC1 += irisSearchMaskC1;
    
    
    float maxSearchR = min(pupilPr * 7.0, processH * 0.45); // empirical ?
    float minSearchR = pupilPr * 2.0;
    irisContourAnalyzer->findIrisEllipse(irisSearchImageC1.getPixels(),
                                          processW, processH, processS1,
                                          pupilPx, pupilPy,
                                          minSearchR, maxSearchR, 0.10);
    
    //--------------------------
	// ESTIMATE IRIS RADIUS; ASSUME CONCENTRIC WITH PUPIL.
    // Warning: Various faulty assumptions here: only circular irises, concentric with pupil, etc. 
	unsigned char *warpedColorPixels = rawImageC3.getPixels();
    irisPr = irisRadius = irisContourAnalyzer->getMedianRadius();
    ofVec3f irisReport  = irisContourAnalyzer->getCircleReport();
    irisPx = irisReport.x;
    irisPy = irisReport.y;
    
    
    //--------------------------
    // UNWARP POLAR TO CARTESIAN
    // Here is where we make the assumption that the iris and the pupil are concentric,
    // but only as far as iris analysis is concerned. 
    unwarper->update ( pupilPx, pupilPy, pupilPr, irisRadius, warpedColorPixels );
    
    unsigned char *unwarpedPixels = unwarper->unwarpedPixels;
    int unwarpedW = unwarper->unwarpedW;
    int unwarpedH = unwarper->unwarpedH;
    unwarpedImageOfxCvC3.setFromPixels(unwarper->unwarpedPixels, unwarpedW, unwarpedH);
    unwarpedImageOfxCvC1.setFromColorImage(unwarpedImageOfxCvC3);
    
    
    //--------------------------
    // histogram equalize.
    unwarpEqer->setActiveROI (unwarpedW, unwarpedH);
    int offsetBytes = 0;
    int heqType = 0;
    float heqAmount = 1.0; 
    
    //IppiSize blackRoi = {unwarpedW, unwarpedH};
    //ippiSet_8u_C1R(0, ipp8uUnwarpedHistEqC1, unwarpedW, blackRoi);
    unsigned char *src = unwarpedImageOfxCvC1.getPixels() + offsetBytes;
    unsigned char *dst = ipp8uUnwarpedHistEqC1 + offsetBytes;
    unwarpEqer->modImageByHistogram(src, dst, heqType, heqAmount);
    unwarpedImageHistEqC1.setFromPixels(ipp8uUnwarpedHistEqC1, unwarpedW, unwarpedH);

}

//--------------------------------------------------------------
void testApp::findFeatures(){
    
    if (bIsTheEyeThere == false){
        
        gfttKeyPoints.clear();
        surfKeyPoints.clear();
        bestSurfKeyPoints.clear();
        keepSurfKeyPoints.clear();
        cartesianSurfFeatureCvPoints.clear();
        cartesianGfttFeatureCvPoints.clear();
        
    } else {
    
        //--------------------------
        // Engage Feature Detectors!
        
        // Find features using GoodFeaturesToTrack. 
        HistEqC1Pix.setFromPixels(ipp8uUnwarpedHistEqC1, featureSearchW, featureSearchH, 1);
        cv::Mat cvImageInWhichToFindFeatures = ofxCv::toCv(HistEqC1Pix);
        gfttKeyPoints.clear();
        goodFeatureTracker->detect(cvImageInWhichToFindFeatures, gfttKeyPoints);
        sort (gfttKeyPoints.begin(), gfttKeyPoints.end(), key_point_x_compare());
        
        // Unwarp the features discovered in polar space. 
        cartesianGfttFeatureCvPoints.clear();
        int nGfttKeyPoints = gfttKeyPoints.size();
        for (int i=0; i<nGfttKeyPoints; i++){
            KeyPoint kp = gfttKeyPoints[i];
            Point2f  pt = kp.pt;
            float    ux  = pt.x;
            float    uy  = pt.y;
            /*
            ofVec2f *warpedVec2f = unwarper->getWarpedXYFromUnwarpedXY (ux,uy);
            float wx = warpedVec2f->x;
            float wy = warpedVec2f->y;
             cartesianGfttFeatureCvPoints.push_back( cv::Point2f(wx,wy) );
             */
            cv::Point2f aCvPoint = unwarper->getPt2fWarpedXYFromUnwarpedXY (ux,uy);
            cartesianGfttFeatureCvPoints.push_back( aCvPoint );
        }
        
        
        // Compute the latest surf points into surfKeyPoints
        surfKeyPoints.clear();
        SFD->detect(cvImageInWhichToFindFeatures, surfKeyPoints);
        
        int nNewSurfKeyPoints = surfKeyPoints.size();
        if (nNewSurfKeyPoints > 0){
            sort (surfKeyPoints.begin(), surfKeyPoints.end(), key_point_response_compare());

            int nFeaturesToAccum = min(nNewSurfKeyPoints, nSurfPointsToKeep);
            int startIndex = nNewSurfKeyPoints - nFeaturesToAccum;
            int endIndex   = nNewSurfKeyPoints;
            float maxAllowableFeatureSize = 0.08 * featureSearchW;
            float maxResponse = surfKeyPoints[endIndex-1].response; 
            
            
            // knock down the buffer, then blur it.
            surfFeatureAccumImage.blurHeavily();
            surfFeatureAccumImage -= 6;
            
            // now add in the spots. For every Surf feature found,
            // Add a blurry spot into the accumulator image. 
            Ipp8u *accumPix = surfFeatureAccumImage.getPixels();
            for (int i=startIndex; i<endIndex; i++){
                KeyPoint kpNew = surfKeyPoints[i];
                
                // To start, discard features that are too large. 
                float newDiam  = kpNew.size;
                float newr = newDiam / 2.0;
                float newr2 = newr*newr;
                
                if (newDiam < maxAllowableFeatureSize){
                    Point2f  ptNew = kpNew.pt;
                    float  kpNewResponse = kpNew.response;
                    int cx  = (int)round(ptNew.x);
                    int cy  = (int)round(ptNew.y);
                    
                    int T = (int) (cy - newr);
                    int B = (int) (cy + newr);
                    int L = (int) (cx - newr);
                    int R = (int) (cx + newr);
                    T = (int)ofClamp(T, 0, featureSearchH-1);
                    B = (int)ofClamp(B, 0, featureSearchH-1);
                    L = (int)ofClamp(L, 0, featureSearchW-1);
                    R = (int)ofClamp(R, 0, featureSearchW-1);
                    
                    float diamFrac = powf(1.0 - min(newDiam, maxAllowableFeatureSize)/maxAllowableFeatureSize, 2.0);
                    float responseFrac = ofMap(kpNewResponse, 0,maxResponse, 0,1);
                    responseFrac = ofClamp(responseFrac, 0,1);
                    responseFrac = powf(responseFrac, 2.0); 
                    
                    for (int y=T; y<B; y++){
                        int featureImgRow = (y * featureSearchW);
                        for (int x=L; x<R; x++){
                            
                            int dx = x - cx;
                            int dy = y - cy;
                            float dh2 = (dx*dx + dy*dy);
                            int index = featureImgRow + x;
                            
                            if (dh2 < newr2){
                                float dh = sqrtf(dh2);
                                
                                // smaller ones should get a bigger add
                                float maxToAdd = 1.0 + (int)(diamFrac * 3.0); // ratio of 3:1 for diam effect
                                maxToAdd *= responseFrac;
                                maxToAdd *= 60;
                                int toAdd = (int) ofMap(dh, 0,newr, maxToAdd,0);
                                accumPix[index] = min(accumPix[index]+toAdd, 255);
                            } 
                        }
                    }
                    
                
                }
            }
            surfFeatureAccumImage.setFromPixels(accumPix, featureSearchW, featureSearchH);
            
            
            // Weight the "response" of (all) the SURF features by the brightness of the corresponding pixel in the accum image.
            // This provides a metric of temporal relevance to each of the features.
            // Then re-sort by the (weighted) response, and select the best ones to keep. 
            
            int maxAccumIndex = (featureSearchW * featureSearchH) -1;
            for (int i=0; i<nNewSurfKeyPoints; i++){
                KeyPoint kpNew = surfKeyPoints[i];
                
                // Once again, discard features that are too large.
                float newDiam  = kpNew.size;
                if (newDiam > maxAllowableFeatureSize){
                    surfKeyPoints[i].response = 0;
                    
                // If it's not too large, weight its "response" by
                //the brightness of its corresponding pixel in the accum.
                } else {
                    Point2f  ptNew = kpNew.pt;
                    int cx  = (int)round(ptNew.x);
                    int cy  = (int)round(ptNew.y);
                    int index = (cy * featureSearchW) + cx;
                    index = (int) ofClamp(index, 0, maxAccumIndex);
                    float brightness = ((float)(accumPix[index])/255.0);
                    if (brightness < 0.1) brightness = 0; // hack threshold
                    kpNew.response *= brightness;
                }
            }
            // Re-sort the surfKeyPoints with the newly-weighted response's. 
            sort (surfKeyPoints.begin(), surfKeyPoints.end(), key_point_response_compare());
            
            // Put the nSurfPointsToKeep best ones, worth keeping, into keepSurfKeyPoints.
            // But, only keep ones, where it is the case that they don't overlap a previously kept one.
            // To ensure that we keep the best ones, work backwards from the top of the sorted list. 
            keepSurfKeyPoints.clear();
            nFeaturesToAccum = min(nNewSurfKeyPoints, nSurfPointsToKeep);
            startIndex   = nNewSurfKeyPoints - 1;
            endIndex     = max(0, nNewSurfKeyPoints - nFeaturesToAccum);
            for (int i=startIndex; i>=endIndex; i--){
                KeyPoint kpNew = surfKeyPoints[i];
                Point2f  ptNew = kpNew.pt;
                float    newx  = ptNew.x;
                float    newy  = ptNew.y;
                
                bool bFilterAlreadyKept = true;
                if (!bFilterAlreadyKept){ // the stupid case
                    keepSurfKeyPoints.push_back(kpNew);
                    
                } else if (bFilterAlreadyKept){
                    int nAlreadyKept = keepSurfKeyPoints.size();
                    if (nAlreadyKept > 0){
                        bool bFoundAmongAlreadyKept = false;
                        for (int k=0; k<nAlreadyKept; k++){
                            KeyPoint kpKept = keepSurfKeyPoints[k];
                            Point2f  ptKept = kpKept.pt;
                            float    keptx = ptKept.x;
                            float    kepty = ptKept.y;
                            float    keptr = kpKept.size;

                            float dx = newx - keptx;
                            float dy = newy - kepty;
                            float dh = sqrtf(dx*dx + dy*dy);
                            if (dh < 6.0){ // pixels
                                bFoundAmongAlreadyKept = true;
                            }
                        }
                        if (!bFoundAmongAlreadyKept){
                            keepSurfKeyPoints.push_back(kpNew);
                        }
                        
                    } else {
                        keepSurfKeyPoints.push_back(kpNew);
                    }
                }
            }
            
            
            
            // Choose the best nSurfPointsToReport (16) of the descending-order sorted keepers.
            // Stash them in another array, then sort that by x position.
            bestSurfKeyPoints.clear();
            int nBest = min(nSurfPointsToReport, (int)keepSurfKeyPoints.size()); //keepSurfKeyPoints.size();//
            for (int i=0; i<nBest; i++){
                KeyPoint kpKept = keepSurfKeyPoints[i];
                bestSurfKeyPoints.push_back(kpKept);           
            }
            sort (bestSurfKeyPoints.begin(), bestSurfKeyPoints.end(), key_point_x_compare());
            
            // Convert unwarped coords to warped coords. 
            cartesianSurfFeatureCvPoints.clear();
            int nBestSurfKeyPoints = bestSurfKeyPoints.size();
            for (int i=0; i<nBestSurfKeyPoints; i++){
                KeyPoint kp = bestSurfKeyPoints[i];
                Point2f  pt = kp.pt;
                float    ux  = pt.x;
                float    uy  = pt.y;
                ofVec2f *warpedVec2f = unwarper->getWarpedXYFromUnwarpedXY (ux,uy);
                float wx = warpedVec2f->x;
                float wy = warpedVec2f->y;
                cartesianSurfFeatureCvPoints.push_back( cv::Point2f(wx,wy) );
            }
            
            //-------------------------------------
            // OFXCV POINT TRACKER
            // wait for half a second before forgetting something
            surfTracker.setPersistence(3);
            gfttTracker.setPersistence(3);
            
            // an object can move up to some number of pixels per frame
            surfTracker.setMaximumDistance(16);
            gfttTracker.setMaximumDistance(16);
            
            // Update the feature trackers!
            surfTracker.track(cartesianSurfFeatureCvPoints);
            gfttTracker.track(cartesianGfttFeatureCvPoints);
            
            
            
        }
    
    
    }
    
}

//--------------------------------------------------------------
void testApp::update(){
    
    
	process();
    findHighlights();
    clobberHighlights(); 
    analyzeIris();
    
    if (bComputeFeatures){
        findFeatures();
    }
    
    addDataToXMLInformaticReport();
    
    if (bUseMovie){
        getNextFrameImage();
    }
}

//--------------------------------------------------------------
void testApp::draw(){
	ofBackground(50,50,50);
    
    if (bShowThings){
        ofPushMatrix();
        
        float drawRatio = 1.0; 
        if (drawRatio > 1){
            drawRatio *= ((float) ofGetWidth() / (float)processW);
        }
        if (bDrawHalfScale){
            drawRatio *= 0.5;
        }
        ofScale(drawRatio, drawRatio, 1.0);

        
        
        glColor3f(1,1,1); 
        histEqImageC3.draw(0,0); 
        pupilContourFinder.draw(0,0);
        highlightContourFinder.draw(0,0);
        
        glColor3f(1,1,1);
        pupilSearchImageDiffC1.draw(processW, 0);
        irisSearchImageC1.draw (0,processH*1);
        irisContourAnalyzer->draw(0,0);
        
        ofFill();
        ofEnableAlphaBlending();
        ofSetColor(255,0,0, 128); 
        ofEllipse(pupilPx, pupilPy, pupilPr*2, pupilPr*2);
        
        ofNoFill();
        glLineWidth(2.0); 
        ofSetColor(0,200,255);
        ofEllipse(irisPx, irisPy, irisPr*2, irisPr*2);
        glLineWidth(1.0);
        
        
        // The unwarp histogram equalizer's histogram gets drawn, cuz it looks cool.
        unwarpEqer->renderHistogram(processW+20, processH, 1024,processH, false, 0);
        

        if (bComputeFeatures){
            
            //--------------------
            
            bool bDrawTrackedSurfFeatures = true; // yep
            if (bDrawTrackedSurfFeatures){
                
                vector<unsigned int> trackedPointLabels = surfTracker.getCurrentLabels();
                int nTrackedPointLabels = trackedPointLabels.size();
                int sCount = 0; 
                for (int i=0; i<nTrackedPointLabels; i++){
                    unsigned int aLabel = trackedPointLabels[i];
                    cv::Point2f pt = surfTracker.getCurrent(aLabel); // a tracked point
                    
                    ofFill();
                    ofSetColor(153,0,255, 80);
                    ofEllipse(pt.x,pt.y, 19,19);
                    ofSetColor(255,255,255);
                    ofEllipse(pt.x,pt.y, 5,5);
                    
                    string msg = ofToString(aLabel);// + ":" + ofToString(surfTracker.getAge(aLabel));
                    ofDrawBitmapString(msg, pt.x+12, pt.y);
                }
                const vector<unsigned int>& deadLabels = surfTracker.getDeadLabels();
                for(int i = 0; i < deadLabels.size(); i++) {
                    int aDeadLabel = deadLabels[i];
                    cv::Point2f pt = surfTracker.getPrevious(aDeadLabel); // a tracked point
                    
                    ofFill();
                    ofSetColor(153,0,255, 80);
                    ofEllipse(pt.x,pt.y, 19,19);
                    ofSetColor(255,255,255);
                    ofEllipse(pt.x,pt.y, 5,5);
                    
                    string msg = ofToString(aDeadLabel);// + ":" + ofToString(surfTracker.getAge(aLabel));
                    ofDrawBitmapString(msg, pt.x+12, pt.y);
                }
                
            
            } else {
            
                // draw cartesian SURF feature points!
                // Doesn't use tracker -- OLD OLD OLD
                int nCartesianSurfPts = cartesianSurfFeatureCvPoints.size();
                for (int i=0; i< nCartesianSurfPts; i++){
                    Point2f pt = cartesianSurfFeatureCvPoints[i];
                    ofFill();
                    ofSetColor(153,0,255, 80);
                    ofEllipse(pt.x,pt.y, 19,19);
                }
                for (int i=0; i< nCartesianSurfPts; i++){
                    Point2f pt = cartesianSurfFeatureCvPoints[i];
                    ofFill();
                    ofSetColor(255,255,255);
                    ofEllipse(pt.x,pt.y, 5,5);
                    //sprintf(handyString, "%d", i);
                    //ofDrawBitmapString(handyString, pt.x+12, pt.y);
                }
            }
            
            //--------------------
            // draw cartesian GFTT feature points!
            
            bool bDrawTrackedGfttFeatures = true;
            if (bDrawTrackedGfttFeatures){
                
                vector<unsigned int> trackedPointLabels = gfttTracker.getCurrentLabels();
                int nTrackedPointLabels = trackedPointLabels.size();
                for (int i=0; i<nTrackedPointLabels; i++){
                    unsigned int aLabel = trackedPointLabels[i];
                    cv::Point2f pt = gfttTracker.getCurrent(aLabel); // a tracked point
                    
                    ofFill();
                    ofSetColor(255,153,0, 80);
                    ofEllipse(pt.x,pt.y, 19,19);
                    ofSetColor(255,255,255);
                    ofEllipse(pt.x,pt.y, 5,5);

                    //string msg = "G:" + ofToString(aLabel);
                    //ofDrawBitmapString(msg, pt.x+12, pt.y);
                }
                
                const vector<unsigned int>& deadLabels = gfttTracker.getDeadLabels();
                for(int i = 0; i < deadLabels.size(); i++) {
                    int aDeadLabel = deadLabels[i];
                    cv::Point2f pt = gfttTracker.getPrevious(aDeadLabel); // a tracked point
                    
                    ofFill();
                    ofSetColor(153,0,255, 80);
                    ofEllipse(pt.x,pt.y, 19,19);
                    ofSetColor(255,255,255);
                    ofEllipse(pt.x,pt.y, 5,5);
                    
                    string msg = ofToString(aDeadLabel);
                    ofDrawBitmapString(msg, pt.x+12, pt.y);
                }
                
                
            } else {
            
                ofFill();
                int nCartesianGfttPts = cartesianGfttFeatureCvPoints.size();
                for (int i=0; i< nCartesianGfttPts; i++){
                    Point2f pt = cartesianGfttFeatureCvPoints[i];
                    ofSetColor(255,153,0, 80);
                    ofEllipse(pt.x,pt.y, 19,19);
                }
                for (int i=0; i< nCartesianGfttPts; i++){
                    Point2f pt = cartesianGfttFeatureCvPoints[i];
                    ofSetColor(255,250,80);
                    ofEllipse(pt.x,pt.y, 5,5);
                }
                
            }
            
            
            // Draw the unwarped iris images. 
            float unwarperDrawX = 0;
            float unwarperDrawY = (processH*2);
            ofSetColor(255,255,255);
            int unwarpedW = unwarper->unwarpedW;
            int unwarpedH = unwarper->unwarpedH;
            unwarper->draw (unwarperDrawX,unwarperDrawY);
            unwarpedImageHistEqC1.draw(unwarperDrawX,             unwarperDrawY+ unwarpedH);
            unwarpedImageOfxCvC1.draw (unwarperDrawX + unwarpedW, unwarperDrawY+ unwarpedH);
            surfFeatureAccumImage.draw(unwarperDrawX + unwarpedW, unwarperDrawY);
            
            //--------------------
            // Draw indicators for the iris features on the unwarped images. 
            ofNoFill();
            ofPushMatrix();
            ofTranslate( unwarperDrawX, unwarperDrawY);
            
            // http://docs.opencv.org/modules/features2d/doc/common_interfaces_of_feature_detectors.html
            ofSetColor(255,200,0);
            int nGfttKeyPoints = gfttKeyPoints.size();
            for (int i=0; i<nGfttKeyPoints; i++){
                KeyPoint kp = gfttKeyPoints[i];
                Point2f pt = kp.pt;
                ofEllipse(pt.x, pt.y, 7, 7);
            }
            
            // draw bestSurfKeyPoints
            int endi = bestSurfKeyPoints.size();
            for (int i=0; i<endi; i++){
                KeyPoint kpBest = bestSurfKeyPoints[i];
                Point2f pt      = kpBest.pt;
                float   ptSize  = kpBest.size;
                ofSetColor(255*(pt.x/1024.0),255,100, 100);
                ofEllipse(pt.x, pt.y, ptSize, ptSize);
            }
            
            for (int i=0; i<endi; i++){
                KeyPoint kpBest = bestSurfKeyPoints[i];
                Point2f pt      = kpBest.pt;
                float   ptSize  = kpBest.size;
                ofSetColor(255*(pt.x/1024.0),255,100, 100);
                ofEllipse(unwarpedW + pt.x, pt.y, ptSize, ptSize);
            }

            ofPopMatrix();
        }
        
        ofPopMatrix();
    }
	
	// DRAW THE FRAME/RECORDING INDICATOR. 
	if (bRecordingXML){
        float indicatorW = 180;
        
        ofFill();
        if (bIsTheEyeThere){
            ofSetColor(0,180,0, 100);
        } else {
            ofSetColor(250,0,0, 100);
        }
		ofRect(0,0, indicatorW,40);
        ofSetColor(255,255,255);
        sprintf(handyString, "RECORDING f# %d", eyeMovie.getCurrentFrame());
        ofDrawBitmapString(handyString, 10,25);
        
        
        ofSetColor(255,255,255, 100);
        ofRect(0,40,indicatorW,20);
        int nFrames = eyeMovie.getTotalNumFrames() -1;
        float indicatorX = ofMap(eyeMovie.getCurrentFrame(), 0, nFrames-1, 0, indicatorW);
        ofSetColor(0,0,0);
        ofRect(indicatorX-1,40,3,20);
	}
	
	
}


//--------------------------------------------------------------
void testApp::histogramNormalizeC3 (){
	// Normalize the histogram of the main color image. 
	IppiSize frameRoi = {processW, processH}; 
	IppiSize entireRoi = {PROCESS_W, PROCESS_H_MAX}; 
	
	// blank out the images. 
	bool bBlackOutImages = false; 
	if (bBlackOutImages){
		Ipp8u black[] = {0,0,0};
		ippiSet_8u_C3R(black, ipp8u_incoming, processS3, entireRoi);
		ippiSet_8u_C1R(0, ipp8u_Y, processS1, entireRoi);
		ippiSet_8u_C1R(0, ipp8u_U, processS1, entireRoi);
		ippiSet_8u_C1R(0, ipp8u_V, processS1, entireRoi);
	}
	
	// copy over color image
	unsigned char *rawScaledPixelsC3 = rawImageC3.getPixels();
	memcpy (ipp8u_incoming, rawScaledPixelsC3, processW*processH*3);
	
	// Here a straightforward histogram equalization on each one of the 
	// color channels will lead to color shifts. A better way of modifying 
	// histograms for color images is to convert the image into another 
	// color space (like YUV), apply the histogram modification to 
	// the luminance channel and then convert the image back to RGB.
	
	//--------------------
	// convert pixel-order ipp8u_incoming (RGB) into planar ipp8u_YUV (a YUV copy).
	// here's the conversion from pixel-order to planar data.
	Ipp8u *pYUV[3] = {ipp8u_Y, ipp8u_U, ipp8u_V};
	ippiRGBToYUV_8u_C3P3R (ipp8u_incoming, processS3, pYUV, processS1, frameRoi);
	
	bool bJustApplyGammaToLuminance = false;
	if (bJustApplyGammaToLuminance){
		
		int ind;
		for (int y=0; y<processH; y++){
			for (int x=0; x<processW; x++){
				ind = y*processS1 + x;
				float val = ofMap(ipp8u_Y[ind], 0,255, 0,1, true); 
				val = powf(val, 0.5); 
				ipp8u_Y2[ind] = (unsigned char) ((int)ofMap(val,0,1, 0,255)); 
			}
		}
		
	} else {
	
		// equalize the histogram of the luminance, into ipp8u_Y2
		// "amountOfHistogramEqualization" is a slider (0...1) set by the app's control panel. 
		bool bUseGolanMethod = true; // it's good now, I promise!
		if (bUseGolanMethod){
			float amountOfHistogramEqualization = 0.9;
			histEqer->setActiveROI (processW, processH);
			histEqer->modImageByHistogram (ipp8u_Y, ipp8u_Y2, 1, amountOfHistogramEqualization);
			
		} else {
			// OpenCV technique; produces less good results. Currently turned off.
			CvRect histEqRect;
			histEqRect.x 		= 0;
			histEqRect.y 		= 0;
			histEqRect.width	= processW;
			histEqRect.height	= processH;
			cvSetImageData ( iplImage_ipp8uY,  ipp8u_Y,  processS3/3); 
			cvSetImageData ( iplImage_ipp8uY2, ipp8u_Y2, processS3/3); 
			cvSetImageROI  ( iplImage_ipp8uY,  histEqRect);
			cvSetImageROI  ( iplImage_ipp8uY2, histEqRect);
			cvEqualizeHist ( iplImage_ipp8uY,  iplImage_ipp8uY2 );
		}
	}
	
	bool bSharpenFrame = false; 
	if (bSharpenFrame){
		// take the opportunity offered by a copy operation to perform sharpening
		// -- on the luminance channel only! me so smart
		IppiSize filteredFrameRoi = {processW-2, processH-2};
		int filterByteOffset = (processS1 * 1) + 1;
		ippiFilterSharpen_8u_C1R(ipp8u_Y2 		+ filterByteOffset, processS1,
								 ipp8u_Y     	+ filterByteOffset, processS1, filteredFrameRoi);
	} else {
		// Otherwise, just copy the luminance data into the correct (pointed-to) location
		ippiCopy_8u_C1R (ipp8u_Y2, processS1, ipp8u_Y,  processS1, frameRoi);
	}
	
	// copy the histeq'd image back into the luminance channel of ipp8u_HSL; 
	// Conversion from planar-order to pixel-order data.
	// Then copy this into the histogram-equalized RGB version
	ippiYUVToRGB_8u_P3C3R((const Ipp8u**)pYUV, processS1, ipp8u_incoming, processS3, frameRoi);
	histEqImageC3.setFromPixels(ipp8u_incoming, processW, processH);
}



//--------------------------------------------------------------
void testApp::keyPressed(int key){
	
	switch (key){
		case 's':
            if (!bUseMovie){
                XML.saveFile("pupilRecording.xml");
                printf("pupilRecording saved to xml.\n");
                bRecordingXML = false;
            }
			break;
			
		case 'r':
            bRecordingXML = !bRecordingXML; 
            if (bRecordingXML){
                lastTagNumber	= XML.addTag("RECORDING");
            }
			break;
            
            
		case 'l':
			loadImage();
			break;
	
        case '0':
            if (bUseMovie){
                eyeMovie.firstFrame();
            }
			break;
        case 'h':
            // draw everything twice as big, if you want!
            bDrawHalfScale = !bDrawHalfScale;
            break;
            
        case OF_KEY_RIGHT:
            if (bUseMovie){
                getNextFrameImage();
            }
			break;
	}

}


//--------------------------------------------------------------
void testApp::keyReleased(int key){
}
//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){
}
//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){
}
//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
}
//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){
}
//--------------------------------------------------------------
void testApp::windowResized(int w, int h){
}
//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){
}
//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 
}