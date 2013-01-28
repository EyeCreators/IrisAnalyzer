#pragma once

#include "ofMain.h"
#include "ipp.h"

#include "ofxCv.h"
#include "ofxOpenCv.h"
#include "ofxXmlSettings.h"
#include "ofxUI.h"


#include "tmVisThresholderC1.h"
#include "PupilContourAnalyzer.h"
#include "threesixtyUnwarp.h"


#define	INVALID_BIOMETRIC -1

#define PROCESS_W		1024
#define	PROCESS_H_MAX	1024 // maximum height, actually. 

using namespace ofxCv;
using namespace cv;



struct key_point_response_compare{
    inline bool operator() (const cv::KeyPoint& struct1, const cv::KeyPoint& struct2){
        return (struct1.response < struct2.response);
    }
};

struct key_point_x_compare {
    inline bool operator() (const cv::KeyPoint& struct1, const cv::KeyPoint& struct2){
        return (struct1.pt.x < struct2.pt.x);
    }
};



class testApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();
		void exit(); 
		void loadImage();
	
        bool bUseMovie;
        bool bSourceSuccessfullyLoaded;
        bool bExitWhenFinished;
        bool bShowThings;
        bool bDrawHalfScale; 
		ofVideoPlayer 	eyeMovie;
    
        

		void getNextFrameImage();
		void loadEyeMovie();
	
		ofxUICanvas *gui;
		void guiEvent(ofxUIEventArgs &e);
		float sliderFloat01;

		void keyPressed  (int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
    
    bool bIsTheEyeThere;
    bool bReportDeadLabels; 
	
		ofxCvColorImage		rawImageC3;
		ofxCvColorImage		histEqImageC3;
		ofxCvGrayscaleImage handyOfxCvC1;
		ofxCvGrayscaleImage pupilSearchMaskC1;
		ofxCvGrayscaleImage pupilSearchImageC1;
		ofxCvGrayscaleImage pupilSearchImageThreshedC1;
        ofxCvGrayscaleImage highlightSearchImageC1;
		tmVisThresholderC1	*histEqer;
		tmVisThresholderC1	*pupilThresholder;
        tmVisThresholderC1  *unwarpEqer;
    
        ofxCvGrayscaleImage pupilSearchImagePrevC1;
        ofxCvGrayscaleImage pupilSearchImageDiffC1;
    
        ofxCvGrayscaleImage irisSearchImageC1;
        ofxCvGrayscaleImage irisSearchMaskC1;
        bool bAllocatedSomeImages;
    
		int					processW;
		int					processH;
		int					processS1;
		int					processS3;
		ofxCvContourFinder	pupilContourFinder;
        ofxCvContourFinder	highlightContourFinder;
		ofPoint				pupilPointPx;
		ofPixels			pupilSearchMaskPix;
        ofPixels			irisSearchMaskPix;
        float               irisRadius;
    
    
        // processing the unwarped iris
        ofxCvColorImage         unwarpedImageOfxCvC3;
        ofxCvGrayscaleImage     unwarpedImageOfxCvC1;
        ofxCvGrayscaleImage     unwarpedImageHistEqC1;
        Ipp8u                   *ipp8uUnwarpedHistEqC1;
        ofPixels                HistEqC1Pix;
    
    
    cv::Mat unwarpDstMat;
    cv::Mat unwarpDst_norm;
    cv::Mat unwarpDst_norm_scaled;
    
    //-------------------------
    GoodFeaturesToTrackDetector *goodFeatureTracker;
    vector<KeyPoint> gfttKeyPoints;
    int featureSearchW; 
    int featureSearchH;
    
    bool bComputeFeatures; 
    vector<KeyPoint> surfKeyPoints;
    vector<KeyPoint> bestSurfKeyPoints;
    vector<KeyPoint> keepSurfKeyPoints;
    cv::SurfFeatureDetector *SFD;
    int nSurfPointsToKeep;
    int nSurfPointsToReport;
    ofxCvGrayscaleImage surfFeatureAccumImage;
    void findFeatures();
    
    vector<cv::Point2f> cartesianSurfFeatureCvPoints;
    vector<cv::Point2f> cartesianGfttFeatureCvPoints;
    int maxXmlSurfFeatures;
    int maxXmlGfttFeatures;
    
    PointTracker surfTracker;
    PointTracker gfttTracker;
    int maxNSurfTrackerLabels; 
    int maxNGfttTrackerLabels;
    
	//-------------------------
		PupilContourAnalyzer *pupilContourAnalyzer;
		PupilContourAnalyzer *irisContourAnalyzer;
		void process ();
		void histogramNormalizeC3 ();
    
        void findHighlights();
        int maxNHighlights = 0;
        ofVec3f getColorAt (float x, float y, ofxCvColorImage colImg);
    
    threesixtyUnwarp *unwarper;
    void clobberHighlights();
    void analyzeIris();
	
		float pupilPx;
		float pupilPy;
		float pupilPr; 
		int pupilThreshold;
    float irisPx;
    float irisPy;
    float irisPr;
    
    float pupilPrevPx;
    float pupilPrevPy;
    
    
		
		bool bVerbose;
		float analysisScaleRatio; 
	
        ofxXmlSettings propertiesXML;
		ofxXmlSettings XML;
        void addDataToXMLInformaticReport();
		bool bRecordingXML;
		int lastTagNumber;
        int currentFrameTag;
        bool bExportedInformaticsFile;
		
	
		char handyString[64];
        
    
        
	
		Ipp8u			*ipp8u_incoming;	// storage for incoming frame
		Ipp8u			*ipp8u_processed;	// the processed image;
		Ipp8u			*ipp8u_Y;			// Y of incoming
		Ipp8u			*ipp8u_U;			// U of incoming
		Ipp8u			*ipp8u_V;			// V of incoming
		Ipp8u			*ipp8u_Y2;			// temp;
		Ipp8u			*ipp8u_HandyC1;		// temp
		IplImage		*iplImage_ipp8uY;
		IplImage		*iplImage_ipp8uY2;
	
};
