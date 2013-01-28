#ifndef _BLOGGIE_APP
#define _BLOGGIE_APP

#include "ofMain.h"
#include "ofxXmlSettings.h"
#include "ofTexture.h"
#include "ofxOpenCv.h"
#include "cv.h"

// For OpenFrameworks 0.0061 MacOSX 10.6.3
class threesixtyUnwarp {
	
	public:
		
		//----------------------------------------
		/* standard openFrameworks app stuff */
        threesixtyUnwarp ();
    
		void setup(int processW, int processH);
        void update ( float newWarpedCx, float newWarpedCy, float innerR, float outerR, unsigned char *inputPixels );
		void draw(float x, float y);
    
        ofVec2f *handyVec2f;
        ofVec2f *getWarpedXYFromUnwarpedXY (float ux, float uy);
        cv::Point2f getPt2fWarpedXYFromUnwarpedXY (float ux, float uy);
		
		//----------------------------------------
		/* Panoramic unwarp stuff */
		char *handyString;
	
		void computeInversePolarTransform();
		 
		bool testMouseInPlayer();

		ofImage unwarpedImage;
		ofxCvColorImage	warpedImageOpenCV;
		ofxCvColorImage unwarpedImageOpenCV;
		ofxCvFloatImage srcxArrayOpenCV; 
		ofxCvFloatImage srcyArrayOpenCV; 
	
		unsigned char *warpedPixels;
		unsigned char *unwarpedPixels;
		
		int   warpedW;
		int   warpedH;
		float unwarpedW;
		float unwarpedH;
		float warpedCx;
		float warpedCy;
		float savedWarpedCx;
		float savedWarpedCy;
		float savedAngularOffset;
		float angularOffset;
	
		float maxR;
		float minR;
		int   interpMethod; 
		float playerScaleFactor;

		unsigned char *blackColor;
		CvScalar	blackOpenCV;
		IplImage	*warpedIplImage;
		IplImage	*unwarpedIplImage;
	
		float *xocvdata;
		float *yocvdata;
	
		float yWarpA; // for parabolic fit for Y unwarping
		float yWarpB;
		float yWarpC;
	
		//-----------------------------------
		/* For the texture-mapped cylinder */
		ofTexture unwarpedTexture;
		int   cylinderRes;
		float *cylinderX;
		float *cylinderY;
		float cylinderWedgeAngle;
		float blurredMouseX;
		float blurredMouseY;
		
};

#endif	

