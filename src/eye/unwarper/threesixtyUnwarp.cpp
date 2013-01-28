#include "threesixtyUnwarp.h"

//--------------------------------------------------------------
threesixtyUnwarp::threesixtyUnwarp(){

	handyString = new char[128];
    blackOpenCV = cvScalarAll(0);

    unwarpedW = 1024;  // target image size
	unwarpedH = 256;
}


//--------------------------------------------------------------
void threesixtyUnwarp::setup (int processW, int processH ){
	
		
    // The "warped" original source video produced by the Bloggie.
	warpedW = processW;
	warpedH = processH;
    warpedCx = warpedW/2;
    warpedCy = warpedH/2;
    angularOffset = 0;
	
	//unwarpedW = 1024;  // target image size
	//unwarpedH = 256;
    
    handyVec2f = new ofVec2f(); 
	
	// Interpolation method: 
	// 0 = CV_INTER_NN, 1 = CV_INTER_LINEAR, 2 = CV_INTER_CUBIC.
	interpMethod = (int) CV_INTER_LINEAR;
	
		
	/*
	// straight rectilinearization
	yWarpA = -0.2047;
	yWarpB =  0.8632;
	yWarpC =  0.3578;
	 */
	yWarpA =   0.1850;
	yWarpB =   0.8184;
	yWarpC =  -0.0028;
    
    maxR = 0.75 * processH /2.0;
    minR = 0.25 * processH /2.0;

    
	int nWarpedBytes = warpedW * warpedH * 3;
	printf("warpedW = %d, warpedH = %d\n", warpedW, warpedH);
	
	warpedImageOpenCV.allocate(warpedW, warpedH);
	warpedPixels = new unsigned char[nWarpedBytes];	
	warpedIplImage = warpedImageOpenCV.getCvImage();
	cvSetImageROI(warpedIplImage, cvRect(0, 0, warpedW, warpedH));
	
	int nUnwarpedPixels = unwarpedW * unwarpedH;
	int nUnwarpedBytes  = unwarpedW * unwarpedH * 3;
	unwarpedImage.allocate(unwarpedW, unwarpedH, OF_IMAGE_COLOR);
	unwarpedPixels = new unsigned char[nUnwarpedBytes];
	unwarpedTexture.allocate(unwarpedW, unwarpedH,GL_RGB);
	
	unwarpedImageOpenCV.allocate(unwarpedW, unwarpedH);
	unwarpedImageOpenCV.setROI(0, 0, unwarpedW, unwarpedH);
	unwarpedIplImage = unwarpedImageOpenCV.getCvImage();
	
	srcxArrayOpenCV.allocate(unwarpedW, unwarpedH);
	srcyArrayOpenCV.allocate(unwarpedW, unwarpedH);
	srcxArrayOpenCV.setROI(0, 0, unwarpedW, unwarpedH);
	srcyArrayOpenCV.setROI(0, 0, unwarpedW, unwarpedH);
	
	xocvdata = (float*) srcxArrayOpenCV.getCvImage()->imageData; 
	yocvdata = (float*) srcyArrayOpenCV.getCvImage()->imageData; 


	computeInversePolarTransform(); 
	
}





//Used for the by hand portion and OpenCV parts of the shootout. 
//For the by hand, use the normal unwarpedW width instead of the step
//For the OpenCV, get the widthStep from the CvImage and use that for quarterstep calculation
//=============================================
void threesixtyUnwarp::computeInversePolarTransform(){

	// we assert that the two arrays have equal dimensions, srcxArray = srcyArray
	float radius, angle;
	float circFactor = 0 - TWO_PI/(float)unwarpedW;
	float difR = maxR-minR;
	int   dstRow, dstIndex;
	
	xocvdata = (float*) srcxArrayOpenCV.getCvImage()->imageData; 
	yocvdata = (float*) srcyArrayOpenCV.getCvImage()->imageData; 
	
	for (int dsty=0; dsty<unwarpedH; dsty++){
		float y = ((float)dsty/(float)unwarpedH);
		float yfrac = y; //yWarpA*y*y + yWarpB*y + yWarpC;
		yfrac = MIN(1.0, MAX(0.0, yfrac)); 

		radius = (yfrac * difR) + minR;
		dstRow = dsty * unwarpedW; 	
		
		for (int dstx=0; dstx<unwarpedW; dstx++){
			dstIndex = dstRow + dstx;
			angle    = ((float)dstx * circFactor) + (DEG_TO_RAD * angularOffset);
			
			xocvdata[dstRow + dstx] = warpedCx + radius*cosf(angle);
			yocvdata[dstRow + dstx] = warpedCy + radius*sinf(angle);
		}
	}
	
	srcxArrayOpenCV.setFromPixels(xocvdata, unwarpedW, unwarpedH);
	srcyArrayOpenCV.setFromPixels(yocvdata, unwarpedW, unwarpedH);
}

//--------------------------------------------------------------
ofVec2f *threesixtyUnwarp::getWarpedXYFromUnwarpedXY (float ux, float uy){
    
    // given ux in 0...unwarpedW and
    // given uy in 0...unwarpedH
    // return wx
    
    ux = ofClamp(ux, 0, unwarpedW-1);
    uy = ofClamp(uy, 0, unwarpedH-1);

    float yfrac = uy / (float)unwarpedH;
    yfrac = MIN(1.0, MAX(0.0, yfrac));

    float difR = maxR-minR;
    float radius = (yfrac * difR) + minR;

    int dstRow = ((int)uy) * unwarpedW;
    int dstIndex = dstRow + ((int)ux);
    
    float circFactor = 0 - TWO_PI/(float)unwarpedW;
    float angle = (ux * circFactor) + (DEG_TO_RAD * angularOffset);
		
    float wx = warpedCx + radius*cosf(angle);
    float wy = warpedCy + radius*sinf(angle);

    handyVec2f->x = wx;
    handyVec2f->y = wy;
    return handyVec2f;
}


//--------------------------------------------------------------
cv::Point2f threesixtyUnwarp::getPt2fWarpedXYFromUnwarpedXY (float ux, float uy){
    
    // given ux in 0...unwarpedW and
    // given uy in 0...unwarpedH
    // return wx
    
    ux = ofClamp(ux, 0, unwarpedW-1);
    uy = ofClamp(uy, 0, unwarpedH-1);
    
    float yfrac = uy / (float)unwarpedH;
    yfrac = MIN(1.0, MAX(0.0, yfrac));
    
    float difR = maxR-minR;
    float radius = (yfrac * difR) + minR;
    
    int dstRow = ((int)uy) * unwarpedW;
    int dstIndex = dstRow + ((int)ux);
    
    float circFactor = 0 - TWO_PI/(float)unwarpedW;
    float angle = (ux * circFactor) + (DEG_TO_RAD * angularOffset);
    
    float wx = warpedCx + radius*cosf(angle);
    float wy = warpedCy + radius*sinf(angle);
    
    cv::Point2f aCvPoint;
    aCvPoint.x = wx;
    aCvPoint.y = wy;
    return aCvPoint;
}


//--------------------------------------------------------------
void threesixtyUnwarp::draw (float x, float y){
	ofSetColor(255,255,255);
	unwarpedImage.draw(x,y);
}


//--------------------------------------------------------------
void threesixtyUnwarp::update ( float newWarpedCx, float newWarpedCy,float innerR, float outerR, unsigned char *inputColorPixels ){
    
    
    minR = innerR;
    maxR = outerR;
    warpedCx = newWarpedCx;
    warpedCy = newWarpedCy;

    computeInversePolarTransform();
    memcpy(warpedPixels, inputColorPixels, warpedW*warpedH*3);
    warpedIplImage->imageData = (char*) warpedPixels; 
    
    cvRemap(warpedIplImage,  unwarpedIplImage, 
            srcxArrayOpenCV.getCvImage(), 
            srcyArrayOpenCV.getCvImage(), 
            interpMethod | CV_WARP_FILL_OUTLIERS, blackOpenCV );
     
    unwarpedPixels = (unsigned char*) unwarpedIplImage->imageData;
    unwarpedImage.setFromPixels(unwarpedPixels, unwarpedW, unwarpedH, OF_IMAGE_COLOR, true);
    unwarpedTexture.loadData(unwarpedPixels, unwarpedW, unwarpedH, GL_RGB);
    

}




