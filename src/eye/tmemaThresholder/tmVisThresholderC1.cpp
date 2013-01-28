#include "tmVisThresholderC1.h"


//=================================================================================
tmVisThresholderC1::tmVisThresholderC1 (int w, int h){
	
	// The current threshold as determined 
	// by the user's choice of algorithm.
	theThreshold 		= 128;
	
	// State storage
	currentMethod		= (ThresholdMethod)0;
	handyString = new char[255];
	sprintf(handyString, "");
	
	// A running history of theThreshold
	theThresholdHistoryLen = 8;
	int threshHistoryStep;
	theThresholdHistory	 = ippiMalloc_8u_C1(theThresholdHistoryLen, 1, &threshHistoryStep); 
	theThreshHistoryCopy = ippiMalloc_8u_C1(theThresholdHistoryLen, 1, &threshHistoryStep);  
	for (int i=0; i<theThresholdHistoryLen; i++){ theThresholdHistory[i] = theThreshold; }
	
	// Pointers to source and destination images.
	// The thresholder owns the destination image.
	buffer_w 			= w;
	buffer_h 			= h;
	buffer_step			= w;
	buffer_roi.width	= buffer_w;
	buffer_roi.height	= buffer_h;
	buffer_n			= buffer_w * buffer_h;
		
	buffer_src 			= NULL;
	buffer_dst			= ippiMalloc_8u_C1(buffer_w, buffer_h, &buffer_step);
	buffer_tmp			= ippiMalloc_8u_C1(buffer_w, buffer_h, &buffer_step); 
	
	// Histogram storage
	histStep32s				= 256;
	nHistogramBuckets 	= 256;
	
	srcHistogram32s 		= ippiMalloc_32s_C1 (nHistogramBuckets, 1, &histStep32s);
	srcHistCopy32s			= ippiMalloc_32s_C1 (nHistogramBuckets, 1, &histStep32s); 
	rawHistogram32s			= ippiMalloc_32s_C1 (nHistogramBuckets, 1, &histStep32s); 
	prvHistogram32s			= ippiMalloc_32s_C1 (nHistogramBuckets, 1, &histStep32s);  
	pLevels32s				= ippiMalloc_32s_C1 (nHistogramBuckets, 1, &histStep32s);
	srcHistCopyA16s			= ippiMalloc_16s_C1 (nHistogramBuckets, 1, &histStep16s); 
	srcHistCopyB16s			= ippiMalloc_16s_C1 (nHistogramBuckets, 1, &histStep16s);
	
	distr					= ippiMalloc_32s_C1 (nHistogramBuckets, 1, &histStep32s);
	transform				= ippiMalloc_32s_C1 (nHistogramBuckets, 1, &histStep32s);
	lutLevels				= ippiMalloc_32s_C1 (nHistogramBuckets, 1, &histStep32s);
	
	for (int i=0; i<nHistogramBuckets; i++){  
		srcHistogram32s[i]		= 0;
		srcHistCopy32s[i]		= 0;
		rawHistogram32s[i]		= 0;
		prvHistogram32s[i]		= 0;
		srcHistCopyA16s[i]		= 0;
		srcHistCopyB16s[i]		= 0;
		
		pLevels32s[i]			= 0;
		distr[i]				= 0;
		transform[i]			= 0;
		lutLevels[i]			= i;
	}
	
	b_threshingCalled = false;
	
	// Special classes used by Otsu and MixtureModelling algorithms.
	otsuC1 = new tmVisThreshC1_OtsuGrayLevelClass (buffer_n, true);
	otsuC2 = new tmVisThreshC1_OtsuGrayLevelClass (buffer_n, false);
	mixMod = new tmVisThreshC1_MixtureModeling2 (buffer_w, buffer_h);
	
	// Additional storage used by Maximum Entropy method.
	normalizedHist32f	= ippiMalloc_32f_C1 (nHistogramBuckets, 1, &histStep32f);
	pT32f 				= ippiMalloc_32f_C1 (nHistogramBuckets, 1, &histStep32f);
	hB32f 				= ippiMalloc_32f_C1 (nHistogramBuckets, 1, &histStep32f);
	hW32f 				= ippiMalloc_32f_C1 (nHistogramBuckets, 1, &histStep32f);
	
}

//=================================================================================
void 	tmVisThresholderC1::setActiveROI (int w, int h){
	buffer_roi.width	= MAX(0, MIN(buffer_w, w));
	buffer_roi.height	= MAX(0, MIN(buffer_h, h));
}
void	tmVisThresholderC1::resetActiveROI (){
	buffer_roi.width	= buffer_w;
	buffer_roi.height	= buffer_h;
}

//=================================================================================
char *tmVisThresholderC1::getCurrentInfo(){
	
	
	switch(currentMethod){
		default:
		case THRESHOLD_METHOD_SIMPLE:
			sprintf(handyString, "Fixed Threshold");
			break;
		case THRESHOLD_METHOD_OTSU:
			sprintf(handyString, "Otsu's Method");
			break;
		case THRESHOLD_METHOD_MAXIMUM_ENTROPY:
			sprintf(handyString, "Maximum Entropy");
			break;
		case THRESHOLD_METHOD_RECURSIVE_ISODATA:
			sprintf(handyString, "Recursive Isodata");
			break;
		case THRESHOLD_METHOD_TRIANGLE_LEFT:
			sprintf(handyString, "Triangle Below Peak");
			break;
		case THRESHOLD_METHOD_TRIANGLE_RIGHT:
			sprintf(handyString, "Triangle Above Peak");
			break;
		case THRESHOLD_METHOD_VALLEY_SEARCH:
			sprintf(handyString, "Valley Search");
			break;
		case THRESHOLD_METHOD_MIXTURE_MODELING:
			sprintf(handyString, "Mixture Modelling");
			break;	
		case THRESHOLD_METHOD_HYBRID:
			sprintf(handyString, "Hybrid");
			break;
	}
	return handyString;
}


//=================================================================================
Ipp32s*	tmVisThresholderC1::computeHistogram (Ipp8u *input, int controlFlags){
	if (input != NULL){
	
	
		// copy the current histogram into the previous histogram.
		for (int i=0; i<nHistogramBuckets; i++){
			prvHistogram32s[i] = srcHistogram32s[i];
		}
		
		
		//---------------
		// Compute the current histogram
		buffer_src = input;
		ippiHistogramEven_8u_C1R (
			buffer_src, buffer_step, buffer_roi, 
			srcHistogram32s, pLevels32s, nHistogramBuckets, 0, 255);
		for (int i=0; i<nHistogramBuckets; i++){ 
			srcHistogram32s[i] = MAX(0,MIN(32767, srcHistogram32s[i]));
		}
		
		//---------------
		// copy (store) the raw histogram
		for (int i=0; i<nHistogramBuckets; i++){ 
			rawHistogram32s[i] = (Ipp16s) srcHistogram32s[i];
		}
		
		//---------------
		// clobber (ramp) the bottom and top end of the histograms, which are causing solarization!
		int rampLen = 3;
		for (int j=0; j<=rampLen; j++){
			float frac = MAX(0,MIN(1, (float)(j)/(float)(rampLen)));
			srcHistogram32s[j]     = (Ipp32s)(frac       * srcHistogram32s[    j]);
			srcHistogram32s[255-j] = (Ipp32s)((1.0-frac) * srcHistogram32s[255-j]);
		}
		

		//---------------
		if (controlFlags & THRESH_CONTROL_HIST_MEDIAN){
			// Copy the histogram into a 16s array for Ipps processing.
			// Make a note of the array's boundary values for later on.
			for (int i=0; i<nHistogramBuckets; i++){ 
				srcHistCopy32s[i] = (Ipp16s) srcHistogram32s[i];
			}
			int valn = srcHistogram32s[nHistogramBuckets-2];
			int val0 = srcHistogram32s[0];
			int medianWidth = 15;
			
			srcHistogram32s[0] = 0; // (this fixes a boundary problem encountered during median).
			ippsFilterMedian_32s (srcHistCopy32s, srcHistogram32s, nHistogramBuckets, medianWidth);
			
			// Restore the boundary values that got clobbered in the median.
			srcHistogram32s[0] = val0;
			srcHistogram32s[nHistogramBuckets-1] = valn;
		}
		
		//---------------
		if (controlFlags & THRESH_CONTROL_HIST_AVERAGE){
			// Blur the histogram binwise to eliminate broad-spectrum noise.
			// Averaging AFTER median filter works better!
			Ipp32s pKernel[3] = {1,1,1};
			
			// average it n times.
			for (int n=0; n<3; n++){
				
				IppiSize histAveragingRoi = {nHistogramBuckets-2, 1};
				for (int i=0; i<nHistogramBuckets; i++){ 
					srcHistCopyA16s[i] = (Ipp16s) srcHistogram32s[i]; 
				}
				ippiFilterRow_16s_C1R (srcHistCopyA16s+1, histStep16s, srcHistCopyB16s+1, histStep16s, histAveragingRoi, pKernel, 3, 1, 3);
				ippiFilterRow_16s_C1R (srcHistCopyB16s+1, histStep16s, srcHistCopyA16s+1, histStep16s, histAveragingRoi, pKernel, 3, 1, 3);
				for (int i=0; i<nHistogramBuckets; i++){ 
					srcHistogram32s[i] = srcHistCopyA16s[i]; 
				}
			}
		}
		
		//---------------
		if (controlFlags & THRESH_CONTROL_HIST_INTEGRATE){
			// evntually use ippiAddWeighted.
			// The histogram is an average of current information and previous information.
			for (int i=0; i<nHistogramBuckets; i++){
				srcHistogram32s[i] = (prvHistogram32s[i] + srcHistogram32s[i])/2;
			}
		}
		
	}
	
	return srcHistogram32s;
}

//=================================================================================
void tmVisThresholderC1::renderHistogram (float x, float y, float w, float h, bool logarithmic, int whichToDraw){
	
	Ipp32s *histToDraw = srcHistogram32s;
	switch(whichToDraw){
		case 0:
			histToDraw = srcHistogram32s;
			break;
		case 1:
			histToDraw = rawHistogram32s;
			break;
		case 2:
			histToDraw = transform;
			break;
		case 3:
			histToDraw = distr;
			break;
	}
	
	
	glPushMatrix();
	glTranslatef(x, y, 0);
	glScalef( (float)w/(float)nHistogramBuckets, 1,1); 
	
	ofNoFill();
	ofSetColor(128,0,0);
	ofRect(0,0,nHistogramBuckets,h);
	
	if (b_threshingCalled){
		glColor3f(1,0,0); 
		glBegin(GL_LINES);
		glVertex2f(theThreshold, 0); 
		glVertex2f(theThreshold, h);
		glEnd();
	}
	
		
	float val;
	if (!logarithmic){
		// linear scaling
		Ipp32s maximumValue; 
		int indexOfMaximum;
		ippsMaxIndx_32s (histToDraw, nHistogramBuckets, &maximumValue, &indexOfMaximum);
		
		glColor3f(0.5,0.5,0.5);
		glBegin(GL_LINES);
		float scale = h/(float)maximumValue;
		for (int i=0; i<nHistogramBuckets; i++){
			val = scale * histToDraw[i];
			glVertex2f(i, 0);
			glVertex2f(i, val);
		}
		glEnd();
		
	} else {
		
		glColor3f(0.5,0.5,0.5);
		glBegin(GL_LINES);
		float scale = h / log(32767.0 + 1.0);
		for (int i=0; i<nHistogramBuckets; i++){
			val = scale * log(1.0+histToDraw[i]);
			glVertex2f(i, 0);
			glVertex2f(i, val);
		}
		glEnd();
	}
	
	
	
	
	glPopMatrix();
}



//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
// Methods for modifying images according to their histogram.

void	tmVisThresholderC1::modImageByHistogram (Ipp8u *input, Ipp8u *buffer_out, int slopeSign, float amount01){
	if (input != NULL){
		buffer_src 	= input;
		
		float constant;
		float intercept;			// histogram 0-intercept
    	float slope;				// histogram slope
    	float z;
    	int maxBinVal = 255;		// maximum value of histogram
    	long ntotal = 0;
		
		// construct histogram, and smooth it.
		int hflags = THRESH_CONTROL_HIST_MEDIAN + THRESH_CONTROL_HIST_AVERAGE + THRESH_CONTROL_HIST_INTEGRATE;
		computeHistogram (buffer_src, hflags);
		
		// construct cumumlative distribution, based on the smoothed histogram
		histdistr (srcHistogram32s, nHistogramBuckets, distr);
		
		// total summation under cumulative distribution omits 
		// 1/2 of lowest and highest occupied bins.
 		ntotal = distr[nHistogramBuckets - 1];
 		
 		// calculate slope 
		switch (slopeSign) {
			case -1:
				slope = -2.0 * ntotal / (float)(maxBinVal * maxBinVal);
				break;
			case 0:
				slope = 0.0;
				break;
			case 1:
				slope = 2.0 * ntotal / (float)(maxBinVal * maxBinVal);
				break;
		}
		
		// compute transform table for histogram equalization or histogram ramp
		switch (slopeSign) {

			case -1:                     
				// negative ramp 
				intercept = ntotal / (float) maxBinVal - slope / 2.0 * maxBinVal;
				constant = intercept / slope;
				for (int i=0; i<=maxBinVal; i++) {
					z = constant * constant + 2.0 / slope * distr[i];
					if (z < 0.0){ z = 0.0; }
					int val = MAX(0, (-constant + 0.5 - sqrt (z)));
					transform[i] = (unsigned char) (val);
				}
				break;

			case 0:
				// histogram equalization 
				constant = ((float) maxBinVal) / (float) ntotal;
				for (int i=0; i<=maxBinVal; i++){
					transform[i] = (unsigned char) (constant * distr[i] + 0.5);
				}
				break;

			case 1: 
				// positive ramp 
				intercept = ntotal / (float) maxBinVal - slope / 2.0 * maxBinVal;
				constant = intercept / slope;
				for (int i=0; i<=maxBinVal; i++) {
					z = constant * constant + 2.0 / slope * distr[i];
					if (z < 0.0){ z = 0.0; }
					int val = MAX(0, (-constant + 0.5 + sqrt (z)));
					transform[i] = (unsigned char) (val);
				}
				break;
		}
		
		// apply amount01, 1=full effect
		for (int i=0; i<=maxBinVal; i++){
			transform[i] = (Ipp32s) round (amount01*(float)transform[i] + (1.0-amount01)*(float)lutLevels[i]);
		}
		
		// transform each pixel value using a LUT.
		// IppiSize roiSize = {buffer_w, buffer_h};
		IppiSize roiSize = {buffer_roi.width, buffer_roi.height};
		ippiLUT_8u_C1R (
			buffer_src, buffer_step, 
			buffer_out, buffer_step, 
			roiSize, transform, lutLevels, nHistogramBuckets);		
	}
}

//----------------------------------------------------------------------
void tmVisThresholderC1::histdistr (Ipp32s *hist, int nBins, Ipp32s *distr) {
	int binLow;        	/* lowest and highest occupied bins */
   int binHigh, bin; 	/* histogram bin incrementor */
  	float total;         /* floating point cumulative distr. */

	// find lowest and highest occupied bins */
	for (binLow=0; binLow < nBins; binLow++)
		if (hist[binLow] != 0)
			break;
	for (binHigh = nBins-1; binHigh >= 0; binHigh--)
		if (hist[binHigh] != 0)
			break;

	// compute cumulative distribution */
	for (bin = 0; bin <= binLow; bin++){
		distr[bin] = 0;
	}

 	total = 0.0;
	for (bin = binLow+1; bin <= binHigh; bin++) {
		total += ((hist[bin-1] + hist[bin]) / 2.0);
		distr[bin] = (long) (total + 0.5);
	}

	for (bin = binHigh+1; bin<nBins; bin++){
	 	distr[bin] = distr[bin-1];
	}
}



//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//
/* The Main thresholding method for the class. 
// It computes the histogram (with settings from controlFlags);
// Obtains an "optimal" threshold (according to the user's choice of method);
// Smoothes/stabilizes this threshold temporally;
// And computes a thresholded image into the buffer_dst.
*/

void tmVisThresholderC1::threshold (Ipp8u *input, Ipp8u *output, ThresholdMethod method, int controlFlags, int offset){
	if ((input != NULL) && (output != NULL)){
	
		buffer_src 		= input;
		currentMethod	= method;
		prvThreshold 	= theThreshold;
		computeHistogram (buffer_src, controlFlags);
		
		// compute an optimal threshold, using some method.
		switch(method){
			default:
			case THRESHOLD_METHOD_SIMPLE:
				theThreshold = 127;
				break;
			case THRESHOLD_METHOD_OTSU:
				theThreshold = getThresholdOtsu (buffer_src);
				break;
			case THRESHOLD_METHOD_MAXIMUM_ENTROPY:
				theThreshold = getThresholdMaximumEntropy (buffer_src);
				break;
			case THRESHOLD_METHOD_RECURSIVE_ISODATA:
				theThreshold = getThresholdIsodata (buffer_src);
				break;
			case THRESHOLD_METHOD_TRIANGLE_LEFT:
				theThreshold = getThresholdTriangle (buffer_src, false);
				break;
			case THRESHOLD_METHOD_TRIANGLE_RIGHT:
				theThreshold = getThresholdTriangle (buffer_src, true);
				break;
			case THRESHOLD_METHOD_VALLEY_SEARCH:
				theThreshold = getThresholdValleySearch (buffer_src, 64);
				break;
			case THRESHOLD_METHOD_MIXTURE_MODELING:
				theThreshold = getThresholdMixtureModeling (buffer_src);
				break;	
			case THRESHOLD_METHOD_HYBRID:
				theThreshold = getThresholdHybrid (buffer_src);
				break;
		}
		
		//-------------------
		theThreshold += offset;
		theThreshold = MIN(255, MAX(0, theThreshold)); 
		
		// Process the threshold to eliminate jitter
		if (controlFlags & THRESH_CONTROL_STABILIZE){
			stabilizeThreshold();
		}
		
		// Actually perform the thresholding.
		ippiThreshold_GTVal_8u_C1R (
			buffer_src, buffer_step,
			buffer_tmp, buffer_step, buffer_roi,
			theThreshold, 255);
		ippiThreshold_LTVal_8u_C1R (
			buffer_tmp, buffer_step,
			output, buffer_step, buffer_roi,
			theThreshold+1, 0);// note the +1
	}
}

//-----------------------------------------------------------------------------------------------------------
Ipp8u* tmVisThresholderC1::threshold (Ipp8u *input, ThresholdMethod method, int controlFlags, int offset){
	
	/*
	// input:			the buffer to process
	// method: 			the method to use (recommend THRESHOLD_METHOD_HYBRID)
	// controlFlags: 	how to process the histogram (recommend THRESH_CONTROL_HIST_AVERAGE + 
																THRESH_CONTROL_HIST_MEDIAN + 
																THRESH_CONTROL_HIST_INTEGRATE + 
																THRESH_CONTROL_STABILIZE);
	// offset:			a user-specified offset in the final thresholding, generally 0.
	*/

	if (input != NULL){
	
		buffer_src 		= input;
		currentMethod	= method;
		prvThreshold 	= theThreshold;
		
		
		computeHistogram (buffer_src, controlFlags);
		
		//-------------------
		// compute an optimal threshold, using some method.
		switch(method){
			default:
			case THRESHOLD_METHOD_SIMPLE:
				theThreshold = 127;
				break;
			case THRESHOLD_METHOD_OTSU:
				theThreshold = getThresholdOtsu (buffer_src);
				break;
			case THRESHOLD_METHOD_MAXIMUM_ENTROPY:
				theThreshold = getThresholdMaximumEntropy (buffer_src);
				break;
			case THRESHOLD_METHOD_RECURSIVE_ISODATA:
				theThreshold = getThresholdIsodata (buffer_src);
				break;
			case THRESHOLD_METHOD_TRIANGLE_LEFT:
				theThreshold = getThresholdTriangle (buffer_src, false);
				break;
			case THRESHOLD_METHOD_TRIANGLE_RIGHT:
				theThreshold = getThresholdTriangle (buffer_src, true);
				break;
			case THRESHOLD_METHOD_VALLEY_SEARCH:
				theThreshold = getThresholdValleySearch (buffer_src, 64);
				break;
			case THRESHOLD_METHOD_MIXTURE_MODELING:
				theThreshold = getThresholdMixtureModeling (buffer_src);
				break;	
			case THRESHOLD_METHOD_HYBRID:
				theThreshold = getThresholdHybrid (buffer_src);
				break;
		}
		
		//-------------------
		// Add in the user-specified offset, and clamp.
		theThreshold += offset;
		theThreshold = MIN(255, MAX(0, theThreshold)); 
		
		//-------------------
		// Process the threshold to eliminate jitter
		if (controlFlags & THRESH_CONTROL_STABILIZE){
			stabilizeThreshold();
		}
		
		//-------------------
		// Actually perform the thresholding.
		
		/*
		// Golan discovered that this IPP double-thresholding method is flawed!
		// It retains gray pixels whose brightness == threshold.
		// Instead the thresholding must be done in two steps.
		ippiThreshold_LTValGTVal_8u_C1R (
			buffer_src, buffer_step,
			buffer_dst, buffer_step, buffer_roi,
			theThreshold, 0,
			theThreshold, 255);
			*/
			
		ippiThreshold_GTVal_8u_C1R (
			buffer_src, buffer_step,
			buffer_tmp, buffer_step, buffer_roi,
			theThreshold, 255);
		
		ippiThreshold_LTVal_8u_C1R (
			buffer_tmp, buffer_step,
			buffer_dst, buffer_step, buffer_roi,
			theThreshold+1, 0);// note the +1
		
	}
	return buffer_dst;
}



//=================================================================================
void tmVisThresholderC1::stabilizeThreshold(){
	// theThreshold is a running average of a running median.
	
	for (int i=(theThresholdHistoryLen-1); i>0; i--){
		theThresholdHistory[i] = theThresholdHistory[i-1];
	}
	theThresholdHistory[0] = theThreshold;
	for (int i=0; i<theThresholdHistoryLen; i++){
		theThreshHistoryCopy[i] = theThresholdHistory[i];
	}
	ippsSortAscend_8u_I(theThreshHistoryCopy, theThresholdHistoryLen);
	theThreshold = theThreshHistoryCopy[(theThresholdHistoryLen/2)-1];
	
	float A = 0.750;
	float B = 1.0-A;
	theThreshold = (int)((A*prvThreshold + B*theThreshold) + 0.5);
}



//=================================================================================
//=================================================================================
//=================================================================================
// A somewhat ad-hoc hybrid method:
// compute the average returned by several of the best other methods;
// discard the one value which is most different from the average;
// use the average of those remaining.
		
int tmVisThresholderC1::getThresholdHybrid (Ipp8u *input){
	if (input != NULL){
		b_threshingCalled = true;
		buffer_src = input;
		
		int to = getThresholdOtsu 				(buffer_src); 
		int te = getThresholdMaximumEntropy 	(buffer_src); 
		int tr = getThresholdIsodata 			(buffer_src);
		int tm = getThresholdMixtureModeling 	(buffer_src);
		
		float average = ((to + te + tr + tm)/4.0);
		int vals[] = {to,te,tr,tm};
		
		float greatestDifference = 0;
		int   idOfTheMostDifferent = -1;
		for (int i=0; i<4; i++){
			float dif = fabs(vals[i] - average);
			if (dif > greatestDifference){
				greatestDifference = dif;
				idOfTheMostDifferent = i;
			}
		}
		float sum=0;
		for (int i=0; i<4; i++){
			if (i != idOfTheMostDifferent){
				sum += vals[i];
			}
		}
		float newAverage = sum/3.0;
		theThreshold = (int)(newAverage + 0.5);
		
	}
	return theThreshold;
}

//-----------------------------------
int tmVisThresholderC1::getMedianOfThreeNumbers (int aa, int bb, int cc){
	int a = aa;
	int b = bb;
	int c = cc;
	int temp;
	
	// If b is less than a, swap a and b
	if (b < a) {
		temp = a;
		a = b;
		b = temp;
	}

    // If c is less than b, swap b and c
	// and see if the new b is less than a
    if (c < b) {
        temp = b;
        b = c;
        c = temp;

        // if (the new) b is less than a, swap a and b
        if (a > b) {
            temp = a;
            a = b;
            b = temp;
        }
    }
    return b;
}




//=================================================================================
//=================================================================================
//=================================================================================
Ipp8u*	tmVisThresholderC1::thresholdUsingValue (Ipp8u *input, int thresh){
	if (input != NULL){
		b_threshingCalled = true;
		buffer_src = input;
		theThreshold = thresh;
		
		ippiThreshold_LTValGTVal_8u_C1R (
			buffer_src, buffer_step,
			buffer_dst, buffer_step, buffer_roi,
			theThreshold, 0,
			theThreshold, 255);
	}
	return buffer_dst;
}	

//=================================================================================
//=================================================================================
//=================================================================================
/**
 *  This algorithm is an implementation of Otsu's thresholding technique 
 *  based on the minimization of inter-class variance [otsu79].
 *
 *  @Article{otsu79,
 *    author =       "N. Otsu",
 *    title =        "A threshold selection method from gray level
 *                    histograms",
 *    journal =      "{IEEE} Trans. Systems, Man and Cybernetics",
 *    year =         "1979",
 *    volume =       "9",
 *    pages =        "62--66",
 *    month =        mar,
 *    keywords =     "threshold selection",
 *    note =         "minimize inter class variance",
 *  }
 *  
 **/

int	tmVisThresholderC1::getThresholdOtsu (Ipp8u *input){
	
	if (input != NULL){
		b_threshingCalled = true;
		buffer_src = input;

		otsuC1->initialize (srcHistogram32s);
		otsuC2->initialize (srcHistogram32s);

		double sigmaMaxd 	= 0;
		double thresholdd 	= 0;
		double sigmad		= 0;
		float fullMu = 	(otsuC1->getOmega()*otsuC1->getMu()) + 
						(otsuC2->getOmega()*otsuC2->getMu());

		for (int i=0; i<255; i++) {
			float dmu1 = otsuC1->getMu() - fullMu;
			float dmu2 = otsuC2->getMu() - fullMu;
			
			sigmad = otsuC1->getOmega()*(dmu1*dmu1) + 
					 otsuC2->getOmega()*(dmu2*dmu2);

			if(sigmad>sigmaMaxd) {
			    sigmaMaxd 	= sigmad;
			    thresholdd 	= otsuC1->getThreshold();
			}
			otsuC1->addToEnd();
			otsuC2->removeFromBeginning();
		}
		theThreshold = (int)(thresholdd + 0.5);// rounding
		
	}
	return theThreshold;
}	



//=================================================================================
//=================================================================================
//=================================================================================
/**
* Automatic thresholding technique based on the entopy of the histogram.
* See: P.K. Sahoo, S. Soltani, K.C. Wong and, Y.C. Chen "A Survey of
* Thresholding Techniques", Computer Vision, Graphics, and Image
* Processing, Vol. 41, pp.233-260, 1988.
*
* @author Jarek Sacha
*/

int tmVisThresholderC1::getThresholdMaximumEntropy (Ipp8u *input){

	if (input != NULL){
		b_threshingCalled = true;
		buffer_src = input;
		theThreshold = maxEntropySplit (srcHistogram32s);
	}
	return theThreshold;
}

//========================================
int tmVisThresholderC1::maxEntropySplit(Ipp32s *hist) {

   	// Normalize histogram, that is make the sum of all bins equal to 1.
	double sum = 0;
	for (int i = 0; i<nHistogramBuckets; ++i) {
		sum += hist[i];
	} if (sum == 0) { 
		return 255;
	}
	
	for (int i=0; i<nHistogramBuckets; i++) {
		normalizedHist32f[i] = hist[i] / sum;
	}

	pT32f[0] = normalizedHist32f[0];
	for (int i=1; i<nHistogramBuckets; i++) {
		pT32f[i] = pT32f[i-1] + normalizedHist32f[i];
	}

	// Entropy for black and white parts of the histogram
	const double epsilon = 5e-324; //Double.MIN_VALUE;
	for (int t=0; t<nHistogramBuckets; t++) {
     
		// Black entropy
		if (pT32f[t] > epsilon) {
			double hhB = 0;
			for (int i = 0; i <= t; i++) {
				if (normalizedHist32f[i] > epsilon) {
					hhB -= normalizedHist32f[i] / pT32f[t] * log(normalizedHist32f[i] / pT32f[t]);
				}
			}
			hB32f[t] = hhB;
		} else {
			hB32f[t] = 0;
		}

		// White  entropy
		double pTW = 1 - pT32f[t];
		if (pTW > epsilon) {
			double hhW = 0;
			for (int i = t+1; i<nHistogramBuckets; ++i) {
				if (normalizedHist32f[i] > epsilon) {
					hhW -= normalizedHist32f[i] / pTW * log(normalizedHist32f[i] / pTW);
				}
			}
			hW32f[t] = hhW;
		} else {
			hW32f[t] = 0;
		}
	}

	// Find histogram index with maximum entropy
	double jMax = hB32f[0] + hW32f[0];
	int tMax = 0;
	for (int t = 1; t < nHistogramBuckets; ++t) {
		double j = hB32f[t] + hW32f[t];
		if (j > jMax) {
			jMax = j;
			tMax = t;
		}
	}

	return tMax;
}




//=================================================================================
//=================================================================================
//=================================================================================
/*
http://www.ph.tn.tudelft.nl/Courses/FIP/frames/fip-Segmenta.html
This iterative technique for choosing a threshold was developed by Ridler and Calvard . 
The histogram is initially segmented into two parts using a starting threshold value such 
as th0 = 2^(B-1), half the maximum dynamic range. The sample mean (mf,0) of the gray values
associated with the foreground pixels and the sample mean (mb,0) of the gray values 
associated with the background pixels are computed. A new threshold value th1 is now computed 
as the average of these two sample means. The process is repeated, based upon the new threshold, 
until the threshold value does not change any more. 
*/

int	tmVisThresholderC1::getThresholdIsodata (Ipp8u *input){

	if (input != NULL){
		b_threshingCalled = true;
		buffer_src = input;
		
		int thresh = 128;
		int tnew = thresh;
		int thr  = 0;
		int sum  = 0;
		int mean1, mean2;
		int ntries = 0;
		do {
			thr = tnew;
			sum = mean1 = mean2 = 0;

			for (int i=0; i<thr; i++){
				mean1 += (srcHistogram32s[i] * i);
				sum   += (srcHistogram32s[i]);
			}     
			if (sum != 0){ mean1 = mean1 / sum;}

			sum = 0;
			for (int i=thr; i<255; i++){
				mean2 += (srcHistogram32s[i] * i);
				sum   += (srcHistogram32s[i]);
			}

			if (sum != 0){ mean2 = mean2 / sum;}
			tnew = (mean1 + mean2) / 2;
			ntries++;

		} while ((tnew != thr) && (ntries < 64));
		theThreshold = tnew;
	}
	return theThreshold;

}



//=================================================================================
//=================================================================================
//=================================================================================
/*
// http://www.ph.tn.tudelft.nl/Courses/FIP/frames/fip-Segmenta.html
A line is constructed between the maximum of the histogram at brightness bmax and 
the lowest value bmin = (p=0)% in the image. The distance d between the line and 
the histogram h[b] is computed for all values of b from b = bmin to b = bmax. 
The brightness value bo where the distance between h[bo] and the line is maximal 
is the threshold value, that is, th = bo. This technique is particularly effective 
when the object pixels produce a weak peak in the histogram.

// NOTE THAT THIS METHOD ASSUMES A PEAK ON THE BRIGHT END.
// we'll have to reverse it otherwise.
*/

int	tmVisThresholderC1::getThresholdTriangle (Ipp8u *input, bool direction){

	if (input != NULL){
		b_threshingCalled = true;
		buffer_src = input;
		int thresh = 128;
		
		// derive the histogram line (a,b),
		// which connects the histogram's tallest point
		// with the histogram's lowest non-zero bucket.
		float ax = 0;
		float ay = 0; 
		float bx = 0; 
		float by = 0; 
		int firstNonzeroBucketIndex = nHistogramBuckets;
		int lastNonzeroBucketIndex  = 0;
		bool bFoundFirstNonzeroBucket = false;
		int largestBucketIndex = 0; 
		int largestBucketValue = 0; 
		for (int i=0; i<nHistogramBuckets; i++){
			if (srcHistogram32s[i] > largestBucketValue){
				largestBucketValue = srcHistogram32s[i];
				largestBucketIndex = i;
			}
			if (srcHistogram32s[i] != 0){
				lastNonzeroBucketIndex = i;
				if (bFoundFirstNonzeroBucket == false){
					bFoundFirstNonzeroBucket = true;
					firstNonzeroBucketIndex = i;
				}
			} 
		}
		
		bx = (float)largestBucketIndex/(float)nHistogramBuckets;
		by = (float)largestBucketValue/32767.0;
		if (direction == false) { // triangle on left side of peak
			ax = (float)firstNonzeroBucketIndex/(float)nHistogramBuckets;
			ay = (float)srcHistogram32s[firstNonzeroBucketIndex]/32767.0;
		} else {
			ax = (float)lastNonzeroBucketIndex/(float)nHistogramBuckets;
			ay = (float)srcHistogram32s[lastNonzeroBucketIndex]/32767.0;
		}
		
		
		// find the histogram point furthest from this line.
		float hx,hy,hd;
		float maxDist = 0;
		for (int i=0; i<nHistogramBuckets; i++){
			// create 2D histogram "coordinates" in (0,1); 
			hx = (float)i/(float)nHistogramBuckets;
			hy = (float)srcHistogram32s[i]/32767.0;
			hd = distanceFromPointToLine(hx,hy, ax,ay,bx,by);
			if (hd > maxDist){
				maxDist = hd;
				thresh = i;
			}
		}
		
		theThreshold = thresh;
	}
	return theThreshold;

}

//----------------------------
// bourke: distance from point to a line
// http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/source.c
float tmVisThresholderC1::distanceFromPointToLine (
	float px, float py,		/* point */
	float ax, float ay,		/* line start */
	float bx, float by){	/* line end */
	
    float abdx = bx-ax;
    float abdy = by-ay;
    float abMagSq = abdx*abdx + abdy*abdy; 
 
   	float U = ( ((px - ax) * abdx) +
      			((py - ay) * abdy) ) / abMagSq;
 
    if ((U < 0.0f) || (U > 1.0f)){
        return 0;   // closest point does not fall within the line segment
    }
 
    float intersectionX = ax + U*abdx;
    float intersectionY = ay + U*abdy;
    float pidx = px - intersectionX;
    float pidy = py - intersectionY;
    float distance = sqrt(pidx*pidx + pidy*pidy);
 
    return distance;
}




//=================================================================================
//=================================================================================
//=================================================================================
/*
Find a valley between two peaks, set the threshold there.
*/

int	tmVisThresholderC1::getThresholdValleySearch (Ipp8u* input, int marginBetweenPeaks){
	
	if (input != NULL){
		b_threshingCalled = true;
		buffer_src = input;
		
		for (int i=0; i<nHistogramBuckets; i++){
			srcHistCopyA16s[i] = (Ipp16s)srcHistogram32s[i];
		}
		
		// obtain the index of the maximum bin
		Ipp16s pMax, pMin; int pIndx;
		ippsMaxIndx_16s (srcHistCopyA16s, nHistogramBuckets, &pMax, &pIndx);
		int 	indexOfMaximum1 = pIndx;
		int 	indexOfMaximum2 = 0;
						
		int 	away = 0;
		Ipp16s	guessMax2a = 0;
		Ipp16s 	guessMax2b = 0;
		int 	guessIndex2a = 0;
		int 	guessIndex2b = 0;
		bool	guess2aExists = false;
		bool	guess2bExists = false;

		// search for a 2nd maximum located at a higher bin than the first.
		away = marginBetweenPeaks + indexOfMaximum1;		
		if (away < nHistogramBuckets){
			ippsMaxIndx_16s (srcHistCopyA16s+away, nHistogramBuckets-away, &guessMax2a, &guessIndex2a);
			guessIndex2a += away;
			guess2aExists = true;
		}
		// search for a 2nd maximum located lower than the first.
		away = indexOfMaximum1 - marginBetweenPeaks;
		if (away > 0){
			ippsMaxIndx_16s (srcHistCopyA16s, away, &guessMax2b, &guessIndex2b);
			guess2bExists = true;
		}
		// take the index of whichever max is larger
		if (guess2aExists && guess2aExists){
			if (guessMax2a > guessMax2b){
				indexOfMaximum2 = guessIndex2a;
			} else {
				indexOfMaximum2 = guessIndex2b;	
			}
		} else if (guess2aExists && !guess2bExists){
			indexOfMaximum2 = guessIndex2a;
		} else if (guess2bExists && !guess2aExists){
			indexOfMaximum2 = guessIndex2b;
		} else {
			indexOfMaximum2 = 0;
		}
		
		// Now search for the minimum between the two maxima.
		// Swap the two maxima, if necessary, so they are sorted.
		if (indexOfMaximum2 < indexOfMaximum1){
			int tmp = indexOfMaximum1;
			indexOfMaximum1 = indexOfMaximum2;
			indexOfMaximum2 = tmp;
		}
		// clamp the maxima appropriately.
		indexOfMaximum2 = MAX(marginBetweenPeaks, indexOfMaximum2);
		indexOfMaximum1 = MIN(nHistogramBuckets - marginBetweenPeaks, indexOfMaximum1);
			
		int start = indexOfMaximum1;
		int nToDo = indexOfMaximum2 - indexOfMaximum1;
		ippsMinIndx_16s(srcHistCopyA16s +start, nToDo, &pMin, &pIndx);
		int indexOfMinimum = pIndx + start;
		indexOfMinimum = MIN(indexOfMinimum, indexOfMaximum2-0);
		indexOfMinimum = MAX(indexOfMinimum, indexOfMaximum1+0);
		
		//--------------
		theThreshold = indexOfMinimum;
	}

	return theThreshold;
}




//=================================================================================
//=================================================================================
//=================================================================================
/*
 * Mixture Modeling algorithm
 * Copyright (c) 2003 by Christopher Mei (christopher.mei@sophia.inria.fr)
 *                    and Maxime Dauphin
 *  This algorithm thresholds the image using a gray-level 
 *  histogram Gaussian characterisation.
 *  
 **/
 
int	tmVisThresholderC1::getThresholdMixtureModeling (Ipp8u *input){
	
	if (input != NULL){
		b_threshingCalled = true;
		buffer_src = input;
		mixMod->index = mixMod->INDEX_MIN-1;
		mixMod->setHistogram (srcHistogram32s);
		
		
		double sigmaMax = 0;
		int tmpThreshold = 0;

		float error = 0;
		float errorMin = 9999999;
		float mu1 = 0;
		float mu2 = 0;
		float variance1 = 0;
		float variance2 = 0;
		
		// Start
		while (mixMod->addToIndex()) {
			error = calculateMMError();

			if(error<errorMin) {
				errorMin = error;
				tmpThreshold = mixMod->getThreshold();
				mu1 = mixMod->getMu1();
				variance1 = mixMod->getVariance1();
				mu2 = mixMod->getMu2();
				variance2 = mixMod->getVariance2();
			}
		}
		mixMod->setIndex(tmpThreshold);
		theThreshold = findMMThreshold ((int)mu1, (int)mu2);
	}
	return theThreshold;
}

//----------------
float tmVisThresholderC1::calculateMMError() {
	float error = 0;
	for (int i=0; i<=(mixMod->MAX_VAL_8u); i++) {
		float dgamma = mixMod->gamma(i) - srcHistogram32s[i];
		error += (dgamma * dgamma);
	}
	return error/(mixMod->MAX_VAL_8u +1);
}

//----------------
int tmVisThresholderC1::findMMThreshold (int mu1, int mu2) {
	float min = 999999999;
	int thr = 0;

	for (int i=mu1; i<mu2; i++) {
	    float val = (float)pow(mixMod->differenceGamma(i),2); 
	    if (min>val) {
			min = val;
			thr = i;
	    }
	}
	return thr;
}
   
