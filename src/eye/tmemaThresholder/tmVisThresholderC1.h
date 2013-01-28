#ifndef _TM_VIS_THRESHOLDER_C1
#define _TM_VIS_THRESHOLDER_C1

/* 
TMEMA DYNAMIC ("OPTIMAL") THRESHOLDING CLASS 
Golan Levin, 16 January 
*/

#include "ofMain.h"
#include "ipp.h"


// support classes for various thresholding methods
#include "tmVisThreshC1_OtsuGrayLevelClass.h"
#include "tmVisThreshC1_MixtureModeling2.h"


#define N_THRESHOLD_METHODS 9
enum ThresholdMethod { 
	THRESHOLD_METHOD_SIMPLE, 
	THRESHOLD_METHOD_OTSU,
	THRESHOLD_METHOD_RECURSIVE_ISODATA,
	THRESHOLD_METHOD_MAXIMUM_ENTROPY,
	THRESHOLD_METHOD_TRIANGLE_LEFT,
	THRESHOLD_METHOD_TRIANGLE_RIGHT,
	THRESHOLD_METHOD_VALLEY_SEARCH,
	THRESHOLD_METHOD_MIXTURE_MODELING,
	THRESHOLD_METHOD_HYBRID
};

enum HistogramFlag { 
	THRESH_CONTROL_HIST_AVERAGE 		= 1,	/* blur the histogram binwise */
	THRESH_CONTROL_HIST_MEDIAN 		= 2, 	/* take a binwise median of the histogram */
	THRESH_CONTROL_HIST_INTEGRATE 	= 4,	/* blur the histogram over time */
	THRESH_CONTROL_STABILIZE		= 8		/* take a running median of the threshold */
};





class tmVisThresholderC1 {
	
	public:
		
		tmVisThresholderC1 (int w, int h);
		
		void		setActiveROI (int w, int h);
		void 		resetActiveROI ();
		void 		threshold 							(Ipp8u *input, Ipp8u *output, ThresholdMethod method, int controlFlags, int offset);
		Ipp8u*	    threshold							(Ipp8u *input, ThresholdMethod method, int controlFlags, int offset);
		Ipp8u*	    thresholdUsingValue 				(Ipp8u *input, int thresh);
		
		int		    getThresholdOtsu 					(Ipp8u *input);
		int		    getThresholdMaximumEntropy 	        (Ipp8u *input); 
		int		    getThresholdIsodata				    (Ipp8u *input);
		int		    getThresholdTriangle				(Ipp8u *input, bool direction);
		int		    getThresholdValleySearch	     	(Ipp8u *input, int marginBetweenPeaks); 
		int		    getThresholdMixtureModeling     	(Ipp8u *input); 
		int 		getThresholdHybrid			    	(Ipp8u *input);
		
		Ipp32s*	    computeHistogram 					(Ipp8u *input, int histFlags);
		void 		renderHistogram 					(float x, float y, float w, float h, bool logarithmic, int whichToDraw);
		
		void		modImageByHistogram			    	(Ipp8u *input, Ipp8u *output, int slopeSign, float amount01); 
		void		histdistr 							(Ipp32s *hist, int nBins, Ipp32s *distr);
		
		int 		theThreshold;
		char		*getCurrentInfo();
		
	private:
		
		Ipp8u 		*buffer_src;
		Ipp8u		*buffer_tmp;
		Ipp8u 		*buffer_dst;
		int 	     buffer_w;
		int 		 buffer_h;
		int 		 buffer_step;
		int 		 buffer_n;		// number of pixels (w*h)
		IppiSize 	 buffer_roi;	// always {w,h};
		
		
		// Used by various histogram-based methods:
		void 			computeSrcHistogram (int histogramOptions);
		Ipp32s		*srcHistogram32s;
		Ipp32s		*srcHistCopy32s;
		Ipp32s		*prvHistogram32s;
		Ipp32s		*rawHistogram32s;
		Ipp32s		*pLevels32s;
		Ipp16s		*srcHistCopyA16s;
		Ipp16s		*srcHistCopyB16s;
		
		// Used by histogram equalizer:
		Ipp32s		*distr; 
		Ipp32s		*transform;
		Ipp32s		*lutLevels;
		
		
		int 			nHistogramBuckets;
		int 			histStep32s;
		int			histStep32f;
		int			histStep16s;
		int			histStep8u;
		
		bool			b_threshingCalled;
		
		// Used by Otsu method:
		tmVisThreshC1_OtsuGrayLevelClass	*otsuC1;
		tmVisThreshC1_OtsuGrayLevelClass	*otsuC2;
		
		// Used by Mixture Modeling method:
		tmVisThreshC1_MixtureModeling2		*mixMod;
		float 		calculateMMError();
		int 			findMMThreshold (int mu1, int mu2);
		
		// Used by Max Entropy method:
		Ipp32f		*normalizedHist32f;
		Ipp32f		*pT32f;
		Ipp32f		*hB32f;
		Ipp32f		*hW32f;
		int maxEntropySplit(Ipp32s *hist);
		
		// Used by Triangle method:
		float distanceFromPointToLine (
			float px, float py,	
			float ax, float ay,
			float bx, float by);
			
		int getMedianOfThreeNumbers (int aa, int bb, int cc);
		
		void 			stabilizeThreshold();
		int			theThresholdHistoryLen;
		Ipp8u*		theThresholdHistory;
		Ipp8u*		theThreshHistoryCopy;
		int			prvThreshold;
		
		ThresholdMethod currentMethod;
		char*			handyString;
};

#endif	



/*

http://www.ph.tn.tudelft.nl/Courses/FIP/frames/fip-Segmenta.html
Background-symmetry algorithm - 
Triangle algorithm - !!!!


    *  P-tile-thresholding
    http://www.icaen.uiowa.edu/~dip/LECTURE/Segmentation1.html#thresh
          o choose a threshold T (based on the image histogram) such that 1/p of the image area has gray values less than T and the rest has gray values larger than T
          o in text segmentation, prior information about the ratio between the sheet area and character area can be used
          o if such a priori information is not available - another property, for example the average width of lines in drawings, etc. can be used - the threshold can be determined to provide the required line width in the segmented image

*/