#ifndef _TM_VIS_THRESHOLDERC1_OTSUGLC
#define _TM_VIS_THRESHOLDERC1_OTSUGLC


class tmVisThreshC1_OtsuGrayLevelClass {
	
	public:
	
		tmVisThreshC1_OtsuGrayLevelClass (int nPixels, bool f);
		void 	initialize (int *srcHistogram); // must be done on each frame
    	
    	float 	getMu();
		float	getOmega();
		int  	getThreshold();
		
		void 	addToEnd();
		void 	removeFromBeginning();
	
	private:
	
		int 	N;
   		int 	index;
    	float 	omega;
    	float 	mu;
		float	*probabilityHistogram;
		bool	first;
};

#endif /* _TM_VIS_THRESHOLDERC1_OTSUGLC */